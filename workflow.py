#!/usr/bin/env python3
"""amplicon-seq-workflow driver: 16S / ITS amplicon reads -> ASVs -> taxonomy -> tree -> phyloseq -> analysis.

One YAML config per project/run (see config/). The driver resolves it per locus (profile, {placeholders},
locus settings) into <results>/<key>.config.json, which the R and shell steps read, and runs the steps:

  reference  ITS: download the UNITE release if missing (get-unite.sh); 16S: check the SILVA files
  dada2      dada2-asv.R: primer check, cutadapt (ITS), filtering, denoising, merging, chimeras, taxonomy
  tree       16S: align-asv-v4.sh (mothur V4 alignment + FastTree); skipped when tree.method is none
  phyloseq   dada2-asv.R --resume=outputs: rebuild the phyloseq objects with the tree
  analyze    analyze-phyloseq.R: barplots, alpha diversity, Venn, NMDS/PERMANOVA, source sharing, distances, betapart

usage:
  workflow.py resolve -c CONFIG [--locus 16S]          print the resolved config
  workflow.py check   -c CONFIG                        check inputs, references and tools
  workflow.py run     -c CONFIG [--locus 16S ITS] [--steps dada2,tree,...] [--threads N] [--nsamples N]
                      [--merge | --no-merge] [--pdf] [--slurm] [--dry-run] [-- extra dada2-asv.R options]
"""
import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")
STEPS = ["reference", "dada2", "tree", "phyloseq", "analyze"]
PLACEHOLDER = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}")


def load_yaml(path):
    """PyYAML when installed, else R's yaml package (R is needed by the workflow anyway)."""
    try:
        import yaml
        with open(path) as fh:
            return yaml.safe_load(fh)
    except ImportError:
        expr = ("cat(jsonlite::toJSON(yaml::read_yaml(commandArgs(TRUE)[1]), auto_unbox = TRUE, null = 'null', "
                "digits = NA))")
        out = subprocess.run(["Rscript", "-e", expr, path], check=True, capture_output=True, text=True).stdout
        return json.loads(out)


def expand(value, variables):
    """Fill {placeholders} in every string of value (nested dicts/lists); repeat so vars can use vars."""
    if isinstance(value, dict):
        return {k: expand(v, variables) for k, v in value.items()}
    if isinstance(value, list):
        return [expand(v, variables) for v in value]
    if not isinstance(value, str):
        return value
    for _ in range(10):
        new = PLACEHOLDER.sub(lambda m: str(variables[m.group(1)]) if m.group(1) in variables else m.group(0), value)
        if new == value:
            break
        value = new
    return value


def merged(base, override):
    """Recursive dict merge, override wins."""
    out = dict(base)
    for k, v in (override or {}).items():
        out[k] = merged(out[k], v) if isinstance(out.get(k), dict) and isinstance(v, dict) else v
    return out


def pick_profile(cfg, name):
    profiles = cfg.get("profiles") or {}
    if name:
        if name not in profiles:
            sys.exit(f"profile '{name}' not in config (profiles: {', '.join(profiles) or 'none'})")
        return name, profiles[name]
    for pname, prof in profiles.items():
        if not prof.get("detect") or os.path.exists(prof["detect"]):
            return pname, prof
    if profiles:
        sys.exit("no profile matches this machine (detect paths: "
                 + ", ".join(str(p.get("detect")) for p in profiles.values()) + "); use --profile")
    return None, {}


def resolve(cfg, locus, profile=None, merge=None):
    """Per-locus settings with every placeholder filled in; this is what the R/shell steps read."""
    if locus not in (cfg.get("locus_settings") or {}):
        sys.exit(f"locus {locus} has no locus_settings in the config")
    pname, prof = pick_profile(cfg, profile)
    variables = {"project": cfg["project"], "run": cfg["run"], "run_id": cfg.get("run_id", cfg["run"]), "locus": locus}
    variables.update(cfg.get("vars") or {})
    variables.update(prof.get("vars") or {})
    variables = {k: expand(v, variables) for k, v in variables.items()}   # vars may refer to vars
    params = merged(cfg.get("dada2") or {}, cfg["locus_settings"][locus])
    if merge is not None:
        params["merge"] = merge
    key = expand(cfg["key"], variables) + ("" if params.get("merge", True) else ".R1")
    res = {
        "project": cfg["project"], "run": cfg["run"], "run_id": variables["run_id"], "locus": locus,
        "key": key, "profile": pname,
        "paths": expand(cfg["paths"], variables),
        "samples": expand(cfg["samples"], variables),
        "params": expand(params, variables),
        "analysis": expand(cfg.get("analysis") or {}, variables),
    }
    return res


def write_resolved(res):
    os.makedirs(res["paths"]["results"], exist_ok=True)
    path = os.path.join(res["paths"]["results"], res["key"] + ".config.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=2)
        fh.write("\n")
    return path


def check(res):
    """Problems that would stop a run, as a list of messages (empty = ready)."""
    problems = []
    p, s, prm = res["paths"], res["samples"], res["params"]
    if not os.path.isfile(p["metadata"]):
        problems.append(f"metadata not found: {p['metadata']}")
    else:
        with open(p["metadata"]) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idcol = header.index(s["id_column"]) if s["id_column"] in header else None
            if idcol is None:
                problems.append(f"metadata has no column {s['id_column']} (columns: {', '.join(header)})")
            else:
                ids = [line.rstrip("\n").split("\t")[idcol] for line in fh if line.strip()]
                rename = s.get("rename") or {}
                names = [re.sub(rename["pattern"], rename["replacement"], i) if rename else i for i in ids]
                missing = [f for n in names for f in (s["r1"], s["r2"])
                           if not os.path.isfile(os.path.join(p["raw"], f.replace("{sample}", n)))]
                if missing:
                    problems.append(f"{len(missing)} fastq files missing in {p['raw']} (first: {missing[0].replace('{sample}', names[0])})"
                                    if len(missing) < 2 * len(names) else f"no fastq files in {p['raw']}")
    tax = prm.get("taxonomy") or {}
    if tax.get("db") and not os.path.isfile(tax["db"]) and not tax.get("unite_doi"):
        problems.append(f"taxonomy reference not found: {tax['db']}")
    tools = ["Rscript"] + (["cutadapt"] if prm.get("run_cutadapt") else [])
    if (prm.get("tree") or {}).get("method") == "mothur-v4":
        tools += [os.environ.get("MOTHUR", "mothur"), os.environ.get("FASTTREE", "FastTreeMP")]
    problems += [f"tool not on PATH: {t}" for t in tools if shutil.which(t) is None]
    return problems


def commands(res, conf_json, steps, args):
    """(step, argv, env) for each step to run."""
    p, prm, key = res["paths"], res["params"], res["key"]
    rscript = ["Rscript"]
    threads = ["--threads=" + str(args.threads)] if args.threads else []
    cmds = []
    for step in steps:
        if step == "reference":
            tax = prm.get("taxonomy") or {}
            if tax.get("unite_doi") and not os.path.isfile(tax["db"]):
                cmds.append((step, ["bash", os.path.join(SCRIPTS, "get-unite.sh"), tax["unite_dir"], str(tax["unite_doi"])], {}))
        elif step == "dada2":
            extra = (["--nsamples=" + str(args.nsamples)] if args.nsamples else []) + list(args.extra)
            cmds.append((step, rscript + [os.path.join(SCRIPTS, "dada2-asv.R"), "--config=" + conf_json] + threads + extra, {}))
        elif step == "tree":
            tree = prm.get("tree") or {}
            if tree.get("method") == "mothur-v4":
                env = {"REFDIR": tree["reference"], "V4START": str(tree["v4_start"]), "V4END": str(tree["v4_end"])}
                cmds.append((step, ["bash", os.path.join(SCRIPTS, "align-asv-v4.sh"),
                                    os.path.join(p["results"], key + ".ASV.fa"), os.path.join(p["results"], "v4align"),
                                    str(args.threads or 8)], env))
        elif step == "phyloseq":
            if (prm.get("tree") or {}).get("method", "none") != "none":
                cmds.append((step, rscript + [os.path.join(SCRIPTS, "dada2-asv.R"), "--config=" + conf_json,
                                              "--resume=outputs"] + threads, {}))
        elif step == "analyze":
            cmds.append((step, rscript + [os.path.join(SCRIPTS, "analyze-phyloseq.R"), "--config=" + conf_json]
                         + (["--pdf"] if args.pdf else []), {}))
    return cmds


def slurm_script(cfg, args, locus):
    sl = cfg.get("slurm") or {}
    argv = ["python3", os.path.abspath(__file__), "run", "-c", os.path.abspath(args.config), "--locus", locus,
            "--steps", ",".join(args.steps)]
    argv += (["--profile", args.profile] if args.profile else []) + (["--nsamples", str(args.nsamples)] if args.nsamples else [])
    argv += ["--merge"] if args.merge is True else ["--no-merge"] if args.merge is False else []
    argv += ["--pdf"] if args.pdf else []
    argv += ["--"] + list(args.extra) if args.extra else []
    lines = ["#!/bin/bash"]
    for opt, k in [("-A", "account"), ("-C", "constraint"), ("-q", "qos"), ("-N", "nodes"), ("-c", "cpus"), ("-t", "time")]:
        if sl.get(k) is not None:
            lines.append(f"#SBATCH {opt} {sl[k]}")
    lines += [f"#SBATCH -J asw-{cfg['run_id'] if 'run_id' in cfg else cfg['run']}-{locus}",
              "#SBATCH -o %x.%j.log", "set -euo pipefail", (sl.get("setup") or "").rstrip(),
              "export THREADS=${SLURM_CPUS_PER_TASK:-4}",
              " ".join(shlex.quote(a) for a in argv) + " --threads $THREADS"]
    return "\n".join(l for l in lines if l) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ["resolve", "check", "run"]:
        sp = sub.add_parser(name)
        sp.add_argument("-c", "--config", required=True)
        sp.add_argument("--locus", nargs="+", help="default: loci in the config")
        sp.add_argument("--profile", help="default: first profile whose detect path exists")
        mg = sp.add_mutually_exclusive_group()
        mg.add_argument("--merge", dest="merge", action="store_true", default=None, help="denoise R1 + R2 and merge")
        mg.add_argument("--no-merge", dest="merge", action="store_false", help="R1 only")
        if name == "run":
            sp.add_argument("--steps", default=",".join(STEPS), type=lambda s: s.split(","),
                            help="comma-separated, default: " + ",".join(STEPS))
            sp.add_argument("--threads", type=int)
            sp.add_argument("--nsamples", type=int, help="test run on the first N samples (dada2 step)")
            sp.add_argument("--pdf", action="store_true", help="analysis figures also as PDF")
            sp.add_argument("--slurm", action="store_true", help="submit one SLURM job per locus")
            sp.add_argument("--dry-run", action="store_true", help="print the commands only")
            sp.add_argument("extra", nargs="*", help="after --: extra dada2-asv.R options, e.g. --resume")
    args = ap.parse_args()

    cfg = load_yaml(args.config)
    loci = args.locus or cfg.get("loci") or list(cfg["locus_settings"])
    if args.cmd == "run":
        bad = [s for s in args.steps if s not in STEPS]
        if bad:
            sys.exit(f"unknown steps: {', '.join(bad)} (steps: {', '.join(STEPS)})")
        args.steps = [s for s in STEPS if s in args.steps]   # always in pipeline order
    status = 0
    for locus in loci:
        res = resolve(cfg, locus, args.profile, args.merge)
        if args.cmd == "resolve":
            print(json.dumps(res, indent=2))
            continue
        if args.cmd == "check":
            problems = check(res)
            print(f"{locus} ({res['key']}, profile {res['profile']}): " + ("ok" if not problems else ""))
            for msg in problems:
                print("  " + msg)
            status |= bool(problems)
            continue
        if args.slurm:
            script = slurm_script(cfg, args, locus)
            path = os.path.join(res["paths"]["work"], f"{res['key']}.sbatch")
            if args.dry_run:
                print(f"# {path}\n{script}")
                continue
            os.makedirs(res["paths"]["work"], exist_ok=True)
            with open(path, "w") as fh:
                fh.write(script)
            subprocess.run(["sbatch", path], check=True, cwd=res["paths"]["work"])
            continue
        conf_json = os.path.join(res["paths"]["results"], res["key"] + ".config.json")
        if not args.dry_run:
            conf_json = write_resolved(res)
        for step, argv, env in commands(res, conf_json, args.steps, args):
            print(f"[{locus}] {step}: " + " ".join(f"{k}={shlex.quote(v)}" for k, v in env.items())
                  + (" " if env else "") + " ".join(shlex.quote(a) for a in argv), flush=True)
            if not args.dry_run:
                rc = subprocess.run(argv, env={**os.environ, **env}).returncode
                if rc != 0:
                    sys.exit(f"[{locus}] step {step} failed (exit {rc})")
    sys.exit(status)


if __name__ == "__main__":
    main()
