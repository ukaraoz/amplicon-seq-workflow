#!/usr/bin/env python3
"""amplicon-seq-workflow driver: 16S / ITS amplicon reads -> ASVs -> taxonomy -> tree -> phyloseq -> analysis.

One YAML config per project/run (see config/). The driver resolves it per locus (profile, {placeholders},
locus settings) into <results>/<key>.config.json, which the R and shell steps read, and runs the steps:

\b
  reference  ITS: download the UNITE release and the guild databases if missing (get-unite.sh, get-guild-dbs.py);
             16S: check the SILVA files
  dada2      dada2-asv.R: primer check, cutadapt (ITS), filtering, denoising, merging, chimeras, taxonomy
  tree       16S: align-asv-v4.sh (mothur V4 alignment + FastTree); skipped when tree.method is none
  phyloseq   dada2-asv.R --resume=outputs: rebuild the phyloseq objects with the tree
  guilds     ITS: fungal-guilds.R: FungalTraits lifestyles and FUNGuild guilds per ASV; skipped without `guilds`
  figures    plot-readtracking.R: read tracking per sample by run -> <results>/figures/
  analyze    analyze-phyloseq.R: barplots, alpha diversity, Venn, NMDS/PERMANOVA, source sharing, distances, betapart

\b
usage:
  workflow.py resolve -c CONFIG [--locus 16S]          print the resolved config
  workflow.py check   -c CONFIG                        check inputs, references and tools
  workflow.py run     -c CONFIG [--locus 16S,ITS] [--steps dada2,tree,...] [--threads N] [--nsamples N]
                      [--merge | --no-merge] [--pdf] [--slurm] [--dry-run] [-- extra dada2-asv.R options]
  workflow.py pull    -c CONFIG [--locus 16S] [--dry-run]    copy results/<key>/ from the remote clone
"""
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from types import SimpleNamespace

import click

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")
STEPS = ["reference", "dada2", "tree", "phyloseq", "guilds", "figures", "analyze"]
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


def resolve(cfg, locus, profile=None, merge=None, repo=None):
    """Per-locus settings with every placeholder filled in; this is what the R/shell steps read.
    repo: the repository folder for {repo} (default: this one; `pull` passes the remote clone's)."""
    if locus not in (cfg.get("locus_settings") or {}):
        sys.exit(f"locus {locus} has no locus_settings in the config")
    pname, prof = pick_profile(cfg, profile)
    # {repo}: this repository (results go to {repo}/results/{key} by default); {dataset}: {project}_{run_id}, the
    # dataset's folder config/{dataset}/ (its metadata/); {key} is added once it is known
    run_id = cfg.get("run_id", cfg["run"])
    variables = {"project": cfg["project"], "run": cfg["run"], "run_id": run_id, "locus": locus,
                 "repo": repo or HERE, "dataset": f"{cfg['project']}_{run_id}"}
    variables.update(cfg.get("vars") or {})
    variables.update(prof.get("vars") or {})
    variables = {k: expand(v, variables) for k, v in variables.items()}   # vars may refer to vars
    params = merged(cfg.get("dada2") or {}, cfg["locus_settings"][locus])
    if merge is not None:
        params["merge"] = merge
    key = expand(cfg["key"], variables) + ("" if params.get("merge", True) else ".R1")
    variables["key"] = key
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


FASTQ_SUFFIX = re.compile(r"(_L[0-9]{3})?_R1(_[0-9]{3})?\.(fastq|fq)(\.gz)?$")


def read_manifest(path, id_column):
    """samples.manifest: tab-separated with a header, one line per sample, columns R1 and R2 (relative paths are
    relative to the manifest's folder) and optionally the sample name (column id_column or "sample"; without one,
    the R1 file name minus _L001_R1_001.fastq.gz). Same rules as dada2-asv.R. Returns [(sample, r1, r2)]."""
    with open(path) as fh:
        header = [h.strip() for h in fh.readline().rstrip("\n").split("\t")]
        rows = [[c.strip() for c in line.rstrip("\n").split("\t")] for line in fh if line.strip()]
    lower = [h.lower() for h in header]
    col = lambda *names: next((lower.index(n.lower()) for n in names if n.lower() in lower), None)
    r1, r2, sid = col("R1"), col("R2"), col(id_column, "sample")
    if r1 is None or r2 is None:
        raise ValueError(f"manifest needs columns R1 and R2: {path} (columns: {', '.join(header)})")
    base = os.path.dirname(os.path.abspath(path))
    full = lambda f: os.path.join(base, os.path.expanduser(f))
    out = []
    for row in rows:
        row += [""] * (len(header) - len(row))
        if not (row[r1] or row[r2]):
            continue
        if not (row[r1] and row[r2]):
            raise ValueError(f"manifest line with only one of R1 / R2: {path}: {row}")
        name = row[sid] if sid is not None else FASTQ_SUFFIX.sub("", os.path.basename(row[r1]))
        out.append((name, full(row[r1]), full(row[r2])))
    return out


def check_manifest(path, id_column, metadata_ids, run_column=None, run_settings=None):
    if not os.path.isfile(path):
        return [f"fastq manifest not found: {path}"]
    try:
        rows = read_manifest(path, id_column)
    except ValueError as e:
        return [str(e)]
    problems = []
    if run_column:
        with open(path) as fh:
            header = [h.strip() for h in fh.readline().rstrip("\n").split("\t")]
            if run_column not in header:
                problems.append(f"manifest has no run column {run_column}")
            else:
                i = header.index(run_column)
                runs = {line.rstrip("\n").split("\t")[i].strip() for line in fh if line.strip()}
                unknown = sorted(set(run_settings or {}) - runs)
                if unknown:
                    problems.append(f"run settings for runs not in the manifest: {', '.join(unknown)}")
    names = [r[0] for r in rows]
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        problems.append(f"duplicated sample names in manifest: {', '.join(dups[:10])}")
    nometa = [n for n in names if n not in metadata_ids]
    if nometa:
        problems.append(f"{len(nometa)} manifest samples not in the metadata (first: {nometa[0]})")
    missing = [f for r in rows for f in r[1:] if not os.path.isfile(f)]
    if missing:
        problems.append(f"{len(missing)} fastq files in the manifest missing (first: {missing[0]})")
    return problems


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
                if s.get("manifest"):
                    problems += check_manifest(s["manifest"], s["id_column"], set(names), s.get("run_column"),
                                               res["params"].get("runs"))
                else:
                    missing = [f for n in names for f in (s["r1"], s["r2"])
                               if not os.path.isfile(os.path.join(p["raw"], f.replace("{sample}", n)))]
                    if missing:
                        problems.append(f"{len(missing)} fastq files missing in {p['raw']} (first: {missing[0].replace('{sample}', names[0])})"
                                        if len(missing) < 2 * len(names) else f"no fastq files in {p['raw']}")
    tax = prm.get("taxonomy") or {}
    if tax.get("db") and not os.path.isfile(tax["db"]) and not tax.get("unite_doi"):
        problems.append(f"taxonomy reference not found: {tax['db']}")
    gu = prm.get("guilds") or {}
    for k in ("fungaltraits", "funguild"):
        if gu and not os.path.isfile(gu.get(k, "")):
            problems.append(f"guild table missing: {gu.get(k)} (the reference step downloads it into guilds.dir)")
    tools = ["Rscript"] + (["cutadapt"] if prm.get("run_cutadapt") else [])
    if (prm.get("tree") or {}).get("method") == "mothur-v4":
        tools += [os.environ.get("MOTHUR", "mothur"), os.environ.get("FASTTREE", "FastTreeMP")]
    problems += [f"tool not on PATH: {t}" for t in tools if shutil.which(t) is None]
    return problems


def commands(res, conf_json, steps, args):
    """(step, argv, env) for each step to run."""
    p, prm, key = res["paths"], res["params"], res["key"]
    rscript = [shutil.which("Rscript") or "Rscript"]   # the same R for every step (tree relabelling too)
    threads = ["--threads=" + str(args.threads)] if args.threads else []
    cmds = []
    for step in steps:
        if step == "reference":
            tax = prm.get("taxonomy") or {}
            if tax.get("unite_doi") and not os.path.isfile(tax["db"]):
                cmds.append((step, ["bash", os.path.join(SCRIPTS, "get-unite.sh"), tax["unite_dir"], str(tax["unite_doi"])], {}))
            gu = prm.get("guilds") or {}
            if gu and not all(os.path.isfile(gu[k]) for k in ("fungaltraits", "funguild")):
                cmds.append((step, [sys.executable, os.path.join(SCRIPTS, "get-guild-dbs.py"), gu["dir"]], {}))
        elif step == "dada2":
            extra = (["--nsamples=" + str(args.nsamples)] if args.nsamples else []) + list(args.extra)
            cmds.append((step, rscript + [os.path.join(SCRIPTS, "dada2-asv.R"), "--config=" + conf_json] + threads + extra, {}))
        elif step == "tree":
            tree = prm.get("tree") or {}
            if tree.get("method") == "mothur-v4":
                env = {"REFDIR": tree["reference"], "V4START": str(tree["v4_start"]), "V4END": str(tree["v4_end"]),
                       "RSCRIPT": rscript[0]}
                cmds.append((step, ["bash", os.path.join(SCRIPTS, "align-asv-v4.sh"),
                                    os.path.join(p["results"], key + ".ASV.fa"), os.path.join(p["results"], "v4align"),
                                    str(args.threads or prm.get("threads") or 8)], env))
        elif step == "phyloseq":
            if (prm.get("tree") or {}).get("method", "none") != "none":
                cmds.append((step, rscript + [os.path.join(SCRIPTS, "dada2-asv.R"), "--config=" + conf_json,
                                              "--resume=outputs"] + threads, {}))
        elif step == "guilds":
            if prm.get("guilds"):
                cmds.append((step, rscript + [os.path.join(SCRIPTS, "fungal-guilds.R"), "--config=" + conf_json], {}))
        elif step == "figures":
            cmds.append((step, rscript + [os.path.join(SCRIPTS, "plot-readtracking.R"), "--config=" + conf_json], {}))
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
    # threads: --threads given here, else the config's dada2 threads, else the job's CPUs
    params = merged(cfg.get("dada2") or {}, (cfg.get("locus_settings") or {}).get(locus))
    if args.threads or not params.get("threads"):
        argv += ["--threads", str(args.threads) if args.threads else "$THREADS"]
    argv += ["--"] + list(args.extra) if args.extra else []
    lines = ["#!/bin/bash"]
    for opt, k in [("-A", "account"), ("-C", "constraint"), ("-q", "qos"), ("-N", "nodes"), ("-c", "cpus"), ("-t", "time")]:
        if sl.get(k) is not None:
            lines.append(f"#SBATCH {opt} {sl[k]}")
    lines += [f"#SBATCH -J asw-{cfg['run_id'] if 'run_id' in cfg else cfg['run']}-{locus}",
              "#SBATCH -o %x.%j.log", "set -euo pipefail", (sl.get("setup") or "").rstrip(),
              "export THREADS=${SLURM_CPUS_PER_TASK:-4}",
              " ".join(a if a == "$THREADS" else shlex.quote(a) for a in argv)]
    return "\n".join(l for l in lines if l) + "\n"


def common_options(f):
    """Options shared by resolve, check and run."""
    for deco in reversed([
        click.option("-c", "--config", required=True, type=click.Path(exists=True, dir_okay=False), help="project YAML config"),
        click.option("--locus", multiple=True, help="repeat or comma-separate (default: loci in the config)"),
        click.option("--profile", help="default: first profile whose detect path exists"),
        click.option("--merge/--no-merge", default=None, help="denoise R1 + R2 and merge / R1 only (default: config)"),
    ]):
        f = deco(f)
    return f


def setup(config, locus):
    cfg = load_yaml(config)
    loci = [l for item in locus for l in item.split(",") if l] or cfg.get("loci") or list(cfg["locus_settings"])
    return cfg, loci


def parse_steps(ctx, param, value):
    steps = [s for s in value.split(",") if s]
    bad = [s for s in steps if s not in STEPS]
    if bad:
        raise click.BadParameter(f"unknown steps: {', '.join(bad)} (steps: {', '.join(STEPS)})")
    return [s for s in STEPS if s in steps]   # always in pipeline order


@click.group(help=__doc__, context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 120})
def cli():
    pass


@cli.command("resolve")
@common_options
def resolve_cmd(config, locus, profile, merge):
    """Print the resolved per-locus config."""
    cfg, loci = setup(config, locus)
    for l in loci:
        click.echo(json.dumps(resolve(cfg, l, profile, merge), indent=2))


@cli.command("check")
@common_options
def check_cmd(config, locus, profile, merge):
    """Check inputs, references and tools."""
    cfg, loci = setup(config, locus)
    status = 0
    for l in loci:
        res = resolve(cfg, l, profile, merge)
        problems = check(res)
        click.echo(f"{l} ({res['key']}, profile {res['profile']}): " + ("ok" if not problems else ""))
        for msg in problems:
            click.echo("  " + msg)
        status |= bool(problems)
    sys.exit(status)


@cli.command("run")
@common_options
@click.option("--steps", default=",".join(STEPS), show_default=True, callback=parse_steps, help="comma-separated steps")
@click.option("--threads", type=int, help="overrides the config's dada2 threads (default there: SLURM CPUs, else 4)")
@click.option("--nsamples", type=int, help="test run on the first N samples (dada2 step)")
@click.option("--pdf", is_flag=True, help="analysis figures also as PDF")
@click.option("--slurm", is_flag=True, help="submit one SLURM job per locus")
@click.option("--dry-run", is_flag=True, help="print the commands only")
@click.argument("extra", nargs=-1, type=click.UNPROCESSED)
def run_cmd(config, locus, profile, merge, steps, threads, nsamples, pdf, slurm, dry_run, extra):
    """Run the workflow steps. After --: extra dada2-asv.R options, e.g. -- --resume"""
    cfg, loci = setup(config, locus)
    args = SimpleNamespace(config=config, profile=profile, merge=merge, steps=steps, threads=threads,
                           nsamples=nsamples, pdf=pdf, extra=extra)
    for l in loci:
        res = resolve(cfg, l, profile, merge)
        if slurm:
            script = slurm_script(cfg, args, l)
            path = os.path.join(res["paths"]["work"], f"{res['key']}.sbatch")
            if dry_run:
                click.echo(f"# {path}\n{script}")
                continue
            os.makedirs(res["paths"]["work"], exist_ok=True)
            with open(path, "w") as fh:
                fh.write(script)
            subprocess.run(["sbatch", path], check=True, cwd=res["paths"]["work"])
            continue
        conf_json = os.path.join(res["paths"]["results"], res["key"] + ".config.json")
        if not dry_run:
            conf_json = write_resolved(res)
        for step, argv, env in commands(res, conf_json, steps, args):
            click.echo(f"[{l}] {step}: " + " ".join(f"{k}={shlex.quote(v)}" for k, v in env.items())
                       + (" " if env else "") + " ".join(shlex.quote(a) for a in argv))
            if not dry_run:
                rc = subprocess.run(argv, env={**os.environ, **env}).returncode
                if rc != 0:
                    sys.exit(f"[{l}] step {step} failed (exit {rc})")


# Result files that `pull` copies (rsync filter rules, relative to results/<key>/). They stay git-ignored; the small
# text outputs (<key>.config.json, readtracking, failed_samples, run.log) are committed in the remote clone instead
# (see .gitignore), and checkpoints (*.RData), the full alignment and earlier-run subfolders stay on the remote.
PULL_INCLUDE = ["{key}.phyloseq.rds", "{key}.phyloseq.filtered.rds", "{key}.ASVtable.txt", "{key}.ASV.fa",
                "{key}.taxonomy.txt", "{key}.guilds.txt", "{key}.guild_abundance.txt",
                "{key}.phyloseq.filtered.guilds.rds", "*.pdf", "v4align/", "v4align/*.ASV.tree", "figures/", "figures/**", "analysis/", "analysis/**"]


def ssh_out(host, command):
    """stdout of a command on host ('' if it fails)."""
    r = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "LogLevel=ERROR", host, command],
                       capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


@cli.command("pull")
@click.option("-c", "--config", required=True, type=click.Path(exists=True, dir_okay=False), help="project YAML config")
@click.option("--locus", multiple=True, help="repeat or comma-separate (default: loci in the config)")
@click.option("--merge/--no-merge", default=None, help="as for run: selects the key (R1-only runs end in .R1)")
@click.option("--dry-run", is_flag=True, help="list what would be copied")
def pull_cmd(config, locus, merge, dry_run):
    """Copy results from the remote clone (config `remote`) into results/<key>/ here.

    \b
    Copies the analysis files (phyloseq objects, ASV table, fasta, taxonomy, tree, PDFs, analysis/); the small text
    outputs come with `git pull`, and checkpoints stay on the remote. Warns when the remote clone has uncommitted
    results or is at a different commit than this one."""
    cfg, loci = setup(config, locus)
    rem = cfg.get("remote") or {}
    if not rem.get("host") or not rem.get("repo"):
        sys.exit("config has no remote.host / remote.repo (the clone the workflow runs from)")
    host, rrepo = rem["host"], rem["repo"].rstrip("/")
    here_head = subprocess.run(["git", "-C", HERE, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    there_head = ssh_out(host, f"git -C {shlex.quote(rrepo)} rev-parse --short HEAD")
    if not there_head:
        sys.exit(f"cannot reach the clone {host}:{rrepo}")
    if there_head != here_head:
        click.echo(f"note: remote clone is at {there_head}, this repository at {here_head} "
                   "(git pull here, or commit/push results there)")
    for l in loci:
        res = resolve(cfg, l, rem.get("profile"), merge, repo=rrepo)
        key, src = res["key"], res["paths"]["results"].rstrip("/")
        dest = os.path.join(HERE, "results", key)
        dirty = ssh_out(host, f"git -C {shlex.quote(rrepo)} status --porcelain -- {shlex.quote(src)}")
        if dirty:
            click.echo(f"[{l}] note: uncommitted result changes in the remote clone ({len(dirty.splitlines())} files): "
                       "commit and push them there, then git pull here")
        rules = [f"--include={p.replace('{key}', key)}" for p in PULL_INCLUDE] + ["--exclude=*"]
        # options that also work with the old rsync 2.6.9 shipped with macOS
        argv = ["rsync", "-a", "-v", "--prune-empty-dirs", "-e", "ssh -o LogLevel=ERROR"] + rules \
            + (["--dry-run"] if dry_run else []) + [f"{host}:{src}/", dest + "/"]
        click.echo(f"[{l}] {key}: {host}:{src} -> {dest}")
        os.makedirs(dest, exist_ok=True)
        if subprocess.run(argv).returncode != 0:
            sys.exit(f"[{l}] rsync failed")


if __name__ == "__main__":
    cli()
