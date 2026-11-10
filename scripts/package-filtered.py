#!/usr/bin/env python3
"""Bundle the filtered phyloseq exports (and the ITS guild results) into <key>.filtered.tar.gz with a README.

usage: package-filtered.py --config=<key>.config.json        (workflow.py `package` step, the last step)

version 1.1

The tarball holds one folder <key>.filtered/ with
  README.txt                       <key>.filtered.README.txt (written by dada2-asv.R), plus a guild section when
                                   guild files are present
  <key>.filtered.otu_table.txt, .tax_table.txt, .sample_data.txt, .refseq.fa, .tree (16S)   from dada2-asv.R
  <key>.readtracking.txt, .metadata_columns.txt, .taxa_filter.txt, .failed_samples.txt    provenance, dada2-asv.R
  <key>.guilds.txt, .guild_summary.txt, .guild_abundance.txt                                 fungal-guilds.R (ITS)
No figures.
Files that do not exist (e.g. no tree for ITS, no guilds for 16S) are left out. Standard library only.
"""
import datetime
import io
import json
import os
import sys
import tarfile

VERSION = "1.1"
EXPORTS = [".filtered.otu_table.txt", ".filtered.tax_table.txt", ".filtered.sample_data.txt", ".filtered.refseq.fa",
           ".filtered.tree"]
PROVENANCE = [".readtracking.txt", ".metadata_columns.txt", ".taxa_filter.txt", ".failed_samples.txt"]
GUILDS = [
    (".guilds.txt", "one row per ASV (all ASVs; in_filtered = in the filtered object): asv_id, asv_seq, taxonomy, "
                    "reads, samples, FungalTraits (FT_*) and FUNGuild (FG_*) assignments and how they were matched"),
    (".guild_summary.txt", "per database / field / value: ASVs, reads, % of fungal reads, samples, over the "
                           "filtered object, with the unassigned share and why"),
    (".guild_abundance.txt", "FUNGuild, wide: one row per sample, fungal_reads, then the relative abundance (share of "
                             "the sample's fungal reads) of each trophic mode (\"trophicMode: ...\") and each guild "
                             "(\"guild: ...\"); each set sums to 1 per sample"),
]


def main():
    cfg_arg = [a for a in sys.argv[1:] if a.startswith("--config=")]
    if len(cfg_arg) != 1:
        sys.exit(__doc__)
    cfg = json.load(open(cfg_arg[0].split("=", 1)[1]))
    key, res = cfg["key"], cfg["paths"]["results"]
    dirs = [res] + [d for d in [cfg["paths"].get("analysis")] if d and d != res]

    def find(name):
        return next((os.path.join(d, name) for d in dirs if os.path.isfile(os.path.join(d, name))), None)

    readme_src = os.path.join(res, key + ".filtered.README.txt")
    if not os.path.isfile(readme_src):
        sys.exit(f"error: {readme_src} missing (written by dada2-asv.R; run the dada2 or phyloseq step)")
    readme = open(readme_src).read().rstrip("\n")

    files = [f for f in (find(key + s) for s in EXPORTS + PROVENANCE) if f]
    rt = find(key + ".readtracking.txt")
    if rt:
        cols = open(rt).readline().rstrip("\n").split("\t")
        tf = [c for c in cols if c.startswith("taxa_filtered_") and not c.endswith("_pct_nonchim")]
        readme += (f"\n\nRead tracking\n  {os.path.basename(rt)}\n      one row per sample that passed filterAndTrim: reads at "
                   f"each step ({', '.join(c for c in cols if c not in ('sample', 'run'))})"
                   + (f"; {tf[0]} = reads kept by the taxa filter (what minfinal applies to), with its % of nonchim"
                      if tf else ""))
    guilds = [(f, desc) for f, desc in ((find(key + s), d) for s, d in GUILDS) if f]
    if guilds:
        readme += ("\n\nFungal guilds (fungal-guilds.R; FungalTraits 1.2 by genus, FUNGuild by the most specific rank)\n"
                   + "\n".join(f"  {os.path.basename(f)}\n      {d}" for f, d in guilds)
                   + "\n  ASV IDs (asv_id) are the same as in the filtered tables.")
    readme += (f"\n\npackaged {datetime.datetime.now().astimezone().strftime('%Y-%m-%d %H:%M %Z')} "
               f"by package-filtered.py {VERSION}\n")
    files += [f for f, _ in guilds]

    out = os.path.join(res, key + ".filtered.tar.gz")
    top = key + ".filtered"
    with tarfile.open(out + ".part", "w:gz") as tar:
        data = readme.encode()
        ti = tarfile.TarInfo(f"{top}/README.txt")
        ti.size, ti.mtime, ti.mode = len(data), int(datetime.datetime.now().timestamp()), 0o644
        tar.addfile(ti, io.BytesIO(data))
        for f in files:
            tar.add(f, arcname=f"{top}/{os.path.basename(f)}")
    os.replace(out + ".part", out)
    print(f"package-filtered.py {VERSION}: {len(files)} files + README.txt -> {out}")
    for f in files:
        print("  " + os.path.basename(f))


if __name__ == "__main__":
    main()
