# results

`workflow.py` writes each run's outputs here, one folder per key (`paths.results: "{repo}/results/{key}"` in the
config), e.g. `results/neon_er_ER18_ER22.16S/` and `results/neon_er_ER18_ER22.ITS.R1/` (R1-only runs get `.R1`).
Each folder holds `<key>.config.json` (the resolved config), the ASV table, taxonomy, ASV fasta, tree, phyloseq
objects, read tracking and logs; community analysis goes to `<key>/analysis/`.

Run the workflow from a clone of this repository on the machine with the data (config `remote`). There:

- the small text outputs are committed: `<key>.config.json`, `readtracking.txt`, `failed_samples.txt`,
  `run.log` (see .gitignore). Commit and push them after a run; `git pull` brings them elsewhere.
- `workflow.py pull -c CONFIG [--locus L]` copies the analysis files into results/<key>/ on another machine:
  phyloseq objects, ASV table, ASV fasta, taxonomy, tree, PDFs and analysis/ (git-ignored).
- checkpoints (`*.RData`), the full V4 alignment and earlier-run subfolders stay in the clone.
