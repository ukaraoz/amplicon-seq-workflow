# results

`workflow.py` writes each run's outputs here, one folder per key (`paths.results: "{repo}/results/{key}"` in the
config), e.g. `results/neon_er_ER18_ER22.16S/` and `results/neon_er_ER18_ER22.ITS.R1/` (R1-only runs get `.R1`).
Each folder holds `<key>.config.json` (the resolved config), the ASV table, taxonomy, ASV fasta, tree, phyloseq
objects, read tracking and logs; community analysis goes to `<key>/analysis/`.

Only this README is tracked: results stay on the machine that ran the workflow (run from a clone of this
repository there).
