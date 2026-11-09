# amplicon-seq-workflow

16S and ITS amplicon sequencing workflow: demultiplexed paired-end reads → DADA2 ASVs → taxonomy → tree →
phyloseq → community analysis. Steps are R and shell scripts; a Python driver (`workflow.py`) runs them from one
YAML config per project/run, locally or as SLURM jobs (NERSC Perlmutter).

## Steps

| step | script | what it does |
|---|---|---|
| `reference` | `scripts/get-unite.sh` | ITS: download the UNITE release if missing; 16S: SILVA DADA2 files must be in place |
| `dada2` | `scripts/dada2-asv.R` | primer check, cutadapt (ITS), filterAndTrim, error models, pseudo-pooled denoising, merging (or R1 only), ASV length window, chimeras, read tracking, taxonomy (assignTaxonomy/addSpecies), ASV table/fasta, phyloseq, failed samples |
| `tree` | `scripts/align-asv-v4.sh` | 16S: mothur alignment to the SILVA V4 reference (MiSeq SOP), screen/filter, FastTree GTR |
| `phyloseq` | `scripts/dada2-asv.R --resume=outputs` | rebuild the phyloseq objects with the tree |
| `analyze` | `scripts/analyze-phyloseq.R` | taxa bar plots, alpha diversity, Venn, NMDS + PERMANOVA (UniFrac or Bray-Curtis/Jaccard), pairwise PERMANOVA, source sharing, distances to a reference treatment, betapart turnover/nestedness |

Main outputs (in `paths.results`, prefix `<key>`): `ASVtable.txt`, `taxonomy.txt`, `ASV.fa`, `readtracking.txt`,
`failed_samples.txt`, `phyloseq.rds` (all ASVs and samples), `phyloseq.filtered.rds` (`remove_taxa` /
`require_ranks` ASVs and failed samples removed; input of the analysis). Analysis figures and tables go to
`paths.analysis`.

## Setup

```bash
conda env create -f environment.yml
conda activate amplicon-seq-workflow
```

References (paths set in the config):
- 16S: SILVA DADA2 training files (`silva_nr99_v144_toGenus_trainset.fa.gz`, `silva_v144_assignSpecies.fa.gz`) from
  https://www.arb-silva.de/current-release/DADA2/1.36.0/SSU; for the tree, mothur's `silva.nr_v138_2.tgz`
  (https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_2.tgz) in the `tree.reference` folder (the V4
  alignment is built from it on first use).
- ITS: UNITE general release, downloaded by the `reference` step (`taxonomy.unite_doi`).

## Usage

```bash
./workflow.py check   -c config/biochar-M007649.yaml               # inputs, references, tools
./workflow.py resolve -c config/biochar-M007649.yaml --locus 16S   # resolved settings (what the scripts read)
./workflow.py run     -c config/biochar-M007649.yaml --dry-run     # print the commands

# everything, both loci, on this machine
./workflow.py run -c config/biochar-M007649.yaml --threads 16
# test on 4 samples, denoising only
./workflow.py run -c config/biochar-M007649.yaml --locus ITS --steps dada2 --nsamples 4 -- --no-taxonomy
# analysis only, with PDFs
./workflow.py run -c config/biochar-M007649.yaml --steps analyze --pdf
# ITS with merged pairs instead of R1 only
./workflow.py run -c config/biochar-M007649.yaml --locus ITS --merge
# SLURM: one job per locus (settings under `slurm` in the config)
./workflow.py run -c config/biochar-M007649.yaml --slurm
```

Extra options after `--` go to `dada2-asv.R` (`--no-taxonomy`, `--resume`, `--resume=outputs`).

## Config

One YAML file per project/run (`config/template.yaml`; worked example `config/biochar-M007649.yaml`):

- `profiles`: per-machine paths (`vars`); the first profile whose `detect` path exists is used, or `--profile`.
- `samples.manifest` (optional): tab-separated file with one line per sample and columns `R1`, `R2` (read paths)
  plus an optional sample name column; otherwise fastq files are `paths.raw/<samples.r1|r2>` per metadata sample.
- `samples.run_column` (optional, with a manifest): the sequencing run of each sample. Runs get their own error
  models and denoising, and `locus_settings.<locus>.runs.<run>` can override `truncLen` / `maxEE` / `trimLeft`.
- `paths`, `samples`: metadata, fastq naming, output folders; `{project}`, `{run}`, `{run_id}`, `{locus}` and any
  `vars` can be used in strings.
- `dada2` + `locus_settings.<locus>`: primers, cutadapt, merge / R1 only, truncLen, trimLeft, maxEE, ASV length window,
  taxonomy references, taxa removed from the filtered object, tree method.
- `analysis`: experimental design and analysis settings. Metadata columns are mapped to the internal names
  `habitat`, `treatment`, `pot`, `is_T0` (baseline flag) and `is_pellet` (source flag), which appear in the tables.
  Treatment levels, labels (and where long labels wrap) and colors; ordination subsets and PERMANOVA formulas.
  Optional sections run only when present: `baseline`, `source`, `venn`, `source_sharing`,
  `reference_distances`, `betapart`.

`workflow.py` writes the resolved per-locus config to `<results>/<key>.config.json`; the R scripts read only that
file, so a run can be repeated with `Rscript scripts/<script>.R --config=<that file>`.

## Notes

- Script versions are in each script header (`version X.Y`); bump on every change.
- PyYAML is optional: without it the driver parses the config with R's `yaml` package.
