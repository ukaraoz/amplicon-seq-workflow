# neon_er_ER18_ER22 metadata

Sample metadata and fastq manifests for `config/neon-ER18-ER22.yaml` (`{dataset}` = `{project}_{run_id}`).
Copies of the files in `/global/cfs/cdirs/m3523/ertraits/ersynthesis/data/neon_itag/` on Perlmutter (identical
md5 when copied, 2026-11-09).

| file | what it is | used by the config |
|---|---|---|
| `NEON_ER_amplicon_metadata_soil.tsv` | one row per 16S/ITS sample (501): `sample`, `SampleSiteCode`, `SamplingCode`, `year`, `replicate`, `run_16S`, `run_ITS`, then the site's 2018 soil physical/chemical measurements | `paths.metadata` |
| `NEON_ER_16S_manifest.tsv`, `NEON_ER_ITS_manifest.tsv` | one row per sample: `sample`, `R1`, `R2` (absolute Perlmutter paths), `SampleSiteCode`, `year`, `replicate`, `run_folder` | `samples.manifest` |
| `NEON_ER_amplicon_metadata.tsv` | the amplicon metadata before the soil join | input of the merge |
| `ER_2018_NEON_soil-physical-chemical-measurements.csv` | 2018 soil physical/chemical measurements, one row per site (438) | input of the merge |
| `merge_soil_metadata.R` | joins the two by `SamplingCode` (3-digit site number) | makes `NEON_ER_amplicon_metadata_soil.tsv` |

The soil columns are 2018 measurements, also on the 2022 samples (`NNN-ER18`, `NNN-ER22` and `NNN-ER18_A/_B` of a site
share the row); their names are as in the soil table. To rebuild the merged file:

    Rscript merge_soil_metadata.R NEON_ER_amplicon_metadata.tsv ER_2018_NEON_soil-physical-chemical-measurements.csv NEON_ER_amplicon_metadata_soil.tsv

The manifests were made from `NEON_ER_amplicon_fastq_lookup.tsv` (site -> fastq paths per year and locus, built by
`map_site_fastqs.py`), which stay in the Perlmutter data folder.
