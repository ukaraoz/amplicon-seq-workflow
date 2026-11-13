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
| `NEON_LeafTraits_master_11.28.2023.csv` | leaf traits per site (438 sites, `SampleSiteCode`): site, area, collection date, vegetation type, coordinates, elevation, d15N, d13C, leaf N / C (wt %), C:N, bare / NPV / PV cover, LMA, leaf water content, dominant plant and its fractional cover and taxonomy | `paths.site_tables` (joined by `samples.site_key`) |
| `NEON_RS_master_updated_11.28.2023.csv` | remote sensing products per site (438 sites, `SampleSiteCode`): the same site columns, Lidar terrain (aspect, curvature, flow accumulation, slope, solar radiation, TPI, TWI, upslope area at 5 / 9 cell windows, elevation, canopy height) and spectral products | `paths.site_tables` (joined by `samples.site_key`) |

The leaf-trait and remote-sensing files are not in git (see `.gitignore`): on the laptop they are symlinks to
`~/Work/EastRiver/synthesis/data/biogeo_remote/`, on Perlmutter copies (same md5) here in the clone and in the data
folder. Both are keyed by `SampleSiteCode` (`NNN-ER18`) and are joined to every sample of the site (2018 and 2022) by
`dada2-asv.R` (`paths.site_tables`, `samples.site_key`): site columns already in the metadata (`Sampling.Area`,
`Latitude`, `Longitude`) are checked against it and not repeated; all other columns are added under their own names
(collection dates are June 2018). `<key>.metadata_columns.txt` in the results lists each column's source.

The soil columns are 2018 measurements, also on the 2022 samples (`NNN-ER18`, `NNN-ER22` and `NNN-ER18_A/_B` of a site
share the row); their names are as in the soil table. To rebuild the merged file:

    Rscript merge_soil_metadata.R NEON_ER_amplicon_metadata.tsv ER_2018_NEON_soil-physical-chemical-measurements.csv NEON_ER_amplicon_metadata_soil.tsv

The manifests were made from `NEON_ER_amplicon_fastq_lookup.tsv` (site -> fastq paths per year and locus, built by
`map_site_fastqs.py`), which stay in the Perlmutter data folder.

`neon_er_ER18_ER22.sample_map.png` maps every sampling location of the filtered 16S and ITS objects on Esri World
Imagery (fill: vegetation type; shape: years sampled; labels: site and sample count), and `neon_er_ER18_ER22.sample_locations.txt` lists those samples (site, vegetation type,
year, coordinates, loci). Made from the results with

    Rscript scripts/sample-map.R --out=config/neon_er_ER18_ER22/metadata --name=neon_er_ER18_ER22 \
      results/neon_er_ER18_ER22.16S results/neon_er_ER18_ER22.ITS.R1

`neon_er_ER18_ER22.Snodgrass.slope_map.png`: the Snodgrass site only, each location labelled with its sampling code,
over terrain slope (USGS 3DEP elevation, slope computed in UTM, 10-degree classes) with 20 m contours (add
`--site=Snodgrass --slope --colour=family` to the command above); points are coloured by plant family within
vegetation type (Meadow greens, Shrubs oranges, Trees purples; the most common families over all samples get their own
shade, the rest "other"), shape = vegetation type, black outline = sampled in 2018 and 2022. Both maps have a scale bar.

The map images carry only a title. Elevation for the slope map: USGS 3D Elevation Program (3DEP). Base-map credit, to cite with the figures: Tiles © Esri - Source: Esri, i-cubed,
USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community. Samples without
coordinates (10, all ER18: 049, 159, 175, 209, 210, 230, 231, 339, 350, 380) are not shown; they are the rows
without Latitude / Longitude in `neon_er_ER18_ER22.sample_locations.txt`.
