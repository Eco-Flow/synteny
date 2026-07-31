## Unreleased

### `Added`

- `main.nf` is now a thin dispatcher between two named workflows under `workflows/`, matching standard nf-core layout: `PAIRWISE_SYNTENY` (`workflows/pairwise_synteny.nf`, the existing default workflow, moved out of `main.nf` unchanged) and `ALGO` (`workflows/algo.nf`, new -- see below), selected via `--mode`
- New `ALGO` subworkflow (`--mode algo`), reconstructing ancestral linkage groups and ancestral gene order across many species from BUSCO single-copy orthologues, following the core method of Maulana et al. 2026 (bioRxiv 2026.07.17.739156): BUSCO -> Syngraph (ALGs + fission/fusion events) -> AGORA (ancestral gene order / CARs) -> a gene-order fragmentation index
  - Structured as three nf-core-style local subworkflows (`subworkflows/local/`): `GENOME_ACQUISITION` (download/validate genomes -- shared with `PAIRWISE_SYNTENY`, not duplicated), `SPECIES_TREE` (optional, see `--iqtree_species_tree` below), and `ANCESTRAL_RECONSTRUCTION` (BUSCO -> Syngraph -> AGORA -> fragmentation index), wired together by the thin `workflows/algo.nf`
  - New modules: `BUSCO`, `BUSCO_FILTER`, `SYNGRAPH`, `AGORA_PREP`, `AGORA`, `FRAGMENTATION_INDEX` (`modules/local/algo/`)
  - New `conf/test_algo.config` smoke-test profile, new `containers/syngraph/` and `containers/agora/` Dockerfiles (local-only, not yet published)
  - AGORA's real input requirements were confirmed against an actual `agora-generic.py` run (its docs describe a flat orthology-groups file, but it actually expects one file per ancestor node of the species tree, each restricted to that node's descendant species, and its own gene-tree loader silently misparses a single combined file): `bin/busco_to_agora.py` and `modules/local/algo/agora.nf` were corrected accordingly
  - The manuscript's fragmentation-index formula was a typeset equation not recoverable as text from the supplied methods; `fragmentation_index.tsv` reports the raw `M`/`A`/`B` components plus a documented placeholder `FI_placeholder` column pending the real equation
- Optional `--iqtree_species_tree` flag for the `ALGO` subworkflow, building the species tree from the same genomes (LONGEST/GFFREAD -> protein extraction -> OrthoFinder -> single-copy supermatrix -> IQ-TREE2 -> rooting) instead of requiring one via `--species_tree`, following the tree-building steps of [Eco-Flow/excon (tree_subsampling branch)](https://github.com/Eco-Flow/excon/tree/tree_subsampling)
  - New modules: `EXTRACT_PROTEINS`, `EXTRACT_SINGLE_COPY`, `ALIGN_SINGLE_COPY`, `CONCAT_SINGLE_COPY`, `ROOT_TREE` (`modules/local/algo/`), plus vendored nf-core `orthofinder` and `iqtree` modules
  - New params: `iqtree_species_tree`, `iqtree_outgroup`, `iqtree_args`, `iqtree_partition_model`, `mafft_args`, `orthofinder_args`

## v4.1.0 - 30.07.26

### `Added`

- Documented the `--score` flag, which is required to produce the `tables` directory of syntenic change summaries, and the `--chromopaint` flag, which is required to produce the painted chromosome figures. Both are off by default since v4.0.1
- `process.resourceLimits` in the base and AWS Batch configs, capping requested resources at 16 CPUs / 128 GB / 48 h (2 CPUs / 6 GB / 6 h in the test profiles)
- `version` and `nextflowVersion` to the manifest. The startup banner now reads the version from `workflow.manifest.version` instead of a hardcoded string, and Nextflow itself will refuse to run the pipeline on anything older than `v26.04.0`
- Renamed `CHANGELOG.nd` to `CHANGELOG.md` so it renders on GitHub

### `Fixed`

- Updated the pipeline to the strict configuration and script syntax required by Nextflow v26:
  - Replaced the `check_max()` config function with `process.resourceLimits`
  - Moved the top level `timestamp` variable into `params.timestamp`, and replaced the `if` statements in `nextflow.config` with conditional expressions
  - Corrected `params.seqkit_*` entries that were nested inside the `params` block, so they are read as `params.seqkit_*` again
  - Removed unused Groovy `import` statements, and moved the banner, the hex channel and the `onComplete` handler inside the entry workflow
  - Added the now mandatory `script:` label to the processes that were missing it, and quoted the names in `env()` output declarations
  - Fixed a `container = ` directive and the `tag "$sample_id"` directives in the `summarise_plots` modules, which referenced a variable those processes do not take
- `SYNTENY` now writes the `species.csv` that `shorten_chromnames.pl` needs to match BED files to species. Without it the chromosome names in the karyotype BEDs were left unshortened while the seqids file was renamed, causing jcvi to fail with a `ZeroDivisionError`
- Corrected the output directory names in the README, which still described the v4.0.0 layout (`Results`/`Figures`/`Data`/`Tables` are now `results`/`figures`/`output_data`/`tables`)

### `Dependencies`

- Nextflow `>= v26` is now required (tested on `v26.04.6`); the pipeline no longer runs on `v25`
- Pinned `nf-schema@2.7.3` and `nf-co2footprint@1.4.0`. `nf-amazon` is deliberately left unpinned so it tracks the version bundled with Nextflow; the previous `nf-amazon@2.5.2` pin failed on Nextflow v26 with `NoClassDefFoundError: StrBuilder` when resolving `s3://` inputs
- Migrated the `co2footprint` config to the nested `trace`/`report`/`summary`/`provenance` blocks used by nf-co2footprint v1.4.0. The old flat keys were being silently ignored, which wrote the reports into the launch directory instead of `results/pipeline_info/co2_emissions`

### `Deprecated`

- Removed the `--max_cpus`, `--max_memory` and `--max_time` parameters, along with their schema entries. Use a custom config (or `-c`) to override `process.resourceLimits` for your infrastructure

## v4.0.1 - 23.04.26

### `Added`

- New way to plot ribbon ortho plots with actual lengths of genomes

### `Fixed`

- Removed equals signs in module contsiner directive

## Example

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`