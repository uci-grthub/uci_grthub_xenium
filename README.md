# KadoK Visium HD (order 0526I-16)

Snakemake workflow for this project's Visium HD data on HPC3 (SLURM). Space
Ranger counts each capture area; the repository also carries the lab's Cell
Ranger and downstream single-cell rules, which this project does not currently
use.

Tissue: six TMAs of 2 mm human skin biopsies, two capture areas per Visium HD
6.5 mm slide, three slides. See `metadata/` for the slide layout documents.

## Layout
- `Snakefile` — the workflow. All rules live here except the includes below.
- `config.yaml` — references, paths, resources, per-sample overrides.
- `rules/` — `scenic.smk` is included; `jbrowse.smk`, `numbat.smk` and
  `pileup_and_phase.smk` are not wired into the Snakefile.
- `src/` — R, Python and shell scripts called by the rules (not `scripts/`).
- `data/` — FASTQs, one directory per sequencing run.
- `metadata/` — CytAssist run folders, H&E images, slide layout spreadsheets.
- `envs/` — conda environment YAMLs and the Space Ranger pixi manifest.
- `output/` — results; `output/logs/` for per-rule logs.

## Running

Snakemake lives in a conda environment, not a module:

```bash
source /opt/apps/mamba/24.3.0/etc/profile.d/conda.sh
conda activate snakemake_dfs
```

Dry run first. It prints the discovered capture areas and warns about anything
it could not resolve — read those before submitting:

```bash
snakemake -n
```

Submit the Visium HD counting:

```bash
snakemake --profile slurm spaceranger_all
```

The `slurm` profile is at `~/.config/snakemake/slurm`; it submits one `sbatch`
job per rule with the partition, account and memory each rule declares. Space
Ranger takes hours per capture area, so run the Snakemake process itself inside
an allocation (or a small `sbatch` wrapper), not on a login node and not in the
background of an interactive session.

`spaceranger_all` stops after counting. The default `all` target additionally
pulls in the Cell Ranger summary and the Harmony integration, which have no
inputs in this project — use `spaceranger_all` unless you have added single-cell
samples. `src/submit_workflow.sh` is stale: it refers to a `cluster.yaml` and a
`snakemake` module that do not exist here.

To rerun one capture area:

```bash
snakemake --profile slurm output/spaceranger/A1_TMA4/outs/web_summary.html
```

## Visium HD

### How capture areas are discovered

There is no sample list to maintain. `rule spaceranger_count` builds one job per
capture area by reading, for each run folder under `spaceranger.cytassist_dir`:

- the `*-run-info.csv`, for the slide ID and the sample loaded in each capture area;
- the CytAssist `.tif` for that capture area;
- FASTQs under `spaceranger.fastq_root`, matched as `<prefix>_S*_L*_R1_*.fastq.gz`
  and searched recursively, so a sample sequenced across several runs contributes
  all of its directories;
- a microscope image from `spaceranger.microscope_image_dir`, if one matches.

A capture area with no FASTQs on disk is skipped with a warning, so the workflow
can be run while transfers are still arriving; the rest appear on the next
invocation with no edit.

Slide and capture area are read by Space Ranger from the CytAssist image
metadata, not passed on the command line. The workflow passes `--slide`, `--area`
and `--override-id` only when both are set in the config for that sample, for the
case where the values entered at the instrument were wrong.

### Run ids

The run id is the sample name entered at the CytAssist, with any character
outside `[a-zA-Z0-9_-]` replaced by `_`, since that is all Space Ranger accepts
for `--id`. If one name was used on more than one capture area, each becomes
`<name>_<slide>_<area>`. Both cases are announced at parse time. Results land in
`output/spaceranger/<run id>/`.

In this project the ids are `A1_TMA2`, `A1_TMA4`, `A1_TMA6`, `D1_TMA1`,
`D1_TMA3`, `D1_TMA5`.

### Configuration

```yaml
spaceranger:
  pixi_manifest: "envs/spaceranger/pixi.toml"
  transcriptome: "/dfs8/commondata/spaceranger/ref/refdata-gex-GRCh38-2020-A"
  probe_set: "/dfs8/commondata/spaceranger/probe/Human_V2/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
  fastq_root: "data"
  cytassist_dir: "metadata/H_E_images_out_0526I_16/CytAssist"
  microscope_image_dir: "metadata/H_E_images_out_0526I_16/TIFF"
  microscope_image_regex: "TMA\\s*(\\d+)$"     # captures {token} from the sample name
  microscope_image_pattern: "TMA {token}.tif"  # glob; {sample} {slide} {area} {token}
  create_bam: false
  localcores: 24
  localmem: 112       # GB
  custom_bin_size:    # extra bin size in um, on top of the standard 2/8/16
  mem_mb: 120000
  runtime: 2880       # minutes
  partition: "standard"
  account: "sbsandme_lab"
  samples: {}
```

`samples:` holds per-capture-area overrides, keyed by run id, by CytAssist sample
name, or by `<slide>_<area>`. Prefer `<slide>_<area>`: it is the only key that
stays unambiguous when a sample name repeats across slides.

```yaml
  samples:
    H1-PWCYQ3Y_A1:
      fastq_prefix: "CoreSampleName"   # when the core did not use the CytAssist name
      image: "metadata/.../section.tif"
      loupe_alignment: "metadata/.../alignment.json"
      slide: "H1-PWCYQ3Y"              # slide + area together add --override-id
      area: "A1"
```

### Outputs

Per capture area, under `output/spaceranger/<run id>/outs/`: `web_summary.html`,
`metrics_summary.csv`, and `binned_outputs/square_008um/filtered_feature_bc_matrix{,.h5}`.
The 2 µm and 16 µm bins are written alongside. Logs in `output/logs/`.

Space Ranger emits one matrix per capture area with no notion of which tissue
core a bin came from. Assigning the TMA cores (patient, normal/pigmented,
timepoint — see `metadata/TMA_slide_image_metadata_summary_v3.pdf`) to bins is a
separate step, via Loupe selections or polygon annotation
(`src/read_polygons_with_names.py`, `src/append_polygon_query_metadata.py`).

## Requirements

- **Space Ranger 3.1.2** — there is **no module** for it on HPC3. The workflow
  calls it through `envs/spaceranger/pixi.toml`, a pixi environment that puts the
  cluster-wide install at `/dfs8/commondata/spaceranger/spaceranger-3.1.2` on
  PATH. Space Ranger is proprietary and is not on conda-forge or bioconda, so the
  environment installs nothing. Check it with:
  ```bash
  pixi run --manifest-path envs/spaceranger/pixi.toml spaceranger --version
  ```
  This manifest is deliberately separate from the workspace `pixi.toml`, whose
  PyPI dependencies fail to solve on this cluster.
- Snakemake 9 with the `cluster-generic` executor plugin (conda env `snakemake_dfs`).
- References under `/dfs8/commondata/spaceranger/`.
- CellRanger 8.0.1 (module) for the single-cell rules only.

## Other rules in this repository

Current and used: `spaceranger_count`, `fastqc`, `multiqc`,
`collect_web_summaries`.

Single-cell rules, working but idle in this project (no Cell Ranger samples
configured): `cellranger_count`, `multi_sample_summary`,
`tenx_harmony_integration`, `tenx_scvi_integration`, `tenx_harmony_notebook`,
`plot_filtering_timeline`, `cellbender`, `scrublet`.

Legacy, and broken as written — do not expect these to run without repair:

- `scvelo` and `pileup_and_phasing` call `src/compute_velocity.py` and
  `src/pileup_and_phase.R`, neither of which is in the repository.
- `velocyto` passes `config["..."]` as literal shell text rather than an
  evaluated value.
- `rules/scenic.smk` consumes `output/loom/<sample>.loom`, which was produced by
  a `loompy` rule that has been removed.
- `rules/jbrowse.smk`, `rules/numbat.smk` and `rules/pileup_and_phase.smk` are
  not included by the Snakefile.

## Troubleshooting

- **A capture area is missing from the dry run.** The warning says which piece
  was not found. "No FASTQs" usually means the sequencing core named them
  differently from the CytAssist sample; set `fastq_prefix` for that
  `<slide>_<area>`.
- **"no microscope image".** The run proceeds on the CytAssist image alone, at
  much lower resolution (3200×3000 against roughly 10000×10000). Fix
  `microscope_image_pattern`/`microscope_image_regex`, or set `image` for that
  sample.
- **Tissue looks misaligned in `web_summary.html`.** Check the slide and capture
  area in the CytAssist run info against the physical layout. The slide
  spreadsheet in `metadata/` numbers the areas in the opposite order from the
  CytAssist run info, so confirm which is correct before trusting an alignment.
- **Re-running a capture area** deletes and recreates its output directory;
  Space Ranger does not resume a partial run here.
- Per-rule SLURM logs are under `logs/<rule>/`; Space Ranger's own output is in
  `output/logs/spaceranger_<run id>.log`.

## Documentation
- `rules/README.md` — downstream analysis rules.
- `docs/prep_geo_data.md` — GEO submission staging.
