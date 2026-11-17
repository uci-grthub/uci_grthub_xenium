# Xenium sample naming

How a Xenium output directory under `data/` becomes the `sample_id` used in file
names and in `obs["sample_id"]`.

## Raw layout

Each region on a slide is delivered as one directory:

```
data/0926X-02_GaneA/<run>/output-XETG00221__<slide>__<region>__<date>__<time>/
```

10x separates the fields with a double underscore (`__`). Example:

```
output-XETG00221__0120633__b16_10724_tumor__20260916__214213
```

## Derivation

The directory name is **not parsed**. `get_xenium_samples()` (`Snakefile:92-118`)
globs every `experiment.xenium` manifest under `xenium.data_root` and reads two
fields from it:

- `slide_id`, for example `0120633`
- `region_name`, for example `b16_10724_tumor`. If this field is missing, the
  directory basename is used instead.

It then joins them with a single underscore (`Snakefile:113`):

```python
sample_id = f"{slide}_{region}" if slide else region
```

A duplicate id raises `ValueError`. The instrument serial (`XETG00221`) and the
timestamp are dropped.

| Raw directory | `slide_id` | `region_name` | `sample_id` |
|---|---|---|---|
| `output-XETG00221__0120633__b16_10724_tumor__20260916__214213` | `0120633` | `b16_10724_tumor` | `0120633_b16_10724_tumor` |

## Propagation

1. `src/xenium_to_anndata.py:46-48` sets `obs["sample_id"]` and
   `uns["sample_id"]`, and prefixes every `obs_name` with `<sample_id>_`.
2. Each sample is written to `output/xenium/anndata/{sample_id}.h5ad`. The
   copy in `annotated/` uses the same name.
3. `src/merge_xenium_anndata.py` concatenates the samples of each slide into
   `output/xenium/anndata/merged/slide_<slide>.h5ad`. It does the same for
   `merged_nocorrection/`. `obs["sample_id"]` is kept as a categorical column
   and serves as the batch key for Harmony.

## Underscore caveat

- A `sample_id` contains only single underscores. `__` appears only in the raw
  directory names.
- Region names contain underscores themselves, so a `sample_id` cannot be split
  back into slide and region on `_`. The only safe split is the first token:
  the slide id is always 7 digits and has no underscore.
- As of 2026-11-17, no code in `src/`, `rules/` or `Snakefile` splits
  `sample_id` on `_`.
- Joining with `__` at `Snakefile:113` instead would keep slide and region
  recoverable. It would also rename every output file and force a full rerun.

## Verified sample ids (2026-11-17)

The `sample_id` categories in `merged/` and `merged_nocorrection/` match the
per-sample h5ad file stems exactly:

| Slide | Samples (cells) |
|---|---|
| `0120633` | `0120633_b16_10724_tumor` (27014), `0120633_b16_21225_5e_3` (112565), `0120633_b16_3e_arn22089` (22256), `0120633_yum17_a39_4e_6623` (443938), `0120633_yumm17_arn22089_6e_10223` (93800), `0120633_yummer_resistance` (86913), `0120633_yummer_tumor` (75431) |
| `0120636` | `0120636_ll37_arn22089` (3126), `0120636_ll37_arn22089_6e` (2445), `0120636_ll37_v` (3338), `0120636_ll37_v_2e` (2507), `0120636_w_arn22089` (2459), `0120636_w_arn22089_6e` (1551), `0120636_w_v_2e` (739), `0120636_w_v_skin` (2032) |
