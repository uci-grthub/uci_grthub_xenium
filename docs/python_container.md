# Fast Python session (Apptainer)

`pixi run -e squidpy-dfs python` is slow to start: the env lives in
`.pixi/envs/squidpy-dfs` on /dfs9, and `import scanpy` touches thousands of
small files, each of them a slow metadata lookup on the parallel filesystem
(more than 3 minutes observed). `envs/squidpy/squidpy.sif` packs the same
locked environment into one squashfs image, so the lookups happen inside the
image instead.

## Build (once, and again whenever `pixi.lock` changes)

```bash
pixi run build-sif          # = sbatch src/submit_build_squidpy_sif.sh
sacct -j <jobid>            # logs in logs/build_squidpy_sif.<jobid>.{out,err}
```

`envs/squidpy/squidpy.def` installs `pixi install --locked -e squidpy-dfs` in a
`ghcr.io/prefix-dev/pixi` build stage and copies only the env into an
`ubuntu:24.04` final stage. The image is git-ignored.

## Use

```bash
pixi run ipy                      # IPython with scanpy, squidpy, anndata
pixi run py my_script.py          # any python args pass through
bash src/sif_exec.sh <cmd> ...    # any command inside the image
```

```python
import scanpy as sc
adata = sc.read_h5ad("output/xenium/anndata/merged/slide_0120636.h5ad")
```

pixi re-splits task arguments and drops their quotes, so
`pixi run py -c "print('x')"` breaks. For `python -c` with quotes, call the
wrapper directly:

```bash
bash src/sif_exec.sh python -c "import scanpy as sc; print(sc.read_h5ad('output/xenium/anndata/merged/slide_0120636.h5ad'))"
```

`src/sif_exec.sh` loads `apptainer/1.4.5` if needed. It passes the HPC3 site
binds (`SINGULARITY_BINDPATH`: `/dfs9,/dfs8,/data/homezvol1,/share`) on as
`APPTAINER_BINDPATH`, and adds `/dfs9` if nothing else binds it. `$HOME` and the
current directory are bound by default. Conda activation scripts run on
start, so variables such as `PROJ_DATA` and `GDAL_DATA` are set.

Measured 2026-11-17: `import scanpy` takes about 10 s in the image and loading
`slide_0120636.h5ad` about 11 s. The same import in `.pixi/envs/squidpy-dfs`
had not finished after 200 s. Building the image takes about 3 min.

## Build notes

- `--authfile` points at an empty file, so stale `ghcr.io` credentials in
  `~/.apptainer/docker-config.json` cannot get the public pull denied.
- The site bind variables are unset during the build. Otherwise `%post` fails
  with `destination /dfs9 doesn't exist in container`.
- Only the `export` lines and conda `activate.d` sources from
  `pixi shell-hook` go into the image. The bash-completion part fails in
  `%environment`.

The pixi tasks run in the empty `default` env, so `pixi run py` never touches
`.pixi/envs/squidpy-dfs`. That env stays in place for Snakemake and Jupyter.
