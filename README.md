# SuperGWAS — K-mer GWAS Pipeline

SuperGWAS is a high-performance k-mer-based genome-wide association study pipeline with three analysis modes:

- **1D:** single k-mer association with a phenotype
- **2D:** k-mer × k-mer epistasis within one genome
- **2W:** two-way host × pathogen k-mer interaction across two genomes

The 1D and 2D analyses produce a Manhattan plot and an interactive 3D Manhattan plot for each reference genome.

---

## Quick start

```bash
# Clone on the cluster. Precompiled x86-64 binaries are included in bin/
git clone https://github.com/githubcbrc/SuperGWAS.git && cd SuperGWAS

# Configure
cp config.sh config.local.sh
vim config.local.sh          # Edit paths for your environment
```

Then choose what to run.

**Preview the 1D + 2D submission** without submitting any jobs:

```bash
./scripts/run_pipeline.sh --mode cols --dry-run
```

**Submit the 1D + 2D pipeline:**

```bash
./scripts/run_pipeline.sh --mode cols
```

**For 1D only**, skip the 2D steps:

```bash
./scripts/run_pipeline.sh --mode cols --no-2d
```

**For host × pathogen analysis**, run the 2W pipeline:

```bash
./scripts/run_pipeline_2w.sh
```

The job scripts call the binaries in `bin/` directly (override with `BIN_DIR=...`). To rebuild the binaries, see [Building](#building).

Monitor submitted jobs with `watch squeue --me`. Stage logs are written to `${RESULTS_ROOT}/logs/<stage>/`.

---

## Configuration

Use `config.sh` as a reference and copy it to `config.local.sh` for your settings. When `config.local.sh` is present, the scripts source it instead of `config.sh`.

```bash
# --- Shared inputs ------------------------------------------------------------

# Sample IDs, one per line. Line i names the sample in position i of every
# k-mer bitstring, so the order must match the matrices exactly.
ACCESSIONS="/path/to/all_accessions.txt"

# Phenotype values: "<sample_id> <value>" per line (tab or space), no header.
# Samples in ACCESSIONS but missing here are excluded from the test.
PHENOTYPE="/path/to/phenotype.txt"

# Directory of reference genomes (*.fasta / *.fa). Results are mapped to each one.
REF_DIR="/path/to/references"

# --- 1D + 2D inputs -----------------------------------------------------------

# K-mer presence/absence matrix, row format: *.tsv.gz shards with lines
# "<kmer>\t<bitstring>", one bit per sample in ACCESSIONS order (1 = present).
ROW_MATRIX_DIR="/path/to/kmer_matrix_rows"

# The same matrix in column format (batch_*/ directories), used by --mode cols.
# Generate it from ROW_MATRIX_DIR with HPC/convert_matrix/submit.sh.
COL_MATRIX_DIR="/path/to/kmer_matrix_cols"

# --- 2W inputs ----------------------------------------------------------------

# Host and pathogen k-mer matrices, "<label> <bitstring>" per line. One bit per
# host accession in HOST_MATRIX, and one bit per pathogen accession in
# PATHOGEN_MATRIX.
HOST_MATRIX="/path/to/host_matrix"
PATHOGEN_MATRIX="/path/to/pathogen_matrix"

# One observation per line, "<host_idx> <pathogen_idx> <value>".
# Indices are 0-based bit positions in the host and pathogen bitstrings.
# Must be a complete grid in host-major order (see "Experimental design" in the
# 2W pipeline section).
PHENOTYPE_2W="/path/to/phenotype_2w.txt"

# --- Output -------------------------------------------------------------------

# All outputs go under here (structure/, kgwas/, mappings/, logs/, ...).
RESULTS_ROOT="/path/to/results"

# Short phenotype name used in output directory and file names.
PHENO_NAME="trichome"
```

See [Configuration parameters](#configuration-parameters) for the full list.

---

## Repository structure

```
SuperGWAS/
├── config.sh                 # Central configuration (copy to config.local.sh)
├── Makefile                  # Builds all binaries into bin/
├── scripts/
│   ├── run_pipeline.sh       # 1D + 2D pipeline (submits all stages)
│   └── run_pipeline_2w.sh    # 2W host × pathogen pipeline
├── HPC/                      # SLURM job scripts, one folder per stage
│   ├── sample_kmers/         # K-mer sampling
│   ├── build_structure/      # filter_kmers + compute_structure_all → U.bin, U2D.bin
│   ├── kgwas/                # 1D GWAS (cols or rows) + merge helpers
│   ├── build_index/          # 1D k-mer lookup index
│   ├── map_merge_filter/     # Map to reference + density filter
│   ├── build_index_2d/       # Small index of 2D candidate k-mers
│   ├── kgwas_2d/             # All-pairs epistasis (array) + merge
│   ├── map_2d_pairs/         # Translate 2D pairs to genome positions
│   ├── plot/                 # 1D PNG + 2D HTML plots per reference
│   ├── compute_structure_2w/ # 2W: dual PCA → U0.txt
│   ├── compute_projections_2w/ # 2W: projections → Th.txt, Tp.txt
│   ├── kgwas_2w/             # 2W: all host × pathogen pairs
│   └── convert_matrix/       # Data prep: row → column format
├── src/                      # C++ sources
├── bin/                      # Precompiled binaries
├── tools/
│   ├── plot_1d.py            # Manhattan plot (matplotlib, PNG)
│   └── plot_2d.py            # 3D epistasis Manhattan (plotly, HTML)
└── docker/
    ├── Dockerfile
    └── build.sh
```

---

## 1D + 2D pipeline

```
sample_kmers → build_structure → kgwas ─┬─ build_index → map_merge_filter ──────────┬─→ plot
                                        └─ build_index_2d → kgwas_2d → map_2d_pairs ┘
```

After `kgwas`, the 1D and 2D branches run in parallel, and plotting waits for both. See [1D + 2D stages](#1d--2d-stages) for what each step does.

Paths in this README are relative to `RESULTS_ROOT`. In them, `{pheno}` stands for `PHENO_NAME` (for example `trichome`), `{ref}` for a reference genome's file name without its extension (for example `TAIR10` for `TAIR10.fasta`), and `{name}` for the batch or shard processed by one `kgwas` task (see step 3).

| Step | Name | Description | Job type | Output |
|------|------|-------------|----------|--------|
| 1 | sample_kmers | Sample k-mers for population structure | Array (shards) | `sampled/part_*.tsv.gz` |
| 2 | build_structure | Merge samples, MAF filter + subsample, PCA for 1D and 2D | Single | `sampled/sampled_kmers.tsv.gz`, `structure/U.bin`, `structure_2d/U2D.bin` |
| 3 | kgwas | 1D association tests (cols or rows mode) | Array (batches) | `kgwas/{pheno}/{name}.tsv`, `{name}.2d.tsv`, `{name}.tsv.count` |
| 4 | build_index | 1D k-mer lookup index | Single | `index/{pheno}/` |
| 5 | map_merge_filter | Map to reference(s), merge contigs into genome-wide coordinates, density filter | Array (refs) | `mappings/{pheno}/{ref}/{pheno}.filtered.tsv` |
| 6 | build_index_2d | Index of 2D candidate k-mers | Single | `index/{pheno}_2d/` |
| 7 | kgwas_2d | All-pairs epistasis, split into `NUM_2D_TASKS` + merge | Array + single | `kgwas_2d/{pheno}_2d.tsv`, `.kmers` |
| 8 | map_2d_pairs | Translate k-mer pairs to genome positions | Array (refs) | `kgwas_2d/{ref}/{pheno}_2d.tsv` |
| 9 | plot | 1D Manhattan + 2D 3D Manhattan | Array (refs) | `plots/{pheno}.png`, `{pheno}_2d.html` |

### 1D + 2D stages

#### Step 1: Sample k-mers (`sample_kmers`)

Draws a random sample of k-mers from the full row matrix. Each array task samples its share of `*.tsv.gz` shards into `sampled/part_<N>.tsv.gz` (plus a `.stats` file). Together the parts hold `RESERVOIR_SIZE` k-mers (default 10 million). These files are intermediate outputs: step 2 (`build_structure`) merges them.

#### Step 2: Build population-structure covariates (`build_structure`)

Turns the sample into population-structure covariates:
- Concatenates the parts into `sampled/sampled_kmers.tsv.gz`. Concatenated gzip files are valid gzip, so nothing is recompressed.
- `filter_kmers` keeps only the phenotyped samples, drops k-mers with a minor allele frequency (MAF) below `MAF_THRESHOLD` (5%), and randomly selects `SUBSAMPLE_SIZE` (200,000) from the remaining k-mers. The result is `structure/X_sample.tsv`.
- `compute_structure_all` runs a principal component analysis (PCA) on that matrix. `U.bin` holds the top `RANK` principal components, the covariates for 1D tests. `U2D.bin` extends them with pairwise products of the components, the covariates for 2D tests.

#### Step 3: Run the 1D association scan (`kgwas`)

Tests every k-mer in the full matrix for association with the phenotype, correcting for population structure with `U.bin`. It can run in one of two modes, depending on the k-mer matrix storage format (see [K-mer matrix formats](#k-mer-matrix-formats)). K-mers whose correlation with the phenotype is below `CORR_THRESHOLD` are skipped before the full test. For each batch or shard it writes three files:

| File | Purpose | Columns or contents |
|------|---------|---------------------|
| `{name}.tsv` | 1D association results | `kmer  effect  pvalue` |
| `{name}.2d.tsv` | Candidates for the 2D analysis | `kmer  bitstring  effect  pvalue` |
| `{name}.tsv.count` | Number of k-mers tested, used for the Bonferroni threshold in the 1D plot | Count |

`{name}` identifies the batch or shard the task processed.

`{name}.tsv` contains every k-mer that passed the correlation pre-filter and has -log10(p) ≥ `LOG10_THRESH` (default 6). This cutoff only keeps the output to a manageable size. Passing it does not by itself establish significance; significance is assessed against the Bonferroni threshold.

`{name}.2d.tsv` contains the k-mers with -log10(p) ≥ `LOG10_THRESH_2D_INPUT` (default 15). Each 2D candidate includes its presence/absence bitstring for the phenotyped samples. The 2D scan uses these bitstrings to test candidate pairs without rereading the full k-mer matrix.

##### Choosing output and candidate cutoffs

The defaults for `LOG10_THRESH` and `LOG10_THRESH_2D_INPUT` were chosen for the development dataset and may not suit other datasets. Each setting controls a different part of the workflow:
- `LOG10_THRESH` sets how many results are written. Lower it to keep more of the weaker associations, or raise it if the `*.tsv` files become too large for the index and mapping steps.
- `LOG10_THRESH_2D_INPUT` sets the number of 2D candidates, N. The 2D step tests every pair of candidates, N(N − 1)/2 pairs, so its cost grows with the square of N.

Keep both cutoffs well below the 1D significance threshold, so that significant k-mers are not dropped.

The 2D cutoff is also supported empirically: 2D signals have come from pairs in which at least one k-mer had a strong marginal signal, even if not strictly significant.

#### Step 4: Index the 1D results (`build_index`)

Builds a lookup index from the 1D association results. Step 5 (`map_merge_filter`) uses this index to look up reference k-mers efficiently.

K-mers are stored in canonical form (the smaller of the k-mer and its reverse complement), so either strand matches. They are split into shards (16 by default) by the low bits of their 2-bit encoding. Each shard has:
- `shard_NNN.umap`: fixed 32-byte records sorted by k-mer (16-byte encoded k-mer, 8-byte effect, 8-byte p-value). The 128-bit encoding limits the k-mer length to 64.
- `shard_NNN.bloom`: a Bloom filter over the shard's k-mers.

`index.meta` records the k-mer length (`K`), the number of shards (`M`) and the number of k-mers. The mapping step loads each shard into a hash table and checks the Bloom filter first, which rejects most reference k-mers without a table lookup.

#### Step 5: Map, filter and merge per reference (`map_merge_filter`)

Runs once per reference genome and does three things in one pass, with no intermediate files:
- **Map:** scans every contig of the reference and looks up each k-mer in the index. It keeps hits with -log10(p) ≥ `LOG10_MIN`.
- **Filter:** applies a density filter within each contig. A hit survives only if it lies in a `WINDOW_SIZE` bp window (default 100) where at least `MIN_DENSITY` (95%) of positions start a hit. This removes isolated hits, keeping regions where consecutive overlapping k-mers all pass.
- **Merge:** combines all contigs into one genome-wide coordinate system. Each contig is offset by the total length of the contigs before it in the FASTA, so a hit at position 1,000 of the second contig gets `pos = length(contig 1) + 1000`.

The output `{pheno}.filtered.tsv` contains `pos  effect  pvalue` in genome-wide coordinates. Rows are not sorted by `pos`, but they do not need to be for plotting: their `pos` values determine their locations.

#### Step 6: Index the 2D candidates (`build_index_2d`)

Collects the 2D candidates from every batch or shard (the `{name}.2d.tsv` files from step 3) into a single index in `index/{pheno}_2d/`, so that step 7 (`kgwas_2d`) can load all of them at once.

The index is built the same way as the 1D index in step 4: canonical k-mers split into shards, each with a `.umap` and a `.bloom` file. In addition, each shard has a `shard_NNN.bits` file holding the candidates' presence/absence bitstrings, which step 7 needs to test pairs.

#### Step 7: Test all candidate pairs (`kgwas_2d`)

Tests every pair of candidate k-mers for an interaction effect on the phenotype, correcting for structure with `U2D.bin`. It writes only pairs with -log10(p) ≥ `LOG10_THRESH_2D` (default 18). As in step 3, this cutoff only limits the output size; it is not a significance threshold. The pairs are split across `NUM_2D_TASKS` array tasks, and a merge job joins the parts into `{pheno}_2d.tsv`. That file contains `idx_a  idx_b  effect  pvalue`, where the indices are line numbers in `{pheno}_2d.kmers`. Pairs are identified by k-mer, not by position: no reference is used until step 8 (`map_2d_pairs`) places them on each genome.

#### Step 8: Map candidate pairs to each reference (`map_2d_pairs`)

Runs once per reference genome. It reads the candidate k-mers from `{pheno}_2d.kmers`, which `kgwas_2d` writes from the step 6 index, and builds a small lookup table keyed by canonical k-mer. The pair indices refer to entries in this candidate list.

The stage then scans the reference and rewrites each mapped pair as `pos_a  pos_b  effect  pvalue`, using the same genome-wide coordinates as step 5 (`map_merge_filter`). Pairs where either k-mer doesn't occur in that reference are dropped.

#### Step 9: Plot the results (`plot`)

Runs once per reference genome. It draws the 1D Manhattan plot from `{pheno}.filtered.tsv` and the interactive 3D plot from the mapped 2D pairs. See [Plotting](#plotting) for requirements.

### 1D + 2D options

Options for `scripts/run_pipeline.sh`:

| Option | Description |
|--------|-------------|
| `--mode cols\|rows` | kgwas input format (default: cols) |
| `--no-2d` | Skip the 2D branch (steps 6–8: `build_index_2d`, `kgwas_2d`, `map_2d_pairs`) |
| `--dry-run` | Show what would be submitted without running |

### Auto-skip

Auto-skip is designed for resuming a failed run. Each step starts only if the previous one succeeded, so a failure stops all later steps. Re-running the same command then skips the completed steps and resumes from the one that failed.

For the array steps `kgwas`, `map_merge_filter` and `map_2d_pairs`, the check passes as soon as one task has written output. If only some of their tasks failed, delete that step's outputs before re-running.

The checks look only for existing outputs, not whether they match the current inputs or configuration. Use a new results directory for a changed analysis.

| Step | Skipped when |
|------|--------------|
| sample_kmers | `sampled/sampled_kmers.tsv.gz` exists |
| build_structure | both `U.bin` and `U2D.bin` exist |
| kgwas | `kgwas/{pheno}/` contains any 1D `.tsv` |
| build_index / build_index_2d | `index.meta` exists in the index directory |
| map_merge_filter | any `*.filtered.tsv` under `mappings/{pheno}/` |
| kgwas_2d | `{pheno}_2d.tsv` and `{pheno}_2d.kmers` exist |
| map_2d_pairs | any per-reference `kgwas_2d/*/{pheno}_2d.tsv` exists |

Plotting always runs.

### Multi-reference support

Steps 5 (`map_merge_filter`), 8 (`map_2d_pairs`) and 9 (`plot`) process every `.fasta` / `.fa` file in `REF_DIR` as an array job, one task per reference:

```
mappings/trichome/
└── TAIR10/
    ├── trichome.filtered.tsv
    └── plots/trichome.png
kgwas_2d/
├── trichome_2d.tsv          # pairs as k-mer indices (reference-independent)
├── trichome_2d.kmers
└── TAIR10/
    ├── trichome_2d.tsv      # pairs as genome positions
    └── trichome_2d.html
```

---

## 2W host × pathogen pipeline

```
compute_structure_2w → compute_projections_2w → kgwas_2w
```

The 2W pipeline looks for host–pathogen k-mer pairs associated with a phenotype measured on host × pathogen combinations, for example disease scores from infecting each host line with each pathogen isolate. It tests every host k-mer against every pathogen k-mer.

| Step | Name | Description | Job type | Output |
|------|------|-------------|----------|--------|
| 1 | compute_structure_2w | Host and pathogen PCA, combined into grid structure covariates | Single | `structure_2w/U0.txt` |
| 2 | compute_projections_2w | Precompute each k-mer's projection onto the structure | Single | `structure_2w/Th.txt`, `Tp.txt` |
| 3 | kgwas_2w | Test all host × pathogen k-mer pairs: joint test with 3 degrees of freedom (df), or 1-df interaction-only | Single (multi-threaded) | `kgwas_2w/kgwas_2w.tsv` |

### Experimental design

The analysis assumes a complete grid: every host is scored against every pathogen, exactly once.

- **Size:** with `nh` hosts and `np` pathogens, `PHENOTYPE_2W` must have `nh × np` lines. The run stops with an error otherwise.
- **Order:** use host-major order: list all pathogens for host 0, then all pathogens for host 1, and continue in the same way (`0 0`, `0 1`, …, `0 np−1`, `1 0`, …).
- **Indices:** host `i` is bit `i` of every `HOST_MATRIX` bitstring, and pathogen `j` is bit `j` of every `PATHOGEN_MATRIX` bitstring. So `nh` and `np` are the bitstring lengths of the two matrices.

**Important:** the pipeline does not check the order. A differently ordered file can run without an error and produce incorrect results.

### 2W stages

#### Step 1: Build grid structure covariates (`compute_structure_2w`)

Builds population-structure covariates for the grid. It runs a PCA on each matrix separately, keeping the top `RANK_HOST` host and `RANK_PATHOGEN` pathogen components. The baseline model is then built by concatenating four blocks of covariates column-wise:
- an intercept
- the host's PC scores, for structure among hosts
- the pathogen's PC scores, for structure among pathogens
- the Kronecker product of the previous two, for structure that depends on the combination, such as related hosts responding alike to related pathogens

#### Step 2: Precompute k-mer projections (`compute_projections_2w`)

Computes, for every host k-mer and every pathogen k-mer, its projection onto `U0`. `Th.txt` has one line per host k-mer and `Tp.txt` one line per pathogen k-mer. The pairwise scan reuses these projections, avoiding the need to recompute them for every pair.

#### Step 3: Test all host × pathogen pairs (`kgwas_2w`)

Tests every host k-mer against every pathogen k-mer with a likelihood-ratio test against the null model. By default it is a joint 3-df test of the host k-mer's effect, the pathogen k-mer's effect and their interaction. With `--interaction-only` it tests the interaction term alone, with 1 df. It writes every pair, with no p-value cutoff, so the output has exactly (host k-mers × pathogen k-mers) lines.

### 2W options

Options for `scripts/run_pipeline_2w.sh`:

| Option | Description |
|--------|-------------|
| `--interaction-only` | Use the 1-df interaction-only test instead of the 3-df joint test |
| `--dry-run` | Show what would be submitted without running |

Steps 1 (`compute_structure_2w`) and 2 (`compute_projections_2w`) are skipped if their outputs already exist (`U0.txt`; `Th.txt` and `Tp.txt`). Delete them to force a rerun. Step 3 (`kgwas_2w`) always runs.

---

## Data preparation

### K-mer matrix formats

Both formats hold the same presence/absence matrix of k-mers × samples, with samples in `ACCESSIONS` order and `1` meaning the k-mer is present.

**Row format** (input to sampling and `kgwas --mode rows`): the matrix split into `*.tsv.gz` shards, one k-mer per line. Each line is the k-mer, a tab, and its bitstring, with one character per sample:

```
ACGTACGTTGCA...	01100101...
TTGACGGATCCA...	00111010...
```

The sequences and bitstrings are abbreviated here (`...`); the gap between them is a single tab.

**Column format** (input to `kgwas --mode cols`): the matrix transposed and split into batches of k-mers, with one file per sample instead of one line per k-mer. This lets `kgwas` read only the phenotyped samples' files. Each batch is a directory:

```
batch_XXX/
├── kmers.txt.gz          # one k-mer per line; line i is k-mer i of the batch
└── bits/
    ├── col_000000.txt.gz # sample 0: a single line of 0/1, character i = k-mer i
    ├── col_000001.txt.gz # sample 1
    └── ...               # one file per sample, numbered by position in ACCESSIONS
```

Files may also be uncompressed (`.txt`).

### Converting row to column format

```bash
./HPC/convert_matrix/submit.sh
```

The script groups consecutive row-format shards into batches. `BATCH_SIZE` sets the number of shards per batch and defaults to 10; set it as an environment variable. Each array task converts one batch and writes it to `COL_MATRIX_DIR/batch_NNN/`. `kgwas --mode cols` then runs one task per batch. Larger batches mean fewer jobs but more memory per job. Batches already converted are skipped, so an interrupted run can be resubmitted.

---

## Input files

| File | Format | Description |
|------|--------|-------------|
| `ACCESSIONS` | One ID per line | Sample order (matches bitstring columns) |
| `PHENOTYPE` | `ID<WS>value`, no header | Phenotype values |
| `REF_DIR/*.fasta`, `*.fa` | FASTA | Reference genome(s) |
| `HOST_MATRIX`, `PATHOGEN_MATRIX` | `label<WS>bitstring` | 2W k-mer matrices, one k-mer per line |
| `PHENOTYPE_2W` | `host_idx<WS>pathogen_idx<WS>value` | One line per host × pathogen combination; indices are 0-based positions in the bitstrings. Complete grid in host-major order ([details](#experimental-design)) |

`<WS>` denotes a tab or space.

---

## Output files

| File | Pipeline | Columns / description |
|------|----------|-----------------------|
| `kgwas/{pheno}/{name}.tsv` | 1D | `kmer  effect  pvalue` |
| `kgwas/{pheno}/{name}.2d.tsv` | 1D → 2D | `kmer  bitstring  effect  pvalue` (2D candidates) |
| `kgwas/{pheno}/{name}.tsv.count` | 1D | K-mers tested per batch (for Bonferroni) |
| `mappings/{pheno}/{ref}/{pheno}.filtered.tsv` | 1D | Mapped, density-filtered hits |
| `mappings/{pheno}/{ref}/plots/{pheno}.png` | 1D | Manhattan plot |
| `kgwas_2d/{pheno}_2d.tsv` | 2D | `idx_a  idx_b  effect  pvalue` |
| `kgwas_2d/{pheno}_2d.kmers` | 2D | K-mer sequence per line (line number = index) |
| `kgwas_2d/{ref}/{pheno}_2d.tsv` | 2D | `pos_a  pos_b  effect  pvalue` |
| `kgwas_2d/{ref}/{pheno}_2d.html` | 2D | Interactive 3D Manhattan |
| `structure_2w/U0.txt` | 2W | Structure covariates, one row per phenotype line |
| `structure_2w/Th.txt`, `Tp.txt` | 2W | Projection of each host or pathogen k-mer onto `U0`, one line per k-mer |
| `kgwas_2w/kgwas_2w.tsv` | 2W | `host_idx  path_idx  host_label  path_label  ratio  delta  LRT  pvalue`, one line per pair, unfiltered |

### Merging 1D results

```bash
./HPC/kgwas/merge.sh PHENO IPATH OPATH     # concat all batches (+ .count totals)
./HPC/kgwas/filter_merge.sh 1e-15          # SLURM job: keep p < 1e-15 and merge
```

`merge.sh` concatenates the batch results in the directory `IPATH` into `OPATH/PHENO.all.tsv`, and their `.count` files into `OPATH/PHENO.all.counts`. `PHENO` is the phenotype name used for the output files. `filter_merge.sh` takes a p-value cutoff (default 1e-15) and submits a SLURM job that filters the batch results in `kgwas/{pheno}/` and merges them into `kgwas/{pheno}/{pheno}.all.tsv`.

---

## Plotting

Step 9 (`plot`) runs `tools/plot_1d.py` and `tools/plot_2d.py` with the cluster's Python environment, not the container. That Python environment must include the following packages:

```bash
pip install --user numpy pandas matplotlib plotly
```

`HPC/plot/plot.sh` first runs `module load python/3.11.0`. Module names differ between clusters: change it to a Python module that exists on yours (`module avail python` lists them), or remove the line and rely on `python3` from `PATH`.
If the packages are missing, the plot step logs `WARNING: 1D/2D plot failed` and the job still exits successfully, so check `${LOGS_DIR}/plot/` if the plot files are missing.

You can also copy the results and plot locally:

```bash
python3 tools/plot_1d.py mappings/trichome/TAIR10/trichome.filtered.tsv \
    --fasta TAIR10.fasta --min-logp 6 \
    --counts kgwas/trichome/*.count \
    --title "trichome - TAIR10" --out trichome.png

python3 tools/plot_2d.py --fasta TAIR10.fasta \
    --tsv kgwas_2d/TAIR10/trichome_2d.tsv \
    --marginal mappings/trichome/TAIR10/trichome.filtered.tsv \
    --title "trichome - TAIR10" -o trichome_2d.html
```

Run either script with `--help` for all options.

---

## Configuration parameters

### 1D

| Parameter | Default | Description |
|-----------|---------|-------------|
| `K` | 51 | K-mer size |
| `RANK` | 3 | PCA rank for population structure |
| `CORR_THRESHOLD` | 0.2 | Correlation pre-filter (0 = disabled) |
| `LOG10_THRESH` | 6 | -log10(p) threshold for kgwas output |
| `LOG10_MIN` | 6 | Min -log10(p) for mapping and plotting |
| `RESERVOIR_SIZE` | 10,000,000 | K-mer sample size |
| `SUBSAMPLE_SIZE` | 200,000 | Matrix subsample size |
| `MAF_THRESHOLD` | 5 | Minor allele frequency threshold (%) |
| `SEED` | 42 | Random seed |
| `SHARDS_PER_JOB` | 50 | Shards per sampling job |
| `WINDOW_SIZE` | 100 | Density filter window (bp) |
| `MIN_DENSITY` | 0.95 | Minimum k-mer density in window |

### 2D

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LOG10_THRESH_2D_INPUT` | 15 | -log10(p) for a k-mer to enter the 2D analysis |
| `LOG10_THRESH_2D` | 18 | -log10(p) threshold for kgwas_2d output |
| `NUM_2D_TASKS` | 16 | Array tasks for kgwas_2d (environment variable) |
| `EFFECT_PERCENTILE` | 90 | 2D plot effect-size percentile (environment variable) |

### 2W

| Parameter | Default | Description |
|-----------|---------|-------------|
| `RANK_HOST` | 3 | PCA rank for host |
| `RANK_PATHOGEN` | 3 | PCA rank for pathogen |
| `INTERACTION_ONLY` | false | Use the 1-df interaction-only test instead of the 3-df joint test |

### Threading

| Parameter | Default | Description |
|-----------|---------|-------------|
| `THREADS` | 48 | General thread count |
| `LOAD_THREADS` | 24 | kgwas data-loading threads |
| `BLAS_THREADS` | 112 | kgwas BLAS threads |
| `EMIT_THREADS` | 112 | kgwas output-writing threads |

---

## Building

The `bin/` directory contains precompiled binaries for convenience. To rebuild them after changing `src/`, use Docker. The image compiles all 14 binaries. Copy the compiled binaries from the container into `bin/`:

```bash
docker build -t supergwas -f docker/Dockerfile .
id=$(docker create supergwas)
docker cp "$id":/usr/local/bin/. bin/
docker rm "$id"
```
