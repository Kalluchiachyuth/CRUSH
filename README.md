# CRUSH — Compartmental Refinement for Ultraprecise Stratification in Hi-C

<p align="center">
  <img src="figures/CRUSH_logo.jpg" alt="CRUSH Logo" width="300"/>
</p>

<p align="center">
  <a href="https://pypi.org/project/CRUSH"><img src="https://img.shields.io/pypi/v/CRUSH.svg" alt="PyPI version"/></a>
  <a href="https://github.com/JRowleyLab/CRUSH/blob/main/LICENSE"><img src="https://img.shields.io/github/license/JRowleyLab/CRUSH" alt="License"/></a>
  <img src="https://img.shields.io/badge/python-3.8%2B-blue" alt="Python 3.8+"/>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey" alt="Platform"/>
</p>

---

CRUSH **(Compartmental Refinement for Ultraprecise Stratification within Hi-C)** is a command-line tool that identifies fine-scale A/B chromatin compartments from Hi-C contact matrices. It has successfully identified compartments in Hi-C, Micro-C, and Single-Cell Hi-C data, and specializes in calling compartments at **high resolutions with significantly lower read depth** than other compartment calling tools.

> **Manuscript in preparation** — JRowleyLab, PI: Jordan Rowley

---

## Table of Contents

- [How It Works](#how-it-works)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Key Parameters](#key-parameters)
- [Test Dataset](#test-dataset)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [Contact](#contact)

---

## How It Works

<p align="center">
  <img src="figures/CRUSHdiag.png" alt="CRUSH workflow diagram" width="650"/>
</p>

At its core, CRUSH asks a simple question for every genomic bin: **does this bin interact more with A-type regions (iA) or B-type regions (iB)?**

The algorithm walks from coarse resolutions down to your target resolution, using each level to refine A/B compartment assignments at the next finer level:

1. **Eigenvector initialization** — Computes principal components of the Hi-C contact matrix (or accepts a user-supplied eigenvector) to define initial A (iA) and B (iB) states.
2. **CRUSH score calculation** — At each resolution, calculates a Genome Interaction (GI) score per bin reflecting how much more it contacts iA regions versus iB regions.
3. **Compartment reclassification** — After each resolution pass, A/B bin assignments are updated based on the new scores, then used to seed the next finer resolution.
4. **Resolution walking with midpoint shifting** — A rolling-window alignment step adjusts finer-resolution scores against the coarser baseline, removing systematic biases between resolution levels.
5. **Statistical filtering** — Applies Benjamini–Hochberg FDR correction and outputs a q-value filtered bedGraph.

**A compartments** → positive CRUSH score (gene-rich, open chromatin, active transcription)  
**B compartments** → negative CRUSH score (gene-poor, closed chromatin, transcriptionally silent)

> Unlike eigenvector-based methods, **you never need to flip CRUSH scores** — A is always positive and B is always negative.

---

## Installation

```bash
pip install CRUSH
```

We recommend setting up a dedicated conda environment:

```bash
conda create -n crush_env python=3.10
conda activate crush_env
conda install -c bioconda bedtools
pip install CRUSH hic-straw cooler numpy scipy pandas statsmodels tqdm
```

### Dependencies

| Tool | Purpose | Install |
|---|---|---|
| Python ≥ 3.8 | Runtime | [python.org](https://www.python.org) |
| bedtools | Genomic intersections | `conda install -c bioconda bedtools` |
| mawk | Fast text processing | `sudo apt install mawk` / `brew install mawk` |
| hic-straw | Read `.hic` files | `pip install hic-straw` |
| cooler | Read `.mcool` files | `pip install cooler` |
| numpy / scipy / pandas | Numerical computing | `pip install numpy scipy pandas` |
| statsmodels | FDR correction | `pip install statsmodels` |
| tqdm | Progress bars | `pip install tqdm` |

### Verify installation

```bash
crush --help
```

---

## Quick Start

### With genome build shortcut (supported builds: `hg19`, `hg38`, `mm10`, `mm9`; res ≥ 500 bp)

```bash
crush \
  -i data.hic \
  -gb hg38 \
  -r 10000 \
  -c 8 \
  -o output_prefix_
```

### With manual reference files (any genome, any resolution)

```bash
crush \
  -i data.hic \
  -g hg38.sizes \
  -a hg38_genes.bed \
  -b hg38.fa \
  -r 10000 \
  -c 8 \
  -o output_prefix_
```

> **Chromosome naming:** CRUSH automatically detects and converts chromosome prefix mismatches between your Hi-C file and reference files (e.g., `chr1` vs `1`). If output is empty or unexpected, verify that your Hi-C file itself uses a consistent naming convention throughout.

---

## Input Files

### Always required

| Flag | Description |
|---|---|
| `-i` | Hi-C file (`.hic` from Juicer or `.mcool` from cooler). Local path or HTTPS URL. |
| `-r` | Target resolution in base pairs (e.g., `10000` for 10 kb). Must exist in your Hi-C file. |

### Reference files — choose one of two paths

**PATH A — genome build shortcut** (res ≥ 500 bp only)

| Flag | Description |
|---|---|
| `-gb` | Genome build shortcut. Supported builds: `hg19`, `hg38`, `mm10`, `mm9`. Auto-downloads chr.sizes, genes.bed, and Bbins.bed from JRowleyLab GitHub. Not available for res < 500 bp because the hosted Bbins.bed was pre-computed at 500 bp — for sub-500 bp analysis supply `-g`, `-a`, and `-b` (FASTA) manually so CRUSH can recompute Bbins at your exact resolution. Explicit `-g`/`-a`/`-b` flags override the auto-download for that specific file. |

**PATH B — manual reference files** (any genome, any resolution)

| Flag | Description |
|---|---|
| `-g` | Chromosome sizes file — two tab-separated columns: `chr_name` and `size` (bp). No header. |
| `-a` | BED file (≥ 3 columns) for A-compartment initialization. Gene annotations work well. ChIP-seq peaks for an active histone mark (e.g., H3K27ac) also work. |
| `-b` | Genome FASTA **or** pre-computed Bbins BED for B-compartment initialization. With FASTA, CRUSH generates Bbins at 500 bp (res ≥ 500 bp) or at the input resolution (res < 500 bp). With BED, the file is used directly as B-compartment seeds. |

### Optional

| Flag | Description |
|---|---|
| `-e` | Pre-computed eigenvector bedGraph (4 columns: chr, start, end, value). Positive = A, Negative = B. Skips automatic eigenvector calculation. |

---

## Output Files

CRUSH produces four output files, each prefixed with whatever you supply via `-o`:

| File | Description |
|---|---|
| `{prefix}CRUSHparamters.txt` | Record of all parameters used. Keep this for reproducibility. |
| `{prefix}mergedCrush_{res}.bedgraph` | **Main output.** CRUSH scores for every bin. Positive = A compartment, Negative = B compartment. Unlike eigenvectors, scores never need to be flipped. |
| `{prefix}mergedqvalue_{res}.bedgraph` | Estimated q-value (BH-corrected) for each bin's score. |
| `{prefix}mergedCrush_{res}_qfiltered_reprocess.bedgraph` | CRUSH scores filtered to bins passing the q-value threshold. Note: this filter can be overly stringent — excellent results are often obtained from the unfiltered `mergedCrush` file. |

All bedGraph files include a UCSC track header for direct loading into genome browsers (IGV, UCSC, WashU).

While running, CRUSH creates a temporary working directory named `CRUSHtmp_[randomnumber]` in your current directory. This is removed automatically when the run completes. To keep it (e.g., for debugging), use `-C 0`. You can also name it yourself with `-f`.

---

## Key Parameters

| Flag | Default | Description |
|---|---|---|
| `-c` | `1` | Number of CPU threads. Set to number of chromosomes or available cores, whichever is smaller. |
| `-gb` | *(none)* | Genome build shortcut (`hg19`, `hg38`, `mm10`, `mm9`). Auto-downloads reference files. res ≥ 500 bp only. |
| `-o` | *(none)* | Output file prefix. |
| `-N` | `NONE` | Normalization: `NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`. |
| `-m` | `2500000` | Coarsest resolution to start walking from. |
| `-Z` | `100000` | Resolution for eigenvector calculation (100 kb recommended). |
| `-w` | `5` | Sliding window size (kb) for score averaging. Set to `1` to disable. Set to `0` for legacy auto-calculation from sequencing depth. |
| `-q` | `0.05` | Q-value threshold for filtered output. Set to `0` to disable filtering. |
| `-s` | `0` | Enable boundary smoothing (`1` = on). |
| `-A` | `0` | Adjust score distribution. **Do not use when comparing samples.** |
| `-C` | `1` | Clean up temp files after run (`0` = keep). |
| `-v` | `0` | Verbose output (`1` = on). |

For the complete parameter reference, see the [User Manual](MANUAL.md).

---

## Test Dataset

A small test dataset covering chromosomes 17–19 of hg19 is provided in `examples/TestData/`:

| File | Description |
|---|---|
| `hg19_c17_18_19_1kb.hic.gz` | Hi-C contact file |
| `hg19_c17_18_19_genes.bed.gz` | Gene annotations for A-state initialization |
| `hg19_c17_18.fa.gz` | Genome FASTA for GC-based B-state initialization |
| `hg19_c17_18.fa.fai` | FASTA index |
| `hg19_c17_18_19.sizes.gz` | Chromosome sizes |
| `Eigen_100kb_c17_18_19.bedgraph.gz` | Pre-computed eigenvector (optional `-e` input) |
| `Bbins_hg19_c17_18_19.bed.gz` | Pre-computed B-bins (alternative to FASTA for `-b`) |

### Run the test

```bash
# Decompress
gunzip examples/TestData/*.gz

# Run with FASTA-based B initialization
crush \
  -i examples/TestData/hg19_c17_18_19_1kb.hic \
  -g examples/TestData/hg19_c17_18_19.sizes \
  -a examples/TestData/hg19_c17_18_19_genes.bed \
  -b examples/TestData/hg19_c17_18.fa \
  -r 10000 \
  -c 4 \
  -o test_
```

Expected output: `test_mergedCrush_10000.bedgraph`, `test_mergedqvalue_10000.bedgraph`, and `test_mergedCrush_10000_qfiltered_reprocess.bedgraph`.

Load `test_mergedCrush_10000.bedgraph` into IGV or the UCSC browser to verify the A/B compartment pattern on chr17–19.

---

## Citation

Manuscript in preparation. If you use CRUSH in your research, please check back for the citation or contact us directly.

---

## Contact

**JRowleyLab** | PI: Jordan Rowley  
For questions, bug reports, or feature requests, please open a [GitHub Issue](https://github.com/JRowleyLab/CRUSH/issues).
