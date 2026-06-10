# CRUSH User Manual

**Compartmental Refinement for Ultraprecise Stratification in Hi-C**  
Version 1.0 | JRowleyLab | PI: Jordan Rowley  
Manuscript in preparation

---

## Table of Contents

1. [Overview](#1-overview)
2. [Background: A/B Compartments](#2-background-ab-compartments)
3. [How CRUSH Works](#3-how-crush-works)
4. [Installation](#4-installation)
5. [Input File Requirements](#5-input-file-requirements)
6. [Running CRUSH](#6-running-crush)
7. [Complete Parameter Reference](#7-complete-parameter-reference)
8. [Output Files](#8-output-files)
9. [Test Dataset Walkthrough](#9-test-dataset-walkthrough)
10. [Tips, Recommendations, and Common Issues](#10-tips-recommendations-and-common-issues)
11. [Comparing Samples](#11-comparing-samples)
12. [FAQ](#12-faq)

---

## 1. Overview

CRUSH is a command-line tool for calling **A/B chromatin compartments** from Hi-C data. It accepts a `.hic` (Juicer) or `.mcool` (cooler) file and produces a bedGraph of CRUSH scores across the genome, where:

- **Positive scores** indicate **A compartments** — open, transcriptionally active chromatin
- **Negative scores** indicate **B compartments** — closed, transcriptionally silent chromatin

CRUSH has successfully identified compartments in **Hi-C, Micro-C, and Single-Cell Hi-C** data, and is specifically designed to call compartments at high resolutions with lower read depth than traditional approaches.

A key advantage over eigenvector-based methods: **CRUSH scores never need to be flipped.** A is always positive and B is always negative, regardless of chromosome or dataset — the algorithm handles orientation internally.

What makes CRUSH distinct is its **multi-resolution refinement strategy**: it starts at a coarse resolution (e.g., 2.5 Mb), builds confidence about compartment identity, then walks down to progressively finer resolutions (e.g., 10 kb or 1 kb), re-evaluating and refining assignments at each step.

---

## 2. Background: A/B Compartments

Chromatin in the nucleus is organized into two broad spatial compartments, first described by Lieberman-Aiden et al. (2009):

**A compartments** are associated with:
- Active gene transcription
- Open (accessible) chromatin
- High GC content
- Enrichment for active histone marks (H3K27ac, H3K4me3)

**B compartments** are associated with:
- Silent or lowly expressed genes
- Closed (compact) chromatin
- Low GC content
- Enrichment for repressive histone marks (H3K27me3, H3K9me2/3)

In a Hi-C contact matrix, A compartment bins preferentially contact other A bins, and B bins preferentially contact other B bins, forming a characteristic "plaid" pattern. The **CRUSH score** (Genome Interaction score) formalizes this preference: a bin's score reflects how much more it contacts A-type regions versus B-type regions.

---

## 3. How CRUSH Works

CRUSH performs the following steps:

### Step 1 — Eigenvector Initialization

If no eigenvector is supplied (`-e`), CRUSH computes principal components of the observed/expected Hi-C contact matrix at the eigenvector resolution (`-Z`, default 100 kb). It evaluates up to 10 principal components per chromosome, flipping signs so that positive values correlate with gene density (iA) and negative values correlate with low GC content (iB). The best PC per chromosome is selected based on combined correlation with gene annotations and GC content.

If you supply an eigenvector (`-e`), this step is skipped entirely.

### Step 2 — CRUSH Score Calculation

For each chromosome at each resolution, CRUSH:

1. Dumps the contact matrix from the Hi-C file.
2. Distance-normalizes contacts by dividing each contact by the expected value at that genomic distance.
3. Z-scores each row of the distance-normalized matrix.
4. Uses the eigenvector-derived iA/iB states to classify each bin's contacts as either A-type or B-type interactions.
5. Computes the CRUSH score for each bin as the weighted difference between its mean A-type interaction Z-score and mean B-type interaction Z-score, normalized by total contact count.

### Step 3 — Resolution Walking with Midpoint Shifting

CRUSH processes resolutions from coarsest to finest. After calculating scores at each resolution, it runs a "midpoint shifter" — a rolling-window alignment step that adjusts finer-resolution scores relative to the coarser-resolution baseline. This removes systematic biases introduced by switching resolution and ensures scores are comparable across the genome.

### Step 4 — Compartment Reclassification

After each resolution pass, CRUSH re-evaluates which bins are A and which are B based on the updated CRUSH scores, then uses these updated bins to seed the next resolution's calculation. This feedback loop progressively improves assignment confidence at finer scales.

### Step 5 — Statistical Testing and Output

At the finest resolution, CRUSH calculates a t-statistic comparing each bin's A-type interaction mean to its B-type interaction mean, converts this to a p-value, and applies Benjamini–Hochberg FDR correction. The final output includes both raw CRUSH scores and q-value filtered scores.

---

## 4. Installation

### Via PyPI (recommended)

```bash
pip install crush-hic
```

> ⚠️ **`pip install` does not install `bedtools` or `mawk`** — these are system tools that must be installed separately before running CRUSH:
>
> ```bash
> # Linux (apt)
> sudo apt install bedtools mawk
>
> # macOS (Homebrew)
> brew install bedtools mawk
>
> # Conda (any platform)
> conda install -c bioconda bedtools mawk
> ```

### Setting up a clean environment

We strongly recommend a dedicated conda environment:

```bash
conda create -n crush_env python=3.10
conda activate crush_env
conda install -c bioconda bedtools mawk
pip install crush-hic
```

### Verifying installation

```bash
crush --help
```

You should see the CRUSH usage banner with all available options.

### Dependency summary

| Tool | Purpose | Install command |
|---|---|---|
| Python ≥ 3.8 | Runtime | [python.org](https://www.python.org) |
| bedtools | Genomic intersections | `conda install -c bioconda bedtools` *(not via pip)* |
| mawk | Fast text processing | `sudo apt install mawk` *(not via pip)* |
| hic-straw | Read `.hic` files | installed automatically with `crush-hic` |
| cooler | Read `.mcool` files | installed automatically with `crush-hic` |
| numpy | Numerical arrays | installed automatically with `crush-hic` |
| scipy | Statistics | installed automatically with `crush-hic` |
| pandas | Data frames | installed automatically with `crush-hic` |
| statsmodels | FDR correction | installed automatically with `crush-hic` |
| tqdm | Progress bars | installed automatically with `crush-hic` |

---

## 5. Input File Requirements

> **Note on chromosome naming:** CRUSH automatically detects and converts chromosome prefix mismatches between your Hi-C file and your reference files (e.g., `chr1` vs `1`). If a mismatch is found, CRUSH converts the reference files to match your Hi-C file and prints a notice. No manual action is needed in most cases. If output is empty or unexpected, verify that your Hi-C file itself uses a consistent naming convention throughout.

### Genome build shortcut (`-gb`) — optional

CRUSH can automatically download the chromosome sizes file, gene annotations, and B-compartment bins for supported genome builds:

```bash
-gb hg38    # or hg19, mm10, mm9
```

**Supported builds:** `hg19`, `hg38`, `mm10`, `mm9`  
**Resolution restriction:** only available for resolutions ≥ 500 bp. The Bbins.bed hosted on GitHub was pre-computed at 500 bp resolution, which is appropriate for analyses at 500 bp and coarser. For finer resolutions, Bbins need to be computed at your exact resolution — supply `-g`, `-a`, and `-b` (FASTA) manually so CRUSH can recalculate GC content at the correct bin size.

Files are downloaded from the JRowleyLab GitHub into the temporary working directory. If you supply `-g`, `-a`, or `-b` explicitly alongside `-gb`, those explicit flags take priority over the auto-download for that specific file.

```bash
# Fully automatic
crush -i data.hic -gb hg38 -r 10000 -c 8

# Override only B-bins
crush -i data.hic -gb hg38 -b myCustomBbins.bed -r 10000 -c 8
```

---

### Hi-C file (`-i`)

CRUSH accepts two formats:

**Juicer `.hic` format:**
```bash
-i /path/to/sample.hic
# or a remote HTTPS URL:
-i https://example.com/sample.hic
```

**Cooler `.mcool` format:**
```bash
-i /path/to/sample.mcool
```

The file must contain the resolution you request with `-r`. CRUSH will report an error if the resolution is not available.

### Chromosome sizes file (`-g`)

A plain-text file with exactly two tab-separated columns: chromosome name and size in base pairs. No header line.

```
chr1    249250621
chr2    243199373
chr3    198022430
```

Standard size files are available from UCSC:
`https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.chrom.sizes`

### A-state initialization file (`-a`)

A BED file (minimum 3 columns: `chr`, `start`, `end`) defining regions to treat as A-compartment seeds. Gene annotations are the standard choice because expressed genes reliably mark active compartments. ChIP-seq peaks for active histone marks (H3K27ac, H3K4me1, H3K4me3) are also effective.

```bash
-a hg38_genes.bed
```

### B-state initialization file (`-b`)

`-b` accepts either a FASTA file or a pre-computed Bbins BED file.

**A FASTA file** — CRUSH uses `bedtools nuc` to compute GC content across the genome and identifies bins with below-average GC content as B-compartment seeds (Bbins). This requires no prior knowledge of B-compartment regions and works at any resolution.

- If your target resolution is ≥ 500 bp, Bbins are computed at **500 bp** resolution internally.
- If your target resolution is < 500 bp, Bbins are computed at **the input resolution**.

```bash
-b /path/to/genome.fa
```

The FASTA must be indexed (`.fai`). If it is not already indexed:
```bash
samtools faidx genome.fa
```

**A pre-computed Bbins BED file** — If you already have a BED file of B-compartment seed regions (e.g., from a previous CRUSH run, lamina-associated domains, or low-GC bins), supply it directly:

```bash
-b known_Bbins.bed
```

This file is used as-is; CRUSH does not recalculate GC content. This is appropriate at any resolution as long as the BED bins are suitable for that resolution.

### Eigenvector file (`-e`) — optional

A pre-computed eigenvector in bedGraph format (4 columns: `chr`, `start`, `end`, `value`). Positive values = A compartment; negative values = B compartment. Supply this to skip automatic eigenvector calculation.

```bash
-e my_eigenvector.bedgraph
```

Using a pre-computed eigenvector is faster and strongly recommended when running CRUSH on multiple samples you intend to compare, as it guarantees consistent sign orientation across all samples.

### Resolution (`-r`)

Target resolution in base pairs. CRUSH automatically detects all available resolutions between your target and the coarsest resolution (`-m`, default 2.5 Mb) and walks through them all.

```bash
-r 10000    # 10 kb
-r 5000     # 5 kb
-r 1000     # 1 kb
```

For Juicer-generated `.hic` files, common available resolutions are: 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 1000.

---

## 6. Running CRUSH

### Minimal command — with genome build shortcut

```bash
crush \
  -i sample.hic \
  -gb hg38 \
  -r 10000
```

### Minimal command — with manual reference files

```bash
crush \
  -i sample.hic \
  -g hg38.chrom.sizes \
  -a hg38_genes.bed \
  -b hg38.fa \
  -r 10000
```

### With parallel processing and output prefix

```bash
# Using genome build shortcut
crush \
  -i sample.hic \
  -gb hg38 \
  -r 10000 \
  -c 8 \
  -o sample1_

# Using manual reference files
crush \
  -i sample.hic \
  -g hg38.chrom.sizes \
  -a hg38_genes.bed \
  -b hg38.fa \
  -r 10000 \
  -c 8 \
  -o sample1_
```

### Using a pre-computed eigenvector

```bash
crush \
  -i sample.hic \
  -g hg38.chrom.sizes \
  -a hg38_genes.bed \
  -b hg38.fa \
  -r 10000 \
  -e eigenvector.bedgraph \
  -c 8 \
  -o sample1_
```

### Using an `.mcool` file

```bash
crush \
  -i sample.mcool \
  -g hg38.chrom.sizes \
  -a hg38_genes.bed \
  -b hg38.fa \
  -r 10000 \
  -c 8
```

### Limiting the resolution range

By default CRUSH starts from 2.5 Mb. To start from a finer coarsest resolution (faster, fewer resolution steps):

```bash
crush \
  -i sample.hic \
  -g hg38.chrom.sizes \
  -a hg38_genes.bed \
  -b hg38.fa \
  -r 10000 \
  -m 1000000 \
  -c 8
```

### Keeping temporary files for debugging

```bash
crush ... -C 0
```

---

## 7. Complete Parameter Reference

### Always required

| Flag | Long form | Description |
|---|---|---|
| `-i` | `--hic` | Hi-C file (`.hic` or `.mcool`). Local path or HTTPS URL. |
| `-r` | `--res` | Target resolution in bp (e.g., `10000`). |
| `-c` | `--cpu` | Number of CPU threads (default: `1`). |

### Reference files — choose one of two paths

**PATH A — genome build shortcut** (supported builds: `hg19`, `hg38`, `mm10`, `mm9`; res ≥ 500 bp only)

| Flag | Long form | Description |
|---|---|---|
| `-gb` | `--genomebuild` | Auto-downloads chr.sizes, genes.bed, and Bbins.bed from JRowleyLab GitHub. Overridden by any explicit `-g`/`-a`/`-b` flags. Not available for res < 500 bp. |

**PATH B — manual reference files** (any genome, any resolution)

| Flag | Long form | Description |
|---|---|---|
| `-g` | `--genomesize` | Chromosome sizes file (2 columns: chr name, size in bp). |
| `-a` | `--initialA` | BED file for A-compartment initialization (e.g., gene annotations). |
| `-b` | `--initialB` | Genome FASTA **or** pre-computed Bbins BED for B-compartment initialization. FASTA: CRUSH generates Bbins at 500 bp (res ≥ 500 bp) or at input resolution (res < 500 bp). BED: used directly as B-compartment seeds at any resolution. |

### Output

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-o` | `--outpre` | *(none)* | Prefix for all output files. |
| `-n` | `--no-merge` | `0` | Set to `1` to keep each resolution as a separate output file rather than merging. |
| `-t` | `--trackline` | `1` | Set to `0` to suppress the UCSC bedGraph track header line. |

### Performance

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-c` | `--cpu` | `1` | Number of parallel CPU threads. Set to number of chromosomes or available cores, whichever is smaller. |
| `-C` | `--cleanup` | `1` | Set to `0` to retain the temporary working directory after completion. |

### Normalization

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-N` | `--norm` | `NONE` | Normalization scheme when reading Hi-C contacts. Options: `NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`. |

### Resolution Control

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-m` | `--maxres` | `2500000` | Coarsest resolution to begin the resolution walking process. CRUSH includes all available resolutions between this value and your target `-r`. |
| `-Z` | `--eigenres` | `100000` | Resolution (bp) for eigenvector calculation. 100 kb is recommended for most genomes. |

### Analysis Behavior

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-e` | `--eigenfile` | *(auto)* | Pre-computed eigenvector bedGraph. If not supplied, CRUSH calculates eigenvectors automatically. |
| `-w` | `--window` | `5` | Sliding window size (in kb) for CRUSH score averaging. Set to `1` to disable. Set to `0` for legacy auto-calculation from sequencing depth. |
| `-S` | `--switch` | `1` | Whether to re-evaluate A/B bin assignments between resolutions. Set to `0` to bypass (not recommended). |
| `-d` | `--distance` | `0` | Minimum genomic distance from diagonal (bins) to include in scoring. |
| `-u` | `--upperlim` | `0` | Maximum genomic distance from diagonal (bins) to include. `0` = no upper limit (use whole chromosome). |
| `-T` | `--threshold` | `10` | Distance-normalized upper cap on individual contact values. Acts as an outlier filter. |
| `-l` | `--lowerthresh` | `0` | Lower threshold for contact value filtering. |
| `-s` | `--smoothing` | `0` | Set to `1` to apply a 3-bin sliding window smoothing at compartment boundaries. |
| `-x` | `--exclbed` | *(none)* | BED file of regions to exclude from scoring (e.g., centromeres, ENCODE blacklist). |
| `-A` | `--adjustment` | `0` | Adjusts CRUSH values to balance A and B score distributions. **Do not use when comparing samples** — this changes the absolute scale of scores and makes them incomparable across datasets. |

### Statistical Options

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-p` | `--pcalculation` | `1` | Set to `0` to skip p-value and q-value calculation. Speeds up runs when statistics are not needed. |
| `-q` | `--qvalue` | `0.05` | Q-value threshold for the filtered output file. Set to `0` to disable filtering and retain all bins. |

### Advanced / Rarely Changed

| Flag | Long form | Default | Description |
|---|---|---|---|
| `-v` | `--verbose` | `0` | Set to `1` for detailed progress messages at each step. |
| `-f` | `--tmpfolder` | `CRUSHtmp_RANDOM` | Custom name for the temporary working directory. Must not already exist. |
| `-E` | `--endZ` | *(auto)* | End Z-score normalization threshold for chromosome-end bins. |
| `-use` | | `o` | `o` = overwrite/recalculate GI tracks. `u` = reuse existing tracks from a prior run (useful when re-merging resolutions without recomputing from scratch). |

---

## 8. Output Files

CRUSH writes all output files to the **directory from which you ran the command**. Temporary files are cleaned up automatically unless `-C 0` is set.

### Primary output files

**`{prefix}CRUSHparamters.txt`**  
A plain-text record of all parameters used in the run. Keep this file alongside your results for reproducibility.

**`{prefix}mergedCrush_{res}.bedgraph`**  
The main output. A 4-column bedGraph of CRUSH scores normalized to a mean absolute value of 100. Positive = A compartment; negative = B compartment. Includes a UCSC track header.

> Unlike eigenvector-based methods, you never need to flip these scores. A is always positive and B is always negative.

**`{prefix}mergedqvalue_{res}.bedgraph`**  
Benjamini–Hochberg corrected q-values per genomic bin at the finest resolution.

**`{prefix}mergedCrush_{res}_qfiltered_reprocess.bedgraph`**  
CRUSH scores retained only for bins passing the q-value threshold (`-q`, default 0.05). Note: this filter can be overly stringent. Excellent compartment results are commonly obtained from the unfiltered `mergedCrush` file, so we recommend examining both.

### Interpreting CRUSH scores

| Score | Interpretation |
|---|---|
| Large positive | Strong A compartment — high confidence active chromatin |
| Small positive | Weak A compartment — check q-value |
| ~0 | Boundary region or ambiguous assignment |
| Small negative | Weak B compartment — check q-value |
| Large negative | Strong B compartment — high confidence inactive chromatin |

### Loading into a genome browser

The output bedGraph files include UCSC track headers and can be loaded directly into:
- **UCSC Genome Browser** — Genomes → Add custom tracks
- **IGV** — File → Load from File
- **WashU Epigenome Browser** — Upload as bedGraph track

The default track configuration uses `viewLimits=-150:150` with a red (A compartment) / black (B compartment) color scheme (`color=204,0,0 altColor=0,0,0`).

---

## 9. Test Dataset Walkthrough

A test dataset for chromosomes 17–19 of hg19 is provided in `examples/TestData/`.

### Step 1: Decompress the test files

```bash
cd examples/TestData
gunzip *.gz
cd ../..
```

### Step 2: Run with FASTA-based B initialization

```bash
crush \
  -i examples/TestData/hg19_c17_18_19_1kb.hic \
  -g examples/TestData/hg19_c17_18_19.sizes \
  -a examples/TestData/hg19_c17_18_19_genes.bed \
  -b examples/TestData/hg19_c17_18.fa \
  -r 10000 \
  -c 4 \
  -o test_fasta_
```

### Step 3: Alternative — Run with pre-computed eigenvector and B-bins

If you want to skip eigenvector calculation entirely:

```bash
crush \
  -i examples/TestData/hg19_c17_18_19_1kb.hic \
  -g examples/TestData/hg19_c17_18_19.sizes \
  -a examples/TestData/hg19_c17_18_19_genes.bed \
  -b examples/TestData/Bbins_hg19_c17_18_19.bed \
  -r 10000 \
  -e examples/TestData/Eigen_100kb_c17_18_19.bedgraph \
  -c 4 \
  -o test_eigen_
```

### Expected output

After a successful run you should see:
```
test_fasta_mergedCrush_10000.bedgraph
test_fasta_mergedqvalue_10000.bedgraph
test_fasta_mergedCrush_10000_qfiltered_reprocess.bedgraph
test_fasta_CRUSHparamters.txt
```

Load `test_fasta_mergedCrush_10000.bedgraph` into IGV or UCSC to verify that chr17–19 show the expected A/B compartment pattern. Positive values (A compartment) should align with gene-dense, GC-rich regions; negative values (B compartment) with gene-poor, GC-low regions.

---

## 10. Tips, Recommendations, and Common Issues

### Choosing resolution

| Hi-C sequencing depth | Recommended target resolution |
|---|---|
| < 500M reads | 25 kb or coarser |
| 500M – 1B reads | 10 kb |
| 1B – 3B reads | 5 kb |
| > 3B reads | 1–2 kb |

Attempting very fine resolution with shallow sequencing will produce sparse, low-confidence compartment calls. CRUSH's q-value filtering will exclude most of these, but the signal will be limited.

### Chromosome name matching

CRUSH automatically detects chromosome prefix mismatches between your Hi-C file and your reference files (e.g., `chr1` vs `1`, or `chr1` vs `CHR1`). When a mismatch is found, CRUSH converts the reference files to match your Hi-C chromosome naming convention and prints a notice. No manual file conversion is required in most cases.

If CRUSH output is empty or unexpected after auto-conversion, the most likely cause is that your Hi-C file itself uses an inconsistent naming convention (e.g., a mix of `chr1` and `1` within the same file). You can verify what chromosome names your Hi-C file uses:

```bash
# Check .hic file chromosomes
python3 -c "import hicstraw; h = hicstraw.HiCFile('sample.hic'); print(h.getChromosomes())"

# Check .mcool file chromosomes
cooler info sample.mcool | grep chromnames
```

### Choosing the coarsest resolution (`-m`)

The default is 2.5 Mb. If compartment calls look unusual at low resolution — for example, entire small chromosomes called as one compartment — try `-m 1000000` to start the walk at 1 Mb. This is especially relevant for smaller chromosomes or non-human genomes with shorter chromosomes.

### CPU threads (`-c`)

Set `-c` to the number of chromosomes in your genome or the number of available CPU cores, whichever is smaller. Each thread processes one chromosome at a time. For human data (24 chromosomes), `-c 24` will maximize parallelism if resources allow.

### Q-value filtering caveat

The q-value filtered output (`_qfiltered_reprocess.bedgraph`) can be overly stringent, particularly for datasets with lower sequencing depth or at very fine resolutions. We recommend always examining the unfiltered `mergedCrush` output as your primary result, and using the filtered file as a conservative complement.

### Normalization

CRUSH reads raw contacts by default (`-N NONE`). If your Hi-C data has known coverage biases, try `KR` normalization. Do not mix normalization settings across samples you intend to compare.

### Excluding problematic regions

Centromeres, telomeres, and assembly gaps can introduce spurious compartment signals. Use `-x` to exclude them:

```bash
# Example: exclude ENCODE blacklist
-x hg38_blacklist.bed
```

### When NOT to use `-A` (adjustment)

The `-A 1` flag adjusts CRUSH scores to balance the A and B score distributions within a single sample. This can improve visualization of a single dataset but **changes the absolute magnitude of scores**. Never use `-A 1` when you plan to compare CRUSH scores across samples.

### Keeping temp files

If a run fails partway through, use `-C 0` and then rerun with `-use u` to reuse any GI tracks already calculated rather than recomputing from scratch.

### Common error messages

**`WARNING: no dumped reads for chrN at XXXXX`**  
The Hi-C file has no data for that chromosome at that resolution. Expected for small scaffolds or unsequenced regions. CRUSH skips these and continues.

**`Error: Input must be a valid .hic or .mcool file`**  
Check that your file path is correct and ends in `.hic` or `.mcool`.

**`Missing one or more required arguments`**  
`-i` and `-r` are always required. For reference files, either supply `-gb` (for a supported build at res ≥ 500 bp) or all three of `-g`, `-a`, `-b`.

**`size file is incorrect format`**  
The chromosome sizes file must have exactly two tab-separated columns with no header line.

**Empty output or all zeros**  
CRUSH auto-corrects chr prefix mismatches, so this is now less common. If it occurs, check that your Hi-C file uses a consistent chromosome naming convention throughout (not a mix of `chr1` and `1` in the same file). Also run with `-v 1` to see where processing stops.

---

## 11. Comparing Samples

CRUSH scores are directly comparable across samples as long as:
- You do **not** use `-A 1` (adjustment changes absolute scale)
- You use the **same normalization** (`-N`) for all samples
- You use the **same eigenvector file** (`-e`) for all samples

### Recommended workflow for differential compartment analysis

1. Calculate a reference eigenvector from one sample (or a merged contact matrix) at the desired resolution.
2. Run CRUSH on every sample with identical parameters, supplying the reference eigenvector via `-e`.
3. Compare the `mergedCrush_{res}.bedgraph` files directly. Bins that shift from negative to positive (or vice versa) represent compartment switching events.

### Extracting the eigenvector from a CRUSH run

The eigenvector file is generated inside the temporary directory during the run. To keep it:

```bash
crush -i sample1.hic ... -C 0 -f sample1_tmp
# After the run:
cp sample1_tmp/EV_full_genome.bedgraph reference_eigenvector.bedgraph

# Use for all subsequent samples:
crush -i sample2.hic ... -e reference_eigenvector.bedgraph
crush -i sample3.hic ... -e reference_eigenvector.bedgraph
```

---

## 12. FAQ

**Q: What is the CRUSH score exactly?**  
A: The CRUSH score for a bin is the weighted difference between that bin's mean interaction Z-score with A-type bins (iA) and its mean interaction Z-score with B-type bins (iB), normalized by the total number of contacts contributing to the calculation. Higher positive = stronger A identity; more negative = stronger B identity.

**Q: Do I need to flip CRUSH scores like I do with eigenvectors?**  
A: No. CRUSH always assigns positive scores to A compartments and negative scores to B compartments. The algorithm handles sign orientation internally using gene density and GC content as reference.

**Q: Can I use CRUSH with organisms other than human?**  
A: Yes. You need an appropriate chromosome sizes file, gene annotation BED (for A initialization), and either a genome FASTA or a B-compartment BED. CRUSH has no human-specific components.

**Q: My output bedGraph is empty. What happened?**  
A: First, check chromosome name consistency across all input files — this is the most common cause. Second, try setting `-q 0` to disable q-value filtering and see if the unfiltered output has data. If it is still empty, run with `-v 1` to see verbose output and identify where the process stops.

**Q: Can I run CRUSH on a single resolution?**  
A: Yes. If your Hi-C file only contains one resolution in the range you specify (e.g., `-r 10000 -m 10000`), CRUSH processes only that resolution. Multi-resolution walking generally produces more confident calls, so it is recommended when multiple resolutions are available.

**Q: The q-filtered output is very sparse. Should I be worried?**  
A: Not necessarily. The q-value filter can be overly stringent, particularly at fine resolutions or with lower sequencing depth. Use the unfiltered `mergedCrush` bedGraph as your primary output and treat the filtered file as a conservative subset.

**Q: How do I load the output into IGV?**  
A: Open IGV, select your genome build, then go to File → Load from File and select the `.bedgraph` output. The UCSC track header ensures the track displays with appropriate color and scale.

**Q: CRUSH is slow. How can I speed it up?**  
A: Increase `-c` (more threads), increase `-m` to reduce the number of resolution steps, increase `-Z` (coarser eigenvector resolution), and set `-p 0` to skip p-value calculation if statistics are not needed.

**Q: Can I use Docker instead of installing locally?**  
A: Yes. A Docker configuration is provided in the `docker/` directory of the repository.

---

*For additional help, please open a [GitHub Issue](https://github.com/JRowleyLab/CRUSH/issues).*

*CRUSH v1.0 | JRowleyLab | Manuscript in preparation*
