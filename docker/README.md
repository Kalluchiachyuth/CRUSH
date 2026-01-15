# CRUSH Docker

This folder provides the official, fully reproducible Docker container for CRUSH
(Chromatin Recursive Ultra-fine Segmentation of Hi-C).

CRUSH identifies fine-scale A/B chromatin compartments from Hi-C and Micro-C data using a recursive multi-resolution refinement strategy.

This Docker container freezes the operating system, Python runtime, and all dependencies to guarantee identical results across computing environments.

## Why Docker?

Running CRUSH natively requires installing:

bedtools, cooler, hic-straw, numpy, scipy, statsmodels, and several system libraries — a process that is error-prone and platform-dependent.

The Docker container packages everything into a single portable scientific environment that runs with one command.

## Quick start

# 1. Install Docker

Download Docker for your operating system:

https://docs.docker.com/get-docker/

# 2. Download CRUSH

git clone https://github.com/Kalluchiachyuth/CRUSH.git
cd CRUSH/docker

# 3. Build the CRUSH image
## This builds the complete CRUSH computational environment.
docker build -t crush:1.0 .   

# 4. Run CRUSH
## Mounting the working directory and run CRUSH.
docker run --rm -it \
  -v $(pwd):/data \
  -w /data \
  crush:1.0 \
  -i sample.hic -g hg19.sizes -r 50000

| Option            | Meaning                                        |
| ----------------- | ---------------------------------------------- |
| `-v $(pwd):/data` | Mount your current folder inside the container |
| `-w /data`        | Start CRUSH inside the mounted folder          |
| `--rm`            | Remove container after execution               |
| `-it`             | Interactive terminal                           |

| Format   | Supported |
| -------- | --------- |
| `.hic`   | Yes       |
| `.mcool` | Yes       |
| BED      | Yes       |

## Example
docker run --rm -it \
  -v $(pwd):/data \
  -w /data \
  crush:1.1 \
  -i HFFc6.hic -g hg19.sizes -r 10000

## Reproducibility Guarantee

This container freezes:

• Linux OS
• Python 3.10
• bedtools
• cooler
• hic-straw
• numpy / scipy / statsmodels / tqdm
• CRUSH source code

ensuring bit-for-bit reproducibility across systems.

## Versioning

Each CRUSH release is distributed as a versioned Docker image:

docker build -t crush:1.1 .
docker build -t crush:1.2 .

Older versions remain fully reproducible.

## Citation

If you use CRUSH, please cite:

(Your paper here)

## License

See LICENSE file.
