#!/usr/bin/env bash

################################################################################
# CRUSH BBINS GENERATOR WRAPPER
# Isolates the GC-content-based B-compartment seed generation from CRUSH.
################################################################################

set -e # Exit on error

# Default parameters
cpu=1
minres=500
verbose=1

# Usage function
usage() {
    echo "Usage: $0 -g <chrom.sizes> -b <genome.fasta> -r <resolution_bp> [-c <cpu_threads>] [-o <output_bed>]"
    echo ""
    echo "Parameters:"
    echo "  -g : Chromosome sizes file (2 columns: chr_name, size)"
    echo "  -b : Genome FASTA file"
    echo "  -r : Resolution in bp (must be >= 500bp)"
    echo "  -c : Number of CPU threads for parallel processing (default: 1)"
    echo "  -o : Output Bbins file name (default: Bbins.bed)"
    exit 1
}

# Parse arguments
while getopts "g:b:r:c:o:h" opt; do
    case "$opt" in
        g) sizefile=$OPTARG ;;
        b) fastafile=$OPTARG ;;
        r) minres=$OPTARG ;;
        c) cpu=$OPTARG ;;
        o) outfile=$OPTARG ;;
        h|*) usage ;;
    esac
done

# Validate required inputs
if [ -z "$sizefile" ] || [ -z "$fastafile" ] || [ -z "$minres" ]; then
    echo "Error: Missing required arguments."
    usage
fi

if [ "$minres" -lt 500 ]; then
    echo "Error: Resolution must be 500bp or greater for this pre-computed logic."
    exit 1
fi

# Set default output name if not provided
outfile=${outfile:-"Bbins.bed"}

# Setup temporary names
resbins="tmp_size_bins.bed"
gcfile="tmp_gc_bins.txt"
gc_g_ga="tmp_gc_g_ga_bins.bed"

# Token-based semaphore for parallel processing
open_sem() {
    local n=$1
    local pipe
    pipe="pipe-$$-$RANDOM"
    mkfifo "$pipe"
    exec 3<>"$pipe"
    rm -f "$pipe"
    local i
    for (( i=0; i<n; i++ )); do
        printf %s 000 >&3
    done
}

run_with_lock() {
    local x
    read -u 3 -n 3 x && ((0==x)) || {
        printf '%s' 000 >&3
        return 1
    }
    (
        "$@"
        printf '%.3d' 0 >&3
    ) &
}

# Clean up traps
cleanup() {
    rm -f "$resbins" "$gcfile" "$gc_g_ga"
    rm -rf gc_tmp
}
trap cleanup EXIT

# 1. Tile the genome into 500bp bins if res >= 500 (CRUSH dynamic floor resolution)
if [ "$minres" -ge 500 ]; then
    cat "$sizefile" | awk -v myres=500 '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > "$resbins"
    [ "$verbose" -gt 0 ] && echo "Genomic bins calculated at 500bp resolution."
else
    cat "$sizefile" | awk -v myres=$minres '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > "$resbins"
    [ "$verbose" -gt 0 ] && echo "Genomic bins calculated at ${minres}bp resolution."
fi

# 2. Process GC content in parallel by chromosome
[ "$verbose" -gt 0 ] && echo "Calculating GC content in parallel using $cpu thread(s)..."
mkdir -p gc_tmp
awk '{print >> "gc_tmp/"$1".bed"; close("gc_tmp/"$1".bed")}' "$resbins"

open_sem $cpu

for chrom in $(cut -f 1 "$sizefile"); do
    run_with_lock bash -c "
        chrom='${chrom}'
        fastafile='${fastafile}'
        if [ -f gc_tmp/\${chrom}.bed ]; then
            bedtools nuc -fi \"\${fastafile}\" -bed gc_tmp/\${chrom}.bed > gc_tmp/\${chrom}.gc 2> gc_tmp/\${chrom}.err
            if [ \$? -ne 0 ]; then
                echo \"Error: bedtools nuc failed for chromosome \${chrom}\" >&2
                cat gc_tmp/\${chrom}.err >&2
            fi
        fi
    "
done

wait

# Check processing completeness
expected_files=$(cut -f 1 "$sizefile" | wc -l)
actual_files=$(ls gc_tmp/*.gc 2>/dev/null | wc -l)
if [ "$actual_files" -lt "$expected_files" ]; then
    echo "Warning: Only obtained GC content for $actual_files chromosomes out of $expected_files expected."
fi

# Combine outputs
cat gc_tmp/*.gc > "$gcfile"
rm -rf gc_tmp

# 3. Filter and parse GC data
cat "$gcfile" | grep -v user | cut -f 1-5 | awk '{print $0"\t"($4+$5)}' | \
    awk '{if ($6 > 0) print $0}' | \
    awk '{print $0"\t"($5/$6)}' > "$gc_g_ga"

# 4. Compute standard deviation-based threshold for B-compartment mapping
gc_thresh=$(cat "$gc_g_ga" | \
    awk '{s+=$7; ss+=$7*$7; linecount+=1} END{print m=s/linecount, sqrt(ss/linecount-m^2)}' | \
    awk '{print $1 - 1*$2}')

if [ "$verbose" -gt 0 ]; then
    echo "Calculated GC threshold: $gc_thresh"
fi

# 5. Extract final Bbins falling below the background threshold
cat "$gc_g_ga" | awk -v thresh=$gc_thresh '{if ($7 < thresh) print $0}' > "$outfile"

echo "Success! Finished generating Bbins file: $outfile"
