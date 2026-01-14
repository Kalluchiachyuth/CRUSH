
################################################################################
#                                                                              #
#                    CRUSH - Chromatin Compartment Analysis                   #
#                      Hi-C Resolution-based Compartment Caller                #
#                                                                              #
################################################################################
#
# CRUSH (Compartment Recognition Using Shifted Hierarchies) analyzes Hi-C data
# to identify A/B chromatin compartments at multiple resolutions.
#
# Key Features:
# - Multi-resolution compartment calling with recursive refinement
# - Eigenvector-based initialization of A/B states
# - Parallel processing for improved performance
# - Statistical significance testing with FDR correction
#
# Background on A/B Compartments:
# - A compartments: Gene-rich, open chromatin, active transcription (positive GI)
# - B compartments: Gene-poor, closed chromatin, silent (negative GI)
# - GI Score: Genome Interaction score measuring preference for A vs B interactions
#
################################################################################




################################################################################
#                         VERSION & DEFAULTS                                 #
################################################################################

#CRUSH Script Version
version=1.0 # A version where EigenVector approach is considered for initialization states and recursive iteration of resolutions is implemented.


################################################################################
#                    PARAMETER INITIALIZATION                                #
################################################################################
# These variables define default values for all CRUSH parameters
# They can be overridden via command-line arguments

#Initial variables
hicpath=0
myres=0
help=0
genome=0
cpu=1
doNotMerge=0
adjustment=0
distance=0
upperlim=0
trackline=1
threshold=10
totreads=0
lowerthresh=0
window=0
verbose=0
switch=1
qthresh=0.05
whichchoice=`echo "o"`
cleanup=1
endZ=0
exclbed=0
crushdir=`echo "CRUSHtmp_""$RANDOM"`
outpre=""
norm=NONE #Default
coarsestres=2500000
pcalculation=1
smoothing=0
eigenres=100000 #Default resolution for Eigen calculations


################################################################################
#                       SYSTEM UTILITIES                                     #
################################################################################

# Trap for cleanup
function shutdown() {
    tput cnorm # Reset cursor
}
trap shutdown EXIT

# Move cursor back
function cursorBack() {
    echo -en "\033[$1D"
}


################################################################################
#                    PARALLEL PROCESSING FRAMEWORK                           #
################################################################################
#
# These semaphore functions implement a token-based system for controlling
# parallel job execution. They ensure we don't exceed the specified CPU limit.
#
# How it works:
# 1. open_sem creates a FIFO pipe with N tokens (where N = number of CPUs)
# 2. run_with_lock consumes a token before starting a job
# 3. When job completes, it returns the token to the pool
# 4. If no tokens available, jobs wait until one is freed
#
################################################################################

# Semaphore functions for parallel processing (used globally)
open_sem() {  # Initialize semaphore with N tokens (N = CPU count)
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

run_with_lock() {  # Execute with semaphore lock (parallel processing)
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
        "$@" &
        local pid=$!
        wait "$pid"
        local exit_code=$?
        printf '%.3d' "$exit_code" >&3
    ) &
}


################################################################################
#                        HELP & USAGE FUNCTIONS                              #
################################################################################

# Display usage
function usage {
    echo -e "\n\nusage : crush -i HIC  -g SIZEFILE -a ABED -b BBED | FASTA -r FINERESOLUTION [-e EIGENVECTORBED] [-cpu CPU] [-w WINDOW] [-h]"
    echo -e "Use option -h|--help for more information"
}

# Help menu function
function help {
    usage;
    echo ""
    echo "================================================================================"
    echo "                         CRUSH - Hi-C Compartment Analysis"
    echo "================================================================================"
    echo ""
    echo "CRUSH analyzes Hi-C data to identify chromatin compartments at multiple"
    echo "resolutions. It will create a temporary folder in your current directory,"
    echo "so please ensure you have write access."
    echo ""
    echo "For full documentation, visit: https://github.com/JRowleyLab/CRUSH/"
    echo ""
    echo "================================================================================"
    echo ""

    echo "GETTING HELP:"
    echo "-------------"
    echo ""
    echo "  -h, --help              Display this help menu"
    echo ""
    echo "================================================================================"
    echo ""

    echo "REQUIRED PARAMETERS:"
    echo "--------------------"
    echo ""
    echo "  -i, --hic FILE          Input Hi-C file (.hic, .cooler, or .mcool format)"
    echo "                          Example: '/path/to/file.hic'"
    echo ""
    echo "  -g, --genomesize FILE   Chromosome size file (2 columns: chr_name, size)"
    echo "                          Example: 'hg19.chrom.sizes'"
    echo ""
    echo "  -a, --initialA FILE     BED file for initializing A compartments"
    echo "                          Example: Gene annotations like 'hg19genes.bed'"
    echo ""
    echo "  -b, --initialB FILE     BED or FASTA file for initializing B compartments"
    echo "                          If FASTA: initializes from GC content"
    echo ""
    echo "  -e, --eigenfile FILE    Eigenfile to initialize A and B states (optional)"
    echo "                          If not specified, CRUSH will calculate eigenvectors"
    echo "                          automatically and select the best principal component"
    echo ""
    echo "  -r, --res VALUE         Desired resolution for analysis"
    echo "                          Example: 10000 (for 10kb resolution)"
    echo ""
    echo "================================================================================"
    echo ""

    echo "OPTIONAL PARAMETERS:"
    echo "--------------------"
    echo ""
    echo "Output Options:"
    echo "  -o, --outpre PREFIX     Prefix for output files"
    echo "                          Default: No prefix"
    echo ""
    echo "  -n, --no-merge [0|1]    Keep each resolution as separate file"
    echo "                          0 = Merge resolutions (default)"
    echo "                          1 = Keep separate files"
    echo ""
    echo "  -t, --trackline [0|1]   Include bedGraph track header"
    echo "                          0 = Disable header"
    echo "                          1 = Include header (default)"
    echo ""
    echo "Performance Options:"
    echo "  -c, --cpu NUMBER        Number of CPU threads to use"
    echo "                          Default: 1"
    echo ""
    echo "  -C, --cleanup [0|1]     Clean up temporary files after completion"
    echo "                          Default: 1 (cleanup enabled)"
    echo ""
    echo "Analysis Options:"
    echo "  -A, --adjustment [0|1]  Adjust CRUSH values based on data distribution"
    echo "                          WARNING: Do NOT use when comparing Hi-C maps"
    echo "                          Default: 0 (disabled)"
    echo ""
    echo "  -d, --distance VALUE    Filter distance from diagonal (in bins)"
    echo "                          Default: 0 (consider all distances)"
    echo ""
    echo "  -u, --upperlim VALUE    Upper limit of distance from diagonal to consider"
    echo "                          Default: 0 (considers whole chromosome)"
    echo ""
    echo "  -T, --threshold VALUE   Distance-normalized threshold for outlier filtering"
    echo "                          Default: 10"
    echo ""
    echo "  -l, --lowerthresh VALUE Lower threshold value"
    echo "                          Default: 0"
    echo ""
    echo "  -w, --window VALUE      Sliding window size for score averaging"
    echo "                          Default: Auto-calculated from sequencing depth"
    echo "                          Set to 1 to disable sliding window"
    echo ""
    echo "  -S, --switch [0|1]      Enable re-initialization of bins between resolutions"
    echo "                          0 = Bypass re-initialization"
    echo "                          1 = Enable (default)"
    echo ""
    echo "  -N, --norm TYPE         Normalization scheme"
    echo "                          Options: NONE (default), VC, VC_SQRT, KR, SCALE"
    echo ""
    echo "  -m, --maxres VALUE      Coarsest resolution to consider"
    echo "                          Default: 2500000 (checks 1000000 to your resolution)"
    echo ""
    echo "  -Z, --eigenres VALUE    Resolution for eigenvector calculations"
    echo "                          Default: 100000 (100kb)"
    echo ""
    echo "  -s, --smoothing [0|1]   Apply smoothing to compartment boundaries"
    echo "                          0 = No smoothing (default)"
    echo "                          1 = Enable smoothing"
    echo ""
    echo "Statistical Options:"
    echo "  -p, --pcalculation [0|1] Calculate p-values"
    echo "                           0 = Skip p-value calculation"
    echo "                           1 = Calculate (default)"
    echo ""
    echo "  -q, --qvalue VALUE       Q-value threshold for significance filtering"
    echo "                           Default: 0.05"
    echo "                           Set to 0 to disable filtering"
    echo ""
    echo "Advanced Options:"
    echo "  -v, --verbose [0|1]     Enable verbose output messages"
    echo "                          Default: 0 (minimal output)"
    echo ""
    echo "  -x, --exclbed FILE      BED file of regions to exclude from analysis"
    echo ""
    echo "  -E, --endZ VALUE        End Z-score value"
    echo ""
    echo "  -u, --use [o|u]         Use existing GI tracks from previous calculations"
    echo "                          'o' = Overwrite/recalculate (default)"
    echo "                          'u' = Use existing (useful for merging resolutions)"
    echo ""
    echo "  -f, --tmpfolder NAME    Custom name for temporary folder"
    echo "                          Default: CRUSHtmp_RANDOM"
    echo "                          (Must not already exist)"
    echo ""
    echo "================================================================================"
    echo ""
    echo "EXAMPLES:"
    echo ""
    echo "  Basic usage:"
    echo "    crush -i data.hic -g hg19.sizes -a genes.bed -b genome.fa -r 10000"
    echo ""
    echo "  With multiple threads and custom output:"
    echo "    crush -i data.hic -g hg19.sizes -a genes.bed -b genome.fa -r 10000 \\"
    echo "          -c 8 -o myproject_"
    echo ""
    echo "  With smoothing and custom resolution range:"
    echo "    crush -i data.hic -g hg19.sizes -a genes.bed -b genome.fa -r 5000 \\"
    echo "          -s 1 -m 1000000 -c 16"
    echo ""
    echo "================================================================================"
}

# Modified URL validation function
validate_url() {
    local url=$1
    if ! command -v curl &> /dev/null; then
        echo "Error: curl is not installed. Please install curl to validate remote files."
        exit 1
    fi

    curl --output /dev/null --silent --head --fail "$url"
    return $?
}

# Parsing command line arguments
while test $# -gt 0; do
    case "$1" in
        -h|--help)
            help
            exit 0
            ;;
        -i|--hic)
            shift  # Move to the next argument which contains the actual path
            if [ $# -eq 0 ]; then
                echo "Error: No file specified after -i|--hic"
                exit 1
            fi
            hicpath="$1"
            # Check if it's a URL or local file
            if [[ "$hicpath" =~ ^https?:// ]]; then
                # Validate URL
                if ! validate_url "$hicpath"; then
                    echo "Error: Unable to access URL: $hicpath"
                    exit 1
                fi
                echo "Using remote file: $hicpath"
            else
                # Check if it's a valid local file
                if [[ ! -f "$hicpath" ]]; then
                    echo "Error: Local file not found: $hicpath"
                    exit 1
                fi
                hicpath=$(readlink -e "$hicpath")
                echo "Using local file: $hicpath"
            fi
            ;;
        -o|--outpre)
            outpre=$2
            ;;
        -r|--res)
            res=$2
            ;;
        -g|--genomesize)
            sizefile=`readlink -e $2`
            ;;
        -a|--initialA)
            genesfile=`readlink -e $2`
            ;;
        -b|--initialB)
            fastafile=`readlink -e $2`
            ;;
        -e|--eigenfile)
            eigenfile=`readlink -e $2`
            shift
            ;;
        -q|--qvalue)
            qthresh=$2
            ;;
        -x|--exclbed)
            exclbed=`readlink -e $2`
            ;;
        -c|--cpu)
            cpu=$2
            ;;
        -n|--no-merge)
            doNotMerge=$2
            ;;
        -A|--adjustment)
            adjustment=$2
            ;;
        -d|--distance)
            distance=$2
            ;;
        -u|--upperlim)
            upperlim=$2
            ;;
        -l|--lowerthresh)
            lowerthresh=$2
            ;;
        -t|--trackline)
            trackline=$2
            ;;
        -T|--threshold)
            threshold=$2
            ;;
        -w|--window)
            window=$2
            ;;
        -v|--verbose)
            verbose=$2
            ;;
        -S|--switch)
            switch=$2
            ;;
        -N|--norm)
            norm=$2
            ;;
        -use)
            whichchoice=$2
            ;;
        -f|--tmpfolder)
            crushdir=$2
            ;;
        -C|--cleanup)
            cleanup=$2
            ;;
        -E|--endZ)
            endZ=$2
            ;;
        -m|--maxres)
            coarsestres=$2
            ;;
        -p|--pcalculation)
            pcalculation=$2
            ;;
        -s|--smoothing)
            smoothing=$2
            ;;
        -Z|--eigenres)
            eigenres=$2
            ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$hicpath" ]]; then
    echo "You must specify a .hic or .mcool file! (local file or HTTPS link)"
    exit 1
fi

# Check file extension and URL format
if [[ "$hicpath" =~ ^https?:// ]]; then
    # For URLs, check if they end with .hic or .mcool
    if [[ ! "$hicpath" =~ \.(hic|mcool)$ ]]; then
        echo "Error: URL must point to a .hic or .mcool file"
        exit 1
    fi
else
    # For local files, check if they exist and have the right extension
    if [[ ! -f "$hicpath" ]] || [[ ! "$hicpath" =~ \.(hic|mcool)$ ]]; then
        echo "Error: Input must be a valid .hic or .mcool file"
        exit 1
    fi
fi

if [ $res == 0 ]; then
    echo "You must specify a resolution!......"
    exit 0
fi

if [ $switch == 0 ]; then
    echo "Warning: Re-evaluation & Re-iteration of A and B Bins has been deactivated. We recommend leaving this parameter alone to achieve the highest possible confidence and resolution."
fi

# Check if the size file is provided
if [ -z "$sizefile" ] || [ -z "$genesfile" ] || [ -z "$fastafile" ]; then
    echo "Missing one or more required arguments: --genomesize, --initialA, or --initialB!..."
    help
    exit 1
fi

if [ $doNotMerge == 1 ]; then
    echo "Warning: doNotMerge set. We recommend merging to achieve the highest possible confidence and resolution. Only set doNotMerge if you plan to merge the individual resolutions later."
fi

if [ $adjustment == 1 ]; then
    echo "Warning: adjusting end compartmental values which may result in some loss of quantitative power. We recommend adjustment only for troubleshooting or when not comparing two samples."
fi

echo "CRUSH_v""$version"" --hic ""$hicpath"" --res ""$res"" --genomesize ""$sizefile"" --initialA ""$genesfile"" --initialB ""$fastafile"" --cpu ""$cpu"" --no-merge ""$doNotMerge"" --adjustment ""$adjustment"" --distance ""$distance"" --upperlim ""$upperlim"" --lowerthresh ""$lowerthresh"" --trackline ""$trackline"" --threshold ""$threshold"" --window ""$window"" --switch ""$switch"" --maxres ""$coarsestres"" --norm ""$norm"" --eigenres ""$eigenres"" --smoothing ""$smoothing" > "$outpre"CRUSHparamters.txt

# Check if the input file is .hic or .mcool
juiceorcool=`echo "$hicpath" | sed 's/\./\t/g' | sed 's/\./\t/g' | awk '{if ($NF == "hic") print 0; else if ($NF == "mcool") print 1; else print 2}'`

if [ $juiceorcool == 2 ]; then
    echo "You must choose a Hi-C file that ends in either .hic (juicer) or .mcool (cooler)"
    exit 0
elif [ $juiceorcool == 0 ]; then
    echo "Identified .hic juicer format. We will use the pythonic hicstraw to extract the data and assume you have it installed, using pip install hic-straw."
else
    echo "Identified .mcool format. We will use cooler to extract the data and assume you have it installed in your path."
fi



################################################################################
#                      CORE PROCESSING FUNCTIONS                             #
################################################################################

# Functions in the order of their precedence


# ------------------------------------------------------------------------------
# HI-C DATA EXTRACTION
# ------------------------------------------------------------------------------
# Extract contact matrices from .hic files using hicstraw (Python)
# ------------------------------------------------------------------------------

# Function to run dumper.py script

# ------------------------------------------------------------------------------
# HI-C DATA EXTRACTION
# ------------------------------------------------------------------------------
# Extract contact matrices from .hic files using hicstraw (Python)
# ------------------------------------------------------------------------------

# Function to run dumper.py script (safe under multithreading)
run_dumper() {
    local readtype=$1  # observed / oe / expected
    local norm=$2      # NONE, VC, VC_SQRT, SCALE, KR
    local hicinnie=$3
    local chrom=$4
    local res=$5
    local dumpoutie=$6
    local dumperrors=$7

    echo "Dumping $readtype reads using $norm normalization at resolution $res for $chrom."

    # Make the script file unique per output (so threads don't collide)
    local script_dir
    script_dir="$(dirname "$dumpoutie")"
    local script_base
    script_base="$(basename "$dumpoutie")"
    local script_path="${script_dir}/run_dumper_${script_base}.py"

    # Ensure directory exists
    mkdir -p "$script_dir"

    # Write a tiny, self-contained Python script *for this chromosome only*
    cat > "$script_path" << EOF
import hicstraw

readtype = "$readtype"
norm = "$norm"
hicfile = "$hicinnie"
chrom = "$chrom"
res = int("$res")

result = hicstraw.straw(readtype, norm, hicfile, chrom, chrom, "BP", res)
for record in result:
    # binX and binY in BP; counts is the contact count
    print("{0}\t{1}\t{2}".format(record.binX, record.binY, record.counts))
EOF

    # Run it
    python3 "$script_path" 2> "$dumperrors" | grep -v WARN > "$dumpoutie"

    # Optional: if you want to keep things tidy, you can uncomment:
    # rm -f "$script_path"
}



# ------------------------------------------------------------------------------
# A/B COMPARTMENT INITIALIZATION
# ------------------------------------------------------------------------------
# Process gene (A compartment) and GC/Eigen (B compartment) files
# to create initial compartment state definitions at each resolution
# ------------------------------------------------------------------------------

# Function to process the initial A (Genes)/ Eigen-A and B (GC)/ Eigen-B files
process_genes_Bbins() {
    local myres=$1
    local mychr=$2
    local genesblock=$3
    local sortedbed=$4
    local pseudoB=$5

    # Adding in the code for re-initializing Bbins for next following resolutions
    if [[ "$myres" == "$maxres" && "$genesblock" -eq 1 ]]; then

        echo "Using initial Genes and Bbins files for $maxres"

        maxres_processed=1 # Set maxres_processed to prevent future initial processing for maxres

        # Processing the genes & Bbins files
        cat $EVAstates | mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        # Bbins
        cat $EVBstates |  mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

    elif [ "$switch" -lt 1 ]; then

        # Processing the genes & Bbins files
        cat $EVAstates | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        # Bbins
        cat $EVBstates | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

    else

        echo "Using combined bins for $myres"

        # Processing the genes & Bbins files for other resolutions
        cat $combined_newAbins | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        # Bbins
        cat $combined_newBbins |  mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

    fi
}



# ------------------------------------------------------------------------------
# EIGENVECTOR ANALYSIS
# ------------------------------------------------------------------------------
# Calculate eigenvectors (principal components) from Hi-C correlation matrix
# PC1 typically separates A/B compartments; we test multiple PCs to find best
# ------------------------------------------------------------------------------

# Function to run the Eigen vector analysis
run_EigenVector() {
    local hic_file="$1"
    local chrom="$2"
    local res="$3"
    local n_components=10
    local eigendumpoutie="Eigen_dump_${chrom}.txt"
    local eigendumperrors="Eigen_errors_${chrom}.log"

    echo "Generating EV's from 1 to $n_components for $chrom at resolution $res."

    file_ext="${hic_file##*.}"
    if [[ "$file_ext" == "hic" ]]; then
        #Calling dumper for HiC files sing hicstraw
        run_dumper oe VC_SQRT "$hic_file" "$chrom" "$res" "$eigendumpoutie" "$eigendumperrors"
    elif [[ "$file_ext" == "cool" || "$file_ext" == "mcool" ]]; then
        echo "Detected cooler format: dumping with cooler CLI..."
        local cooler_uri="$hic_file"
        if [[ "$file_ext" == "mcool" ]]; then
            cooler_uri="${hic_file}::/resolutions/${res}"
        fi
        #Dumping using Cooler functionality
        cooler dump "$cooler_uri" -r "$chrom" --join | cut -f 2,5,7 > "$eigendumpoutie"
    else
        echo "Unsupported Hi-C file type: $hic_file"
        exit 1
    fi

    # Generate Python script to run eigenvector analysis
    cat << EOF > run_EigenVector.py
import numpy as np
import pandas as pd
from scipy.linalg import eigh
import sys

def save_eigenvector_to_file(eigenvector, bin_starts, res, chrom, output_file):
    with open(output_file, 'w') as f:
        for bin_start, value in zip(bin_starts, eigenvector):
            f.write(f"{chrom}\\t{bin_start}\\t{bin_start + res}\\t{value}\\n")

def process_hic_to_eigen(dump_file, res, chrom, n_components=10):
    dump_df = pd.read_csv(dump_file, delim_whitespace=True, header=None, names=['binX', 'binY', 'counts'])
    dump_df = dump_df.apply(pd.to_numeric, errors='coerce').dropna()

    binX, binY, counts = dump_df['binX'].values, dump_df['binY'].values, dump_df['counts'].values
    n_bins = max(max(binX), max(binY)) // res + 1  # Adjust the size dynamically based on binX and binY
    matrix = np.zeros((n_bins, n_bins))

    for x, y, c in zip(binX, binY, counts):
        if x // res < n_bins and y // res < n_bins:  # Add bounds check to prevent out of bounds errors
            matrix[x // res, y // res] = c
            if x != y:
                matrix[y // res, x // res] = c

    mask = np.sum(matrix, axis=0) > 0
    matrix = matrix[mask][:, mask]
    pearson_matrix = np.corrcoef(matrix)

    eigenvalues, eigenvectors = eigh(pearson_matrix)
    idx = eigenvalues.argsort()[::-1]
    eigenvectors = eigenvectors[:, idx]
    eigenvectors = eigenvectors[:, :n_components]

    for i in range(n_components):
        save_eigenvector_to_file(eigenvectors[:, i], np.arange(0, n_bins * res, res)[mask], res, chrom, f"PC{i+1}_{chrom}_{res}.bedgraph")

if __name__ == "__main__":
    dump_file = sys.argv[1]
    res = int(sys.argv[2])
    chrom = sys.argv[3]

    process_hic_to_eigen(dump_file, res, chrom)
EOF

    # Execute the Python script with the provided arguments
    python run_EigenVector.py "$eigendumpoutie" "$res" "$chrom"
}


# ------------------------------------------------------------------------------
# EIGENVECTOR ORIENTATION
# ------------------------------------------------------------------------------
# Orient eigenvector so positive values = A compartments (gene-rich)
# and negative values = B compartments (gene-poor)
# ------------------------------------------------------------------------------

# Function to check and flip sign of PC based on gene content using bedtools
check_and_flip_sign() {
    local pc_file="$1"
    local genes_file="$2"
    local output_file="sign_check_output.txt"

    echo "Checking and potentially flipping signs for $pc_file based on gene content from $genes_file..."

    # Check overlap and count positive and negative correlations with gene annotations
    pos_neg=$(awk 'BEGIN {FS=OFS="\t"} {if ($4 > 0) print $0, "positive"; else if ($4 < 0) print $0, "negative"}' "$pc_file" |
        bedtools intersect -u -a stdin -b "$genes_file" |
        sort -k1,1 -V -k5b,5b |
        bedtools groupby -i stdin -g 1,5 -c 1 -o count |
        bedtools groupby -i stdin -g 1 -c 3 -o collapse |
        sed 's/,/\t/g' |
        awk 'BEGIN {FS=OFS="\t"} {if ($2 > $3) print "-"; else print "+"}'
    )

    echo "Overlap check result for $pc_file: $pos_neg" >> "$output_file"

    # Flip signs if needed
    if [ "$pos_neg" == "-" ]; then
        echo "Flipping signs for $pc_file..." >> "$output_file"
        awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, -1 * $4}' "$pc_file" > temp && mv temp "$pc_file"
    else
        echo "No flipping needed for $pc_file." >> "$output_file"
    fi

    echo "Check and flip sign process completed for $pc_file. See $output_file for details."
}

# Function to process all PCs for a specific chromosome and resolution
process_all_pcs_for_chromosome() {
    local chrom="$1"
    local res="$2"
    local genes_file="$3"

    # Initialize or clear the output file
    echo "Starting sign check and flip for all PCs for chromosome $chrom at resolution $res." > "sign_check_output.txt"

    # Loop through all PC files for the given chromosome and resolution
    for pc_file in PC*_${chrom}_${res}.bedgraph; do
        echo "Processing $pc_file..."
        check_and_flip_sign "$pc_file" "$genes_file"
    done

    echo "All PCs for chromosome $chrom at resolution $res have been processed. See sign_check_output.txt for details."
}


# ------------------------------------------------------------------------------
# PRINCIPAL COMPONENT SELECTION
# ------------------------------------------------------------------------------
# Find PC with highest correlation to gene density minus GC content
# This identifies the component that best separates A/B compartments
# ------------------------------------------------------------------------------

# Function to find the best PC among the 10 PCs. Still working on this to make it better. Adapt from the working versions later
find_best_pc() {
    local chrom="$1"
    local res="$2"
    local genes_file="$3"
    local gc_file="$4"
    local chrom_size="$5"
    local best_pc=""
    local best_corr=-100
    local log_file="correlation_log_${chrom}_${res}.txt"

    echo "Calculating best PC for chromosome $chrom at resolution $res..." > "$log_file"

    # Get total number of bins
    local bin_size=$res
    local num_bins=$((chrom_size / bin_size))

    # Function to calculate density
    calculate_density() {
        local feature_file="$1"
        awk -v bin_size="$bin_size" -v chrom="$chrom" '
        BEGIN { OFS="\t"; for (i = 0; i < ENVIRON["num_bins"]; i++) density[i] = 0 }
        $1 == chrom && $2 >= 0 && $3 >= 0 {
            bin_start = int($2 / bin_size)
            bin_end = int($3 / bin_size)
            for (i = bin_start; i <= bin_end; i++) {
                if (i < ENVIRON["num_bins"]) density[i]++
            }
        }
        END { for (i = 0; i < ENVIRON["num_bins"]; i++) print density[i] }
        ' "$feature_file"
    }

    export num_bins=$num_bins
    local gene_density=($(calculate_density "$genes_file"))
    local gc_density=($(calculate_density "$gc_file"))

    # Iterate over each PC file
    for pc_file in PC*_"${chrom}"_"${res}".bedgraph; do
        if [[ ! -s "$pc_file" ]]; then
            echo "Skipping empty or missing PC file: $pc_file" >> "$log_file"
            continue
        fi

        # Extract PC values
        local pc_values=($(awk -v bin_size="$bin_size" -v num_bins="$num_bins" '
        BEGIN { for (i = 0; i < num_bins; i++) pc[i] = 0 }
        { bin = int($2 / bin_size); if (bin < num_bins) pc[bin] = $4 }
        END { for (i = 0; i < num_bins; i++) print pc[i] }' "$pc_file"))

        if [[ ${#pc_values[@]} -eq 0 ]]; then
            echo "No valid data in $pc_file" >> "$log_file"
            continue
        fi

        # Convert arrays to strings for Python
        local gene_density_str=$(IFS=,; echo "${gene_density[*]}")
        local gc_density_str=$(IFS=,; echo "${gc_density[*]}")
        local pc_values_str=$(IFS=,; echo "${pc_values[*]}")

        # Calculate correlations
        read -r gene_corr gc_corr <<< $(python3 - <<EOF
import numpy as np
from scipy.stats import pearsonr

try:
    pc_values = np.array([${pc_values_str}])
    gene_density = np.array([${gene_density_str}])
    gc_density = np.array([${gc_density_str}])

    gene_corr, _ = pearsonr(pc_values, gene_density)
    gc_corr, _ = pearsonr(pc_values, gc_density)

    print(gene_corr, gc_corr)
except Exception as e:
    print(0, 0)
EOF
)

        # Handle invalid correlations
        if [[ -z "$gene_corr" || -z "$gc_corr" || "$gene_corr" == "0" && "$gc_corr" == "0" ]]; then
            echo "Skipping $pc_file: Invalid correlations (gene_corr=$gene_corr, gc_corr=$gc_corr)" >> "$log_file"
            continue
        fi

        # Calculate total correlation
        local total_corr=$(echo "$gene_corr - $gc_corr" | bc -l)

        # Update best PC
        if (( $(echo "$total_corr > $best_corr" | bc -l) )); then
            best_pc="$pc_file"
            best_corr="$total_corr"
            echo "New best PC: $pc_file with total correlation: $total_corr" >> "$log_file"
        fi
    done

    # Copy the best PC to a new file
    if [[ -n "$best_pc" ]]; then
        cp "$best_pc" "best_Eigen_${chrom}_${res}.bedgraph"
        echo "Best PC for chromosome $chrom: $best_pc with correlation: $best_corr"
    else
        echo "No best PC identified for chromosome $chrom at resolution $res." >> "$log_file"
    fi
}



# ------------------------------------------------------------------------------
# EIGENVECTOR PIPELINE WRAPPER
# ------------------------------------------------------------------------------
# Master function that coordinates: EV generation → sign flipping → PC selection
# ------------------------------------------------------------------------------

# Function to run the Eigen related commands
generateEV() {
    local hic_file="$1"
    local chrom="$2"
    local res="$3"
    local genes_file="$4"
    local gc_file="$5"
    local chrom_size="$6"

    echo "Starting EV generation for chromosome $chrom at resolution $res..."

    # Run Python script to generate PCs
    run_EigenVector "$hic_file" "$chrom" "$eigenres"

    echo "PCs generated for $chrom at resolution $eigenres."

    # Check and flip sign of PC1 if necessary
    echo "Flipping the signs for $chrom PC's at resolution $eigenres..."

    process_all_pcs_for_chromosome "$chrom" "$eigenres" "$genes_file"

    # Find the best PC based on correlation
    echo "Finding best PC for $chrom at resolution $eigenres..."

    find_best_pc "$chrom" "$eigenres" "$genes_file" "$gc_file" "$chrom_size"

    echo "Best PC identification completed for $chrom at resolution $eigenres."
}


# ------------------------------------------------------------------------------
# Concatenate best eigenvectors from all chromosomes into genome-wide file
# ------------------------------------------------------------------------------

concatenate_best_pcs() {
    echo "Concatenating the best eigenvectors for all chromosomes into a full genome file."
    cat best_Eigen_*.bedgraph > EV_full_genome.bedgraph
    echo "Concatenation completed. Full genome eigenvector file: EV_full_genome.bedgraph"
}


# ------------------------------------------------------------------------------
# SMOOTHING FUNCTION
# ------------------------------------------------------------------------------
# Apply sliding window smoothing to reduce noise at compartment boundaries
# Uses shifted left/right averaging
# ------------------------------------------------------------------------------

# Function to perform shifting and smoothing on bedgraph data
run_smoothing() {
    local innieS=$1  # Input file path
    local outieS=$2  # Output file name

    # Step 1: Sort the input file and save to finalshifter
    cat "$innieS" | sort -k 1,1 -V -k 2bn,2b --stable > "$outieS"

    # Step 2: Prepare the right-shifted version of the file (skip the first line)
    cat "$outieS" | tail -n +2 > tmp1_shiftedright

    # Step 3: Prepare the left-shifted version of the file (duplicate first line)
    cat "$outieS" | awk '{if (NR == 1) print $0"\n"$0; else print $0}' > tmp2_shiftedleft

    # Step 4: Paste the original, right-shifted, and left-shifted files together
    # Calculate the average of the fourth columns from the three files and output the result
    paste "$outieS" tmp1_shiftedright tmp2_shiftedleft | awk '{print $1"\t"$2"\t"$3"\t"(($4 + $8 + $12)/3)}' > shifted_outie

    # Step 5: Move the output back to newoutie
    mv shifted_outie "$innieS"
}


# ------------------------------------------------------------------------------
# GI SCORE CALCULATION
# ------------------------------------------------------------------------------
# Calculate Genome Interaction (GI) scores for each bin:
#   GI = (A_interactions - B_interactions) / (A_interactions + B_interactions)
# Positive GI → prefers A compartments (likely A)
# Negative GI → prefers B compartments (likely B)
# ------------------------------------------------------------------------------

# Function to generate A and B compartmental states
process_oppocheck_statement() {
    local innie_genes="$1"
    local innie_gc="$2"
    local outie_oppo="$3"
    local chr="$4"
    local res="$5"
    local input_dir="$6"
    local output_dir="$7"

    # Paths for Zfile and RowZ
    local zfile="${input_dir}/zerofile"
    local rowZ="${input_dir}/RowZs"

    # Check if files exist
    if [ ! -f "$zfile" ] || [ ! -f "$rowZ" ]; then
        echo "Required file missing: $zfile or $rowZ"
        return 1  # Exit if files are missing
    fi

    # Process innie_genes, rowZ, and zfile
    mawk 'NR==FNR { c1[$1] = $2; next } { if ($1 in c1) print $2 "\t" $3; next } { print $0 }' "$innie_genes" "$rowZ" "$zfile" | awk '{ b[$1] += $2; c[$1]++; d[$1] += ($2**2) } END { for (i in b) { print i "\t" (b[i])*(c[i]) "\t" c[i] "\t" b[i]/c[$1] "\t" ((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } }' > "${output_dir}/A2removedfile" &

    # Process innie_gc, rowZ, and zfile
    mawk 'NR==FNR { c1[$1] = $2; next } { if ($1 in c1) print $2 "\t" $3; next } { print $0 }' "$innie_gc" "$rowZ" "$zfile" | awk '{ b[$1] += $2; c[$1]++; d[$1] += ($2**2) } END { for (i in b) { print i "\t" (b[i])*(c[i]) "\t" c[i] "\t" b[i]/c[$1] "\t" ((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } }' > "${output_dir}/B2removedfile" &

    wait

    # Calculate the sum of A and B bins
    myABsum=$(cat "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | mawk '{ sum += $3 } END { print sum }')

    wait

    if [ "$verbose" -gt 0 ]; then
        end=$(date +%s)
        seconds=$(echo "$end - $start" | bc)
        echo "It took $seconds seconds. Comparing A vs B scores for $chr chromosome."
        start=$(date +%s)
    fi

    # Calculate the output
    mawk -v var="$res" -v myABsum="$myABsum" '{ if ($1 in b) b[$1] -= $2; else b[$1] = $2 } END { for (i in b) { if (b[i] != 0) print i "\t" (b[i])/myABsum } }' "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | sort -k 1bn,1b --stable | mawk -v mchr="$chr" -v myres="$res" -v window="$newwindow" -v OFS="\t" 'BEGIN { slide = 1 } { mod = NR % window; if (NR <= window) { print mchr, int($1), $2; count++ } else { sum -= array[mod] } sum += $2; array[mod] = $2;} (NR % slide) == 0 { print mchr, int($1-((myres*window)/2)), sum/count }' | mawk -v myres="$res" '{ print $1 "\t" $2-int(myres/2)+myres "\t" $2+(int(myres/2))+myres "\t" $3 }' | mawk '{ if ($2 > 0) print $0 }' > "$outie_oppo" 2> "${output_dir}/errorfile"

    # Process exclusion regions if exclbed file is provided
    if [ "$exclbed" != 0 ]; then
        exclbins="${output_dir}/exclbins"
        excloutie="${output_dir}/exclout_${res}_${chr}"

        # Processing the exclusion regions
        cat "$exclbed" | mawk -v myres="$res" -v mychr="$chr" '{ if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres }' | awk -v myres="$res" '{ for (i = $2; i <= $3; i += myres) b[i] += 1 } END { for (j in b) print j "\t" b[j] }' > "$exclbins"

        wait

        # Exclude regions from the output
        mawk 'NR==FNR { c1[$1] = $2; next } { if ($2 in c1); else print $0 }' "$exclbins" "$outie_oppo" > "$excloutie"

        wait

        mv "$excloutie" "$outie_oppo"
    fi

    # Adjust compartmental values if needed
    if [[ "$oppocheck" -gt 0 && "$adjustment" -gt 0 ]]; then

        amed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1) print $4}' $innie_genes $outie_oppo | mawk '{sum +=$1} END {print sum/NR}'`
        bmed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $4}' $innie_genes $outie_oppo | mawk '{sum += $1} END {print sum/NR}'`

        wait

        cat $outie_oppo | mawk -v varA=$amed -v varB=$bmed '{print $1"\t"$2"\t"$3"\t"$4-(varA-(varB*-1))}' > $outie2

        wait

        mv $outie2 $outie_oppo

    fi

    if [ "$myres" -ge "$endZ" ]; then

        myendmean=`cat $outie_oppo | mawk '{sum += $4; cnum++} END {print sum/cnum}'`
        echo "$myendmean"
        tmpendZ="${output_dir}/tmpendnorm"
        cat $outie_oppo | mawk -v mymean=$myendmean '{print $1"\t"$2"\t"$3"\t"($4-mymean)}' > $tmpendZ

        wait

        mv $tmpendZ $outie_oppo
    fi
}


# ------------------------------------------------------------------------------
# RESOLUTION SHIFTER
# ------------------------------------------------------------------------------
# Adjust finer resolution scores based on coarser resolution trends
# This implements the hierarchical refinement approach
# ------------------------------------------------------------------------------

#Shifter function
run_midshifter() {
    local finer_res=$1
    local coarser_res=$2
    local high_innie=$3
    local coarse_innie=$4
    local shifter_outie=$5

    echo "Running shifter with Finer Resolution: $finer_res and Coarser Resolution: $coarser_res"
    echo "High-res file: $high_innie, Coarse-res file: $coarse_innie"

    # Generate Python script for shifter
    cat << EOF > hdrescompare.py
import numpy as np
from tqdm import tqdm

# File paths and parameters
sizefile = "$sizefile"
res = int("$finer_res")
coarseres = int("$coarser_res")
infile1 = "$high_innie"
infile2 = "$coarse_innie"
outiewrite = "$shifter_outie"

# Read chromosome names and sizes
chrnames = []
chrsizes = []
with open(sizefile, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        chrnames.append(lines[0])
        chrsizes.append(int(lines[1]))

# Initialize matrices
chrmat = [[] for _ in chrnames]
smat = [[] for _ in chrnames]
emat = [[] for _ in chrnames]
scoremat = [[] for _ in chrnames]
coscoremat = [[] for _ in chrnames]

# Populate base matrices
for i in tqdm(range(len(chrnames))):
    for j in range(0, chrsizes[i] + coarseres + coarseres, res):
        chrmat[i].append(chrnames[i])
        smat[i].append(j)
        emat[i].append(j + res)
        scoremat[i].append(0)
        coscoremat[i].append(0)

# Read high-resolution scores
with open(infile1, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        try:
            myindex = chrnames.index(lines[0])
            mypos = int(float(lines[1])/res)
            scoremat[myindex][mypos] = float(lines[3])
        except (ValueError, IndexError):
            continue

# Rolling window size (number of high-resolution bins to consider)
window_size = 3 * (coarseres // res)

# Rolling average calculation
rolling_avg_scores = [[] for _ in chrnames]
for c in range(len(scoremat)):
    for i in range(len(scoremat[c])):
        window_start = max(0, i - window_size // 2)
        window_end = min(len(scoremat[c]), i + window_size // 2)
        window_scores = scoremat[c][window_start:window_end]
        rolling_avg_scores[c].append(np.mean(window_scores))

# Initialize shifted_scores with zeros
shifted_scores = [[0 for _ in range(len(scoremat[c]))] for c in range(len(scoremat))]

# Calculate shifted scores
for c in range(len(scoremat)):
    for i in range(len(scoremat[c])):
        mymean = rolling_avg_scores[c][i]
        coarscore = coscoremat[c][i] if i < len(coscoremat[c]) else 0  # Handle index if out of range
        shiftval = mymean - coarscore
        newval = scoremat[c][i] - shiftval
        shifted_scores[c][i] = newval

# Write output
with open(outiewrite, "w") as ofile:
    for w in range(len(shifted_scores)):
        for p in range(len(shifted_scores[w])):
            ofile.write(f"{chrmat[w][p]}\t{smat[w][p]}\t{emat[w][p]}\t{shifted_scores[w][p]}\n")
EOF

    # Execute the Python script
    python hdrescompare.py

    echo "After Shifter: Executed with Coarser Resolution: $coarser_res, Finer Resolution: $finer_res"
    echo "Shifter executed and output written to $shifter_outie"

}

# Shifter function



# ------------------------------------------------------------------------------
# COMPARTMENT RECLASSIFICATION
# ------------------------------------------------------------------------------
# Separate bins into A (positive GI) and B (negative GI) based on current scores
# If switch=1, re-initialize A/B states for next resolution iteration
# ------------------------------------------------------------------------------

# Separating out A and B regions from the outiefull and intersecting with the genes and Bbins to improve compartments identification
separating_ABbins() {
    local separatinginnie=$1

AbinsR=`echo "AbinsfromGI_""$myres"".txt"`
BbinsR=`echo "BbinsfromGI_""$myres"".txt"`

newBbins=`echo "newBbinsGI_""$myres"".txt"`
newAbins=`echo "newAbinsGI_""$myres"".txt"`

non_newBbins=`echo "non_newBbinsGI_""$myres"".txt"`
non_newAbins=`echo "non_newAbinsGI_""$myres"".txt"`

combined_newBbins=`echo "combined_newBbinsGI_""$myres"".txt"`
combined_newAbins=`echo "combined_newAbinsGI_""$myres"".txt"`

cat $separatinginnie | awk '{if ($4 > 0) print $0}' > $AbinsR

cat $separatinginnie | awk '{if ($4 < 0) print $0}' > $BbinsR

if [ "$switch" -eq 1 ]; then

    echo "Re-initialization of Bins is switched ""ON"". Re-evaluating both A and B bins and re-initializing for next resolution."

    cat $Bbins | intersectBed -u -a stdin -b $BbinsR > $newBbins

    awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newBbins $Bbins > $non_newBbins

    cat $newBbins $non_newBbins > $combined_newBbins

    cat $genesfile | cut -f 1-3 | intersectBed -u -a stdin -b $AbinsR > $newAbins

    awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newAbins $genesfile > $non_newAbins

    cat $newAbins $non_newAbins > $combined_newAbins

else

    echo "Generating the output without re-initialization."

fi

}


# ------------------------------------------------------------------------------
# STATISTICAL CORRECTION
# ------------------------------------------------------------------------------
# Apply Benjamini-Hochberg FDR correction for multiple testing
# Converts p-values to q-values (FDR-adjusted)
# ------------------------------------------------------------------------------

# Benjamini-Hochberg Function
BHcorrection() {
    local BHinnie=$1
    echo "Generating Benjamini-Hochberg Corrected Values"

cat << EOF > BHcorrection.py
from statsmodels.stats.multitest import multipletests
inputfile="$BHinnie"
outfile="bhcorrected"
mystarts=[]
mychr=[]
myends=[]
mypvals=[]
with open(inputfile, 'r') as innie:
        for line in innie:
                sline = line.strip()
                lines = sline.split("\t")
                mychr.append(str(lines[0]))
                mystarts.append(str(lines[1]))
                myends.append(str(lines[2]))
                mypvals.append(float(lines[3]))
hyptest, FDRlist, _, alpha_corrected =multipletests(mypvals, alpha=0.5, method='fdr_bh', is_sorted=False, returnsorted=False)
with open(outfile, 'w') as outie:
        for w in range(0,len(FDRlist)):
                outie.write(str(mychr[w]) + "\t" + str(mystarts[w]) + "\t" + str(myends[w]) + "\t" + str(FDRlist[w]) + "\n")

EOF

python BHcorrection.py

}


# ------------------------------------------------------------------------------
# P-VALUE CALCULATION
# ------------------------------------------------------------------------------
# Calculate statistical significance of GI scores using t-test
# Tests whether A and B interaction patterns are significantly different
# ------------------------------------------------------------------------------

calculate_pvalues() {
    local output_dir=$1
    local chr=$2
    local res=$3
    local outie_pvalues=$4

    awk 'NR==FNR { c[$1] = $3; m[$1] = $4; sd[$1] = $5; next } { if (($1 in m) && ($3 != 1) && (c[$1] != 1)) print $1 "\t" (m[$1]-$4) / (((sd[$1]**2)/c[$1]) + (($5**2)/$3))**.5 "\t" (c[$1] + $3 -2) }' "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" > "${output_dir}/ttable"

    zfileforpy="${output_dir}/ttable"
    pfileforpy="${output_dir}/ptable"

    # Generate p-values using Python
    cat << EOF > zlookup.py
import scipy.stats
inputfile = "$zfileforpy"
outfile = "$pfileforpy"
mybins = []
mypvals = []
with open(inputfile, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        mybins.append(str(lines[0]))
        mytstat = float(lines[1])
        mydf = int(lines[2])
        mypvals.append(scipy.stats.norm.sf(abs(mytstat)))
with open(outfile, 'w') as outie:
    for w in range(len(mybins)):
        outie.write(str(mybins[w]) + "\t" + str(mypvals[w]) + "\n")
EOF
    python zlookup.py

    mawk -v mchr="$chr" -v myres="$res" '{ print mchr "\t" $1 "\t" $1+myres "\t" $2 }' "${output_dir}/ptable" > "$outie_pvalues"
    wait
}


################################################################################
#                    MAIN RESOLUTION PROCESSING LOOP                         #
################################################################################

# Define arrays to store filenames
declare -a maxres_processed=0


# ------------------------------------------------------------------------------
# RESOLUTION PROCESSOR
# ------------------------------------------------------------------------------
# Main function to process a single resolution:
# 1. Extract Hi-C data at specified resolution
# 2. Calculate GI scores using current A/B state definitions
# 3. Generate bedGraph output files
# ------------------------------------------------------------------------------

# Function to process a given resolution and generate output
process_resolution() {
    local myres=$1  # Resolution
    local genesblock=$2  # Block for gene regions

    if [ $window == 0 ]; then
        # Calculate total bins for the chromosome size file and set the default window size
        totbins=$(awk -v var=$myres '{print $1"\t"($2/var)**2}' $sizefile | awk '{sum += $2} END {print sum}')
        newwindow=1
    else
        # Set the new window size
        newwindow=$(mawk -v win=$window -v res=$myres 'BEGIN{print int(win/(res/1000))+1}')
    fi

    # Define output filenames based on resolution
    outiefull="Crush_${myres}.bedgraph"

    echo -e "\nNow dumping reads and calculating initial CRUSH for all chromosomes at $myres"

    # Define a task function to process each chromosome
    task() {
        mychr=$(awk -v var=$myiter 'NR==var {print $1}' $sizefile)
        mysize=$(awk -v var=$myiter 'NR==var {print $2}' $sizefile)
        rowbins=$(awk -v myres=$myres -v var=$myiter 'NR==var {print int($2/myres)+1}' $sizefile)

        # Define temporary filenames
        outie="Crush_${myres}_${mychr}_tmp"
        outie2="Crush_${myres}_${mychr}_tmp2"
        tmpfiles="ABtmpfiles_${mychr}_${myres}"

        # Create temporary directory if not exists
        [ ! -d "$tmpfiles" ] && mkdir "$tmpfiles"

        mydist=$((distance + 0))
        myupper=$((upperlim + 0))

        if [ "$verbose" -gt 0 ]; then
            echo "Getting genic bins."
        fi

        # Define more temporary filenames
        sortedbed="${tmpfiles}/genebins"
        pseudoB="${tmpfiles}/pseudoB"
        dumped="${tmpfiles}/dumped_${myres}_${mychr}_tmp"

        # Initialize states and bins
        process_genes_Bbins $myres $mychr $genesblock $sortedbed $pseudoB

        if [ "$verbose" -gt 0 ]; then
            echo "Dumping ${mychr} chromosomal reads."
        fi

        # Define filenames for observed and expected reads
        newofile="${tmpfiles}/observeds"
        newefile="${tmpfiles}/expecteds"
        dumperrors="dumperrors_${myres}"

        # Use straw or cooler for dumping reads based on file format
        if [ $juiceorcool -eq 0 ]; then
            # Dump reads using straw
            run_dumper 'observed' $norm $hicpath $mychr $myres $dumped $dumperrors
        else

            echo "Dumping using cooler..."

            # Dump reads using cooler
            cooler dump $hicpath::/resolutions/$myres -r $mychr --join | cut -f 2,5,7 > $dumped 2> $dumperrors
        fi

        # Sanity check: did dump succeed?
        if [ ! -s "$dumped" ]; then
            echo "WARNING: no dumped reads for $mychr at $myres (file: $dumped). Check $dumperrors" >&2
            return
        fi

        wait

        # Filter reads based on upper limit if specified
        if [ "$upperlim" -gt 0 ]; then
            awk -v var=$myres -v distfilt=$mydist -v upper=$upperlim '{if (($2-$1 >= distfilt) && ($2-$1 <= upper)) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' $dumped | \
            awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile }' &
        else
            awk -v var=$myres -v distfilt=$mydist '{if ($2-$1 >= distfilt) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' $dumped | \
            awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile }' &
        fi

        wait

        if [ "$verbose" -gt 0 ]; then
            echo "Done dumping."
            echo "Now getting distance normalized values for ${mychr} chromosome."
            start=$(date +%s)
        fi

        # Define filenames for error log, temporary mean, and symm file
        errorfile="errorlog"
        tmpmean="${tmpfiles}/tmpmean"
        symmfile="${tmpfiles}/symmonfile"

        # Calculate rows-Mean using expected and observed files
        mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 != $3) print $2"\t"$3"\t"($4+1)/(c1[$1]+1)}' $newefile $newofile | mawk -v thresh=$threshold '{if ($3 < thresh) print $1"\t"$2"\t"$3"\n"$2"\t"$1"\t"$3; else print $1"\t"$2"\t"thresh"\n"$2"\t"$1"\t"thresh}' | awk -v symmfile=$symmfile -v tmpmean=$tmpmean -v rowbins=$rowbins '{m[$1]+=$3;d[$1]+=($3**2); print $0 >> symmfile} END {for (i in m) print i"\t"m[i]/rowbins"\t"((d[i]/rowbins-((m[i]/rowbins)**2))**.5) >> tmpmean }'  2> $tmpfiles/$errorfile

        wait

        if [ "$verbose" -gt 0 ]; then
            end=$(date +%s)
            seconds=$(echo "$end - $start" | bc)
            start=$(date +%s)
        fi

        # Check if there are enough bins overlapping with bed file
        if [ "$myres" -gt 50000 ]; then
            oppocheck=$(mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) ; else if ($2 in c1); else print $0}' $sortedbed $symmfile | wc -l)
        else
            oppocheck=2
        fi

        # Perform real calculation
        if [ "$oppocheck" -lt 1 ] && [ "$verbose" -gt 0 ]; then
            echo "WARNING: All bins on chromosome ${mychr} overlap with your bed file at ${myres} resolution. This should not affect the scores if using multiple resolutions. If using a single resolution, then maybe not use this resolution?"
        fi

        if [ "$verbose" -gt 0 ]; then
            end=$(date +%s)
            seconds=$(echo "$end - $start" | bc)
            start=$(date +%s)
        fi

        # Calculate Z-scores
        RowZs=`echo "RowZs"`
        mawk 'NR==FNR { m[$1] = $2; d[$1] = $3; next} {if (d[$2] > 0) print $1 "\t" $2 "\t" ($3 - m[$2]) / (d[$2])}' $tmpmean $symmfile > $tmpfiles/$RowZs

        # Create a zero's file
        zfile="zerofile"
        mawk -v var=$mychr '{if ($1 == var) print $0}' $sizefile | mawk -v wlkr=$myres '{for (i=0;i<=$2;i+=wlkr) print i"\t"i"\t0.1"}' > $tmpfiles/$zfile

        # Run the process_oppocheck_statement function
        process_oppocheck_statement $sortedbed $pseudoB $outie $mychr $myres "ABtmpfiles_${mychr}_${myres}" "ABtmpfiles_${mychr}_${myres}"

    }

    # Use the globally defined semaphore functions
    open_sem $cpu  # Initialize semaphore with N tokens (N = CPU count)

    if [ "$whichchoice" == "o" ]; then
        for (( myiter=1; myiter<=$totchroms; myiter++ )); do
            run_with_lock task $myiter  # Execute with semaphore lock (parallel processing)
        done
    fi

    wait

    # Merging CRUSH files
    echo "Resolution: $myres for Resolution output"
    echo -ne "Merging individual chromosome files\033[0K\r"

    cat Crush_${myres}_*_tmp | grep -v -i nan | \
    awk '
    {
        key = $1"\t"$2"\t"$3;  # Create a unique key for each bin based on chr, start, and end
        sum[key] += $4;        # Sum the scores for each unique bin
        count[key]++;          # Count the occurrences for each bin
    }
    END {
        for (key in sum) {
            if (count[key] > 1) {
                avg = sum[key] / count[key];  # Calculate average for duplicated bins
                print key "\t" avg;
            } else {
                print key "\t" sum[key];      # Print unique bin as is
            }
        }
    }' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull

    wait

    rm Crush_${myres}_*_tmp

    #Setting up for shifter

    if [ "$myres" == "$maxres" ]; then
        coarse_infile=`echo "GI_""$myres"".bedgraph"`
        cat $outiefull > $coarse_infile
        wait
    fi

    res_values=`echo "res_values"`
    echo "$myres" >> $res_values

    # Re-evaluating both A and B bins for re-initialization for next resolution

    if [ "$countres" -gt 0 ]; then
        # Creating a python output file
        python_output=`echo "shifter.bedgraph"`

        if [ "$coarse_res" -eq "$maxres" ] && [ "$coarse_res" != "$high_res" ]; then

            # Shifter Input files for the specified resolutions
            high_res_infile=`echo "$outiefull"`
            coarse_res_infile=`echo "$coarse_infile"`

        else

            # Shifter Input files for the specified resolutions
            high_res_infile=`echo "$outiefull"`
            coarse_res_infile=`echo "$python_output"`

        fi

        # Conditional checking for multiples of resolutions for shifter
        # If coarseres is divisible by highres, use these resolutions for shifter as coarse and high
        if [[ "$(( $coarse_res % $high_res ))" -eq 0 ]]; then

            crtmp=$coarse_res

        else

            if [[ "$(( $crtmp % $high_res ))"  -eq 0 ]]; then
                coarse_res=$crtmp
            fi

        fi

        # Executing the Shifter Python Script
        run_midshifter $high_res $coarse_res $high_res_infile $coarse_res_infile $python_output

        #Adding this to check python output
        echo "Resolution: "$myres" for shifter output"

        #Separating A/B bins using Shifter python output
        separating_ABbins $python_output

        # Update coarse and high resolutions
        coarse_res=$high_res

    else

        #Separating A/B bins using Resolution based Output
        separating_ABbins $outiefull

    fi

    wait

    rm Crush*_tmp 2> smallerrors
    rm ACrush*_tmp 2> smallerrors
    rm BCrush*_tmp 2> smallerrors

    countres=$((countres+1))
}

reprocess_resolutions_with_shifter() {
    local current_res=$1
    local end_index=$2
    local reprocess=$3

    countres_reprocess=0
    coarse_res_reprocess=$maxres
    crtmp_reprocess=$maxres


    for (( j=0; j<=$end_index; j++ )); do

        local prev_res=${res_array[$j]}

        echo "Reprocessing for resolution: $prev_res"

        highres_reprocess=$prev_res

        local outiefull_reprocess=`echo "Crush_Re_""$prev_res"".bedgraph"`

        # Loop through each chromosome
        local chromosomes=($(cut -f1 $sizefile))

        for mychr in "${chromosomes[@]}"; do
            local sortedbed_reprocess="genebins_reprocess_${mychr}_${prev_res}"
            local pseudoB_reprocess="pseudoB_reprocess_${mychr}_${prev_res}"
            local outie_reprocess="Crush_reprocess_${prev_res}_${mychr}_tmp"

            reprocesstmpfiles=`echo "ABreprocesstmpfiles_""$mychr""_""$prev_res"`

            #Creating tmp directory
            [ ! -d "$reprocesstmpfiles" ] && mkdir "$reprocesstmpfiles"

            echo "Using combined bins for $prev_res and the current value of resolution in this loop is $current_res"

            # Processing the genes & Bbins files for other resolutions
            cat $combined_newAbins | mawk -v myres=$prev_res -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$prev_res '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed_reprocess

            # Bbins
            cat $combined_newBbins |  mawk -v myres=$prev_res -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$prev_res '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB_reprocess

            echo "Files: $sortedbed_reprocess, $pseudoB_reprocess"

            # Assuming the existence of process_oppocheck_statement function with correct parameters
            process_oppocheck_statement $sortedbed_reprocess $pseudoB_reprocess $outie_reprocess $mychr $prev_res "ABtmpfiles_${mychr}_${prev_res}" "ABreprocesstmpfiles_${mychr}_${prev_res}"

            wait

            if [ $prev_res -eq $minres ] && [ $pcalculation -gt 0 ]; then
                pscores_reprocess="ttest_reprocess_${prev_res}_${mychr}_tmp"
                calculate_pvalues "ABreprocesstmpfiles_${mychr}_${prev_res}" $mychr $prev_res $pscores_reprocess
            fi
        done

        # Merge the reprocessed CRUSH files
        echo -ne "Merging individual chromosome files\033[0K\r"

        cat Crush_reprocess_${prev_res}_*_tmp | grep -v -i nan | \
        awk '
        {
            key = $1"\t"$2"\t"$3;  # Create a unique key for each bin based on chr, start, and end
            sum[key] += $4;        # Sum the scores for each unique bin
            count[key]++;          # Count the occurrences for each bin
        }
        END {
            for (key in sum) {
                if (count[key] > 1) {
                    avg = sum[key] / count[key];  # Calculate average for duplicated bins
                    print key "\t" avg;
                } else {
                    print key "\t" sum[key];      # Print unique bin as is
                }
            }
        }' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull_reprocess

        wait

        rm Crush_reprocess_${prev_res}_*_tmp

        if [ $prev_res -eq $minres ] && [ $pcalculation -gt 0 ]; then

            outiefullPval_reprocess="pvalues_Re_${prev_res}.bedgraph"
            echo "Resolution: $prev_res for Resolution output and minres is : $minres"
            cat ttest_reprocess_${prev_res}_*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefullPval_reprocess
            wait
        fi

        if [ "$prev_res" == "$maxres" ]; then

            coarse_infile_reprocess=`echo "GI_Re_""$prev_res"".bedgraph"`
            cat $outiefull_reprocess > $coarse_infile_reprocess

            wait
        fi

        prev_values=`echo "prev_values"`
        echo "$prev_res" >> $prev_values

        # Re-evaluating both A and B bins for re-initialization for next resolution

        echo "$j"

        if [ "$j" -gt 0 ]; then

            if [ "$coarse_res_reprocess" -eq "$maxres" ] && [ "$coarse_res_reprocess" != "$highres_reprocess" ]; then
                # Shifter Input files for the specified resolutions
                highres_infile_reprocess=`echo "$outiefull_reprocess"`
                coarse_res_infile_reprocess=`echo "$coarse_infile_reprocess"`
            else
                highres_infile_reprocess=`echo "$outiefull_reprocess"`
                coarse_res_infile_reprocess=`echo "$python_output"`

            fi

            # Conditional checking for multiples of resolutions for shifter
            # If coarseres is divisible by highres, use these resolutions for shifter as coarse and high
            if [[ "$(( $coarse_res_reprocess % $highres_reprocess ))" -eq 0 ]]; then
                crtmp_reprocess=$coarse_res_reprocess
            else
                if [[ "$(( $crtmp_reprocess % $highres_reprocess ))" -eq 0 ]]; then
                coarse_res_reprocess=$crtmp_reprocess
                fi
            fi

            if [ -f "$highres_infile_reprocess" ] && [ -f "$coarse_res_infile_reprocess" ]; then
                echo "Files exist: $highres_infile_reprocess, $coarse_res_infile_reprocess"
            else
                echo "One or both files do not exist: $highres_infile_reprocess, $coarse_res_infile_reprocess"
            fi

            echo "Running shifter for $prev_res"

            # Executing the Shifter Python Script
            run_midshifter $highres_reprocess $coarse_res_reprocess $highres_infile_reprocess $coarse_res_infile_reprocess $python_output

            echo "Shifter run completed for $prev_res"

            #Adding this to check python output
            echo "Resolution: "$prev_res" for shifter output"

            #Separating A/B bins using Shifter python output
            separating_ABbins $python_output

            # Update coarse and high resolutions
            coarse_res_reprocess=$highres_reprocess
        else
            #Separating A/B bins using Resolution based Output
            separating_ABbins $outiefull_reprocess
        fi
    done
}

#The show unfolds from here:

# Creating a directory for temporary file
echo "Creating and moving into temporary directory"
mkdir $crushdir
cd $crushdir

# Reading Resolutions
if [ $juiceorcool -eq 0 ]
then
cat << EOF > listres.py
import hicstraw
hic = hicstraw.HiCFile("$hicpath")
totres=hic.getResolutions()
keepres=[]
for i in totres:
	if (int(i) >= int("$res")) and (int(i) <= int("$coarsestres")):
		keepres.append(str(i))
print(",".join(keepres))

EOF

res=`python listres.py`

else
reslist=`cooler ls $hicpath | sed 's/\//\t/g' | awk -v var=$res -v cres=$coarsestres '{if (($NF >= var) && ($NF <= cres)) print $NF}' | sort -k 1bnr,1b --stable | awk '{if (NR == 1) printf "%s", $1; else printf ",%s", $1}' | awk '{print $0}'`
res=$reslist
fi


# Check for fasta
isfasta=`head -1 $fastafile | awk '{if ($1 ~ /^>/) print 1; else print 0}'`

# Display the file paths
echo "initial-A: $genesfile"
echo "sizes-file: $sizefile"

# Check sizefile
numcolumns=`mawk '{if (NF != 2) print $0}' $sizefile | wc -l`
if [ "$numcolumns" -gt 0 ]
then
echo "size file is incorrect format. It should be two columns with chromosome name and size"
exit 0
fi

totchroms=`cat $sizefile | wc -l`
maxres=`echo "$res" | sed "s/,/ /g" | mawk '{max=$1;for(i=2;i<=NF;i++){if($i > max) max = $i} print max}'`
minres=`echo "$res" | sed "s/,/ /g" | mawk '{min=$1;for(i=2;i<=NF;i++){if($i < min) min = $i} print min}'`

if [ $isfasta == "1" ]; then
    echo "Fasta-file: $fastafile"
else
    echo "initial-B: $fastafile"
fi

resbins=`echo "size_bins"".bed"`
gcfile=`echo "gc_bins"".txt"`
gc_g_ga=`echo "gc_g_ga_bins"".bed"`
Bbins=`echo "Bbins"".bed"`

if [ $minres -ge 500 ]; then
    cat $sizefile | awk -v myres=500 '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
    echo "$resbins calculated at 500bp resolution"
else
    cat $sizefile | awk -v myres=$minres '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
    echo "$resbins calculated at ""$minres"" resolution"
fi

wait

if [ $isfasta == "1" ]; then
    echo "Now, calculating the gc content in parallel by chromosome:"

    # Split resbins by chromosome
    mkdir -p gc_tmp
    awk '{print >> "gc_tmp/"$1".bed"; close("gc_tmp/"$1".bed")}' $resbins

    # Process each chromosome in parallel
    open_sem $cpu  # Initialize semaphore with N tokens (N = CPU count)

    for chrom in $(cut -f 1 $sizefile); do
        run_with_lock bash -c "  # Execute with semaphore lock (parallel processing)
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

    wait  # Wait for all bedtools jobs to complete

    # Check if we got results for all chromosomes
    expected_files=$(cut -f 1 $sizefile | wc -l)
    actual_files=$(ls gc_tmp/*.gc 2>/dev/null | wc -l)

    if [ "$actual_files" -lt "$expected_files" ]; then
        echo "Warning: Only got GC content for $actual_files chromosomes out of $expected_files expected"
    fi

    # Merge results
    cat gc_tmp/*.gc > $gcfile

    # Clean up temporary files
    rm -rf gc_tmp

    wait

    # Created a %G/%G+%A file in 7th column
    cat $gcfile | grep -v user | cut -f 1-5 | awk '{print $0"\t"($4+$5)}' | awk '{if ($6 > 0) print $0}' | awk '{print $0"\t"($5/$6)}' > $gc_g_ga
    wait

    # Calculating mean and SD for the whole genome and then subtracting 1SD from the mean
    gc_thresh=`cat $gc_g_ga | awk '{s+=$7; ss+=$7*$7; linecount+=1} END{print m=s/linecount, sqrt(ss/linecount-m^2)}' | awk '{print $0}' | awk '{print $1 - 1*$2}'`

    if [ $verbose -gt 0 ]; then
        echo "gc threshold is ""$gc_thresh"
    fi

    # Generating the bins by considering the bins below
    cat $gc_g_ga | awk -v thresh=$gc_thresh '{if ($7 < thresh) print $0}' > $Bbins
    wait
    echo "Finished generating Bbins file for all chromosomes "

else
    cat $fastafile > $Bbins
fi

if [ $endZ == 0 ]; then
    endZ=`echo "$minres"`
fi

# Eigen Vector Block After Reslist
EVAstates="EVAstates.bed"
EVBstates="EVBstates.bed"

process_eigenvector_file() {
    local eigenfile="$1"
    local output_positive="$2"
    local output_negative="$3"

    # Process positive eigenvalues
    awk '{if ($4 > 0) print $0}' "$eigenfile" > "$output_positive"
    if [ $? -ne 0 ]; then
        echo "Error creating $output_positive from $eigenfile. Please check the awk command."
        exit 1
    fi

    # Process negative eigenvalues
    awk '{if ($4 < 0) print $0}' "$eigenfile" > "$output_negative"
    if [ $? -ne 0 ]; then
        echo "Error creating $output_negative from $eigenfile. Please check the awk command."
        exit 1
    fi

    # Verify that both files are created and not empty
    if [ ! -s "$output_positive" ]; then
        echo "Error: $output_positive is empty or not created."
        exit 1
    fi
    if [ ! -s "$output_negative" ]; then
        echo "Error: $output_negative is empty or not created."
        exit 1
    fi

    echo "Eigen Vector files created successfully: $output_positive, $output_negative"
}

# Handle user-provided eigenfile
if [ -n "$eigenfile" ]; then
    echo "Using user-provided Eigen Vector file: $eigenfile"
    process_eigenvector_file "$eigenfile" "$EVAstates" "$EVBstates"
else
    echo "No user-provided Eigen Vector file. Preprocessing binned files for all chromosomes at resolution $eigenres."

    gcbinned="gcbinned_all.bed"

    # Generate binned GC content file for all chromosomes
    awk -v res=$eigenres -v chr_size_file=$sizefile '
    BEGIN { OFS="\t"; while ((getline < chr_size_file) > 0) sizes[$1] = $2 }
    $1 in sizes && $2 >= 0 && $3 >= 0 {
        start = int($2 / res) * res
        end = int($3 / res) * res
        if (end > sizes[$1]) end = sizes[$1]
        if (start < sizes[$1]) print $1, start, end
    }' "$Bbins" > "$gcbinned"

    # Verify the files are created and not empty
    if [[ ! -s $gcbinned ]]; then
        echo "Error: One or both of the preprocessed binned files are empty. Exiting..."
        exit 1
    fi

    echo "Binned files generated for all chromosomes: $genesbinned, $gcbinned."

    # Default behavior when no eigenfile is provided
    # Parallelize eigenvector generation using semaphore (same pattern as chromosome processing)
    open_sem $cpu  # Initialize semaphore with N tokens (N = CPU count)

    for chrom in $(cut -f 1 $sizefile); do
        run_with_lock bash -c "  # Execute with semaphore lock (parallel processing)
            chrom='$chrom'
            sizefile='$sizefile'
            gcbinned='$gcbinned'
            hicpath='$hicpath'
            eigenres='$eigenres'
            genesfile='$genesfile'

            chrom_size=\$(awk -v chr=\"\$chrom\" '\$1 == chr {print \$2}' \"\$sizefile\")
            if [ -z \"\$chrom_size\" ]; then
                echo \"Error: Chromosome size for \$chrom not found in \$sizefile. Skipping...\"
                exit 0
            fi

            # Filter out binned data for the current chromosome
            gcbinned_chr=\"gcbinned_\${chrom}.bed\"

            awk -v chr=\$chrom '\$1 == chr' \"\$gcbinned\" > \"\$gcbinned_chr\"

            if [[ ! -s \$gcbinned_chr ]]; then
                echo \"Error: Binned files for chromosome \$chrom are empty. Skipping...\"
                exit 0
            fi

            generateEV \"\$hicpath\" \"\$chrom\" \"\$eigenres\" \"\$genesfile\" \"\$gcbinned_chr\" \"\$chrom_size\"
        "
    done

    wait  # Wait for all eigenvector generation jobs to complete

    # Verify we got eigenvector files for all chromosomes
    echo "Checking eigenvector generation completion..."
    missing_evs=0
    for chrom in $(cut -f 1 $sizefile); do
        if [ ! -f "best_Eigen_${chrom}.bedgraph" ]; then
            echo "Warning: Missing eigenvector file for chromosome ${chrom}"
            missing_evs=$((missing_evs + 1))
        fi
    done

    if [ $missing_evs -gt 0 ]; then
        echo "Warning: $missing_evs chromosomes are missing eigenvector files"
        echo "Continuing anyway, but results may be incomplete..."
    fi

    # Concatenate best eigenvectors for full genome BEDGraph
    concatenate_best_pcs

    process_eigenvector_file "EV_full_genome.bedgraph" "$EVAstates" "$EVBstates"
fi



# Initialize variables
countres=0
coarse_res=$maxres
crtmp=$maxres
IFS=',' read -r -a res_array <<< "$res"
maxres_processed=0

# Loop through each resolution
for (( i=0; i<${#res_array[@]}; i++ )); do
    myres=${res_array[$i]}
    echo "Current resolution: $myres"
    echo "Array content: ${res_array[@]}"
    high_res=$myres

    # Determine if genesblock should be used
    genesblock=0
    if [ "$myres" == "$maxres" ] && [ "$maxres_processed" -eq 0 ]; then
        genesblock=1
    fi

    echo "Processing with myres = $myres, genesblock = $genesblock"

    # Call the function for the current resolution
    process_resolution $myres $genesblock

    if [ $i -gt 0 ]; then
        reprocess_resolutions_with_shifter $myres $i 1

        #Final Processing of bedgraphs

        finaloutie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`
        Crush_todelete=`echo "Crush_todelete_""$myres""_reprocess"`

        finalpoutie=`echo "$outpre""mergedqvalue_""$myres"".bedgraph"`
        pval_todelete=`echo "pval_todelete_""$myres""_reprocess"`

        cat Crush_Re*.bedgraph | mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4/int($3-$2)}' | mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | mawk '{c1[$1] += $2; c2[$1]++} END {for (i in c1) print i"\t"c1[i]/c2[i]}' | sed 's/:/\t/g' > $finaloutie 2> smallerrors

        if [ $myres -eq $minres ]; then
            cat pvalues_Re*.bedgraph | mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | awk '{c1[$1] += log($2)/log(10); c2[$1]++} END {for (i in c1) print i"\t"10**(c1[i]/c2[i])}' | sed 's/:/\t/g' | sort -k 1,1 -V -k 2bn,2b --stable > $finalpoutie 2> smallerrors
        fi

        echo "Running BH correction"
        BHcorrection $finaloutie
        wait
        mv bhcorrected $finalpoutie
        wait

        GIave=`cat $finaloutie | mawk '{if ($4 < 0) sum+=($4*-1); else sum+=($4)} END {print sum/NR}' `

        if [ $(bc <<< "$qthresh > 0") -eq 1 ]; then

            finaloutiefilt=`echo "$outpre""mergedCrush_""$myres""_qfiltered_reprocess.bedgraph"`

            mawk -v fdr=$qthresh -v var=$GIave 'NR==FNR {a[$1":"$2":"$3] = $4; next} {if (a[$1":"$2":"$3] <= fdr) print $1"\t"$2"\t"$3"\t"$4/(var/100)}'  $finalpoutie $finaloutie> $finaloutiefilt

        fi

        if [ $trackline -eq 0 ]; then

            cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' > $Crush_todelete
            cat $finalpoutie | mawk '{print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
        else
            cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$0; else print $0}' > $Crush_todelete
            cat $finalpoutie | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$1"\t"$2"\t"$3"\t"$4; else print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
            wait
        fi

        if [ $trackline -gt 0 ] && [ $(bc <<< "$qthresh > 0") -eq 1 ]; then
            filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`
            cat $finaloutiefilt | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20,20 autoScale off\n"$0; else print $0}' > $filt_todelete

            wait

            mv $filt_todelete $finaloutiefilt
            mv $Crush_todelete $finaloutie
            mv $pval_todelete $finalpoutie

            wait
        fi

        ## Refixing the extra bins
        # Creating a size bed file
        sizeBed=`echo "Chrsizes.bed"`
        finaloutie2=`echo "mergedCrush2_""$myres"".bedgraph"`

        cat $sizefile | awk '{print $1"\t""1""\t"$2}' > $sizeBed

        cat $finaloutie | grep -v track | intersectBed -wa -a stdin -wb -b $sizeBed | awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | sort -k 1,1 -V -k 2bn,2b --stable | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' > $finaloutie2
        mv $finaloutie2 $finaloutie

        # Apply smoothing if required
        if [ "$smoothing" -gt 0 ]; then
            echo "Smoothing Function is turned on"
            smoother="smoother.bedgraph"
            run_smoothing $newoutie $smoother
            rm $smoother tmp2_shiftedleft tmp1_shiftedright  # Clean up smoothing-related temporary files
        fi

        if [ $(bc <<< "$qthresh > 0") -eq 1 ] && [ $pcalculation -eq 1 ] && [ $myres -eq $minres ]; then
            filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`

            # Perform intersection with sizeBed and format for visualization
            cat $finaloutiefilt | grep -v track | intersectBed -wa -a stdin -wb -b $sizeBed | \
                awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | \
                mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' > $filt_todelete

            mv $filt_todelete $finaloutiefilt   # Move the temporary filtered file to the final filtered output
            mv $finalpoutie ../                 # Move the p-value output to the parent directory
        fi

        # Move final outputs to the parent directory
        mv $finaloutie ../

        if [ $(bc <<< "$qthresh > 0") -eq 1 ] && [ $pcalculation -eq 1 ] && [ $myres -eq $minres ]; then
            mv $finaloutiefilt ../
        fi
    fi
done

if [ $doNotMerge -eq 1 ]; then
    echo "Finished with each resolution listed, but not merging...."
    exit 0
fi

countres=0
for totres in $(echo $res | sed "s/,/\n/g" | sort -k 1bn,1b --stable | sed "s/\n/ /g" ); do
    echo "$totres" >> GIrestmpfile
    countres=$((countres+1))
done

if [ $countres -le 1 ]; then
    echo "Only one resolution listed, so nothing to merge. Check Crush_""$totres"".bedgraph for the single resolution GI score"

    rm AbinsfromGI_*.txt 2> smallerrors
    rm BbinsfromGI_*.txt 2> smallerrors
    rm newAbinsGI_*.txt 2> smallerrors
    rm newBbinsGI_*.txt 2> smallerrors
    rm non_newAbinsGI_*.txt 2> smallerrors
    rm non_newBbinsGI_*.txt 2> smallerrors
    rm combined_newAbinsGI_*.txt 2> smallerrors
    rm combined_newBbinsGI_*.txt 2> smallerrors
    rm smallerrors
    exit 0
fi

rm GIrestmpfile

wait

cd ../

if [ "$cleanup" -gt 0 ]; then
    echo "cleaning up"
    rm -r $crushdir
fi

echo "Finished! Check the output."
echo -e "Note: Please keep in mind that we are using resolution walking.\nIf the coarsest resolution doesn't match the compartment pattern at that resolution, please consider re-running with the -m parameter set to start the walking with smaller bins."

