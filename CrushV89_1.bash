#!/usr/bin/env bash
################################################################################
#                                                                              #
#                    CRUSH - Chromatin Compartment Analysis                   #
#         Compartmental Refinement for Ultraprecise Stratification in Hi-C    #
#                                                                              #
################################################################################
#
# CRUSH analyzes Hi-C data to identify A/B chromatin compartments at multiple resolutions.
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
# - CRUSH Score: Genome Interaction score measuring preference for A vs B interactions
#
################################################################################




################################################################################
#                         VERSION & DEFAULTS                                 #
################################################################################

#CRUSH Script Version
version=1.0


################################################################################
#                    PARAMETER INITIALIZATION                                #
################################################################################

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
window=5
verbose=0
switch=1
qthresh=0.05
whichchoice=`echo "o"`
cleanup=1
endZ=0
exclbed=0
crushdir=`echo "CRUSHtmp_""$RANDOM"`
outpre=""
norm=NONE
coarsestres=2500000
pcalculation=1
smoothing=0
eigenres=100000
genomebuild=0


################################################################################
#                       SYSTEM UTILITIES                                     #
################################################################################

function shutdown() {
    tput cnorm
}
trap shutdown EXIT

function cursorBack() {
    echo -en "\033[$1D"
}


################################################################################
#                    PARALLEL PROCESSING FRAMEWORK                           #
################################################################################
#
# Token-based semaphore system for controlling parallel job execution.
#
################################################################################

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
        echo "WARNING: a worker exited with code $x — continuing." >&2
        printf '%s' 000 >&3
        return 1
    }
    (
        "$@"
        local exit_code=$?
        if [[ $exit_code -ne 0 ]]; then
            echo "WARNING: worker exited with code $exit_code" >&2
        fi
        printf '%.3d' 0 >&3
    ) &
}


################################################################################
#                        HELP & USAGE FUNCTIONS                              #
################################################################################

function usage {
    echo -e "\n\nusage : crush -i HIC -r RESOLUTION { -gb GENOMEBUILD | -g SIZEFILE -a ABED -b BBED|FASTA } [-c CPU] [-h] [--advanced-help]"
    echo -e "Use option -h|--help for basic help, --advanced-help for all options"
}

function help {
    usage;
    echo ""
    echo "================================================================================"
    echo "                         CRUSH - Hi-C Compartment Analysis"
    echo "         Compartmental Refinement for Ultraprecise Stratification in Hi-C"
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
    echo "  --advanced-help         Display all options including advanced parameters"
    echo ""
    echo "================================================================================"
    echo ""

    echo "ALWAYS REQUIRED:"
    echo "----------------"
    echo "  -i, --hic FILE          Input Hi-C file (.hic or .mcool format)"
    echo "  -r, --res VALUE         Resolution in bp  (e.g. 10000 for 10kb)"
    echo "  -c, --cpu NUMBER        CPU threads  (default: 1)"
    echo ""
    echo "================================================================================"
    echo ""
    echo "REFERENCE FILES — choose ONE path:"
    echo ""
    echo "  PATH A  -gb  (supported builds: hg19, hg38, mm10, mm9 | res >= 500bp)"
    echo "  -----------------------------------------------------------------------"
    echo "  -gb, --genomebuild BUILD"
    echo "      Auto-downloads chrom.sizes, genes.bed, and Bbins.bed."
    echo "      -g, -a, and -b are NOT needed."
    echo "      Example:  crush -i data.hic -gb hg38 -r 10000 -c 8"
    echo ""
    echo "  PATH B  -g -a -b  (any genome | any resolution)"
    echo "  -----------------------------------------------------------------------"
    echo "  -g, --genomesize FILE   Chromosome sizes file (chr_name  size)"
    echo "  -a, --initialA   FILE   Gene annotations BED  (initialises A compartments)"
    echo "  -b, --initialB   FILE   Bbins BED or genome FASTA"
    echo "                          Use FASTA when res < 500bp (GC recalculated)"
    echo "      Example:  crush -i data.hic -g hg19.sizes -a genes.bed -b genome.fa -r 10000 -c 8"
    echo ""
    echo "================================================================================"
    echo ""
    echo "  Run --advanced-help for output control, analysis tuning, and all other options."
    echo ""
    echo "================================================================================"
}

function advanced_help {
    usage;
    echo ""
    echo "================================================================================"
    echo "                         CRUSH - Hi-C Compartment Analysis"
    echo "         Compartmental Refinement for Ultraprecise Stratification in Hi-C"
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
    echo "  -h, --help              Display basic help menu"
    echo "  --advanced-help         Display this full help menu"
    echo ""
    echo "================================================================================"
    echo ""

    echo "REQUIRED PARAMETERS:"
    echo "--------------------"
    echo ""
    echo "  -i, --hic FILE          Input Hi-C file (.hic, .cooler, or .mcool format)"
    echo "                          Example: '/path/to/file.hic'"
    echo ""
    echo "  -r, --res VALUE         Desired resolution for analysis"
    echo "                          Example: 10000 (for 10kb resolution)"
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
    echo "================================================================================"
    echo ""

    echo "OUTPUT OPTIONS:"
    echo "---------------"
    echo ""
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
    echo "================================================================================"
    echo ""

    echo "PERFORMANCE OPTIONS:"
    echo "--------------------"
    echo ""
    echo "  -c, --cpu NUMBER        Number of CPU threads to use"
    echo "                          Default: 1"
    echo ""
    echo "  -C, --cleanup [0|1]     Clean up temporary files after completion"
    echo "                          Default: 1 (cleanup enabled)"
    echo ""
    echo "================================================================================"
    echo ""

    echo "ANALYSIS OPTIONS:"
    echo "-----------------"
    echo ""
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
    echo "  -w, --window VALUE      Sliding window size for score averaging (in kb)"
    echo "                          Default: 5"
    echo "                          Set to 1 to disable smoothing entirely"
    echo "                          Set to 0 for legacy auto-calculation from depth"
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
    echo "  -e, --eigenfile FILE    Supply your own eigenfile to initialize A and B states"
    echo "                          If not provided, CRUSH calculates eigenvectors"
    echo "                          automatically and selects the best principal component"
    echo ""
    echo "  -s, --smoothing [0|1]   Apply smoothing to compartment boundaries"
    echo "                          0 = No smoothing (default)"
    echo "                          1 = Enable smoothing"
    echo ""
    echo "================================================================================"
    echo ""

    echo "STATISTICAL OPTIONS:"
    echo "--------------------"
    echo ""
    echo "  -p, --pcalculation [0|1] Calculate p-values"
    echo "                           0 = Skip p-value calculation"
    echo "                           1 = Calculate (default)"
    echo ""
    echo "  -q, --qvalue VALUE       Q-value threshold for significance filtering"
    echo "                           Default: 0.05"
    echo "                           Set to 0 to disable filtering"
    echo ""
    echo "================================================================================"
    echo ""

    echo "ADVANCED OPTIONS:"
    echo "-----------------"
    echo ""
    echo "  -v, --verbose [0|1]     Enable verbose output messages"
    echo "                          Default: 0 (minimal output)"
    echo ""
    echo "  -x, --exclbed FILE      BED file of regions to exclude from analysis"
    echo ""
    echo "  -E, --endZ VALUE        End Z-score value"
    echo ""
    echo "  -use [o|u]              Use existing GI tracks from previous calculations"
    echo "                          'o' = Overwrite/recalculate (default)"
    echo "                          'u' = Use existing (useful for merging resolutions)"
    echo ""
    echo "  -f, --tmpfolder NAME    Custom name for temporary folder"
    echo "                          Default: CRUSHtmp_RANDOM"
    echo "                          (Must not already exist)"
    echo ""
    echo "  -gb, --genomebuild BUILD  Genome build: hg19, hg38, mm10, or mm9"
    echo "                            Auto-fetches chrom.sizes, genes.bed, and Bbins.bed"
    echo "                            into the crush working folder. -g, -a, -b optional."
    echo "                            Only available for resolution >= 500bp — Bbins.bed"
    echo "                            is pre-computed at 500bp resolution."
    echo "                            For res < 500bp (any genome): supply -g, -a, and"
    echo "                            -b FASTA manually so GC content is recalculated."
    echo "                            For builds outside hg19/hg38/mm10/mm9: same."
    echo ""
    echo "================================================================================"
    echo ""
    echo "EXAMPLES:"
    echo ""
    echo "  With --genomebuild (automatic reference files):"
    echo "    crush -i data.hic -gb hg38 -r 10000"
    echo "    crush -i data.hic -gb hg38 -r 10000 -c 8 -o myproject_"
    echo ""
    echo "  With manual reference files:"
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

validate_url() {
    local url=$1
    if ! command -v curl &> /dev/null; then
        echo "Error: curl is not installed. Please install curl to validate remote files."
        exit 1
    fi
    curl --output /dev/null --silent --head --fail "$url"
    return $?
}


################################################################################
#                     ARGUMENT PARSING                                       #
################################################################################

while test $# -gt 0; do
    case "$1" in
        -h|--help)
            help
            exit 0
            ;;
        --advanced-help)
            advanced_help
            exit 0
            ;;
        -i|--hic)
            shift
            if [ $# -eq 0 ]; then
                echo "Error: No file specified after -i|--hic"
                exit 1
            fi
            hicpath="$1"
            if [[ "$hicpath" =~ ^https?:// ]]; then
                if ! validate_url "$hicpath"; then
                    echo "Error: Unable to access URL: $hicpath"
                    exit 1
                fi
                echo "Using remote file: $hicpath"
            else
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
        -gb|--genomebuild)
            genomebuild=$2
            shift
            ;;
    esac
    shift
done

CRUSH_SUPPORTED_BUILDS=("hg19" "hg38" "mm10" "mm9")
genomebuild_supported=0

if [ "$genomebuild" != "0" ]; then

    # Universal rule: below 500bp requires the user to supply FASTA + sizes + genes.
    # Pre-computed Bbins.bed are at 500bp and are too coarse for finer resolutions.
    if [ "$res" != "0" ] && [ "$res" -lt 500 ] 2>/dev/null; then
        echo ""
        echo "Note: --genomebuild auto-fetch is unavailable for resolutions finer than 500bp."
        echo "      Pre-computed Bbins.bed files are at 500bp resolution, which is too coarse"
        echo "      for ${res}bp analysis. GC content must be recalculated at that resolution."
        echo "      Please supply all three reference files manually:"
        echo "        -g  chromosome sizes file"
        echo "        -a  gene annotations BED (initialA)"
        echo "        -b  genome FASTA file (for GC-content-based B-compartment seeds)"
        echo ""
    else
        for _b in "${CRUSH_SUPPORTED_BUILDS[@]}"; do
            if [ "$_b" == "$genomebuild" ]; then
                genomebuild_supported=1
                break
            fi
        done

        if [ "$genomebuild_supported" -eq 0 ]; then
            echo "Note: '${genomebuild}' does not have pre-built reference files."
            echo "      Pre-built references are available for: hg19, hg38, mm10, mm9."
            echo "      Please supply -g, -a, and -b manually for this genome."
        else
            # Download reference files directly into the crush working folder so the
            # algorithm can access them without any path resolution issues.
            CRUSH_REFCACHE="$(pwd)/${crushdir}/refs"
            mkdir -p "${CRUSH_REFCACHE}"
            JROWLEY_BASE="https://raw.githubusercontent.com/Kalluchiachyuth/CRUSH/main/genomes/${genomebuild}"

            _fetch_ref() {
                local url="$1"
                local dest="$2"
                local label="$3"
                if [ -s "$dest" ]; then
                    echo "Using existing ${label}: ${dest}"
                    return 0
                fi
                if ! command -v curl &> /dev/null; then
                    echo "Error: curl is required for --genomebuild but is not installed."
                    exit 1
                fi
                echo "Fetching ${label} for ${genomebuild}..."
                if ! curl -fsSL "$url" -o "$dest"; then
                    rm -f "$dest"
                    echo "Error: Could not fetch ${label} for build '${genomebuild}'."
                    echo "  URL tried: ${url}"
                    exit 1
                fi
                echo "Downloaded ${label} to ${dest}"
            }

            if [ -z "$sizefile" ]; then
                _fetch_ref \
                    "${JROWLEY_BASE}/${genomebuild}.chrom.sizes" \
                    "${CRUSH_REFCACHE}/${genomebuild}.chrom.sizes" \
                    "genome sizes"
                sizefile="${CRUSH_REFCACHE}/${genomebuild}.chrom.sizes"
            fi

            if [ -z "$genesfile" ]; then
                _fetch_ref \
                    "${JROWLEY_BASE}/${genomebuild}_genes.bed" \
                    "${CRUSH_REFCACHE}/${genomebuild}_genes.bed" \
                    "gene annotations (initialA)"
                genesfile="${CRUSH_REFCACHE}/${genomebuild}_genes.bed"
            fi

            if [ -z "$fastafile" ]; then
                _fetch_ref \
                    "${JROWLEY_BASE}/${genomebuild}_Bbins.bed" \
                    "${CRUSH_REFCACHE}/${genomebuild}_Bbins.bed" \
                    "B-compartment seeds (Bbins)"
                fastafile="${CRUSH_REFCACHE}/${genomebuild}_Bbins.bed"
            fi
        fi
    fi
fi

if [[ -z "$hicpath" ]]; then
    echo "You must specify a .hic or .mcool file! (local file or HTTPS link)"
    exit 1
fi

if [[ "$hicpath" =~ ^https?:// ]]; then
    if [[ ! "$hicpath" =~ \.(hic|mcool)$ ]]; then
        echo "Error: URL must point to a .hic or .mcool file"
        exit 1
    fi
else
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

if [ "$genomebuild_supported" -eq 0 ] && { [ -z "$sizefile" ] || [ -z "$genesfile" ] || [ -z "$fastafile" ]; }; then
    echo "Missing required arguments: supply either -gb GENOMEBUILD (hg19/hg38/mm10/mm9)"
    echo "or all three of --genomesize, --initialA, and --initialB manually."
    help
    exit 1
fi

if [ $doNotMerge == 1 ]; then
    echo "Warning: doNotMerge set. We recommend merging to achieve the highest possible confidence and resolution. Only set doNotMerge if you plan to merge the individual resolutions later."
fi

if [ $adjustment == 1 ]; then
    echo "Warning: adjusting end compartmental values which may result in some loss of quantitative power. We recommend adjustment only for troubleshooting or when not comparing two samples."
fi

echo "CRUSH_v""$version"" --hic ""$hicpath"" --res ""$res"" --genomesize ""$sizefile"" --initialA ""$genesfile"" --initialB ""$fastafile"" --cpu ""$cpu"" --no-merge ""$doNotMerge"" --adjustment ""$adjustment"" --distance ""$distance"" --upperlim ""$upperlim"" --lowerthresh ""$lowerthresh"" --trackline ""$trackline"" --threshold ""$threshold"" --window ""$window"" --switch ""$switch"" --maxres ""$coarsestres"" --norm ""$norm"" --eigenres ""$eigenres"" --smoothing ""$smoothing"" --genomebuild ""$genomebuild" > "$outpre"CRUSHparamters.txt

juiceorcool=`echo "$hicpath" | sed 's/\./\t/g' | sed 's/\./\t/g' | awk '{if ($NF == "hic") print 0; else if ($NF == "mcool") print 1; else print 2}'`

if [ $juiceorcool == 2 ]; then
    echo "You must choose a Hi-C file that ends in either .hic (juicer) or .mcool (cooler)"
    exit 0
elif [ $juiceorcool == 0 ]; then
    echo "Identified .hic juicer format. We will use the pythonic hicstraw to extract the data and assume you have it installed, using pip install hic-straw."
else
    echo "Identified .mcool format. We will use cooler to extract the data and assume you have it installed in your path."
fi

# ------------------------------------------------------------------------------
# CHR PREFIX CHECKER
# GitHub reference files are hosted without a chr prefix (plain '1', '2', 'X').
# This block detects the exact prefix the Hi-C file uses (chr / Chr / CHR / none)
# and the prefix in the sizefile, then:
#   --genomebuild runs : auto-converts the downloaded files to match the Hi-C
#   manual -g -a -b    : exits with a clear fix message if there is a mismatch
# ------------------------------------------------------------------------------

# Extract the exact leading alpha prefix from the first chromosome in the sizefile.
# e.g. "chr1" -> "chr", "CHR1" -> "CHR", "1" -> ""
sizechrprefix=$(awk 'NR==1{
    if (tolower(substr($1,1,3)) == "chr") print substr($1,1,3)
    else print ""
}' "$sizefile")

# Extract the exact prefix from the Hi-C file chromosome names.
if [ $juiceorcool -eq 0 ]; then
    hicchrprefix=$(python3 - "$hicpath" <<'PEOF'
import hicstraw, sys, re
try:
    hic = hicstraw.HiCFile(sys.argv[1])
    chroms = [c.name for c in hic.getChromosomes() if c.name.lower() != "all"]
    if chroms:
        m = re.match(r'^(chr)', chroms[0], re.IGNORECASE)
        print(chroms[0][:3] if m else '')
    else:
        print('')
except Exception:
    print('?')
PEOF
)
else
    firstres=$(cooler ls "$hicpath" 2>/dev/null | head -1)
    hicchrprefix=$(cooler info "${firstres}" 2>/dev/null | python3 -c "
import json, sys, re
try:
    info = json.load(sys.stdin)
    chroms = info.get('chromnames', [])
    if chroms:
        m = re.match(r'^(chr)', chroms[0], re.IGNORECASE)
        print(chroms[0][:3] if m else '')
    else:
        print('')
except Exception:
    print('?')
")
fi

if [ "$hicchrprefix" != "?" ] && [ "$sizechrprefix" != "$hicchrprefix" ]; then

    if [ "$genomebuild_supported" -eq 1 ]; then
        # Auto-convert: strip any existing chr prefix from col 1 then add the
        # exact prefix the Hi-C file uses (handles chr->CHR, none->chr, etc.)
        echo ""
        echo "Note: Hi-C chromosomes use '${hicchrprefix:-no prefix}'; reference files"
        echo "      use '${sizechrprefix:-no prefix}'. Auto-converting reference files..."

        _convert_chr() {
            local infile="$1"
            local outfile="$2"
            awk -v p="$hicchrprefix" '
            BEGIN { FS="\t"; OFS="\t" }
            {
                if (tolower(substr($1,1,3)) == "chr") $1 = substr($1,4)
                if (p != "") $1 = p $1
                print
            }' "$infile" > "$outfile"
        }

        _convert_chr "$sizefile"  "${CRUSH_REFCACHE}/converted.chrom.sizes"
        _convert_chr "$genesfile" "${CRUSH_REFCACHE}/converted_genes.bed"
        _convert_chr "$fastafile" "${CRUSH_REFCACHE}/converted_Bbins.bed"

        sizefile="${CRUSH_REFCACHE}/converted.chrom.sizes"
        genesfile="${CRUSH_REFCACHE}/converted_genes.bed"
        fastafile="${CRUSH_REFCACHE}/converted_Bbins.bed"

        echo "      Done. Running with '${hicchrprefix:-no prefix}' chromosome names."
        echo ""

    else
        # Manual files: tell the user exactly what to fix.
        echo ""
        echo "ERROR: Chromosome naming mismatch between your Hi-C file and your reference files."
        if [ -n "$hicchrprefix" ] && [ -z "$sizechrprefix" ]; then
            echo "  Hi-C file uses '${hicchrprefix}' prefix  (e.g. '${hicchrprefix}1')"
            echo "  Your -g / -a / -b files use no prefix   (e.g. '1')"
            echo "  Fix: add '${hicchrprefix}' to column 1 of your -g, -a, and -b files."
            echo "       e.g.  awk 'BEGIN{OFS=\"\t\"} {\$1=\"${hicchrprefix}\"\$1; print}' sizes.txt"
        elif [ -z "$hicchrprefix" ] && [ -n "$sizechrprefix" ]; then
            echo "  Hi-C file uses no prefix               (e.g. '1', '2', 'X')"
            echo "  Your -g / -a / -b files use '${sizechrprefix}' prefix  (e.g. '${sizechrprefix}1')"
            echo "  Fix: strip the prefix from column 1 of your -g, -a, and -b files."
            echo "       e.g.  awk 'BEGIN{OFS=\"\t\"} {sub(/^${sizechrprefix}/,\"\",\$1); print}' sizes.txt"
        else
            echo "  Hi-C file uses '${hicchrprefix}' prefix, reference files use '${sizechrprefix}'."
            echo "  Fix: ensure chromosome names in -g, -a, and -b match your Hi-C file exactly."
        fi
        echo ""
        exit 1
    fi
fi


################################################################################
#                      CORE PROCESSING FUNCTIONS                             #
################################################################################


# ------------------------------------------------------------------------------
# HI-C DATA EXTRACTION
# ------------------------------------------------------------------------------

run_dumper() {
    local readtype=$1
    local norm=$2
    local hicinnie=$3
    local chrom=$4
    local res=$5
    local dumpoutie=$6
    local dumperrors=$7

    echo "Dumping $readtype reads using $norm normalization at resolution $res for $chrom."

    local script_dir
    script_dir="$(dirname "$dumpoutie")"
    local script_base
    script_base="$(basename "$dumpoutie")"
    local script_path="${script_dir}/run_dumper_${script_base}.py"

    mkdir -p "$script_dir"

    cat > "$script_path" << EOF
import hicstraw

readtype = "$readtype"
norm = "$norm"
hicfile = "$hicinnie"
chrom = "$chrom"
res = int("$res")

result = hicstraw.straw(readtype, norm, hicfile, chrom, chrom, "BP", res)
for record in result:
    print("{0}\t{1}\t{2}".format(record.binX, record.binY, record.counts))
EOF

    python3 "$script_path" 2> "$dumperrors" | grep -v WARN > "$dumpoutie"
}


# ------------------------------------------------------------------------------
# A/B COMPARTMENT INITIALIZATION
# ------------------------------------------------------------------------------

process_genes_Bbins() {
    local myres=$1
    local mychr=$2
    local genesblock=$3
    local sortedbed=$4
    local pseudoB=$5

    if [[ "$myres" == "$maxres" && "$genesblock" -eq 1 ]]; then

        echo "Using initial Genes and Bbins files for $maxres"
        maxres_processed=1

        cat $EVAstates | mawk -v myres=$maxres -v mychr=$mychr \
            '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
            awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
            > $sortedbed

        cat $EVBstates | mawk -v myres=$maxres -v mychr=$mychr \
            '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
            awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
            > $pseudoB

    elif [ "$switch" -lt 1 ]; then

        cat $EVAstates | mawk -v myres=$myres -v mychr=$mychr \
            '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
            awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
            > $sortedbed

        # BUG 24 FIX: was "$sortedbed" — must be "$pseudoB"
        cat $EVBstates | mawk -v myres=$myres -v mychr=$mychr \
            '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
            awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
            > $pseudoB

    else

        echo "Using combined bins for $myres"

        cat $combined_newAbins | mawk -v myres=$myres -v mychr=$mychr \
            '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
            awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
            > $sortedbed

        cat $combined_newBbins | mawk -v myres=$myres -v mychr=$mychr \
            '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
            awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
            > $pseudoB

    fi
}


# ------------------------------------------------------------------------------
# EIGENVECTOR ANALYSIS
# ------------------------------------------------------------------------------

run_EigenVector() {
    local hic_file="$1"
    local chrom="$2"
    local res="$3"
    local n_components=10
    local eigendumpoutie="Eigen_dump_${chrom}_${res}.txt"
    local eigendumperrors="Eigen_errors_${chrom}_${res}.log"
    local ev_script="run_EigenVector_${chrom}_${res}.py"

    echo "Generating EV's from 1 to $n_components for $chrom at resolution $res."

    local file_ext="${hic_file##*.}"
    if [[ "$file_ext" == "hic" ]]; then
        run_dumper oe VC_SQRT "$hic_file" "$chrom" "$res" "$eigendumpoutie" "$eigendumperrors"
    elif [[ "$file_ext" == "cool" || "$file_ext" == "mcool" ]]; then
        echo "Detected cooler format: dumping with cooler CLI..."
        local cooler_uri="$hic_file"
        if [[ "$file_ext" == "mcool" ]]; then
            cooler_uri="${hic_file}::/resolutions/${res}"
        fi
        cooler dump "$cooler_uri" -r "$chrom" --join | cut -f 2,5,7 > "$eigendumpoutie" 2>"$eigendumperrors"
    else
        echo "Unsupported Hi-C file type: $hic_file" >&2
        return 1
    fi

    if [[ ! -s "$eigendumpoutie" ]]; then
        echo "NOTE: No Hi-C data for $chrom at resolution $res — skipping." >&2
        return 1
    fi

    cat > "$ev_script" << PYEOF
import numpy as np
import pandas as pd
from scipy.linalg import eigh
import sys

def save_eigenvector_to_file(eigenvector, bin_starts, res, chrom, output_file):
    with open(output_file, 'w') as f:
        for bin_start, value in zip(bin_starts, eigenvector):
            f.write(f"{chrom}\t{bin_start}\t{bin_start + res}\t{value}\n")

def process_hic_to_eigen(dump_file, res, chrom, n_components=10):
    dump_df = pd.read_csv(dump_file, sep='\t', header=None, names=['binX', 'binY', 'counts'])
    dump_df = dump_df.apply(pd.to_numeric, errors='coerce').dropna()
    if dump_df.empty:
        print(f"No valid data in {dump_file}", file=sys.stderr)
        sys.exit(1)

    binX   = dump_df['binX'].values.astype(int)
    binY   = dump_df['binY'].values.astype(int)
    counts = dump_df['counts'].values

    n_bins = max(binX.max(), binY.max()) // res + 1
    matrix = np.zeros((n_bins, n_bins))
    for x, y, c in zip(binX, binY, counts):
        xi, yi = x // res, y // res
        if xi < n_bins and yi < n_bins:
            matrix[xi, yi] = c
            if x != y:
                matrix[yi, xi] = c

    mask = np.sum(matrix, axis=0) > 0
    if mask.sum() < 2:
        print(f"Too few non-zero bins for {chrom}", file=sys.stderr)
        sys.exit(1)

    matrix = matrix[mask][:, mask]
    pearson_matrix = np.corrcoef(matrix)
    pearson_matrix = np.nan_to_num(pearson_matrix)

    eigenvalues, eigenvectors = eigh(pearson_matrix)
    idx = eigenvalues.argsort()[::-1]
    eigenvectors = eigenvectors[:, idx]
    n_components = min(n_components, eigenvectors.shape[1])
    eigenvectors = eigenvectors[:, :n_components]

    bin_starts = np.arange(0, n_bins * res, res)[mask]
    for i in range(n_components):
        out = f"PC{i+1}_{chrom}_{res}.bedgraph"
        save_eigenvector_to_file(eigenvectors[:, i], bin_starts, res, chrom, out)

if __name__ == "__main__":
    dump_file = sys.argv[1]
    res       = int(sys.argv[2])
    chrom     = sys.argv[3]
    process_hic_to_eigen(dump_file, res, chrom)
PYEOF

    python3 "$ev_script" "$eigendumpoutie" "$res" "$chrom"
    local py_exit=$?
    if [[ $py_exit -ne 0 ]]; then
        echo "ERROR: EigenVector Python script failed for $chrom (exit $py_exit)" >&2
        return 1
    fi
}


# ------------------------------------------------------------------------------
# EIGENVECTOR ORIENTATION
# ------------------------------------------------------------------------------

check_and_flip_sign() {
    local pc_file="$1"
    local genes_file="$2"
    local output_file="sign_check_${pc_file}.txt"

    echo "Checking and potentially flipping signs for $pc_file based on gene content from $genes_file..."

    local pos_neg
    pos_neg=$(
        awk 'BEGIN{FS=OFS="\t"} {if ($4>0) print $0,"positive"; else if ($4<0) print $0,"negative"}' "$pc_file" |
        bedtools intersect -u -a stdin -b "$genes_file" |
        sort -k1,1 -V -k5b,5b |
        bedtools groupby -i stdin -g 1,5 -c 1 -o count |
        bedtools groupby -i stdin -g 1   -c 3 -o collapse |
        sed 's/,/\t/g' |
        awk 'BEGIN{FS=OFS="\t"} {if (NF>=3 && $2+0 > $3+0) print "-"; else print "+"}'
    )

    echo "Overlap check result for $pc_file: $pos_neg" > "$output_file"

    if [[ "$pos_neg" == "-" ]]; then
        echo "Flipping signs for $pc_file." >> "$output_file"
        awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,-1*$4}' "$pc_file" > "${pc_file}.tmp" \
            && mv "${pc_file}.tmp" "$pc_file"
    else
        echo "No flipping needed for $pc_file." >> "$output_file"
    fi

    echo "Check and flip sign process completed for $pc_file."
}

process_all_pcs_for_chromosome() {
    local chrom="$1"
    local res="$2"
    local genes_file="$3"

    for pc_file in PC*_${chrom}_${res}.bedgraph; do
        echo "Processing $pc_file..."
        check_and_flip_sign "$pc_file" "$genes_file"
    done

    echo "All PCs for chromosome $chrom at resolution $res have been processed."
}


# ------------------------------------------------------------------------------
# PRINCIPAL COMPONENT SELECTION
# ------------------------------------------------------------------------------

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

    local bin_size=$res
    local num_bins=$((chrom_size / bin_size))

    calculate_density() {
        local feature_file="$1"
        awk -v bin_size="$bin_size" -v chrom="$chrom" \
        'BEGIN { OFS="\t"; for (i = 0; i < ENVIRON["num_bins"]; i++) density[i] = 0 }
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

    for pc_file in PC*_"${chrom}"_"${res}".bedgraph; do
        if [[ ! -s "$pc_file" ]]; then
            echo "Skipping empty or missing PC file: $pc_file" >> "$log_file"
            continue
        fi

        local pc_values=($(awk -v bin_size="$bin_size" -v num_bins="$num_bins" '
        BEGIN { for (i = 0; i < num_bins; i++) pc[i] = 0 }
        { bin = int($2 / bin_size); if (bin < num_bins) pc[bin] = $4 }
        END { for (i = 0; i < num_bins; i++) print pc[i] }' "$pc_file"))

        if [[ ${#pc_values[@]} -eq 0 ]]; then
            echo "No valid data in $pc_file" >> "$log_file"
            continue
        fi

        local gene_density_str=$(IFS=,; echo "${gene_density[*]}")
        local gc_density_str=$(IFS=,; echo "${gc_density[*]}")
        local pc_values_str=$(IFS=,; echo "${pc_values[*]}")

        read -r gene_corr gc_corr <<< $(python3 - <<EOF
import numpy as np
from scipy.stats import pearsonr

try:
    pc_values    = np.array([${pc_values_str}])
    gene_density = np.array([${gene_density_str}])
    gc_density   = np.array([${gc_density_str}])

    gene_corr, _ = pearsonr(pc_values, gene_density)
    gc_corr,   _ = pearsonr(pc_values, gc_density)

    print(gene_corr, gc_corr)
except Exception as e:
    print(0, 0)
EOF
)

        if [[ -z "$gene_corr" || -z "$gc_corr" || "$gene_corr" == "0" && "$gc_corr" == "0" ]]; then
            echo "Skipping $pc_file: Invalid correlations (gene_corr=$gene_corr, gc_corr=$gc_corr)" >> "$log_file"
            continue
        fi

        local total_corr=$(echo "$gene_corr - $gc_corr" | bc -l)

        if (( $(echo "$total_corr > $best_corr" | bc -l) )); then
            best_pc="$pc_file"
            best_corr="$total_corr"
            echo "New best PC: $pc_file with total correlation: $total_corr" >> "$log_file"
        fi
    done

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

generateEV() {
    local hic_file="$1"
    local chrom="$2"
    local res="$3"
    local genes_file="$4"
    local gc_file="$5"
    local chrom_size="$6"

    echo "Starting EV generation for chromosome $chrom at resolution $res..."

    run_EigenVector "$hic_file" "$chrom" "$res" || {
        echo "NOTE: Skipping chromosome $chrom — no Hi-C data found at resolution $res." >&2
        return 1
    }

    echo "PCs generated for $chrom at resolution $res."

    echo "Flipping the signs for $chrom PC's at resolution $res..."
    process_all_pcs_for_chromosome "$chrom" "$res" "$genes_file"

    echo "Finding best PC for $chrom at resolution $res..."
    find_best_pc "$chrom" "$res" "$genes_file" "$gc_file" "$chrom_size"

    echo "Best PC identification completed for $chrom at resolution $res."
}


# ------------------------------------------------------------------------------
# CONCATENATE BEST EIGENVECTORS
# ------------------------------------------------------------------------------

concatenate_best_pcs() {
    echo "Concatenating the best eigenvectors for all chromosomes into a full genome file."
    local files=( best_Eigen_*_${eigenres}.bedgraph )
    if [[ ! -e "${files[0]}" ]]; then
        echo "ERROR: No best_Eigen_* files found to concatenate." >&2
        return 1
    fi
    cat "${files[@]}" > EV_full_genome.bedgraph
    echo "Concatenation completed. Full genome eigenvector file: EV_full_genome.bedgraph"
}


# ------------------------------------------------------------------------------
# SMOOTHING
# ------------------------------------------------------------------------------

run_smoothing() {
    local innieS=$1
    local outieS=$2

    cat "$innieS" | sort -k 1,1 -V -k 2bn,2b --stable > "$outieS"
    cat "$outieS" | tail -n +2 > tmp1_shiftedright
    cat "$outieS" | awk '{if (NR == 1) print $0"\n"$0; else print $0}' > tmp2_shiftedleft
    paste "$outieS" tmp1_shiftedright tmp2_shiftedleft | \
        awk '{print $1"\t"$2"\t"$3"\t"(($4 + $8 + $12)/3)}' > shifted_outie
    mv shifted_outie "$innieS"
}


# ------------------------------------------------------------------------------
# CRUSH SCORE CALCULATION
# ------------------------------------------------------------------------------

process_oppocheck_statement() {
    local innie_genes="$1"
    local innie_gc="$2"
    local outie_oppo="$3"
    local chr="$4"
    local res="$5"
    local input_dir="$6"
    local output_dir="$7"

    local zfile="${input_dir}/zerofile"
    local rowZ="${input_dir}/RowZs"

    if [ ! -f "$zfile" ] || [ ! -f "$rowZ" ]; then
        echo "Required file missing: $zfile or $rowZ"
        return 1
    fi

    mawk 'NR==FNR { c1[$1] = $2; next } { if ($1 in c1) print $2 "\t" $3; next } { print $0 }' \
        "$innie_genes" "$rowZ" "$zfile" | \
        awk '{ b[$1] += $2; c[$1]++; d[$1] += ($2**2) } END { for (i in b) { print i "\t" (b[i])*(c[i]) "\t" c[i] "\t" b[i]/c[$1] "\t" ((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } }' \
        > "${output_dir}/A2removedfile" &

    mawk 'NR==FNR { c1[$1] = $2; next } { if ($1 in c1) print $2 "\t" $3; next } { print $0 }' \
        "$innie_gc" "$rowZ" "$zfile" | \
        awk '{ b[$1] += $2; c[$1]++; d[$1] += ($2**2) } END { for (i in b) { print i "\t" (b[i])*(c[i]) "\t" c[i] "\t" b[i]/c[$1] "\t" ((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } }' \
        > "${output_dir}/B2removedfile" &

    wait

    myABsum=$(cat "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | mawk '{ sum += $3 } END { print sum }')

    wait

    if [ "$verbose" -gt 0 ]; then
        end=$(date +%s)
        seconds=$(echo "$end - $start" | bc)
        echo "It took $seconds seconds. Comparing A vs B scores for $chr chromosome."
        start=$(date +%s)
    fi

    mawk -v var="$res" -v myABsum="$myABsum" \
        '{ if ($1 in b) b[$1] -= $2; else b[$1] = $2 } END { for (i in b) { if (b[i] != 0) print i "\t" (b[i])/myABsum } }' \
        "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | \
        sort -k 1bn,1b --stable | \
        mawk -v mchr="$chr" -v myres="$res" -v window="$newwindow" -v OFS="\t" \
        'BEGIN { slide = 1 } { mod = NR % window; if (NR <= window) { print mchr, int($1), $2; count++ } else { sum -= array[mod] } sum += $2; array[mod] = $2; } (NR % slide) == 0 { print mchr, int($1-((myres*window)/2)), sum/count }' | \
        mawk -v myres="$res" '{ print $1 "\t" $2-int(myres/2)+myres "\t" $2+(int(myres/2))+myres "\t" $3 }' | \
        mawk '{ if ($2 > 0) print $0 }' > "$outie_oppo" 2> "${output_dir}/errorfile"

    if [ "$exclbed" != 0 ]; then
        exclbins="${output_dir}/exclbins"
        excloutie="${output_dir}/exclout_${res}_${chr}"

        cat "$exclbed" | mawk -v myres="$res" -v mychr="$chr" \
            '{ if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres }' | \
            awk -v myres="$res" '{ for (i = $2; i <= $3; i += myres) b[i] += 1 } END { for (j in b) print j "\t" b[j] }' \
            > "$exclbins"

        wait

        mawk 'NR==FNR { c1[$1] = $2; next } { if ($2 in c1); else print $0 }' \
            "$exclbins" "$outie_oppo" > "$excloutie"

        wait

        mv "$excloutie" "$outie_oppo"
    fi

    if [[ "$oppocheck" -gt 0 && "$adjustment" -gt 0 ]]; then

        amed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1) print $4}' $innie_genes $outie_oppo | mawk '{sum +=$1} END {print sum/NR}'`
        bmed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $4}' $innie_genes $outie_oppo | mawk '{sum += $1} END {print sum/NR}'`

        wait

        cat $outie_oppo | mawk -v varA=$amed -v varB=$bmed \
            '{print $1"\t"$2"\t"$3"\t"$4-(varA-(varB*-1))}' > $outie2

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

run_midshifter() {
    local finer_res=$1
    local coarser_res=$2
    local high_innie=$3
    local coarse_innie=$4
    local shifter_outie=$5

    echo "Running shifter with Finer Resolution: $finer_res and Coarser Resolution: $coarser_res"
    echo "High-res file: $high_innie, Coarse-res file: $coarse_innie"

    cat << EOF > hdrescompare.py
import numpy as np
from tqdm import tqdm

sizefile = "$sizefile"
res = int("$finer_res")
coarseres = int("$coarser_res")
infile1 = "$high_innie"
infile2 = "$coarse_innie"
outiewrite = "$shifter_outie"

chrnames = []
chrsizes = []
with open(sizefile, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        chrnames.append(lines[0])
        chrsizes.append(int(lines[1]))

chrmat = [[] for _ in chrnames]
smat = [[] for _ in chrnames]
emat = [[] for _ in chrnames]
scoremat = [[] for _ in chrnames]
coscoremat = [[] for _ in chrnames]

for i in tqdm(range(len(chrnames))):
    for j in range(0, chrsizes[i] + coarseres + coarseres, res):
        chrmat[i].append(chrnames[i])
        smat[i].append(j)
        emat[i].append(j + res)
        scoremat[i].append(0)
        coscoremat[i].append(0)

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

window_size = 3 * (coarseres // res)

rolling_avg_scores = [[] for _ in chrnames]
for c in range(len(scoremat)):
    for i in range(len(scoremat[c])):
        window_start = max(0, i - window_size // 2)
        window_end = min(len(scoremat[c]), i + window_size // 2)
        window_scores = scoremat[c][window_start:window_end]
        rolling_avg_scores[c].append(np.mean(window_scores))

shifted_scores = [[0 for _ in range(len(scoremat[c]))] for c in range(len(scoremat))]

for c in range(len(scoremat)):
    for i in range(len(scoremat[c])):
        mymean = rolling_avg_scores[c][i]
        coarscore = coscoremat[c][i] if i < len(coscoremat[c]) else 0
        shiftval = mymean - coarscore
        newval = scoremat[c][i] - shiftval
        shifted_scores[c][i] = newval

with open(outiewrite, "w") as ofile:
    for w in range(len(shifted_scores)):
        for p in range(len(shifted_scores[w])):
            ofile.write(f"{chrmat[w][p]}\t{smat[w][p]}\t{emat[w][p]}\t{shifted_scores[w][p]}\n")
EOF

    python hdrescompare.py

    echo "After Shifter: Executed with Coarser Resolution: $coarser_res, Finer Resolution: $finer_res"
    echo "Shifter executed and output written to $shifter_outie"
}


# ------------------------------------------------------------------------------
# COMPARTMENT RECLASSIFICATION
# ------------------------------------------------------------------------------

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

        echo "Re-initialization of Bins is switched ON. Re-evaluating both A and B bins and re-initializing for next resolution."

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


calculate_pvalues() {
    local output_dir=$1
    local chr=$2
    local res=$3
    local outie_pvalues=$4

    local zscript="zlookup_${chr}_${res}.py"
    local zfileforpy="${output_dir}/ttable_${chr}_${res}"
    local pfileforpy="${output_dir}/ptable_${chr}_${res}"

    awk 'NR==FNR { c[$1] = $3; m[$1] = $4; sd[$1] = $5; next } { if (($1 in m) && ($3 != 1) && (c[$1] != 1)) print $1 "\t" (m[$1]-$4) / (((sd[$1]**2)/c[$1]) + (($5**2)/$3))**.5 "\t" (c[$1] + $3 -2) }' \
        "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" > "$zfileforpy"

    cat > "$zscript" << PYEOF
import scipy.stats
inputfile = "$zfileforpy"
outfile   = "$pfileforpy"
mybins  = []
mypvals = []
with open(inputfile, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        mybins.append(str(lines[0]))
        mytstat = float(lines[1])
        mydf    = int(lines[2])
        mypvals.append(scipy.stats.norm.sf(abs(mytstat)))
with open(outfile, 'w') as outie:
    for w in range(len(mybins)):
        outie.write(str(mybins[w]) + "\t" + str(mypvals[w]) + "\n")
PYEOF

    python3 "$zscript"

    mawk -v mchr="$chr" -v myres="$res" \
        '{ print mchr "\t" $1 "\t" $1+myres "\t" $2 }' "$pfileforpy" > "$outie_pvalues"

    wait
}


################################################################################
#                    MAIN RESOLUTION PROCESSING LOOP                         #
################################################################################

declare -a maxres_processed=0


# ------------------------------------------------------------------------------
# RESOLUTION PROCESSOR
# ------------------------------------------------------------------------------

process_resolution() {
    local myres=$1
    local genesblock=$2

    if [ $window == 0 ]; then
        totbins=$(awk -v var=$myres '{print $1"\t"($2/var)**2}' $sizefile | awk '{sum += $2} END {print sum}')
        newwindow=1
    else
        newwindow=$(mawk -v win=$window -v res=$myres 'BEGIN{print int(win/(res/1000))+1}')
    fi

    outiefull="Crush_${myres}.bedgraph"

    echo -e "\nNow dumping reads and calculating initial CRUSH for all chromosomes at $myres"

    task() {
        mychr=$(awk -v var=$myiter 'NR==var {print $1}' $sizefile)
        mysize=$(awk -v var=$myiter 'NR==var {print $2}' $sizefile)
        rowbins=$(awk -v myres=$myres -v var=$myiter 'NR==var {print int($2/myres)+1}' $sizefile)

        outie="Crush_${myres}_${mychr}_tmp"
        outie2="Crush_${myres}_${mychr}_tmp2"
        tmpfiles="ABtmpfiles_${mychr}_${myres}"

        [ ! -d "$tmpfiles" ] && mkdir "$tmpfiles"

        mydist=$((distance + 0))
        myupper=$((upperlim + 0))

        if [ "$verbose" -gt 0 ]; then
            echo "Getting genic bins."
        fi

        sortedbed="${tmpfiles}/genebins"
        pseudoB="${tmpfiles}/pseudoB"
        dumped="${tmpfiles}/dumped_${myres}_${mychr}_tmp"

        process_genes_Bbins $myres $mychr $genesblock $sortedbed $pseudoB

        if [ "$verbose" -gt 0 ]; then
            echo "Dumping ${mychr} chromosomal reads."
        fi

        newofile="${tmpfiles}/observeds"
        newefile="${tmpfiles}/expecteds"
        dumperrors="dumperrors_${myres}"

        > "$newofile"
        > "$newefile"

        if [ $juiceorcool -eq 0 ]; then
            run_dumper 'observed' $norm $hicpath $mychr $myres $dumped $dumperrors
        else
            echo "Dumping using cooler..."
            cooler dump $hicpath::/resolutions/$myres -r $mychr --join | cut -f 2,5,7 > $dumped 2> $dumperrors
        fi

        if [ ! -s "$dumped" ]; then
            echo "WARNING: no dumped reads for $mychr at $myres (file: $dumped). Check $dumperrors" >&2
            return
        fi

        wait

        if [ "$upperlim" -gt 0 ]; then
            awk -v var=$myres -v distfilt=$mydist -v upper=$upperlim \
                '{if (($2-$1 >= distfilt) && ($2-$1 <= upper)) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' $dumped | \
            awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile \
                '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile }' &
        else
            awk -v var=$myres -v distfilt=$mydist \
                '{if ($2-$1 >= distfilt) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' $dumped | \
            awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile \
                '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile }' &
        fi

        wait

        if [ "$verbose" -gt 0 ]; then
            echo "Done dumping."
            echo "Now getting distance normalized values for ${mychr} chromosome."
            start=$(date +%s)
        fi

        errorfile="errorlog"
        tmpmean="${tmpfiles}/tmpmean"
        symmfile="${tmpfiles}/symmonfile"

        mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 != $3) print $2"\t"$3"\t"($4+1)/(c1[$1]+1)}' \
            $newefile $newofile | \
            mawk -v thresh=$threshold \
            '{if ($3 < thresh) print $1"\t"$2"\t"$3"\n"$2"\t"$1"\t"$3; else print $1"\t"$2"\t"thresh"\n"$2"\t"$1"\t"thresh}' | \
            awk -v symmfile=$symmfile -v tmpmean=$tmpmean -v rowbins=$rowbins \
            '{m[$1]+=$3;d[$1]+=($3**2); print $0 >> symmfile} END {for (i in m) print i"\t"m[i]/rowbins"\t"((d[i]/rowbins-((m[i]/rowbins)**2))**.5) >> tmpmean }' \
            2> $tmpfiles/$errorfile

        wait

        if [ "$verbose" -gt 0 ]; then
            end=$(date +%s)
            seconds=$(echo "$end - $start" | bc)
            start=$(date +%s)
        fi

        if [ "$myres" -gt 50000 ]; then
            oppocheck=$(mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) ; else if ($2 in c1); else print $0}' $sortedbed $symmfile | wc -l)
        else
            oppocheck=2
        fi

        if [ "$oppocheck" -lt 1 ] && [ "$verbose" -gt 0 ]; then
            echo "WARNING: All bins on chromosome ${mychr} overlap with your bed file at ${myres} resolution. This should not affect the scores if using multiple resolutions. If using a single resolution, then maybe not use this resolution?"
        fi

        if [ "$verbose" -gt 0 ]; then
            end=$(date +%s)
            seconds=$(echo "$end - $start" | bc)
            start=$(date +%s)
        fi

        RowZs=`echo "RowZs"`
        mawk 'NR==FNR { m[$1] = $2; d[$1] = $3; next} {if (d[$2] > 0) print $1 "\t" $2 "\t" ($3 - m[$2]) / (d[$2])}' \
            $tmpmean $symmfile > $tmpfiles/$RowZs

        zfile="zerofile"
        mawk -v var=$mychr '{if ($1 == var) print $0}' $sizefile | \
            mawk -v wlkr=$myres '{for (i=0;i<=$2;i+=wlkr) print i"\t"i"\t0.1"}' > $tmpfiles/$zfile

        process_oppocheck_statement $sortedbed $pseudoB $outie $mychr $myres \
            "ABtmpfiles_${mychr}_${myres}" "ABtmpfiles_${mychr}_${myres}"
    }

    open_sem $cpu

    if [ "$whichchoice" == "o" ]; then
        for (( myiter=1; myiter<=$totchroms; myiter++ )); do
            run_with_lock task $myiter
        done
    fi

    wait

    echo "Resolution: $myres for Resolution output"
    echo -ne "Merging individual chromosome files\033[0K\r"

    cat Crush_${myres}_*_tmp | grep -v -i nan | \
    awk '
    {
        key = $1"\t"$2"\t"$3;
        sum[key] += $4;
        count[key]++;
    }
    END {
        for (key in sum) {
            if (count[key] > 1) {
                avg = sum[key] / count[key];
                print key "\t" avg;
            } else {
                print key "\t" sum[key];
            }
        }
    }' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull

    wait

    rm Crush_${myres}_*_tmp

    if [ "$myres" == "$maxres" ]; then
        coarse_infile=`echo "GI_""$myres"".bedgraph"`
        cat $outiefull > $coarse_infile
        wait
    fi

    res_values=`echo "res_values"`
    echo "$myres" >> $res_values

    if [ "$countres" -gt 0 ]; then
        python_output=`echo "shifter.bedgraph"`

        if [ "$coarse_res" -eq "$maxres" ] && [ "$coarse_res" != "$high_res" ]; then
            high_res_infile=`echo "$outiefull"`
            coarse_res_infile=`echo "$coarse_infile"`
        else
            high_res_infile=`echo "$outiefull"`
            coarse_res_infile=`echo "$python_output"`
        fi

        if [[ "$(( $coarse_res % $high_res ))" -eq 0 ]]; then
            crtmp=$coarse_res
        else
            if [[ "$(( $crtmp % $high_res ))" -eq 0 ]]; then
                coarse_res=$crtmp
            fi
        fi

        run_midshifter $high_res $coarse_res $high_res_infile $coarse_res_infile $python_output

        echo "Resolution: "$myres" for shifter output"

        separating_ABbins $python_output

        coarse_res=$high_res

    else

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

        local chromosomes=($(cut -f1 $sizefile))

        for mychr in "${chromosomes[@]}"; do
            local sortedbed_reprocess="genebins_reprocess_${mychr}_${prev_res}"
            local pseudoB_reprocess="pseudoB_reprocess_${mychr}_${prev_res}"
            local outie_reprocess="Crush_reprocess_${prev_res}_${mychr}_tmp"

            reprocesstmpfiles=`echo "ABreprocesstmpfiles_""$mychr""_""$prev_res"`
            [ ! -d "$reprocesstmpfiles" ] && mkdir "$reprocesstmpfiles"

            echo "Using combined bins for $prev_res and the current value of resolution in this loop is $current_res"

            cat $combined_newAbins | mawk -v myres=$prev_res -v mychr=$mychr \
                '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
                awk -v myres=$prev_res '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
                > $sortedbed_reprocess

            cat $combined_newBbins | mawk -v myres=$prev_res -v mychr=$mychr \
                '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | \
                awk -v myres=$prev_res '{for (i=$2;i<=$3;i+=myres) b[i]+=1} END {for (j in b) print j"\t"b[j]}' \
                > $pseudoB_reprocess

            echo "Files: $sortedbed_reprocess, $pseudoB_reprocess"

            process_oppocheck_statement $sortedbed_reprocess $pseudoB_reprocess $outie_reprocess \
                $mychr $prev_res "ABtmpfiles_${mychr}_${prev_res}" "ABreprocesstmpfiles_${mychr}_${prev_res}"

            wait

            if [ $prev_res -eq $minres ] && [ $pcalculation -gt 0 ]; then
                pscores_reprocess="ttest_reprocess_${prev_res}_${mychr}_tmp"
                calculate_pvalues "ABreprocesstmpfiles_${mychr}_${prev_res}" $mychr $prev_res $pscores_reprocess
            fi
        done

        echo -ne "Merging individual chromosome files\033[0K\r"

        cat Crush_reprocess_${prev_res}_*_tmp | grep -v -i nan | \
        awk '
        {
            key = $1"\t"$2"\t"$3;
            sum[key] += $4;
            count[key]++;
        }
        END {
            for (key in sum) {
                if (count[key] > 1) {
                    avg = sum[key] / count[key];
                    print key "\t" avg;
                } else {
                    print key "\t" sum[key];
                }
            }
        }' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull_reprocess

        wait

        rm Crush_reprocess_${prev_res}_*_tmp

        if [ $prev_res -eq $minres ] && [ $pcalculation -gt 0 ]; then
            outiefullPval_reprocess="pvalues_Re_${prev_res}.bedgraph"
            echo "Resolution: $prev_res for Resolution output and minres is : $minres"
            cat ttest_reprocess_${prev_res}_*_tmp | grep -v -i nan | \
                sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefullPval_reprocess
            wait
        fi

        if [ "$prev_res" == "$maxres" ]; then
            coarse_infile_reprocess=`echo "GI_Re_""$prev_res"".bedgraph"`
            cat $outiefull_reprocess > $coarse_infile_reprocess
            wait
        fi

        prev_values=`echo "prev_values"`
        echo "$prev_res" >> $prev_values

        echo "$j"

        if [ "$j" -gt 0 ]; then

            if [ "$coarse_res_reprocess" -eq "$maxres" ] && [ "$coarse_res_reprocess" != "$highres_reprocess" ]; then
                highres_infile_reprocess=`echo "$outiefull_reprocess"`
                coarse_res_infile_reprocess=`echo "$coarse_infile_reprocess"`
            else
                highres_infile_reprocess=`echo "$outiefull_reprocess"`
                coarse_res_infile_reprocess=`echo "$python_output"`
            fi

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

            run_midshifter $highres_reprocess $coarse_res_reprocess \
                $highres_infile_reprocess $coarse_res_infile_reprocess $python_output

            echo "Shifter run completed for $prev_res"
            echo "Resolution: "$prev_res" for shifter output"

            separating_ABbins $python_output

            coarse_res_reprocess=$highres_reprocess
        else
            separating_ABbins $outiefull_reprocess
        fi
    done
}


################################################################################
#                          MAIN EXECUTION                                    #
################################################################################

echo "Creating and moving into temporary directory"
mkdir -p $crushdir
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
reslist=`cooler ls $hicpath | sed 's/\//\t/g' | \
    awk -v var=$res -v cres=$coarsestres '{if (($NF >= var) && ($NF <= cres)) print $NF}' | \
    sort -k 1bnr,1b --stable | \
    awk '{if (NR == 1) printf "%s", $1; else printf ",%s", $1}' | \
    awk '{print $0}'`
res=$reslist
fi

isfasta=`head -1 $fastafile | awk '{if ($1 ~ /^>/) print 1; else print 0}'`

echo "initial-A: $genesfile"
echo "sizes-file: $sizefile"

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

    mkdir -p gc_tmp
    awk '{print >> "gc_tmp/"$1".bed"; close("gc_tmp/"$1".bed")}' $resbins

    open_sem $cpu

    for chrom in $(cut -f 1 $sizefile); do
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

    expected_files=$(cut -f 1 $sizefile | wc -l)
    actual_files=$(ls gc_tmp/*.gc 2>/dev/null | wc -l)

    if [ "$actual_files" -lt "$expected_files" ]; then
        echo "Warning: Only got GC content for $actual_files chromosomes out of $expected_files expected"
    fi

    cat gc_tmp/*.gc > $gcfile
    rm -rf gc_tmp

    wait

    cat $gcfile | grep -v user | cut -f 1-5 | awk '{print $0"\t"($4+$5)}' | \
        awk '{if ($6 > 0) print $0}' | awk '{print $0"\t"($5/$6)}' > $gc_g_ga
    wait

    gc_thresh=`cat $gc_g_ga | \
        awk '{s+=$7; ss+=$7*$7; linecount+=1} END{print m=s/linecount, sqrt(ss/linecount-m^2)}' | \
        awk '{print $0}' | awk '{print $1 - 1*$2}'`

    if [ $verbose -gt 0 ]; then
        echo "gc threshold is ""$gc_thresh"
    fi

    cat $gc_g_ga | awk -v thresh=$gc_thresh '{if ($7 < thresh) print $0}' > $Bbins
    wait
    echo "Finished generating Bbins file for all chromosomes"

else
    cat $fastafile > $Bbins
fi

if [ $endZ == 0 ]; then
    endZ=`echo "$minres"`
fi

# ------------------------------------------------------------------------------
# EIGENVECTOR BLOCK
# ------------------------------------------------------------------------------

EVAstates="EVAstates.bed"
EVBstates="EVBstates.bed"

process_eigenvector_file() {
    local eigenfile="$1"
    local output_positive="$2"
    local output_negative="$3"

    awk '{if ($4 > 0) print $0}' "$eigenfile" > "$output_positive"
    if [ $? -ne 0 ]; then
        echo "Error creating $output_positive from $eigenfile."
        exit 1
    fi

    awk '{if ($4 < 0) print $0}' "$eigenfile" > "$output_negative"
    if [ $? -ne 0 ]; then
        echo "Error creating $output_negative from $eigenfile."
        exit 1
    fi

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

if [ -n "$eigenfile" ]; then
    echo "Using user-provided Eigen Vector file: $eigenfile"
    process_eigenvector_file "$eigenfile" "$EVAstates" "$EVBstates"
else
    echo "No user-provided Eigen Vector file. Calculating eigenvectors for all chromosomes at resolution $eigenres."

    gcbinned="gcbinned_all.bed"

    awk -v res=$eigenres -v chr_size_file=$sizefile '
    BEGIN { OFS="\t"; while ((getline < chr_size_file) > 0) sizes[$1] = $2 }
    $1 in sizes && $2 >= 0 && $3 >= 0 {
        start = int($2 / res) * res
        end = int($3 / res) * res
        if (end > sizes[$1]) end = sizes[$1]
        if (start < sizes[$1]) print $1, start, end
    }' "$Bbins" > "$gcbinned"

    if [[ ! -s $gcbinned ]]; then
        echo "Error: gcbinned file is empty. Check Bbins at eigenres $eigenres." >&2
        exit 1
    fi

    echo "Binned GC file generated: $gcbinned"


    export -f generateEV
    export -f run_EigenVector
    export -f run_dumper
    export -f process_all_pcs_for_chromosome
    export -f check_and_flip_sign
    export -f find_best_pc

    open_sem $cpu

    for chrom in $(cut -f1 "$sizefile"); do
        run_with_lock bash -c "
            chrom='${chrom}'
            sizefile='${sizefile}'
            gcbinned='${gcbinned}'
            hicpath='${hicpath}'
            eigenres='${eigenres}'
            genesfile='${genesfile}'

            chrom_size=\$(awk -v chr=\"\${chrom}\" '\$1==chr{print \$2}' \"\${sizefile}\")
            if [[ -z \"\$chrom_size\" ]]; then
                echo \"ERROR: size not found for \${chrom} — skipping.\" >&2
                exit 1
            fi

            gcbinned_chr=\"gcbinned_\${chrom}.bed\"
            awk -v chr=\"\${chrom}\" '\$1==chr' \"\${gcbinned}\" > \"\${gcbinned_chr}\"
            if [[ ! -s \"\${gcbinned_chr}\" ]]; then
                echo \"WARNING: No GC bins for \${chrom} — skipping.\" >&2
                exit 1
            fi

            generateEV \"\${hicpath}\" \"\${chrom}\" \"\${eigenres}\" \"\${genesfile}\" \"\${gcbinned_chr}\" \"\${chrom_size}\"
        "
    done

    wait

    # Check all chromosomes produced an EV file
    echo "Checking EV generation completeness..."
    missing_evs=0
    for chrom in $(cut -f1 "$sizefile"); do
        if [[ ! -s "best_Eigen_${chrom}_${eigenres}.bedgraph" ]]; then
            echo "NOTE: No EV generated for chromosome $chrom (no data or skipped)." >&2
            (( missing_evs++ ))
        fi
    done
    if (( missing_evs > 0 )); then
        echo "NOTE: $missing_evs chromosome(s) skipped (no Hi-C data). Continuing with remaining chromosomes." >&2
    fi

    concatenate_best_pcs
    process_eigenvector_file "EV_full_genome.bedgraph" "$EVAstates" "$EVBstates"
fi


# ------------------------------------------------------------------------------
# MAIN RESOLUTION LOOP
# ------------------------------------------------------------------------------

countres=0
coarse_res=$maxres
crtmp=$maxres
IFS=',' read -r -a res_array <<< "$res"
maxres_processed=0

for (( i=0; i<${#res_array[@]}; i++ )); do
    myres=${res_array[$i]}
    echo "Current resolution: $myres"
    echo "Array content: ${res_array[@]}"
    high_res=$myres

    genesblock=0
    if [ "$myres" == "$maxres" ] && [ "$maxres_processed" -eq 0 ]; then
        genesblock=1
    fi

    echo "Processing with myres = $myres, genesblock = $genesblock"

    process_resolution $myres $genesblock

    if [ $i -gt 0 ]; then
        reprocess_resolutions_with_shifter $myres $i 1

        finaloutie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`
        Crush_todelete=`echo "Crush_todelete_""$myres""_reprocess"`

        finalpoutie=`echo "$outpre""mergedqvalue_""$myres"".bedgraph"`
        pval_todelete=`echo "pval_todelete_""$myres""_reprocess"`

        cat Crush_Re*.bedgraph | \
            mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4/int($3-$2)}' | \
            mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | \
            mawk '{c1[$1] += $2; c2[$1]++} END {for (i in c1) print i"\t"c1[i]/c2[i]}' | \
            sed 's/:/\t/g' > $finaloutie 2> smallerrors

        if [ $myres -eq $minres ]; then
            cat pvalues_Re*.bedgraph | \
                mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | \
                mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | \
                awk '{c1[$1] += log($2)/log(10); c2[$1]++} END {for (i in c1) print i"\t"10**(c1[i]/c2[i])}' | \
                sed 's/:/\t/g' | sort -k 1,1 -V -k 2bn,2b --stable > $finalpoutie 2> smallerrors
        fi

        echo "Running BH correction"
        BHcorrection $finaloutie
        wait
        mv bhcorrected $finalpoutie
        wait

        GIave=`cat $finaloutie | mawk '{if ($4 < 0) sum+=($4*-1); else sum+=($4)} END {print sum/NR}'`

        if [ $(bc <<< "$qthresh > 0") -eq 1 ]; then
            finaloutiefilt=`echo "$outpre""mergedCrush_""$myres""_qfiltered_reprocess.bedgraph"`
            mawk -v fdr=$qthresh -v var=$GIave \
                'NR==FNR {a[$1":"$2":"$3] = $4; next} {if (a[$1":"$2":"$3] <= fdr) print $1"\t"$2"\t"$3"\t"$4/(var/100)}' \
                $finalpoutie $finaloutie > $finaloutiefilt
        fi

        if [ $trackline -eq 0 ]; then
            cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' > $Crush_todelete
            cat $finalpoutie | mawk '{print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
        else
            cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' | \
                mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$0; else print $0}' > $Crush_todelete
            cat $finalpoutie | \
                mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$1"\t"$2"\t"$3"\t"$4; else print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
            wait
        fi

        if [ $trackline -gt 0 ] && [ $(bc <<< "$qthresh > 0") -eq 1 ]; then
            filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`
            cat $finaloutiefilt | \
                mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20,20 autoScale off\n"$0; else print $0}' > $filt_todelete

            wait

            mv $filt_todelete $finaloutiefilt
            mv $Crush_todelete $finaloutie
            mv $pval_todelete $finalpoutie

            wait
        fi

        sizeBed=`echo "Chrsizes.bed"`
        finaloutie2=`echo "mergedCrush2_""$myres"".bedgraph"`

        cat $sizefile | awk '{print $1"\t""1""\t"$2}' > $sizeBed

        cat $finaloutie | grep -v track | \
            intersectBed -wa -a stdin -wb -b $sizeBed | \
            awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | \
            sort -k 1,1 -V -k 2bn,2b --stable | \
            mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' > $finaloutie2
        mv $finaloutie2 $finaloutie

        if [ "$smoothing" -gt 0 ]; then
            echo "Smoothing Function is turned on"
            smoother="smoother.bedgraph"
            run_smoothing $newoutie $smoother
            rm $smoother tmp2_shiftedleft tmp1_shiftedright
        fi

        if [ $(bc <<< "$qthresh > 0") -eq 1 ] && [ $pcalculation -eq 1 ] && [ $myres -eq $minres ]; then
            filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`

            cat $finaloutiefilt | grep -v track | \
                intersectBed -wa -a stdin -wb -b $sizeBed | \
                awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | \
                mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' > $filt_todelete

            mv $filt_todelete $finaloutiefilt
            mv $finalpoutie ../
        fi

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
for totres in $(echo $res | sed "s/,/\n/g" | sort -k 1bn,1b --stable | sed "s/\n/ /g"); do
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