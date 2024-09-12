#!/bin/bash
version=0.9 #A version where recursive iteration of resolutions is implemented

# Initial variables
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
keeptracks=0
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
coarsestres=2500000
reshift=1
#creates and chooses a random spinner 

function shutdown() {
  tput cnorm # reset cursor
}
trap shutdown EXIT

function cursorBack() {
  echo -en "\033[$1D"
}

function spinner() {
  # make sure we use non-unicode character type locale 
  # (that way it works for any locale as long as the font supports the characters)
  local LC_CTYPE=C

  local pid=$1 # Process Id of the previous running command

  case $(($RANDOM % 14)) in
  0)
    local spin='⠁⠂⠄⡀⢀⠠⠐⠈'
    local charwidth=3
    ;;
  1)
    local spin='-\|/'
    local charwidth=1
    ;;
  2)
    local spin="▁▂▃▄▅▆▇█▇▆▅▄▃▂▁"
    local charwidth=3
    ;;
  3)
    local spin="▉▊▋▌▍▎▏▎▍▌▋▊▉"
    local charwidth=3
    ;;
  4)
    local spin='←↖↑↗→↘↓↙'
    local charwidth=3
    ;;
  5)
    local spin='▖▘▝▗'
    local charwidth=3
    ;;
  6)
    local spin='┤┘┴└├┌┬┐'
    local charwidth=3
    ;;
  7)
    local spin='◢◣◤◥'
    local charwidth=3
    ;;
  8)
    local spin='◰◳◲◱'
    local charwidth=3
    ;;
  9)
    local spin='◴◷◶◵'
    local charwidth=3
    ;;
  10)
    local spin='◐◓◑◒'
    local charwidth=3
    ;;
  11)
    local spin='⣾⣽⣻⢿⡿⣟⣯⣷'
    local charwidth=3
    ;;
  12)
    local spin='CCCRRRUUUSSSHHH'
    local charwidth=1
    ;;
  13)
    local spin='RRROOOWWWLLLEEEYYYLLLAAABBB'
    local charwidth=1
    ;;
  esac

  local i=0
  tput civis # cursor invisible
  while kill -0 $pid 2>/dev/null; do
    local i=$(((i + $charwidth) % ${#spin}))
    printf "%s" "${spin:$i:$charwidth}"

    cursorBack 1
    sleep .1
  done
  tput cnorm
  wait $pid # capture exit code
  return $?
}


function usage {
    echo -e "\n\nusage : crush -i HIC  -g SIZEFILE -a ABED -b BBED | FASTA -r FINERESOLUTION [-cpu CPU] [-w WINDOW] [-h]"
    echo -e "Use option -h|--help for more information"
}


# Help menu function
function help {
  usage;
  echo

  echo "CRUSH will create a temprory folder in your current directory, so make sure you have write access to your current directory."
  echo "See https://github.com/JRowleyLab/CRUSH/ for full documentation on CRUSH"
  echo "-----------------------------------------"
  echo "OPTIONS:"
  
  echo "-h|--help                :  Display this help menu"
  echo " "
  echo "--------------------------REQUIRED PARAMETERS------------------"
  echo "-i|--hic                     :  Input .hic/.cooler/.mcool file by specifying the path. (e.g., '/path/to/ file.hic or file.mcool or file.cooler)"
  echo "-g|--genomesize              :  Specify path to a chromosome size file with two columns corresponding to chromosome and size respectively."
  echo "-a|--initialA                :  Specify path to a bed file with the regions for initializing A. For example, gene annotations e.g. hg19genes.bed."
  echo "-b|--initialB                :  Specify path to either a fasta file or to a bed file for initializing B. If you specify a fasta file, we will calculate intialB from gc content." 
  echo "-r|--res                     :  Resolution desired."
    echo "-e|--eigenfile           :  Specify path to an eigenfile to initialize A and B states. If you don't specify a file, CRUSH will calculate EigenVectors by default and pick the best PC possible." 

  echo "---------------" 

  echo " "
  echo "--------------------------OPTIONAL PARAMETERS------------------"
  
  echo "-o|--outpre              :  Set this if you want to specify a prefix for the output files"
  echo "-c|--cpu                 :  Set the value for cpu number of threads to use. Default is 1."
  echo "-n|--no-merge            :  Set this option to 1 to keep each resolution as a separate output file. Default is to merge in a way that provides maximum resolution."
  echo "-A|--adjustment          :  Set this option to 1 to include a adjustment of CRUSH values at the end. This adjustment shifts values based on any internal skewing of the data. Do not set this if using CRUSH to compare between two Hi-C maps."
  echo "-d|--distance            :  Using this option will filter out the distance next to the diagonal. Default is 0 which considers everything."
  echo "-u|--upperlim            :  The upperlimit of the distance away from the diagonal to consider. Default is 0 so that it considers the whole chromosome."
  echo "-t|--trackline           :  Set to 0 if you want to disable printing a bedgraph trackline header."
  echo "-T|--threshold           :  Distance normalized threshold to filter out extreme outliers."
  echo "-k|--keeptracks          :  Set to 1 in order to keep separate A and B tracks for the probability of interacting with bed file vs other regions."
  echo "-l|--lowerthresh         :  Set the value for lowerthresh."
  echo "-w|--window              :  Set to perform a sliding window average of the scores at individual resolutions. Default is to calculate the appropriate window based on sequencing depth. Set to 1 to remove sliding window."
  echo "-v|--verbose             :  Set to 1 to enable verbose mode showing extensive messages."
  echo "-S|--switch              :  Set this option to 0 for bypassing re-initialization. Default is 1."  
  echo "-C|--cleanup             :  Set the value for cleanup"
  echo "-E|--endZ                :  Set the value for endZ"
  echo "-x|--exclbed             :  Set the value for exclbed"
  echo "-q|--qvalue              :  Set the qvalue threshold. default 0.05. Set to 0 to not perform qvalue filtering. The qvalues will be reported as a separate track regardelss."
  echo "-u|--use                 :  Whether to use of overwrite existing GI tracks previously calculated at individual resolutions. Set this option to u to use previous calculations. Default is to recalculate. This option is useful for merging resolutions."     
  echo "-f|--tmpfolder           :  Set this if you want to name the temporary folder youself. Make sure it doesn't already exist in your current working directory. Default is to name it CRUSHtmp with a randomnumber."
  echo "-m|--maxres              :  Set this to the coarsest resolution you want to consider. Default is to check every resolution between 1000000 and your desired resolution to inform each other." 
  echo "-R|--reshift             :  Set this to option to 0 for bypassing end-shifting. Default is 1." 

  # Add more options here if needed
}


# Parsing command line arguments
while test $# -gt 0
do
        case "$1" in
                -h| --help)
                help
                exit 0
                ;;
                -i|--hic)
                hicpath=`readlink -e $2`
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
		-k|--keeptracks)
		keeptracks=$2
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
        -R|--reshift)
        reshift=$2
        ;;
        esac
        shift
done

# Checking if arguments are provided

if [ $hicpath == 0 ]
then
echo "You must specify a .hic file!......""$menu"
exit 0
fi
if [ $res == 0 ]
then
echo "You must specify a resolution!......""$menu"
exit 0
fi
if [ $switch == 0 ]
then
echo "Warning: Re-evaluation & Re-interation of A and B Bins has been deactivated. We recommend leaving this parameter alone to achieve highest possible confidence and resolution."
fi

# Check if the size file is provided
if [ -z "$sizefile" ] || [ -z "$genesfile" ] || [ -z "$fastafile" ]; then
    echo "Missing one or more required arguments: --genomesize, --initialA, or --initialB!..."
    help
    exit 1
fi

if [ $doNotMerge == 1 ]
then
echo "Warning: doNotMerge set. We recommend merging to achieve highest possible confidence and resolution. Only set doNotMerge if you plan to merge the individual resolutions later."
fi

if [ $adjustment == 1 ]
then
echo "Warning: adjusting end compartmental values which may result in some loss of quantitative power. We recommend adjustment only for troubleshooting or when not comparing two samples."
fi

echo "CRUSH_v""$version"" --hic ""$hicpath"" --res ""$res"" --genomesize ""$sizefile"" --initialA ""$genesfile"" --initialB ""$fastafile"" --cpu ""$cpu"" --no-merge ""$doNotMerge"" --adjustment ""$adjustment"" --distance ""$distance"" --upperlim ""$upperlim"" --lowerthresh ""$lowerthresh"" --trackline ""$trackline"" --threshold ""$threshold"" --keeptracks ""$keeptracks"" --window ""$window"" --switch ""$switch"" --maxres ""$coarsestres" > "$outpre"CRUSHparamters.txt

# Check for .hic or .mcool
juiceorcool=`echo "$hicpath" | sed 's/\./\t/g' | sed 's/\./\t/g' | awk '{if ($NF == "hic") print 0; else if ($NF == "mcool") print 1; else print 2}'`

if [ $juiceorcool == 2 ]
then
echo "You must choose a Hi-C file that ends in either .hic (juicer) or .mcool (cooler)"
exit 0
elif [ $juiceorcool == 0 ]
then
echo "Identified .hic juicer format. We will use the pythonic hicstraw to extract the data and assume you have it installed, using pip install hic-straw."
else
echo "Identified .mcool format. We will use cooler to extract the data and assume you have it installed in your path."
fi

###################################################################
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

process_oppocheck_statement() {
  local innie_genes="$1"
  local innie_gc="$2"
  local outie_oppo="$3"
  local outie_pvalues="$4"
  local chr="$5"
  local res="$6"
  local input_dir="$7"
  local output_dir="$8"

  # Paths for Zfile and RowZ
  local zfile="${input_dir}/zerofile"
  local rowZ="${input_dir}/RowZs"

  # Debugging output to check paths
  echo "Zfile path: $zfile"
  echo "RowZ path: $rowZ"
  echo "Output directory: $output_dir"

  # Check if files exist
  if [ ! -f "$zfile" ] || [ ! -f "$rowZ" ]; then
    echo "Required file missing: $zfile or $rowZ"
    return 1  # Exit if files are missing
  fi

numbinsA=`wc -l $innie_genes | awk '{print $1}'`
numbinsB=`wc -l $innie_gc | awk '{print $1}'`

mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) print $2"\t"$3; next} {print $0}' $innie_genes $rowZ $zfile | awk '{b[$1]+=$2;c[$1]++;d[$1]+=($2**2)} END { for (i in b) { print i"\t"(b[i])*(c[i])"\t"c[i]"\t"b[i]/c[$1]"\t"((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } } ' > "${output_dir}/A2removedfile" &
mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) print $2"\t"$3; next} {print $0}' $innie_gc $rowZ $zfile | awk '{b[$1]+=$2;c[$1]++;d[$1]+=($2**2)} END { for (i in b) { print i"\t"(b[i])*(c[i])"\t"c[i]"\t"b[i]/c[$1]"\t"((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } } ' > "${output_dir}/B2removedfile" &

wait

myABsum=`cat "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | mawk '{sum+=$3} END {print sum}'`

wait

if [ "$verbose" -gt 0 ]; then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "It took ""$seconds"" seconds. Comparing A vs B scores for ""$mychr"" chromosome."
start=$(date +%s)
fi

# Calculating the output
mawk -v var=$res -v myABsum=$myABsum '{if ($1 in b) b[$1]-=$2; else b[$1]=$2} END { for (i in b) {  if (b[i] != 0) print i"\t"(b[i])/myABsum } } ' "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | sort -k 1bn,1b --stable | mawk -v mchr=$chr -v myres=$res -v window=$newwindow -v OFS="\t" 'BEGIN{slide=1} {mod=NR%window; if(NR<=window){print mchr, int($1), $2; count++}else{sum-=array[mod]}sum+=$2;array[mod]=$2;} (NR%slide)==0{print mchr,int($1-((myres*window)/2)),sum/count}' | mawk -v myres=$res '{print $1"\t"$2-int(myres/2)+myres"\t"$2+(int(myres/2))+myres"\t"$3}' | mawk '{if ($2 > 0) print $0}' > $outie_oppo 2> "${output_dir}/errorfile"

# Calculating pvalues
awk 'NR==FNR { c[$1] = $3; m[$1] = $4; sd[$1]=$5; next} {if (($1 in m) && ($3 != 1) && (c[$1] != 1)) print $1"\t"(m[$1]-$4)/((((sd[$1]**2)/c[$1]) + (($5**2)/$3))**.5)"\t"(c[$1] + $3 -2) }' "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" > "${output_dir}/ttable"

zfileforpy="${output_dir}/ttable"
pfileforpy="${output_dir}/ptable"

cat << EOF > zlookup.py
import scipy.stats
inputfile = "$zfileforpy"
outfile = "$pfileforpy"
mybins=[]
mypvals=[]
with open(inputfile, 'r') as innie:
        for line in innie:
                sline = line.strip()
                lines = sline.split("\t")
                mybins.append(str(lines[0]))
                mytstat = float(lines[1])
                mydf = int(lines[2])
                mypvals.append(scipy.stats.norm.sf(abs(mytstat)))
with open(outfile, 'w') as outie:
        for w in range(0,len(mybins)):
                outie.write(str(mybins[w]) + "\t" + str(mypvals[w]) + "\n")
EOF
python zlookup.py

mawk -v mchr=$chr -v myres=$res '{print mchr"\t"$1"\t"$1+myres"\t"$2}' "${output_dir}/ptable" > $outie_pvalues
wait

#If exclbed file is given as an input, and the option is turned on, Crush works through this
if [ "$exclbed" != 0 ]
then
exclbins="${output_dir}/exclbins"
excloutie="${output_dir}/exclout_${res}_${chr}"

#Processing the exclusion regions
cat $exclbed | mawk -v myres=$res -v mychr=$chr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$res '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $exclbins

wait

#Excluding regions from the output 
mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $0}' $exclbins $outie_oppo > $excloutie

wait

mv $excloutie $outie_oppo
fi

if [ $keeptracks -eq 1 ]
then
rawcomp="${output_dir}/rawcompvectors"

if [ $reprocess == 1 ]
then 
cat "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | sort -k 1bn,1b --stable | groupBy -i stdin -g 1 -c 2 -o collapse | sed 's/,/\t/g' | awk -v var=$prev_res -v mchr=$mychr '{print mchr"\t"$1"\t"$1+myres"\t"$2"\t"$3}' | awk '{print $0"\t"$4"\t"$5}' > $rawcomp
cat $rawcomp | cut -f 1-3,6 > $Aoutie_reprocess
cat $rawcomp | cut -f 1-3,7 > $Boutie_reprocess

else
cat "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | sort -k 1bn,1b --stable | groupBy -i stdin -g 1 -c 2 -o collapse | sed 's/,/\t/g' | awk -v var=$myres -v mchr=$mychr '{print mchr"\t"$1"\t"$1+myres"\t"$2"\t"$3}' | awk '{print $0"\t"$4"\t"$5}' > $rawcomp
cat $rawcomp | cut -f 1-3,6 > $Aoutie
cat $rawcomp | cut -f 1-3,7 > $Boutie

fi

fi

if [[ "$oppocheck" -gt 0 && "$adjustment" -gt 0 ]]
then

amed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1) print $4}' $innie_genes $outie_oppo | mawk '{sum +=$1} END {print sum/NR}'`
bmed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $4}' $innie_genes $outie_oppo | mawk '{sum += $1} END {print sum/NR}'`

wait

#Changed this from outie to outie_oppo
#Here, Crush is taking the amed and bmed values and doing the operations such as multiplying varB with -1 and then substracting these values from varA
cat $outie_oppo | mawk -v varA=$amed -v varB=$bmed '{print $1"\t"$2"\t"$3"\t"$4-(varA-(varB*-1))}' > $outie2

wait

mv $outie2 $outie_oppo
fi

if [ "$res" -ge "$endZ" ]
then

#Changed this from outie to outie_oppo
myendmean=`cat $outie_oppo | mawk '{sum += $4; cnum++} END {print sum/cnum}'`

echo "$myendmean"

tmpendZ="${output_dir}/tmpendnorm"

cat $outie_oppo | mawk -v mymean=$myendmean '{print $1"\t"$2"\t"$3"\t"($4-mymean)}' > $tmpendZ

wait

mv $tmpendZ $outie_oppo
fi

echo "You got this working till here"

}

# Calculations start from here------

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

if [ $isfasta == "1" ]
then
echo "Fasta-file: $fastafile"
else
echo "initial-B: $fastafile"
fi

resbins=`echo "size_bins"".bed"`
gcfile=`echo "gc_bins"".txt"`
gc_g_ga=`echo "gc_g_ga_bins"".bed"`
Bbins=`echo "Bbins"".bed"`

if [ $minres -ge 500 ]
then
cat $sizefile | awk -v myres=500 '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
echo "$resbins calculated at 500bp resolution"

else

cat $sizefile | awk -v myres=$minres '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
echo "$resbins calculated at ""$minres"" resolution"
fi

wait

if [ $isfasta == "1" ]
then
echo "Now, calculating the gc content:"
bedtools nuc -fi $fastafile -bed $resbins > $gcfile
wait

# Created a %G/%G+%A file in 7th column
cat $gcfile | grep -v user | cut -f 1-5 | awk '{print $0"\t"($4+$5)}' | awk '{if ($6 > 0) print $0}' | awk '{print $0"\t"($5/$6)}' > $gc_g_ga
wait

# Calculating mean and SD for the whole genome and then subtracting 2SD from the mean
gc_thresh=`cat $gc_g_ga | awk '{s+=$7; ss+=$7*$7; linecount+=1} END{print m=s/linecount, sqrt(ss/linecount-m^2)}' | awk '{print $0}' | awk '{print $1 - 1*$2}'`

if [ $verbose -gt 0 ]
then
echo "gc threshold is ""$gc_thresh"
fi

# Generating the bins by considering the bins below
cat $gc_g_ga | awk -v thresh=$gc_thresh '{if ($7 < thresh) print $0}' > $Bbins
wait
echo "Finished generating Bbins file for all chromosomes "

else
cat $fastafile > $Bbins
fi

if [ $endZ == 0 ]
then
endZ=`echo "$minres"`
fi

# Here is where I am running Eigen Block after reslist

EVAstates="EVAstates.bed"
EVBstates="EVBstates.bed"

if [ -n "$eigenfile" ]; then
    echo "Using user-provided Eigen Vector file: $eigenfile"
    
    # Correct the awk commands to ensure they process the file correctly
    awk '{if ($4 > 0) print $0}' "$eigenfile" > "$EVAstates"
    if [ $? -ne 0 ]; then
        echo "Error creating EVAstates from eigenfile. Please check the awk command."
        exit 1
    fi

    awk '{if ($4 < 0) print $0}' "$eigenfile" > "$EVBstates"
    if [ $? -ne 0 ]; then
        echo "Error creating EVBstates from eigenfile. Please check the awk command."
        exit 1
    fi

    # Check if files are created
    if [ ! -s "$EVAstates" ]; then
        echo "Error: EVAstates is empty or not created."
        exit 1
    fi

    if [ ! -s "$EVBstates" ]; then
        echo "Error: EVBstates is empty or not created."
        exit 1
    fi

    echo "Eigen Vector files created successfully: $EVAstates, $EVBstates"
else

    echo "No user-provided Eigen Vector file. Generating Eigen Vectors using run_EigenVector script."

    # Default behavior when no eigenfile is provided
    for chrom in $(cut -f 1 $sizefile); do
        
        genesbinned="genesbinned.bed"
        gcbinned="gcbinned.bed"

        cat $genesfile | mawk -v myres=100000 -v mychr=$chrom '{if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres}' > $genesbinned
        cat $Bbins | mawk -v myres=100000 -v mychr=$chrom '{if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres}' > $gcbinned

        generateEV "$hicpath" "$chrom" 100000 "$genesbinned" "$gcbinned" 
    
    done

    # Concatenate best eigenvectors for full genome BEDGraph
    concatenate_best_pcs
fi


#Shifter function
run_shifter() {
    local finer_res=$1
    local coarser_res=$2
    local high_innie=$3
    local coarse_innie=$4
    local shifter_outie=$5
    
    echo "Running shifter with Finer Resolution: $finer_res and Coarser Resolution: $coarser_res"
    echo "High-res file: $high_innie, Coarse-res file: $coarse_innie"

    # Generate intervals
    awk -v var="$coarser_res" '
        {
            for (i=0; i<=$2; i+=var) {   # Loop from 0 to the size of the chromosome in steps of "var"
                start = i-(var*2)        # Calculate the start of the interval (2 bins to the left)
                end = i+(var*3)          # Calculate the end of the interval (3 bins to the right)
                if (start < 0) start = 0  # Ensure the start is not negative
                if (end > $2) end = $2   # Ensure the end does not exceed the chromosome length
                print $1 "\t" start "\t" end    # Print the chromosome name, start, and end of the interval
            }
        }
    ' "$sizefile" > temp_intervals.bed
    
    # Main pipeline (1. Intersect with high_innie and calculate mean 2. Intersect with coarse_innie and calculate first difference 3. Intersect with high_innie again and calculate final difference)
    intersectBed -wa -a temp_intervals.bed -wb -b $high_innie | groupBy -i stdin -g 1,2,3 -c 7 -o mean | intersectBed -wa -a $coarse_innie -wb -b stdin | awk '{print $1"\t"$2"\t"$3"\t"$8-$4}' | intersectBed -wa -a $high_innie -wb -b stdin | groupBy -i stdin -g 1,2,3,4 -c 8 -o mean | awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' > $shifter_outie

    #Cleanup
    rm temp_intervals.bed

    echo "Shifter executed and output written to $shifter_outie"
}

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
        #cat $genesfile | mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed
        cat $EVAstates | mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        # Bbins
        #cat $Bbins |  mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB
        cat $EVBstates |  mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

    elif [ "$switch" -lt 1 ]; then

        # Processing the genes & Bbins files
        #cat $genesfile | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed
        cat $EVAstates | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        # Bbins
        #cat $Bbins |  mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB
        cat $EVBstates | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

    else  

        echo "Using combined bins for $myres"

        # Processing the genes & Bbins files for other resolutions 
        cat $combined_newAbins | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        # Bbins
        cat $combined_newBbins |  mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

    fi

}

#Separating out A and B regions from the outiefull and intersecting with the genes and Bbins to improve compartments identification
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

if [ "$switch" -eq 1 ]
then

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

#Benjamini-Hochberg Function
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

# Define arrays to store filenames
declare -a maxres_processed=0

process_resolution() {
    local myres=$1 # Pass the resolution as a parameter
    local genesblock=$2

    if [ $window == 0 ]
    then
        # Changed the below line from mawk to awk as mawk was throwing an error while working through script line by line
        totbins=`cat $sizefile | awk -v var=$myres '{print $1"\t"($2/var)**2}' | awk '{sum += $2} END {print sum}'`

        newwindow=1
    else
        newwindow=`echo "$window"`

        # Calculating window sizes for the resolutions 
        newwindow=$(mawk -v win=$window -v res=$myres 'BEGIN{print int(win/(res/1000))+1}')
    fi

    outiefull=`echo "Crush_""$myres"".bedgraph"`
    outiefullPval=`echo "pvalues_""$myres"".bedgraph"`
    Aoutiefull=`echo "ACrush_""$myres"".bedgraph"`
    Boutiefull=`echo "BCrush_""$myres"".bedgraph"`

    echo -e "\nNow dumping reads and calculating initial CRUSH for all chromosomes at ""$myres"

    task(){

    mychr=`cat $sizefile | head -n $myiter | tail -1 | cut -f 1`
    mysize=`cat $sizefile | head -n $myiter | tail -1 | cut -f 2`
    rowbins=`cat $sizefile | head -n $myiter | tail -1 | mawk -v myres=$myres '{print int($2/myres)+1}'`
    outie=`echo "Crush_""$myres""_""$mychr""_tmp"`
    pscores=`echo "ttest_""$myres""_""$mychr""_tmp"`
    outie2=`echo "Crush_""$myres""_""$mychr""_tmp2"`
    tmpfiles=`echo "ABtmpfiles_""$mychr""_""$myres"`
    Aoutie=`echo "ACrush_""$myres""_""$mychr""_tmp"`
    Boutie=`echo "BCrush_""$myres""_""$mychr""_tmp"`

    #Atmp=`echo "$tmpfiles""/Atmp"`
    #Btmp=`echo "$tmpfiles""/Btmp"`

    #Atmppseudo=`echo "$tmpfiles""/Atmpp_""$mychr""_""$myres"`
    #Btmppseudo=`echo "$tmpfiles""/Btmpp_""$mychr""_""$myres"`

    #Afile=`echo "tmpA"`
    #Bfile=`echo "tmpB"`
    #Bpfile=`echo "tmpBp"`
    #Apfile=`echo "tmpAp"`

    #A2file=`echo "tmpA2"`
    #B2file=`echo "tmpB2"`
    #Ap2file=`echo "tmpAp2"`
    #Bp2file=`echo "tmpBp2"`

    [ ! -d "$tmpfiles" ] && mkdir "$tmpfiles"

    efile=`echo "expecteds"`
    ofile=`echo "observeds"`
    mydist=$((distance + 0))
    myupper=$((upperlim + 0))

    if [ "$verbose" -gt 0 ]
    then
        echo "Getting genic bins."

    fi
    
    sortedbed=`echo "$tmpfiles""/genebins"`
    pseudoB=`echo "$tmpfiles""/pseudoB"`
    dumped=`echo "$tmpfiles""/dumped_""$myres""_""$mychr""_tmp"`

    #Generating Intializing states
    #Adding in the code for re-initializing Bbins for next following resolutions
    process_genes_Bbins $myres $mychr $genesblock $sortedbed $pseudoB 

    if [ "$verbose" -gt 0 ]
    then
        echo "Dumping ""$mychr"" chromosomal reads."
    fi

        newofile=`echo "$tmpfiles""/observeds"`
        newefile=`echo "$tmpfiles""/expecteds"`
        errofile=`echo "errorlog"`

    # Conditional argument for straw vs cooler usage
    dumperrors=`echo "dumperrors_""$myres"`

if [ $juiceorcool -eq 0 ]
then

# Dumping ofile and efile using no normalization schemes (Uncommented from previous version) using straw
cat << EOF > dumper.py
import hicstraw
result = hicstraw.straw('observed', 'NONE', "$hicpath", "$mychr", "$mychr", "BP", int("$myres"))
for i in range(len(result)):
	print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts)) 

EOF

python dumper.py | grep -v WARN > $dumped 2> $dumperrors

else
echo "Dumping using cooler..."

# Dumping reads for cooler files
cooler dump $hicpath::/resolutions/$myres -r $mychr --join | cut -f 2,5,7 > $dumped 2> $dumperrors

fi

chromcheck=`cat dumperrors_$myres | grep -e KeyError -e name | wc -l`

if [ $chromcheck -gt 0 ]
then
echo "WARNING: A chromosome in your size file was not found in the .mcool file. Make sure the names match up. For example, check for chr1 vs. 1."
fi

wait

if [ "$upperlim" -gt 0 ]
then
cat $dumped | awk -v var=$myres -v distfilt=$mydist -v upper=$upperlim '{if (($2-$1 > distfilt) && ($2-$1 <= upper)) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' | awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile } ' &
else
cat $dumped | awk -v var=$myres -v distfilt=$mydist '{if ($2-$1 > distfilt) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' | awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile } ' &
fi

wait 

if [ "$verbose" -gt 0 ]
then
echo "Done dumping."
fi

if [ "$verbose" -gt 0 ]
then
echo "Now getting distance normalized values for ""$mychr"" chromosome."
start=$(date +%s)
fi


errorfile=`echo "errorlog"`
tmpmean=`echo "$tmpfiles""/tmpmean"`
symmfile=`echo "$tmpfiles""/symmonfile"`
#onfile=`echo "$tmpfiles""/distnorms"`

# Calculating rows-Mean using expected and observed files
mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 != $3) print $2"\t"$3"\t"($4+1)/(c1[$1]+1)}' $newefile $newofile | mawk -v thresh=$threshold '{if ($3 < thresh) print $1"\t"$2"\t"$3"\n"$2"\t"$1"\t"$3; else print $1"\t"$2"\t"thresh"\n"$2"\t"$1"\t"thresh}' | awk -v symmfile=$symmfile -v tmpmean=$tmpmean -v rowbins=$rowbins '{m[$1]+=$3;d[$1]+=($3**2); print $0 >> symmfile} END {for (i in m) print i"\t"m[i]/rowbins"\t"((d[i]/rowbins-((m[i]/rowbins)**2))**.5) >> tmpmean }'  2> $tmpfiles/$errorfile

wait

#Generating on file to revive keep tracks. Can be commented if not needed
#awk 'NR % 2 == 1' $symmfile > $onfile

if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)

start=$(date +%s)
fi

if [ "$myres" -gt 50000 ]
then
oppocheck=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) ; else if ($2 in c1); else print $0}' $sortedbed $symmfile | wc -l`
else
oppocheck=2
fi

# Perform real calculation
if [ "$oppocheck" -lt 1 ] && [ "$verbose" -gt 0 ]
then
echo "WARNING: All bins on chromosome ""$mychr"" overlap with your bed file at ""$myres"" resolution"". This should not affect the scores if using multiple resolutions. If using a single resolution, then maybe not use this resolution?."
fi


if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
start=$(date +%s)
fi

RowZs=`echo "RowZs"`
myzscore=`echo "Zscorefile"`

# Calculating Z-scores
mawk 'NR==FNR { m[$1]=$2; d[$1]=$3; next} {if (d[$2] > 0) print $1"\t"$2"\t"($3-m[$2])/(d[$2])}' $tmpmean $symmfile > $tmpfiles/$RowZs

#Creating a zero's file 
zfile=`echo "zerofile"`
cat $sizefile | mawk -v var=$mychr '{if ($1 == var) print $0}' | mawk -v wlkr=$myres '{for (i=0;i<=$2;i+=wlkr) print i"\t"i"\t0.1"}' > $tmpfiles/$zfile

A2removedfile=`echo "A2removed"`
B2removedfile=`echo "B2removed"`
ttable=`echo "ttable"`
ptable=`echo "ptable"`

#Running the Process_oppocheck function
process_oppocheck_statement $sortedbed $pseudoB $outie $pscores $mychr $myres "ABtmpfiles_${mychr}_${myres}" "ABtmpfiles_${mychr}_${myres}"

}

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    ( "$@"; )
    
    printf '%.3d' $? >&3
    )& spinner $!
}

open_sem $cpu 

if [ "$whichchoice" == "o" ]
then

for (( myiter=1; myiter<=$totchroms; myiter++ ))
do
    run_with_lock task $myiter
done

wait

# Merging CRUSH files

#Adding this to check python output
echo "Resolution: "$myres" for Resolution output"

echo -ne "Merging individual chromosome files\033[0K\r"

cat Crush*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull

cat ttest_*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefullPval

wait 

rm ttest_*_tmp

wait

#Setting up for shifter

if [ "$myres" == "$maxres" ]
then
coarse_infile=`echo "GI_""$myres"".bedgraph"`
cat $outiefull > $coarse_infile

wait
fi

res_values=`echo "res_values"`
echo "$myres" >> $res_values

# Re-evaluating both A and B bins for re-initialization for next resolution

# Re-evaluating both A and B bins for re-initialization for next resolution
    if [ "$countres" -gt 0 ]; then
        python_output="shifter.bedgraph"

        high_res_infile="$outiefull"
        coarse_res_infile="$coarse_infile"

        # Executing the Shifter Python Script
        run_shifter $high_res $coarse_res $high_res_infile $coarse_res_infile $python_output
        
        cat $python_output > $high_res_infile
        
        echo "Resolution: $myres for shifter output"
        separating_ABbins $high_res_infile

        coarse_res=$high_res
        coarse_infile=$high_res_infile
    else
        separating_ABbins $outiefull
    fi

fi

if [ $keeptracks -eq 1 ] && [ "$whichchoice" == "o" ]
then
cat ACrush*_tmp | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $Aoutiefull &
cat BCrush*_tmp | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $Boutiefull &
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

        local outiefull_reprocess=`echo "Crush_reprocess_""$prev_res"".bedgraph"`
        local outiefullPval_reprocess=`echo "pvalues_reprocess_""$prev_res"".bedgraph"`

        # Loop through each chromosome
        local chromosomes=($(cut -f1 $sizefile))

        for mychr in "${chromosomes[@]}"; do
            local sortedbed_reprocess="genebins_reprocess_${mychr}_${prev_res}"
            local pseudoB_reprocess="pseudoB_reprocess_${mychr}_${prev_res}"
            local outie_reprocess="Crush_reprocess_${prev_res}_${mychr}_tmp"
            local pscores_reprocess="ttest_reprocess_${prev_res}_${mychr}_tmp"
            local Aoutie_reprocess="ACrush_reprocess_${prev_res}_${mychr}_tmp"
            local Boutie_reprocess="BCrush_reprocess_${prev_res}_${mychr}_tmp"

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
            process_oppocheck_statement $sortedbed_reprocess $pseudoB_reprocess $outie_reprocess $pscores_reprocess $mychr $prev_res "ABtmpfiles_${mychr}_${prev_res}" "ABreprocesstmpfiles_${mychr}_${prev_res}"

            wait 
        done    

        # Merging CRUSH files

        echo -ne "Merging individual chromosome files\033[0K\r"
        cat Crush_reprocess*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull_reprocess

        #Adding this to check python output
        echo "Resolution: "$prev_res" for Resolution output"

        cat ttest_*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefullPval_reprocess

        #Adding this to check python output
        echo "Resolution: "$prev_res" for Resolution p-val output"

        wait 

        if [ "$prev_res" == "$maxres" ]; then

            coarse_infile_reprocess=`echo "GI_reprocess_""$prev_res"".bedgraph"`
            cat $outiefull_reprocess > $coarse_infile_reprocess

            wait
        fi

        prev_values=`echo "prev_values"`
        echo "$prev_res" >> $prev_values

        # Re-evaluating both A and B bins for re-initialization for next resolution

        echo "$j"

        if [ "$j" -gt 0 ]; then

            highres_infile_reprocess=$outiefull_reprocess
            coarse_res_infile_reprocess=$coarse_infile_reprocess

            # Ensure both highres and coarse resolution files exist
            if [ -f "$highres_infile_reprocess" ] && [ -f "$coarse_res_infile_reprocess" ]; then
                echo "Files exist: $highres_infile_reprocess, $coarse_res_infile_reprocess"
            else
                echo "One or both files do not exist: $highres_infile_reprocess, $coarse_res_infile_reprocess"
            fi

            echo "Running shifter for $prev_res"

            # Run the shifter function
            run_shifter $highres_reprocess $coarse_res_reprocess $highres_infile_reprocess $coarse_res_infile_reprocess $python_output

            cat $python_output > $highres_infile_reprocess

            echo "Shifter run completed for $prev_res"
            echo "Resolution: $prev_res for shifter output"

            # Separate A/B bins using the shifter output
            separating_ABbins $highres_infile_reprocess

            coarse_res_reprocess=$highres_reprocess
            coarse_infile_reprocess=$highres_infile_reprocess
        else
            # Separate A/B bins using the resolution-based output
            separating_ABbins $outiefull_reprocess
        fi
    done
}



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
        finalpoutie=`echo "$outpre""mergedqvalue_""$myres"".bedgraph"`

        Crush_todelete=`echo "Crush_todelete_""$myres""_reprocess"`
        pval_todelete=`echo "pval_todelete_""$myres""_reprocess"`

        cat Crush_reprocess*.bedgraph | mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4/int($3-$2)}' | mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | mawk '{c1[$1] += $2; c2[$1]++} END {for (i in c1) print i"\t"c1[i]/c2[i]}' | sed 's/:/\t/g' > $finaloutie 2> smallerrors
        cat pvalues_reprocess*.bedgraph | mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | awk '{c1[$1] += log($2)/log(10); c2[$1]++} END {for (i in c1) print i"\t"10**(c1[i]/c2[i])}' | sed 's/:/\t/g' | sort -k 1,1 -V -k 2bn,2b --stable > $finalpoutie 2> smallerrors

        echo "Running BH correction"
        BHcorrection $finaloutie
        wait
        mv bhcorrected $finalpoutie
        wait

        GIave=`cat $finaloutie | mawk '{if ($4 < 0) sum+=($4*-1); else sum+=($4)} END {print sum/NR}' `

        if [ $(bc <<< "$qthresh > 0") -eq 1 ]
        then

        finaloutiefilt=`echo "$outpre""mergedCrush_""$myres""_qfiltered_reprocess.bedgraph"`
        filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`

        mawk -v fdr=$qthresh -v var=$GIave 'NR==FNR {a[$1":"$2":"$3] = $4; next} {if (a[$1":"$2":"$3] <= fdr) print $1"\t"$2"\t"$3"\t"$4/(var/100)}'  $finalpoutie $finaloutie> $finaloutiefilt

        fi

        if [ $trackline -eq 0 ]
        then

        cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' > $Crush_todelete
        cat $finalpoutie | mawk '{print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete

        else

        cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$0; else print $0}' > $Crush_todelete
        cat $finalpoutie | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$1"\t"$2"\t"$3"\t"$4; else print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
        wait
        fi

        if [ $trackline -gt 0 ] && [ $(bc <<< "$qthresh > 0") -eq 1 ]
        then

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

        cat $finaloutie | grep -v track | intersectBed -wa -a stdin -wb -b $sizeBed | awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-100:100 autoScale off\n"$0; else print $0}'> $finaloutie2
        mv $finaloutie2 $finaloutie

        if [ $(bc <<< "$qthresh > 0") -eq 1 ]
        then

        cat $finaloutiefilt | grep -v track | intersectBed -wa -a stdin -wb -b $sizeBed | awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-100:100 autoScale off\n"$0; else print $0}'> $filt_todelete
        wait
        mv $filt_todelete $finaloutiefilt
        fi

        mv $finaloutie ../
        mv $finalpoutie ../

        if [ $(bc <<< "$qthresh > 0") -eq 1 ]
        then

        mv $finaloutiefilt ../
        fi

    fi
done


if [ $doNotMerge -eq 1 ]
then
echo "Finished with each resolution listed, but not merging...."
exit 0
fi

countres=0
for totres in $(echo $res | sed "s/,/\n/g" | sort -k 1bn,1b --stable | sed "s/\n/ /g" )
do
echo "$totres" >> GIrestmpfile
countres=$((countres+1))
done

if [ $countres -le 1 ]
then
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


Afinalout_reprocess=`echo "mergedCrush_A_""$minres""_reprocess.bedgraph"`
Bfinalout_reprocess=`echo "mergedCrush_B_""$minres""_reprocess.bedgraph"`

if [ $keeptracks -eq 1 ]
then

cat ACrush_reprocess_*.bedgraph | mawk -v minr=$minres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$minres '{for (i=$2;i<$3;i++) print $1"\t"int(i)*mr"\t"(int(i)+1)*mr"\t"$4}' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable | groupBy -i stdin -g 1,2,3 -c 4 -o sum > $Afinalout_reprocess 2> smallerrors&
cat BCrush_reprocess_*.bedgraph | mawk -v minr=$minres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$minres '{for (i=$2;i<$3;i++) print $1"\t"int(i)*mr"\t"(int(i)+1)*mr"\t"$4}' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable | groupBy -i stdin -g 1,2,3 -c 4 -o sum > $Bfinalout_reprocess 2> smallerrors &

fi
wait

cd ../

if [ "$cleanup" -gt 0 ]
then
echo "cleaning up"
rm -r $crushdir
fi

if [ "$reshift" -gt 0 ]; then

    echo "Going into Final Shifter"

    # Final shifter
    for (( i=0; i<${#res_array[@]}; i++ )); do

        myres=${res_array[$i]}

        if [ $i -eq 0 ]; then
            echo "almost done..."
            oldres=`echo "$myres"`
            oldinnie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`
        else
            newres=`echo "$myres"`
            newinnie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`
            finalshifter=`echo "finalshifter.bedgraph"`
            newoutie=`echo "mergedCrush2_""$myres"".bedgraph"`

            # Run the shifter function
            run_shifter $newres $oldres $newinnie $oldinnie $newoutie

            # Updating the resolution
            oldres=`echo "$myres"`
            oldinnie=`echo "mergedCrush2_""$myres"".bedgraph"`

            wait

            mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-100:100 autoScale off\n"$0; else print $0}' $newoutie > $newinnie    
            wait

        fi

    done
    #rm mergedCrush2_*.bedgraph
fi

echo "Finished! Check the output."
echo -e "Note: Please keep in mind that we are using resolution walking.\nIf the coarsest resolution doesn't match the compartment pattern at that resolution, please consider re-running with the -m parameter set to start the walking with smaller bins."

#Pearson Correlation:
echo "Generating Pearson Correlation values for the bedgraphs."

# Create the PearsonCorrelation.py script
cat << EOF > PearsonCorrelation.py
import pandas as pd
from scipy import stats
import argparse

def average_bins(df, target_resolution):
    df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
    df['End'] = pd.to_numeric(df['End'], errors='coerce')
    df = df.dropna(subset=['Start', 'End'])
    df['Start'] = df['Start'].astype(int)
    df['End'] = df['End'].astype(int)
    df['bin'] = (df['Start'] // target_resolution) * target_resolution
    df_avg = df.groupby(['Chromosome', 'bin']).agg({'Score': 'mean'}).reset_index()
    df_avg['End'] = df_avg['bin'] + target_resolution
    return df_avg[['Chromosome', 'bin', 'End', 'Score']]

def calculate_pearson(control_file, sample_file, resolution):
    cols = ["Chromosome", "Start", "End", "Score"]
    control_df = pd.read_csv(control_file, sep="\t", header=None, names=cols, dtype={'Chromosome': str}, comment='t')
    control_avg = average_bins(control_df, resolution)
    
    sample_df = pd.read_csv(sample_file, sep="\t", header=None, names=cols, dtype={'Chromosome': str}, comment='t')
    sample_avg = average_bins(sample_df, resolution)
    
    control_avg["Start"] = control_avg["bin"].astype(str)
    sample_avg["Start"] = sample_avg["bin"].astype(str)
    control_avg["Chromosome"] = control_avg["Chromosome"].astype(str)
    sample_avg["Chromosome"] = sample_avg["Chromosome"].astype(str)
    control_avg["lame"] = control_avg["Chromosome"] + "_" + control_avg["Start"]
    sample_avg["lame"] = sample_avg["Chromosome"] + "_" + sample_avg["Start"]
    
    merged_df = control_avg.merge(sample_avg, on="lame")
    pearson_corr = stats.pearsonr(merged_df["Score_x"], merged_df["Score_y"])
    
    print(f"{sample_file} Pearson correlation: {pearson_corr[0]}")

def main():
    parser = argparse.ArgumentParser(description='Calculate Pearson correlation between bedgraph files at different resolutions.', add_help=False)
    parser.add_argument("-c", "--control", dest="ctl_file", help="Control eigenvector file. 4 column bedgraph.", required=True)
    parser.add_argument("-s", "--sample", dest="sample_file", help="Sample eigenvector file. 4 column bedgraph.", required=True)
    parser.add_argument("-r", "--resolution", dest="resolution", help="Resolution of the control eigenvector file.", type=int, required=True)
    
    try:
        args = parser.parse_args()
    except argparse.ArgumentError:
        print("Invalid arguments provided. Please provide control, sample, and resolution.")
        return

    calculate_pearson(args.ctl_file, args.sample_file, args.resolution)

if __name__ == "__main__":
    main()
EOF

# Function to determine the control and sample files
determine_files_and_run() {
    local dir=$1
    local prefix=$2
    
    # Find all bedgraph files with the specific pattern
    bedgraph_files=$(find "$dir" -name "${prefix}mergedCrush_*.bedgraph" ! -name "*qfiltered*" ! -name "*reprocess*")
    
    # Extract resolutions and sort them
    resolutions=($(for f in $bedgraph_files; do basename "$f" | sed -E 's/.*_([0-9]+)\.bedgraph/\1/'; done | sort -nr))
    
    if [ ${#resolutions[@]} -lt 2 ]; then
        echo "Not enough resolutions to perform comparison. Exiting."
        return
    fi
    
    # Use the second highest resolution as the control
    control_resolution=${resolutions[0]}
    control_file=$(find "$dir" -name "${prefix}mergedCrush_${control_resolution}.bedgraph" | head -n 1)
    
    # Sample files: all files with lower resolutions than control
    sample_files=()
    for res in "${resolutions[@]:1}"; do
        sample_file=$(find "$dir" -name "${prefix}mergedCrush_${res}.bedgraph" | head -n 1)
        sample_files+=("$sample_file")
    done
    
    # Run the Python script for each sample file
    echo "Calculating Pearson correlation with control: $control_file"
    for sample_file in "${sample_files[@]}"; do
        python PearsonCorrelation.py -c "$control_file" -s "$sample_file" -r "$control_resolution" 2>/dev/null
    done
}

# Example usage
dir="./"  # Directory containing the bedgraph files
prefix="$outpre"  # Prefix of the bedgraph files

determine_files_and_run "$dir" "$prefix"

echo "Finished! Check the output."

##########################################

