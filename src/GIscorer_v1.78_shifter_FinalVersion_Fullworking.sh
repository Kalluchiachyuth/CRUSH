#!/bin/bash
version=1.78_PythonShifter

hicpath=0
myres=0
help=0
strawpath=0
sizefile=0
genesfile=0
fastafile=0
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
norm=0
whichchoice=`echo "o"`
cleanup=1
endZ=0
exclbed=0
menu=`printf "\nThis tool will create and use several temporary files in your current directory, so we recommend running it in an empty directory where you have write access.\noptions: \n-hic\t\tinput .hic file\n-straw\t\tpath to straw\n-res\t\tresolution desired or comma separated list of resolutions to obtain consensus\n-sizes\t\tchromosome size file with two columns corresponding to chromosome and size respectively\n-genes\t\tgene annotations in bed format\n-fasta\t\tfastafile for generating B-bins\n-cpu\t\tnumber of threads to use; default is 1\n\n-doNotMerge\t\tset this option to 1 to keep each resolution as a separate output file. Default is to merge in a way that provides maximum resolution.\n\n-adjustment\t\tSet this option to 1 to include a adjustment of GIvalues at the end. This adjustment shifts values based on any internal skewing of the data. Do not set this if using the GIscorer to compare between two Hi-C maps.\n\n-distance\t\tThe distance next to the diagonal to be filtered out in Mb. Default is 0 which considers everything.\n\n-trackline\t\tSet to 0 if you want to disable printing a bedgraph trackline header.\n\n-threshold\t\tDistance normalized threshold to filter out extreme outliers.\n\n-upperlim\t\tThe upperlimit of the distance away from the diagonal to consider. Default is 0 so that it considers the whole chromosome.\n\n-keeptracks\t\tSet to 1 in order to keep separate A and B tracks for the probability of interacting with bed file vs other regions.\n\n-switch\t\tSwitch on to initialize re-iteration. Default is 1.Set this to 0 for bypassing re-initialization. \n\n-norm\t\tNormalization scheme for dumping reads. Default is 0. Set this to 1 for VC_SQRT. \n\n-window\t\tSet to perform a sliding window average of the scores at individual resolutions. Default is to calculate the appropriate window based on sequencing depth. Set to 1 to remove sliding window.\n\n-v\t\tSet to 1 to run in verbose mode, extensive messages.\n\n-use\t\t Whether to use of overwrite existing GI tracks previously calculated at individual resolutions. Set this option to u to use previous calculations. Default is to recalculate. This option is useful for merging resolutions."`
while test $# -gt 0
do
        case "$1" in
                -res)
                res=$2
                ;;
                -hic)
                hicpath=$2
                ;;
                -straw)
                strawpath=$2
                ;;
                -sizes)
                sizefile=$2
		;;
		-genes)
		genesfile=$2
                ;;
		-fasta)
		fastafile=$2
		;;
		-exclbed)
		exclbed=$2
		;;
		-cpu)
		cpu=$2
		;;
		-doNotMerge)
		doNotMerge=$2
		;;
		-adjustment)
                adjustment=$2
                ;;
		-distance)
		distance=$2
		;;
		-upperlim)
		upperlim=$2
		;;
		-lowerthresh)
		lowerthresh=$2
		;;
		-trackline)
		trackline=$2
		;;
                -threshold)
		threshold=$2
		;;
		-keeptracks)
		keeptracks=$2
		;;
		-window)
		window=$2
		;;
		-v)
		verbose=$2
		;;
		-switch)
		switch=$2
        ;;
        -norm)
		norm=$2
		;;
		-use)
		whichchoice=$2
		;;
		-clean)
		cleanup=$2
		;;
		-endnorm)
		endZ=$2
		;;
		--help) echo "argument $1"
                help=1
                echo "$menu"
                exit 0
                ;;
        esac
        shift
done
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
if [ $strawpath == 0 ]
then
echo "You must specify a path to juicers straw function!......""$menu"
exit 0
fi
if [ $sizefile == 0 ]
then
echo "You must specify a path to a chromosome size file!......""$menu"
exit 0
fi
if [ $genesfile == 0 ]
then
echo "You must specify a path to a bed file with gene annotations!......""$menu"
exit 0
fi
if [ $fastafile == 0 ]
then
echo "You must specify a path to the organism fasta file!......""$menu"
exit 0
fi
if [ $switch == 1 ]
then
echo "Warning: Re-evaluation & Re-interation of A and B Bins is set. We recommend having this parameter to achieve highest possible confidence and resolution."
fi
if [ $norm == 0 ]
then
echo "Warning: Matrix Normalizing scheme is set to NONE. Set this to 1 for VC_SQRT"
fi
if [ $doNotMerge == 1 ]
then
echo "Warning: doNotMerge set. We recommend merging to achieve highest possible confidence and resolution."
fi
if [ $exclbed == 0 ]
then
echo "Warning: we highly recommend including a bed file to exclude known problematic / repetitive regions within the genome. This can typically be obtained from ENCODE: ENCSR636HFF, or generated yourself using github.com/Boyle-Lab/Blacklist."
fi
if [ $adjustment == 1 ]
then
echo "Warning: adjusting end compartmental values which may result in some loss of quantitative power. We recommend adjustment only if not comparing two samples."
fi
echo "GIscorer_v""$version"" -straw ""$strawpath"" -hic ""$hicpath"" -res ""$res"" -sizes ""$sizefile"" -genes ""$genesfile"" -fasta ""$fastafile"" -cpu ""$cpu"" -doNotMerge ""$doNotMerge"" -adjustment ""$adjustment"" -distance ""$distance"" -upperlim ""$upperlim"" -lowerthresh ""$lowerthresh"" -trackline ""$trackline"" -threshold ""$threshold"" -keeptracks ""$keeptracks"" -window ""$window"" -switch ""$switch"" -norm ""$norm" > GIparamters.txt
totchroms=`cat $sizefile | wc -l`
maxres=`echo "$res" | sed "s/,/ /g" | mawk '{max=$1;for(i=2;i<=NF;i++){if($i > max) max = $i} print max}'`
minres=`echo "$res" | sed "s/,/ /g" | mawk '{min=$1;for(i=2;i<=NF;i++){if($i < min) min = $i} print min}'`
if [ $endZ == 0 ]
then
endZ=`echo "$minres"`
fi

#check sizefile
numcolumns=`mawk '{if (NF != 2) print $0}' $sizefile | wc -l`
if [ "$numcolumns" -gt 0 ]
then
echo "size file is incorrect format. It should be two columns with chromosome name and size"
exit 0
fi

echo "Now calculating Bbins for the entire genome using fastafile and sizefile at 500bp"

resbins=`echo "size_500bp_bins"".bed"`
gcfile=`echo "gc_500bp_bins"".txt"`
gc_g_ga=`echo "gc_g_ga_500bp_bins"".bed"`
Bbins=`echo "Bbins_500bp"".bed"`

#Binning the size file into 500bp bins:
cat $sizefile | awk -v myres=500 '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
wait
echo "$resbins calculated at 500bp resolution"

#Now, calculating the gc content:
bedtools nuc -fi $fastafile -bed $resbins > $gcfile 
wait
echo "$gcfile generated at 500bp resolution"

#Created a %G/%G+%A file in 7th column 
cat $gcfile | grep -v user | cut -f 1-5 | awk '{print $0"\t"($4+$5)}' | awk '{if ($6 > 0) print $0}' | awk '{print $0"\t"($5/$6)}' > $gc_g_ga

wait

#Calculating mean and SD for the whole genome and then subtracting 2SD from the mean 
gc_thresh=`cat $gc_g_ga | awk '{s+=$7; ss+=$7*$7; linecount+=1} END{print m=s/linecount, sqrt(ss/linecount-m^2)}' | awk '{print $0}' | awk '{print $1 - 2*$2}'`

echo "$gc_thresh"

#Generating the bins by considering the bins below 
cat $gc_g_ga | awk -v thresh=$gc_thresh '{if ($7 < thresh) print $0}'  > $Bbins &

wait

echo "Finished generating Bbins file for all chromosomes at 500bp resolution"

#Defining the high res and coarse res for shifter script

#resolutions=`echo "resolutions"`
#echo $res | sed "s/,/ /g" > $resolutions

#IFS=',' read -ra resolutions <<< "$res"
#resolutions=($(echo $res | sed "s/,/ /g")) #Converting to array

countres=0
coarse_res=$maxres
crtmp=$maxres
#coarse_res=${resolutions[0]}

for myres in $(echo $res | sed "s/,/ /g")
#for myres in "${!resolutions[@]}"
do
high_res=$myres
if [ $window == 0 ]
then

#Changed the below line from mawk to awk as mawk was throwing an error while working through script line by line
totbins=`cat $sizefile | awk -v var=$myres '{print $1"\t"($2/var)**2}' | awk '{sum += $2} END {print sum}'`

newwindow=1
else
newwindow=`echo "$window"`

#Calculating window sizes for the resolutions 
newwindow=$(mawk -v win=$window -v res=$myres 'BEGIN{print int(win/(res/1000))+1}')
fi
echo "window size for ""$myres"" resolution is ""$newwindow"

outiefull=`echo "GIscore_""$myres"".bedgraph"`
Aoutiefull=`echo "AGIscore_""$myres"".bedgraph"`
Boutiefull=`echo "BGIscore_""$myres"".bedgraph"`

outcheck=0
outcheck=`ls | grep $outiefull`
if [ "$outiefull" == "$outcheck" ]
then
echo "$outiefull"" already exists.... exiting."
exit 0
fi

echo "Now dumping reads and calculating z-scores for all chromosomes at ""$myres"

task(){

mychr=`cat $sizefile | head -n $myiter | tail -1 | cut -f 1`
mysize=`cat $sizefile | head -n $myiter | tail -1 | cut -f 2`
rowbins=`cat $sizefile | head -n $myiter | tail -1 | mawk -v myres=$myres '{print int($2/myres)+1}'`
outie=`echo "GIscore_""$myres""_""$mychr""_tmp"`
outie2=`echo "GIscore_""$myres""_""$mychr""_tmp2"`
tmpfiles=`echo "ABtmpfiles_""$mychr""_""$myres"`
Aoutie=`echo "AGIscore_""$myres""_""$mychr""_tmp"`
Boutie=`echo "BGIscore_""$myres""_""$mychr""_tmp"`

mkdir $tmpfiles


efile=`echo "expecteds"`
ofile=`echo "observeds"`
mydist=$((distance * 1000000))
myupper=$((upperlim * 1000000))
if [ "$verbose" -gt 0 ]
then
echo "Getting genic bins."
fi
sortedbed=`echo "$tmpfiles""/genebins"`
pseudoB=`echo "$tmpfiles""/pseudoB"`
pseudoA=`echo "$tmpfiles""/pseudoA"`
dumped=`echo "$tmpfiles""/dumped_""$myres""_""$mychr""_tmp"`

#Adding in the code for re-initializing Bbins for next following resolutions
if [ "$myres" == "$maxres" ] 
then

#Processing the genes & Bbins files
cat $genesfile | mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

#Bbins
cat $Bbins |  mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

else

#Processing the genes & Bbins files for other resolutions 
cat $combined_newAbins | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

#Bbins
cat $combined_newBbins |  mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

fi

#Adding in the code for bypassing re-initializing Bbins for next following resolutions
if [ "$switch" -lt 1 ]
then

#Processing the genes & Bbins files
cat $genesfile | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

#Bbins
cat $Bbins |  mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

fi

if [ "$verbose" -gt 0 ]
then
echo "Dumping ""$mychr"" chromosomal reads with straw."
fi
newofile=`echo "$tmpfiles""/observeds"`

#Adding VCoutie which has information about vc normalizaion for each chromosome
VCoutie=`echo "$tmpfiles""/VCoutie_""$myres""_""$mychr""_tmp"`

newofile_vcfilt=`echo "$tmpfiles""/observeds_vcfilt"`
newofile2=`echo "$tmpfiles""/observeds_final"`
newefile2=`echo "$tmpfiles""/expecteds_final"`
VCfinal=`echo "$tmpfiles""/VCfinal_""$myres"".bed"`

touch  $newofile
touch  $newefile
errofile=`echo "errorlog"`

if [ "$norm" -eq 0 ]
then

#Dumping the reads using straw at multiple resolutions on multiple chromosomes using VC_SQRT
$strawpath NONE $hicpath $mychr $mychr BP $myres | grep -v WARN > $dumped &

wait

else

echo "Dumping the reads using VCSQRT"

#Dumping the reads using straw at multiple resolutions on multiple chromosomes
$strawpath VC_SQRT $hicpath $mychr $mychr BP $myres | grep -v WARN > $dumped &

wait

fi

#Calculating observed 
cat $dumped | mawk -v var=$myres -v distfilt=$mydist '{if ($2-$1 > distfilt) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' > $newofile &

wait

#Calculating VC using straw at highest resolution:
cat $dumped | awk '{b[$2]+=$3; b[$1]+=$3; c+=$3*2} END { for (i in b) { print i"\t"b[i]"\t"c } }'  | awk '{a[$1] = $2; b[$1] = $3; count+=1} END {for (i in a) {print i"\t"a[i]"\t"b[i]"\t"count}}' | awk '{print $1"\t"$2"\t"$3/$4}' | awk '{print $1"\t"$2/$3}' > $VCoutie &
wait

#Generating multiple files for processing of VC normalization
vczero_expanded=`echo "$tmpfiles""/vczero_expanded_""$mychr""_""$myres"`
vcnonzero=`echo "$tmpfiles""/vcnonzero_""$mychr""_""$myres"`
vcremove_expanded=`echo "$tmpfiles""/vcremove_expanded_""$mychr""_""$myres"`
vczeroremovecomb=`echo "$tmpfiles""/vczeroremovecomb_""$mychr""_""$myres"`
vckeep=`echo "$tmpfiles""/vckeep_""$mychr""_""$myres"`


#Considering both the non-zero as well as zero values from VC file and then log transforming the numeric values from column 4:
cat $VCoutie | awk '{if ($2 == 0) print $0}' | awk -v res=$myres -v chr=$mychr '{print chr"\t"$1-res"\t"$1+2*res"\t"$2}' | awk '{if ($2 >= 0) print $0}'  > $vczero_expanded &

wait

cat $VCoutie | awk '{if ($2 != 0) print $0}' | awk -v chr=$mychr '{print chr"\t"$0"\t"log($2)}'  > $vcnonzero &
wait

#Calculating the mean and SD of vc0filtered log values and further generating a cut-off values of mean - 2*SD for filtering out the repeating values
vcthresh=`cat $vcnonzero | awk '{s+=$4; ss+=$4*$4; linecount+=1} END{print m=s/linecount , sqrt(ss/linecount-m^2)}' | awk '{print $1"\t"$2}' | awk '{print ($1 - 2*$2)}'`

echo "$vcthresh"

#Now filtering out the vc non-zero values and considering the criteria of less than and greater than threshold. Also, moving the coordinates 1 bin to the left and 1 bin to the right with vcoutie.
cat $vcnonzero | awk -v thr=$vcthresh -v  res=$myres '{if ($4 < thr) print $1"\t"$2-res"\t"$2+2*res}' | awk '{if ($2 >= 0) print $0}' > $vcremove_expanded

cat $vcnonzero | awk -v thr=$vcthresh -v res=$myres '{if ($4 >= thr) print $1"\t"$2"\t"$2+res}' > $vckeep

cat $vczero_expanded $vcremove_expanded | sort -k 1b,1b --stable > $vczeroremovecomb

cat $vckeep | intersectBed -v -a stdin -b $vczeroremovecomb > $VCfinal &

wait

#Generating a newofile2 which has all the bins with vc values from VCoutie.
awk 'NR==FNR {c1[$2] = $2; next} {if ($2 in c1) print $0}' $VCfinal $newofile > $newofile_vcfilt

#Re-calulating the observed & expected using vc filtered file.
cat $newofile_vcfilt | mawk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile2 -v newefile=$newefile2 '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile } ' &

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

onfile=`echo "$tmpfiles""/distnorms"`
tmpmean=`echo "$tmpfiles""/tmpmean"`
symmfile=`echo "$tmpfiles""/symmonfile"`
touch  $tmpmean
touch  $symmfile

if [ "$myupper" -eq 0 ]
then
distlim=100000000000000000
else
distlim=`echo "$myupper"`
fi

mawk -v distlim=$distlim 'NR==FNR { c1[$1] = $2; next} {if (($2 != $3) && ($3-$2 < distlim)) print $2"\t"$3"\t"($4+1)/(c1[$1]+1)}' $newefile2 $newofile2 | mawk -v thresh=$threshold '{if ($3 < thresh) print $1"\t"$2"\t"$3"\n"$2"\t"$1"\t"$3; else print $1"\t"$2"\t"thresh"\n"$2"\t"$1"\t"thresh}' | awk -v symmfile=$symmfile -v tmpmean=$tmpmean -v rowbins=$rowbins '{m[$1]+=$3;d[$1]+=($3**2); print $0 >> symmfile} END {for (i in m) print i"\t"m[i]/rowbins"\t"((d[i]/rowbins-((m[i]/rowbins)**2))**.5) >> tmpmean }'  2> $tmpfiles/$errorfile

wait

if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "Done normalizing. It took ""$seconds" "seconds. Now cleaning up a bit for ""$mychr"" chromosome."
start=$(date +%s)
fi

Afile=`echo "tmpA"`
Bfile=`echo "tmpB"`
Bpfile=`echo "tmpBp"`
Apfile=`echo "tmpAp"`

if [ "$myres" -gt 50000 ]
then
oppocheck=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) ; else if ($2 in c1); else print $0}' $sortedbed $symmfile | wc -l`
else
oppocheck=2
fi

if [ "$oppocheck" -lt 1 ]
then
echo "WARNING: All bins on chromosome ""$mychr"" overlap with your bed file at ""$myres"" resolution"". This should not affect the scores if using multiple resolutions. If using a single resolution, then maybe not use this resolution?."

else

myzscore=`echo "Zscorefile"`
if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "Got means and stdevs, It took ""$seconds"" seconds. Now applying z-scores values for ""$mychr"" chromosome."
start=$(date +%s)
fi
Atmp=`echo "$tmpfiles""/Atmp_""$mychr""_""$myres"`
Btmp=`echo "$tmpfiles""/Btmp_""$mychr""_""$myres"`
Atmppseudo=`echo "$tmpfiles""/Atmpp_""$mychr""_""$myres"`
Btmppseudo=`echo "$tmpfiles""/Btmpp_""$mychr""_""$myres"`

touch  $Atmp
touch  $Btmp

mawk 'NR==FNR { m[$1]=$2; d[$1]=$3; next} {if (d[$2] > 0) print $1"\t"$2"\t"($3-m[$2])/(d[$2])}' $tmpmean $symmfile > $tmpfiles/RowZs
mawk -v Atmp=$Atmp -v Btmp=$Btmp 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) print $2"\t"$3 >> Atmp; else print $2"\t"$3 >> Btmp}' $sortedbed $tmpfiles/RowZs

wait

mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) print $2"\t"$3}' $pseudoB $tmpfiles/RowZs > $Btmppseudo

wait

if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "Done z-scoring. It took ""$seconds"" seconds. Now getting B values for ""$mychr"" chromosome."
start=$(date +%s)
fi

wait

if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "It took ""$seconds"" seconds. Now getting A values for ""$mychr"" chromosome."
start=$(date +%s)
fi

cat $Atmp | mawk '{b[$1]+=$2;c[$1]++} END { for (i in b) { print i"\t"(b[i])*(c[i])"\t"c[i] } } ' > $tmpfiles/$Afile 2> $tmpfiles/$errorfile

cat $Btmppseudo | mawk '{b[$1]+=$2;c[$1]++} END { for (i in b) { print i"\t"(b[i])*(c[i])"\t"c[i] } } ' > $tmpfiles/$Bpfile 2> $tmpfiles/$errorfile

wait 

myABsum=`cat $tmpfiles/$Afile $tmpfiles/$Bpfile | mawk '{sum+=$3} END {print sum}'`

if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "It took ""$seconds"" seconds. Accounting for zeros for ""$mychr"" chromosome."
start=$(date +%s)
fi

zfile=`echo "zerofile"`
cat $sizefile | mawk -v var=$mychr '{if ($1 == var) print $0}' | mawk -v wlkr=$myres '{for (i=0;i<=$2;i+=wlkr) print i"\t0"}' > $tmpfiles/$zfile


A2file=`echo "tmpA2"`
B2file=`echo "tmpB2"`
Ap2file=`echo "tmpAp2"`
Bp2file=`echo "tmpBp2"`

mawk '{b[$1]+=$2} END { for (i in b) { print i,b[i] } } ' $tmpfiles/$zfile $tmpfiles/$Afile > $tmpfiles/$A2file 2> $tmpfiles/$errorfile

mawk '{b[$1]+=$2} END { for (i in b) { print i,b[i] } } ' $tmpfiles/$zfile $tmpfiles/$Bpfile > $tmpfiles/$Bp2file 2> $tmpfiles/$errorfile

#Removing Zero' from the Abins & Bbins 
A2removedfile=`echo "A2removed"`
B2removedfile=`echo "B2removed"`
mawk 'NR==FNR { c[$1] = $2; next} {if (($2 == 0) && (c[$1] ==0)); else print $1"\t"c[$1]}' $tmpfiles/$A2file $tmpfiles/$Bp2file > $tmpfiles/$A2removedfile
mawk 'NR==FNR { c[$1] = $2; next} {if (($2 == 0) && (c[$1] ==0)); else print $1"\t"$2}' $tmpfiles/$A2file $tmpfiles/$Bp2file > $tmpfiles/$B2removedfile
 
wait

if [ "$verbose" -gt 0 ]
then
end=$(date +%s)
seconds=$(echo "$end - $start" | bc)
echo "It took ""$seconds"" seconds. Comparing A vs B scores for ""$mychr"" chromosome."
start=$(date +%s)

fi

wait

#Calculating the output
mawk -v var=$myres -v myABsum=$myABsum '{if ($1 in b) b[$1]-=$2; else b[$1]=$2; if ($1 in a); else a[$1]=$2} END { for (i in b) { print i"\t"(b[i])/myABsum } } ' $tmpfiles/$A2removedfile $tmpfiles/$B2removedfile | sort -k 1bn,1b --stable | mawk -v mchr=$mychr -v myres=$myres -v window=$newwindow -v OFS="\t" 'BEGIN{slide=1} {mod=NR%window; if(NR<=window){count++}else{sum-=array[mod]}sum+=$2;array[mod]=$2;} (NR%slide)==0{print mchr,int($1-((myres*window)/2)),sum/count}' | mawk -v myres=$myres '{print $1"\t"$2-int(myres/2)+myres"\t"$2+(int(myres/2))+myres"\t"$3}' | mawk '{if ($2 > 0) print $0}' > $outie 2> $tmpfiles/$errorfile

if [ "$exclbed" != 0 ]
then
exclbins=`echo "$tmpfiles""/exclbins"`
excloutie=`echo "$tmpfiles""/exclout_""$myres""_""$mychr"`

cat $exclbed | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $exclbins
mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $0}' $exclbins $outie > $excloutie
wait
mv $excloutie $outie
fi

if [ $keeptracks -eq 1 ]
then
rawcomp=`echo "rawcompvectors"`

fi
fi
if [ "$oppocheck" -gt 0 ] && [ "$adjustment" -gt 0 ]
then
amed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1) print $4}' $sortedbed $outie | mawk '{sum +=$1} END {print sum/NR}'`
bmed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $4}' $sortedbed $outie | mawk '{sum += $1} END {print sum/NR}'`

cat $outie | mawk -v varA=$amed -v varB=$bmed '{print $1"\t"$2"\t"$3"\t"$4-(varA-(varB*-1))}' > $outie2
wait
mv $outie2 $outie
fi

if [ "$myres" -ge "$endZ" ]
then
echo "Applying end ending score normalization"
myendmean=`cat $outie | mawk '{sum += $4; cnum++} END {print sum/cnum}'`
echo "$myendmean"
tmpendZ=`echo $tmpfiles"/tmpendnorm"`
cat $outie | mawk -v mymean=$myendmean '{print $1"\t"$2"\t"$3"\t"($4-mymean)}' > $tmpendZ
wait
mv $tmpendZ $outie
fi


if [ "$cleanup" -gt 0 ]
then
rm -r $tmpfiles
fi

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
    )&
}


open_sem $cpu
if [ "$whichchoice" == "o" ]
then

for (( myiter=1; myiter<=$totchroms; myiter++ ))
do
    run_with_lock task $myiter
done

wait

echo "Merging individual chromosome files"
cat GIscore*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull

wait

if [ "$myres" == "$maxres" ]
then
coarse_infile=`echo "GI_""$myres"".bedgraph"`
cat $outiefull > $coarse_infile

wait
fi
res_values=`echo "res_values"`
echo "$myres" >> $res_values
#Re-evaluating both A and B bins for re-initialization for next resolution


if [ "$countres" -gt 0 ]
then

# Print current coarse and high resolutions
echo "Coarse resolution: ""$coarse_res"", High resolution: ""$high_res"


#Creating a python output file 
python_output=`echo "shifter.bedgraph"`

if [ "$coarse_res" -eq "$maxres" ] && [ "$coarse_res" != "$high_res" ]
then

# Shifter Input files for the specified resolutions
high_res_infile=`echo "$outiefull"`
coarse_res_infile=`echo "$coarse_infile"`

else
# Shifter Input files for the specified resolutions
high_res_infile=`echo "$outiefull"` 
coarse_res_infile=`echo "$python_output"`

fi

# Conditional checking for multiples of resolutions for shifter 
#If coarseres is divisible by highres, use these resolutions for shifter as coarse and high
if [[ "$(( $coarse_res % $high_res ))" -eq 0 ]]
then
echo "Coarse resolution: ""$coarse_res"", High resolution: ""$high_res"

crtmp=$coarse_res

else




if [[ "$(( $crtmp % $high_res ))"  -eq 0 ]]
then

echo $crtmp

coarse_res=$crtmp

echo $coarse_res

echo "Coarse resolution: ""$crtmp"", High resolution: ""$high_res"

fi

fi

echo "Executing Shifter using Coarse resolution: ""$coarse_res"" and High resolution:  ""$high_res"


#Executing the python script
python /Zulu/achyuth/Projects/GIscores_vs_Eigen/HDDv77/shifter_Test/shifter_integration/hdrescompare.py -r "$high_res" -c "$coarse_res" -1 "$high_res_infile" -2 "$coarse_res_infile" -s "$sizefile" -o "$python_output"

wait


#Seperating out A and B regions from the outiefull and intersecting them with genes and Bbins to improve compartments identification.

AbinsR=`echo "AbinsfromGI_""$myres"".txt"`
BbinsR=`echo "BbinsfromGI_""$myres"".txt"`

newBbins=`echo "newBbinsGI_""$myres"".txt"`
newAbins=`echo "newAbinsGI_""$myres"".txt"`

non_newBbins=`echo "non_newBbinsGI_""$myres"".txt"`
non_newAbins=`echo "non_newAbinsGI_""$myres"".txt"`

combined_newBbins=`echo "combined_newBbinsGI_""$myres"".txt"`
combined_newAbins=`echo "combined_newAbinsGI_""$myres"".txt"`

cat $python_output | awk '{if ($4 > 0) print $0}' > $AbinsR

cat $python_output | awk '{if ($4 < 0) print $0}' > $BbinsR

if [ "$switch" -eq 1 ]
then

echo "Re-initialization of Bins is switched ""ON"". Re-evaluating both A and B bins and re-initializing for next resolution."

cat $Bbins | intersectBed -u -a stdin -b $BbinsR > $newBbins

echo "Printing out the file with bins on the chromosomes not present in ""$newBbins"" at ""$myres""."

awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newBbins $Bbins > $non_newBbins

cat $newBbins $non_newBbins > $combined_newBbins 

cat $genesfile | cut -f 1-3 | intersectBed -u -a stdin -b $AbinsR > $newAbins

awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newAbins $genesfile > $non_newAbins

cat $newAbins $non_newAbins > $combined_newAbins 

else

echo "Generating the output without re-initialization."

fi

# Update coarse and high resolutions
coarse_res=$high_res


else

#Seperating out A and B regions from the outiefull and intersecting them with genes and Bbins to improve compartments identification.

AbinsR=`echo "AbinsfromGI_""$myres"".txt"`
BbinsR=`echo "BbinsfromGI_""$myres"".txt"`

newBbins=`echo "newBbinsGI_""$myres"".txt"`
newAbins=`echo "newAbinsGI_""$myres"".txt"`

non_newBbins=`echo "non_newBbinsGI_""$myres"".txt"`
non_newAbins=`echo "non_newAbinsGI_""$myres"".txt"`

combined_newBbins=`echo "combined_newBbinsGI_""$myres"".txt"`
combined_newAbins=`echo "combined_newAbinsGI_""$myres"".txt"`

cat $outiefull | awk '{if ($4 > 0) print $0}' > $AbinsR

cat $outiefull | awk '{if ($4 < 0) print $0}' > $BbinsR

if [ "$switch" -eq 1 ]
then

echo "Re-initialization of Bins is switched ""ON"". Re-evaluating both A and B bins and re-initializing for next resolution."

cat $Bbins | intersectBed -u -a stdin -b $BbinsR > $newBbins

echo "Printing out the file with bins on the chromosomes not present in ""$newBbins"" at ""$myres""."

awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newBbins $Bbins > $non_newBbins

cat $newBbins $non_newBbins > $combined_newBbins 

cat $genesfile | cut -f 1-3 | intersectBed -u -a stdin -b $AbinsR > $newAbins

awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newAbins $genesfile > $non_newAbins

cat $newAbins $non_newAbins > $combined_newAbins 

else

echo "Generating the output without re-initialization."

fi

fi

fi

if [ $keeptracks -eq 1 ] && [ "$whichchoice" == "o" ]
then
cat AGIscore*_tmp | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $Aoutiefull &
cat BGIscore*_tmp | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $Boutiefull &
fi
wait
rm GIscore*_tmp 2> smallerrors
rm AGIscore*_tmp 2> smallerrors
rm BGIscore*_tmp 2> smallerrors

countres=$((countres+1))
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
echo "Only one resolution listed, so nothing to merge. Check GIscore_""$totres"".bedgraph for the single resolution GI score"

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
finalout=`echo "mergedGIscore_""$minres"".bedgraph"`
Afinalout=`echo "mergedGIscore_A_""$minres"".bedgraph"`
Bfinalout=`echo "mergedGIscore_B_""$minres"".bedgraph"`
cat GIscore_*.bedgraph | mawk -v minr=$minres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4/int($3-$2)}' | mawk -v mr=$minres '{for (i=$2;i<$3;i++) print $1"\t"int(i)*mr"\t"(int(i)+1)*mr"\t"$4}' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable | groupBy -i stdin -g 1,2,3 -c 4 -o mean > $finalout 2> smallerrors
wait
GIave=`cat $finalout | mawk '{if ($4 < 0) print "1\t2\t3\t"$4*-1; else print "1\t2\t3\t"$4}' | groupBy -i stdin -g 1 -c 4 -o mean | cut -f 2`

if [ $trackline -eq 0 ]
then
cat $finalout | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' > GIscores_todelete
else
cat $finalout | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20,20 autoScale off\n"$0; else print $0}' > GIscores_todelete
fi
wait
mv GIscores_todelete $finalout
if [ $keeptracks -eq 1 ]
then
cat AGIscore_*.bedgraph | mawk -v minr=$minres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$minres '{for (i=$2;i<$3;i++) print $1"\t"int(i)*mr"\t"(int(i)+1)*mr"\t"$4}' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable | groupBy -i stdin -g 1,2,3 -c 4 -o sum > $Afinalout 2> smallerrors&
cat BGIscore_*.bedgraph | mawk -v minr=$minres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$minres '{for (i=$2;i<$3;i++) print $1"\t"int(i)*mr"\t"(int(i)+1)*mr"\t"$4}' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable | groupBy -i stdin -g 1,2,3 -c 4 -o sum > $Bfinalout 2> smallerrors &
fi
wait

#rm GIscore_*.bedgraph 2> smallerrors
rm AGIscore_*.bedgraph 2> smallerrors
rm BGIscore_*.bedgraph 2> smallerrors
rm AbinsfromGI_*.txt 2> smallerrors
rm BbinsfromGI_*.txt 2> smallerrors
rm newAbinsGI_*.txt 2> smallerrors
rm newBbinsGI_*.txt 2> smallerrors
rm non_newAbinsGI_*.txt 2> smallerrors
rm non_newBbinsGI_*.txt 2> smallerrors
rm combined_newAbinsGI_*.txt 2> smallerrors
rm combined_newBbinsGI_*.txt 2> smallerrors
rm smallerrors
echo "Finished! Check ""$finalout"