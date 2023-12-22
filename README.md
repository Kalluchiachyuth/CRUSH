![](/examples/figures/CRUSH_logo.jpg)
# CRUSH


### What is CRUSH?
CRUSH (**C**ompartmental **R**efinement for **U**ltraprecise **S**tratification within **H**i-C) is a tool that can identify fine-scale compartments in chromatin conformation matrices. It has successfully identified compartments in Hi-C, Micro-C, and Single-Cell Hi-C. CRUSH specializes in identifying fine-scale compartments at high resolutions with significantly lower read depth than other compartment calling tools.  
#### Visual Algorithmic Explanation:
![CRUSH](/examples/figures/CRUSHdiag.png)

### How to install CRUSH

_Dependencies_
```
BASH environment 
bedtools (intersectBed)
cooler (if using .mcool files)
python3
hicstraw (if working with .hic files)
numpy
scipy
statsmodels
tqdm
```

After these dependencies are installed, simply clone this repository. 

To run CRUSH, you may want to make the tool directly executable. e.g. sudo chmod +x CRUSH/CRUSH_v1.0


### How to use CRUSH

_Required files:_

**HIC**

This can be either a .hic file (from juicer) or a .mcool file (from cooler). 

**SIZEFILE**

This is a two-column tab-delimited file with chromosome names and sizes.

**ABED**

This file sets the initialization A states and can be any bed file with at least 3 tab-separated columns. We typically use genes, but use your imagination. For example, ChIP-seq peaks for an active mark should work well also.

**BBED**

This file sets the initialization B states. It can be a bed file or it can be a fasta file. It should simply correspond to something that correlates with the inactive compartment. If you specify a fasta file, we will use gc content to calculate. This can slow it down a bit, but works well.

**PLEASE ENSURE THAT THE CHROMOSOME NAMES MATCH BETWEEN ALL OF YOUR FILES.**



### Example usage

CRUSH_v1.0 -i Myfile.hic -g hg38.sizes -a hg38_genes.bed -b hg38.fasta -r 1000 -cpu 23 -o output_prefix


CRUSH has other optional parameters. We do not recommend that you change these unless you are confident in what they do. 

CRUSH_v1.0 --help

```
usage : CRUSH_vX.X -i HIC  -g SIZEFILE -a ABED -b BBED | FASTA -r FINERESOLUTION [-cpu CPU] [-w WINDOW] [-h]
Use option -h|--help for more information

CRUSH will create a temporary folder in your current directory, so make sure you have write access to your current directory.

-----------------------------------------
OPTIONS:
-h|--help                :  Display this help menu

--------------------------REQUIRED PARAMETERS------------------
-i|--hic                 :  Input .hic/.cooler/.mcool file by specifying the path. (e.g., '/path/to/ file.hic or file.mcool or file.cooler)
-g|--genomesize              :  Specify path to a chromosome size file with two columns corresponding to chromosome and size respectively.
-a|--initialA                :  Specify path to a bed file with the regions for initializing A. For example, gene annotations e.g. hg19genes.bed.
-b|--initialB                :  Specify path to either a fasta file or to a bed file for initializing B. If you specify a fasta file, we will calculate intialB from gc content.
-r|--res                 :  Resolution desired.
---------------

--------------------------OPTIONAL PARAMETERS------------------
-o|--outpre              :  Set this if you want to specify a prefix for the output files
-c|--cpu               :  Set the value for cpu number of threads to use. Default is 1.
-n|--no-merge            :  Set this option to 1 to keep each resolution as a separate output file. Default is to merge in a way that provides maximum resolution.
-A|--adjustment          :  Set this option to 1 to include a adjustment of CRUSH values at the end. This adjustment shifts values based on any internal skewing of the data. Do not set this if using CRUSH to compare between two Hi-C maps.
-d|--distance            :  Using this option will filter out the distance next to the diagonal. Default is 0 which considers everything.
-u|--upperlim            :  The upperlimit of the distance away from the diagonal to consider. Default is 0 so that it considers the whole chromosome.
-t|--trackline           :  Set to 0 if you want to disable printing a bedgraph trackline header.
-T|--threshold           :  Distance normalized threshold to filter out extreme outliers.
-k|--keeptracks          :  Set to 1 in order to keep separate A and B tracks for the probability of interacting with bed file vs other regions.
-l|--lowerthresh         :  Set the value for lowerthresh.
-w|--window              :  Set to perform a sliding window average of the scores at individual resolutions. Default is to calculate the appropriate window based on sequencing depth. Set to 1 to remove sliding window.
-v|--verbose             :  Set to 1 to enable verbose mode showing extensive messages.
-S|--switch              :  Set this option to 0 for bypassing re-initialization. Default is 1.
-C|--cleanup             :  Set the value to 0 to keep the temporary files. Default is 1.
-E|--endZ                :  Set the value for endZ. 
-x|--exclbed             :  Set the value for exclbed.
-q|--qvalue              :  Set the qvalue threshold. default 0.05. Set to 0 to not perform qvalue filtering. The qvalues will be reported as a separate track regardless.
-u|--use                 :  Whether to use of overwrite existing GI tracks previously calculated at individual resolutions. Set this option to u to use previous calculations. Default is to recalculate. This option is useful for merging resolutions.
-f|--tmpfolder         :  Set this if you want to name the temporary folder yourself. Make sure it doesn't already exist in your current working directory. Default is to name it CRUSHtmp with a random number.
-m|--maxres             : Set this to the coarsest resolution you want to consider. Default is to check every resolution present in the .hic or .mcool file between 2500000 and your desired resolution to inform each other.
```

### Output Files

CRUSH's main output is 4 files, each of which starts with the prefix that you specified with the -o option. The file name endings are:

CRUSHparameters.txt: Contatins a simple record of the parameters and files that you used.

mergedCrush_resolution.bedgraph (resolution is replaced with whatever you specified with the -r option): This is the main output file containing the scores for A (positive) and B (negative). Note that unlike eigenvector, you do not need to flip these calls, because A is always positive and B is always negative.

mergedqvalue_resolution.bedgraph: This bedgraph track contains an estimated q-value for each bin's score.

mergedCrush_resolution_qfiltered.bedgraph: This bedgraph track contains only scores that meet the qvalue threshold. However, the thresholding currently seems overly stringent, and we've obtained excellent results without this filter.  

While running, CRUSH will also create a temporary directory. The default is to name it CRUSHtmp_[randomnumber], but you can specify a name for this directory using the -f option. After completing, CRUSH will remove this directory by default, however, if you want to keep all the temporary files, you can use -C 0.
