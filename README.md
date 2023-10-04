![](/examples/figures/logo.png)
# CRUSH


### What is Crush?
CRUSH (**C**ompartmental **R**efinement for **U**ltraprecise **S**tratification within **H**i-C) is a tool that can identify fine-scale compartments in chromatin conformation matrices. It has successfully identified compartments in Hi-C, Micro-C, and Single-Cell Hi-C. CRUSH specializes in identifying fine-scale compartments at high resolutions with significantly lower read depth than other compartment calling tools.  
#### Visual Algorithmic Explanation:
![](/examples/figures/vis_abstract.png)

### How to install CRUSH
Dependencies  
Any other info (Java version?)  

### How to use CRUSH
```
-hic		input .hic file
-straw		path to straw
-res		resolution desired or comma separated list of resolutions to obtain consensus
-sizes		chromosome size file with two columns corresponding to chromosome and size respectively
-genes		gene annotations in bed format
-fasta		fastafile for generating B-bins
-cpu		number of threads to use; default is 1

-doNotMerge		set this option to 1 to keep each resolution as a separate output file. Default is to merge in a way that provides maximum resolution.

-adjustment		Set this option to 1 to include a adjustment of GIvalues at the end. This adjustment shifts values based on any internal skewing of the data. Do not set this if using the GIscorer to compare between two Hi-C maps.

-distance		The distance next to the diagonal to be filtered out in Mb. Default is 0 which considers everything.

-trackline		Set to 0 if you want to disable printing a bedgraph trackline header.

-threshold		Distance normalized threshold to filter out extreme outliers.

-upperlim		The upperlimit of the distance away from the diagonal to consider. Default is 0 so that it considers the whole chromosome.

-keeptracks		Set to 1 in order to keep separate A and B tracks for the probability of interacting with bed file vs other regions.

-switch		Switch on to initialize re-iteration. Default is 1.Set this to 0 for bypassing re-initialization. 

-norm		Normalization scheme for dumping reads. Default is 0. Set this to 1 for VC_SQRT. 

-window		Set to perform a sliding window average of the scores at individual resolutions. Default is to calculate the appropriate window based on sequencing depth. Set to 1 to remove sliding window.

-v		Set to 1 to run in verbose mode, extensive messages.

-use		 Whether to use of overwrite existing GI tracks previously calculated at individual resolutions. Set this option to u to use previous calculations. Default is to recalculate. This option is useful for merging resolutions.
```

### Author Contributions / How to Cite
Achyuth: King of the Pirates  
Jordan: Our Fearless Leader  

