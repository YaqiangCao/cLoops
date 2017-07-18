# Fine and differential enriched loops calling for protein-centric chromatin conformations with cLoops

# not available yet 

## Introduction
By taking the mapped paired-end tags from ChIA-PET or HiChIP as 2D points, the problem for calling loops is converted to draw significant clusters from sparse points with noise. After classifying the detected clusters into self-ligation and inter-ligation clusters, the significances of the inter-ligation clusters are estimated using permuted local backgrounds. We implemented the approach in the “cLoops (see/ChIA-PET/HiChIP loops)” package. Although without the peak calling step, the anchors determined by cLoops shows a high overlap with the peaks. By comparing to peaks based loop calling tools, we show that cLoops can detect more interactions with better ranked p-values, better supported by Hi-C data, sharper anchors, higher enrichment for TF motifs, work well both for sharp and broad peak like ChIA-PET data.

If you find cLoops is useful, please cite our paper:    
*** Fine and differential enriched loops calling for protein-centric chromatin conformations with cLoops ***

--------
## Install
[scipy](https://www.scipy.org/),[numpy](http://www.numpy.org/), [seaborn](https://seaborn.pydata.org/), [pandas](http://pandas.pydata.org/),[HTSeq](https://github.com/simon-anders/htseq) and [joblib](https://pythonhosted.org/joblib/) are required. If you have problems for installing scipy, please refer to [Anaconda](https://docs.continuum.io/anaconda/) or [SAGE](http://www.sagemath.org/).
```
tar xvzf cLoops.tar.gz
cd cLoops
python setup.py install    
```

or just

```
pip install cLoops
```
Please refer to [here](https://docs.python.org/2/install/index.html) to install cLoops to customized path.

--------
## Usage
Run ***cLoops -h*** to see all options. Key parameters are ***eps***, ***minPts*** and ***twice***. ***minPts*** defines at least how many PETs are required for a candidate loop, 5 is good for ChIA-PET and 20 is good enough for HiChIP. ***eps*** defines the distance requried for two PETs being neighbors, for sharp peaks like ChIA-PET data, cLoops can auto estimate it from the data as 2 fold of the fragment size. For broad peaks like ChIA-PET data, empirical experience is 2000, for HiChIP, 2000 worth a first trial. For practically usage, using the PETs in the smallest chromosome except chrY and chrM, then run a series of ***eps***, choose the smallest ***eps*** that can get well seperated inter-ligation and self-ligation PETs distance distributions. The ***twice*** model means first run with the input or auto estimated ***eps***, then using the self-ligation clusters re-estimate a ***eps*** and the distance cutoff for self-ligation PETs, cLoops then runs the clustering again with the new ***eps*** and distance cutoff, the final candidate loops are combined with the first round and second round clustering. The mode is good for not a ideal ***eps***, HiChIP and broad peaks like ChIA-PET.

--------
### Input  
Mapped PETs in [BEDPE format](http://bedtools.readthedocs.io/en/latest/content/general-usage.html), compressed files with gzip are also accepected, first 6 columns as following are necessary: chrom1,start1,end1,chrom2,start2,end2.

--------
### Output
The main output is a loop file and a PDF file for the plot of self-ligation and inter-ligation PETs distance distributions.
For the .loop file, columns and explaination are as follwing:

column | name | explaination
------ | ---- | ------------
0th | loopId | Id for a loop, like chr1-chr1-1
1th | ES | Enrichment score for the loop, caculated by observed PETs number divided by the mean PETs number of nearby permutated  regions
2th | FDR | false discovery rate for the loop, caculated as the number of permutated regions that there are more observed PETs than the region  
3th | binomal\_p-value | binomal test p-value for the loop
4th | distance | distance (bp) between the centers of the anchors for the loop
5th | hypergeometric\_local\_FDR | FDR for the hypergeometric test p-value compared to permutated regions
6th | hypergeometric\_p-value | hypergeometric test p-value for the loop
7th | iva | genomic coordinates for the left anchor, for example, chr13:50943050-50973634
8th | ivb | genomic coordinates for the right anchor
9th | poisson_p-value | poisson test p-value for the loop
10th | ra | observed PETs number for the left anchor
11th | rab | observed PETs number linking the left and right anchors
12th | rb | observed PETs number for the right anchor
13th | poisson\_p-value\_corrected | Bonferroni corrected poisson p-value according to number of loops for each chromosome
14th | binomal\_p-value\_corrected | Bonferroni corrected binomal p-value according to number of loops for each chromosome
15th | hypergeometric\_p-value\_corrected | Bonferroni corrected hypergeometric p-value according to number of loops for each chromosome
16th | significant | 1 or 0, 1 means we think the loop is significant compared to permutated regions. For ChIA-PET data, significant requiring ES >=1.0, FDR <=0.05, hypergeometric\_local\_FDR <=0.05 and all uncorrected p-values <= 1e-5; For HiChIP data, significant requiring ES >= 2.0, FDR <=0.01, corrected poisson and binomal p-values <=0.01. You can ignore this and customize your cutoffs.

--------
## Examples
1. ChIA-PET data
We provide a test data from GM12878 CTCF ChIA-PET ([GSM1872886](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1872886)), just the chromosome 21 mapped to hg38. Run the command as following then you will get the result if cLoops is successfuly installed. The ***eps*** is auto estimated and default ***minPts*** is 5,**-w** option will generate tracks and loops for visualization in [washU browser](http://epigenomegateway.wustl.edu/browser/).
```
cLoops -f GSM1872886_GSM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz -o chiapet -w 1
```    

1. HiChIP data
We provide a test data of from mESC cohesin HiChIP ([GSE80820](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80820),merged biological and technological replicates for 25 million cells), just the chromosome 19 mapped to mm10. Run the command as following.
```
cLoops -f GSE80820_mESC_cohesin_HiChIP_chr19_mm10.bedpe.gz -o hichip -eps 2000 -minPts 20 -twice 1 -w 1 -hichip 1
```

--------
## Questions & Answers  
Please address questions and bugs to Yaqiang Cao (caoyaqiang@picb.ac.cn) or Xingwei Chen (chenxingwei@picb.ac.cn) or Daosheng Ai (aidaosheng@picb.ac.cn), using the subject as "cLoops: questions about" to escape misjudged as spams.  
