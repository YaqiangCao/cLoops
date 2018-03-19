# cLoops: a clustering based loops calling method for ChIA-PET, HiChIP and Hi-C
![](https://github.com/YaqiangCao/cLoops/raw/master/cLoops.png)

## Introduction
Chromosome conformation capture (3C) derived high-throughput sequencing methods such as ChIA-PET,HiChIP and Hi-C provide genome-wide view of chromatin organization. Fine scale loops formed by interactions of regulatory elements spanning hundreds kilobases can be detected from these data. Here we introduce cLoops ('see loops'),a common loops calling tool for ChIA-PET, HiChIP and high-resolution Hi-C data. Paired-end tags (PETs) are first classified as self-ligation and inter-ligation clusters using an optimized unsupervisied clustering algorithm called cDBSCAN. The significances of the inter-ligation clusters are then estimated using permutated local background. Both steps are data type independent, and thus enable cLoops to be applicable to even new genome-wide interaction mapping to be developed in the future.

Basic workflow of cLoops is as following:
![](https://github.com/YaqiangCao/cLoops/raw/master/workflow.png)

If you find cLoops is useful, please cite our paper:    
### cLoops: a clustering based loops calling method for ChIA-PET, HiChIP and Hi-C ###

--------
## Install
[scipy](https://www.scipy.org/),[numpy](http://www.numpy.org/), [seaborn](https://seaborn.pydata.org/), [pandas](http://pandas.pydata.org/) and [joblib](https://pythonhosted.org/joblib/) are required. If you have problems for installing scipy, please refer to [Anaconda](https://docs.continuum.io/anaconda/) or [SAGE](http://www.sagemath.org/).
```
wget https://github.com/YaqiangCao/cLoops/archive/0.5.tar.gz -O cLoops.tar.gz
tar xvzf cLoops.tar.gz
cd cLoops-0.5
python setup.py install    
```

To test whether cLoops is successfully installed:
```
cd examples
sh run.sh
```


Please refer to [here](https://docs.python.org/2/install/index.html) to install cLoops to customized path.

--------
## Usage
Run ***cLoops -h*** to see all options. Key parameters are ***eps*** and ***minPts*** . ***minPts*** defines at least how many PETs are required for a candidate loop, ***eps*** defines the distance requried for two PETs being neighbors. For ChIA-PET data with sharp peaks, cLoops can auto estimate ***eps*** from the data as 2 fold of the fragment size, and ***minPts***=5 is empirical good. For ChIA-PET data with broad peaks (like H3K4me1), empirical experience is set ***eps*** to 2000,5000. For HiChIP, set a series ***eps***=2000,4000,6000,8000,10000 & ***minPts***=20 worth a first trial, if deeply sequenced, increase ***minPts*** to 30 or 50. For practically usage, using the PETs in the smallest chromosome except chrY and chrM, then run a series of ***eps***, choose the smallest ***eps*** that can get well seperated inter-ligation and self-ligation PETs distance distributions. Or just apply a series of ***eps***,all rounds clustering result will be combined. Actually, there will be more significant loops given more ***eps***, but it will take longer time.
Since version 0.8, cLoops added a new parameter ***--mode(-m)***, which is the pre-set parameters for different types of data. -m 0 accepts user settings; -m 1 equals -eps 500,1000,2000 -minPts 5 for sharp peak like ChIA-PET data; -m 2 equals -eps 1000,2000,5000 -minPts 5 for broad peak like ChIA-PET data; -m 3 equals -eps 2000,5000,10000 -minPts 30 -hic 1 for HiChIP and Hi-C data. 

--------
### Input  
Mapped PETs in [BEDPE format](http://bedtools.readthedocs.io/en/latest/content/general-usage.html), compressed files with gzip are also accepected, following columns are necessary: chrom1 (1st),start1 (2),end1 (3),chrom2 (4),start2 (5),end2 (6),strand1 (9),strand2 (10). For the column of name or score, "." is accepcted.

--------
### Output
The main output is a loop file and a PDF file or PDFs for the plot of self-ligation and inter-ligation PETs distance distributions.
For the .loop file, columns and explaination are as follwing:

column | name | explaination
------ | ---- | ------------
0th | loopId | Id for a loop, like chr1-chr1-1
1th | ES | Enrichment score for the loop, caculated by observed PETs number divided by the mean PETs number of nearby permutated  regions
2th | FDR | false discovery rate for the loop, caculated as the number of permutated regions that there are more observed PETs than the region  
3th | binomal\_p-value | binomal test p-value for the loop
4th | distance | distance (bp) between the centers of the anchors for the loop
5th | hypergeometric\_p-value | hypergeometric test p-value for the loop
6th | iva | genomic coordinates for the left anchor, for example, chr13:50943050-50973634
7th | ivb | genomic coordinates for the right anchor
8th | poisson_p-value | poisson test p-value for the loop
9th | ra | observed PETs number for the left anchor
10th | rab | observed PETs number linking the left and right anchors
11th | rb | observed PETs number for the right anchor
12th | poisson\_p-value\_corrected | Bonferroni corrected poisson p-value according to number of loops for each chromosome
13th | binomal\_p-value\_corrected | Bonferroni corrected binomal p-value according to number of loops for each chromosome
14th | hypergeometric\_p-value\_corrected | Bonferroni corrected hypergeometric p-value according to number of loops for each chromosome
15th | significant | 1 or 0, 1 means we think the loop is significant compared to permutated regions. For ChIA-PET data, significant requiring ES >=1.0, FDR <=0.05, hypergeometric\_local\_FDR <=0.05 and all uncorrected p-values <= 1e-5; For HiChIP and high resolution Hi-C data, significant requiring ES >= 2.0, FDR <=0.05, both binomal\_p-value and poisson\_p-value <=1e-5. You can ignore this and customize your cutoffs by visualization in the [Juicebox](https://github.com/theaidenlab/juicebox) or [washU](http://epigenomegateway.wustl.edu/) to determine your cutoffs.

--------
## Examples
All following examples source data, result and log file can be found in the [examples](https://github.com/YaqiangCao/cLoops_supplementaryData/tree/master/examples).

### 1. ChIA-PET data    
We provide a test data from GM12878 CTCF ChIA-PET ([GSM1872886](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1872886)), just the chromosome 21 mapped to hg38. Run the command as following then you will get the result if cLoops is successfuly installed. The ***eps*** is auto estimated and default ***minPts*** is 5,**-w** option will generate loops for visualization in [washU browser](http://epigenomegateway.wustl.edu/browser/),**-j** option will generate loops for visualization in [Juicebox](https://github.com/theaidenlab/juicebox) .
```
wget https://github.com/YaqiangCao/cLoops/blob/master/examples/GSM1872886_GM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz
cLoops -f GSM1872886_GM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz -o chiapet -w 1 -j 1
```      
For ChIA-PET data with sharp peak, like the CTCF here, you will get the inter-ligation and self-ligation PETs distance distribution like following, the two kinds of PETs well seperated:
![](https://github.com/YaqiangCao/cLoops/raw/master/pngs/chiapet_disCutoff.png)

If your experimental data doesn't look like this by auto estimated ***eps***, which could be true for some ChIA-PET data with broad peak (like H3K27ac), please use the small chromosome (chr21 in human and chr19 in mouse) run a series of ***eps***, then chose the smallest one that generate the well seperated distance distribution to run cLoops, or just using the series. 

We recommend washU to visualize the loops, by the script jd2washU we can convert the cLoops temp files to washU long range track, and [bedtools](http://bedtools.readthedocs.io/en/latest/),[bgzip & tabix](http://www.htslib.org/doc/tabix.html) are needed in the command enviroment. 
```
jd2washU -d chiapet -o chiapet       
``` 

With other ChIP-seq data, you can get following plot:
![](https://github.com/YaqiangCao/cLoops/raw/master/pngs/chiapet_washU.png)

### 2. HiChIP data   
We provide test data of GM12878 cohesin HiChIP two biological replicates, just the chromosome 21 mapped to hg38. Run the command as following to call merged loops. ***-s*** option is used to keep working directory and temp files, which could be used by scripts of deLoops, jd2washU (BEDTOOLS needed), jd2juice (Juicer needed), jd2fingerprint and jd2saturation. ***-hic*** option means using cutoffs design for Hi-C like data, see above. 
```
wget https://github.com/YaqiangCao/cLoops_supplementaryData/blob/master/examples/GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio1.bedpe.gz 
wget https://github.com/YaqiangCao/cLoops_supplementaryData/blob/master/examples/GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio2.bedpe.gz 
cLoops -f GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio1.bedpe.gz,GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio2.bedpe.gz -o hichip -eps 1000,2000,4000,6000,8000,10000 -minPts 50 -s 1 -hic 1 -w 1 -j 1
```    
Then use jd2juice to convert cLoops temp files to hic file for juicebox:
```
jd2juice -d hichip -o hichip -org hg38 
``` 
With the adjustment of resolution, color range and how to show the loops, then you can get following visualization:
![](https://github.com/YaqiangCao/cLoops/raw/master/pngs/hichip_juicebox_example.png) 


### 3. Hi-C data   
We provide test data from GM12878 Hi-C, just the chromosome 21 mapped to hg38. Run the the command as following to call loops.
```
wget https://github.com/YaqiangCao/cLoops_supplementaryData/blob/master/examples/GSM1551552_GM12878_HiC_chr21_hg38.bedpe.gz 
cLoops -f GSM1551552_GM12878_HiC_chr21_hg38.bedpe.gz -o hic -w 1 -j 1 -s 1 -eps 2000,4000,6000,8000,10000 -minPts 30 -s 1 -hic 1
```   

### 4. Fingerprint plot for data qualities comparasion of loops calling 
Run following and you will get [a PDF plot](https://github.com/YaqiangCao/cLoops_supplementaryData/blob/master/examples/compare_fingerprint.pdf), the far from the random line, the better for the data used to call loops by cLoops.
```
jd2fingerprint -d chiapet,hichip,hic -plot 1 -o compare -bs 2000
```
--------
## Other data  
In theory cLoops could be applied to more 3D genomic data as long as there are enriched clusters in the heatmap, however, parameters and significance cutoffs should be tuned. We're now trying to make cLoops work for [GRID-seq](https://www.nature.com/articles/nbt.3968) and [Capture HiC](https://www.nature.com/articles/ng.3286). If you have designed a new sequencing based 3D genomic method and want to try cLoops, please contact caoyaqiang@picb.ac.cn first.


--------
## Questions & Answers  
Please address questions and bugs to Yaqiang Cao (caoyaqiang@picb.ac.cn) or Xingwei Chen (chenxingwei@picb.ac.cn) or Daosheng Ai (aidaosheng@picb.ac.cn), using the subject as "cLoops: questions about" to escape misjudged as spams.  
