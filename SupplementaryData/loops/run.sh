chroms="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr10,chr21,chr22,chrX,chrY,chrM"

#chiapet
cLoops -f HeLa_CTCF.bedpe.gz -o HeLa_CTCF -p 1 -c $chroms -j 1 -w 1 -minPts 3
cLoops -f GM12878_CTCF.bedpe.gz -o GM12878_CTCF -p 1 -c $chroms -j 1 -w 1
cLoops -f GM12878_RAD21.rmdup.bedpe.gz -o GM12878_RAD21 -p 1 -c $chroms -w 1 -j 1
cLoops -f RAD21.rmdup.bedpe.gz -o RAD21 -p 1 -c $chroms -w 1 -j 1
cLoops -f H3K27ac.rmdup.bedpe.gz -o H3K27ac -p 1 -c $chroms -w 1 -j 1 -eps 2000,5000
cLoops -f H3K4me1.rmdup.bedpe.gz -o H3K4me1 -p 1 -c $chroms -w 1 -j 1 -eps 2000,5000
cLoops -f H3K4me3.rmdup.bedpe.gz -o H3K4me3 -p 1 -c $chroms -w 1 -j 1 -eps 2000,5000
cLoops -f Naive_SMC1_rep1.rmdup.bedpe.gz,Naive_SMC1_rep2.rmdup.bedpe.gz -o Naive_SMC1 -p 1 -c $chroms -w 1 -j 1 -eps 1000 
cLoops -f Primed_SMC1_rep1.rmdup.bedpe.gz,Primed_SMC1_rep2.rmdup.bedpe.gz -o Primed_SMC1 -p 1 -c $chroms -w 1 -j 1 -eps 1000 

#hichip 
cLoops -f cohesin_GM.bedpe.gz -o cohesin_GM -p 1 -c $chroms -eps 2000,4000,6000,8000,10000 -minPts 30 -hic 1 -j 1 -w 1
cLoops -f oct4_mESC_25m.bedpe.gz -o oct4_mESC -p 1 -c $chroms -eps 2000,4000,6000,8000,10000 -minPts 50 -hic 1 -w 1 -j 1
cLoops -f cohesin_mESC_all.bedpe.gz -o cohesin_mESC -p 1 -c $chroms -eps 2000,4000,6000,8000,10000 -minPts 50 -hic 1 -w 1 -j 1

#hic 
cLoops -f GM12878_GSM1551552_Cis.bedpe.gz -o GM12878_hic -eps 2000,4000,6000,8000,10000 -minPts 30 -c $chroms -p 1 -j 1 -hic 1 -w 1 
cLoops -f K562_GSM1551619_Cis.bedpe.gz -o K562_hic -eps 2000,4000,6000,8000,10000 -minPts 30 -c $chroms -p 1 -j 1 -hic 1 -w 1 
