cLoops -f GSM1872886_GM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz -o chiapet -w 1 -j 1 -s 1
#cLoops -f GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio1.bedpe.gz,GSE80820_GM12878_cohesin_HiChIP_chr21_hg38_bio2.bedpe.gz -o hichip -eps 1000,2000,4000,6000,8000,10000 -minPts 50 -s 1 -w 1 -j 1 -hic 1
#jd2juice -d hichip -o cohesin_chr21 -org hg38
#jd2washU -d hichip -o cohesin_chr21
