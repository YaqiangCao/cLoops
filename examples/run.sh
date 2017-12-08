cLoops -f GSM1872886_GM12878_CTCF_ChIA-PET_chr21_hg38.bedpe.gz -o chiapet -w 1 -j 1 -s 1
jd2juice -d chiapet -o CTCF_chr21 -org hg38
jd2washU -d chiapet -o CTCF_chr21
jd2fingerprint -d chiapet -plot 1 -o compare -bs 2000
rm -fvr chiapet/
