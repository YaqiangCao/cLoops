library( Sushi )
bg = read.table("GM12878_CTCF.bedgraph")
bg = data.frame(bg) 
chrom = "chr21"
chromstart = 29195000 
chromend = 29198000
#chromstart = 29025000 
#chromend = 29420000
plotBedgraph(bg,chrom,chromstart,chromend,transparency=.80,color=SushiColors(2)(2)[2])
labelgenome(chrom,chromstart,chromend,n=4,scale="Mb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
