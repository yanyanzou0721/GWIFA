library(circlize)
args <- commandArgs(trailingOnly = TRUE)
input_xls = args[1] ## input file
interaction_threshold = as.integer(args[2]) ## int
out_name=argv[3]

af = read.table(input_xls,header = 1)  ## interaction information between the focal amplified region and whole Genome
af = subset(af,count>interaction_threshold)

af$chrom1 = paste0("chr",af$chrom1)
af$chrom2 = paste0("chr",af$chrom2)
af_1 = af[,c("chrom1","start1","end1","count")]
af_2 = af[,c("chrom2","start2","end2")]
names(af_1)=c("chr","start","end","value")
names(af_2)=c("chr","start","end")

pdf(paste0(out_name,".pdf"))
par(mar = c(1, 1, 1, 1))
circos.par("start.degree" = 90)
circos.par("track.height" = 0.1)
colors <- colorRampPalette(c("skyblue","blue"))(max(af_1["value"]))
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", 1:22),"chrX"), plotType = c("ideogram", "labels"), ideogram.height = 0.03)
circos.genomicLink(af_1,af_2,col=colors[c(as.vector(af_1$value))],lwd=0.5)
title(out_name)
circos.clear()
dev.off()
