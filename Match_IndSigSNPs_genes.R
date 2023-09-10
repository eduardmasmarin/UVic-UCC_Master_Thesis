# load the data (adjust the path if necessary!)

genes=read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135670/genes.txt", header = T)
head(genes)
IndSigSNPs=read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135670/IndSigSNPs.txt", header = T)
head(IndSigSNPs)

# make a shorter version of the IndSigSNPs dataframe

IndSigSNPs.short=IndSigSNPs[,c("GenomicLocus","chr","pos")]
head(IndSigSNPs.short)


# merge IndSigSNPs and genes using the genomic locus

genes.IndSigSNPs=merge(genes, IndSigSNPs.short, by= "GenomicLocus", all=TRUE)


# check if the chr of genes and IndSigSNPs match

genes.IndSigSNPs$match.chr=genes.IndSigSNPs$chr.x %in% genes.IndSigSNPs$chr.y
table(genes.IndSigSNPs$match.chr) # 146/146 TRUE


# check if the IndSigSNPs and the genes overlap

genes.IndSigSNPs$overlap=0
genes.IndSigSNPs[genes.IndSigSNPs$start < genes.IndSigSNPs$pos & genes.IndSigSNPs$end > genes.IndSigSNPs$pos,"overlap"]=1
table(genes.IndSigSNPs$overlap) # 141/146 overlap


# for genes that do not overlap, check the distance from the locus

genes.IndSigSNPs$delta.start=genes.IndSigSNPs$end- genes.IndSigSNPs$pos
genes.IndSigSNPs$delta.end=genes.IndSigSNPs$start- genes.IndSigSNPs$pos
genes.IndSigSNPs$out.range=TRUE
genes.IndSigSNPs[genes.IndSigSNPs$overlap==1, "out.range"]=FALSE
table(genes.IndSigSNPs$out.range)
genes.IndSigSNPs[abs(genes.IndSigSNPs$delta.start) <= 10000, "out.range"]=FALSE
table(genes.IndSigSNPs$out.range) # (141 + 3)/146 overlap or are upstream of < 10000 bp
genes.IndSigSNPs[abs(genes.IndSigSNPs$delta.end) <= 10000, "out.range"]=FALSE
table(genes.IndSigSNPs$out.range) # (141 + 3 + 2)/146 overlap or are upstream or downstream of < 10000 bp


# classify genes

genes.IndSigSNPs$match.type=0
genes.IndSigSNPs[genes.IndSigSNPs$overlap==1, "match.type"]="overlap"
genes.IndSigSNPs[genes.IndSigSNPs$match.type==0 & abs(genes.IndSigSNPs$delta.start) <=10000, "match.type"]="upstream_10kb"
genes.IndSigSNPs[genes.IndSigSNPs$match.type==0 & abs(genes.IndSigSNPs$delta.end) <=10000, "match.type"]="downstream_10kb"
table(genes.IndSigSNPs$match.type) # final classification

