# load the data (adjust the path if necessary!)

genes=read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_job136401/genes_136401.txt", header = T)
head(genes)
loci=read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_job136401/GenomicRiskLoci.txt", header = T)
head(loci)

# make a shorter version of the loci dataframe

loci.short=loci[,c("GenomicLocus","chr","start","end")]
head(loci.short)


# merge loci and genes using the genomic locus

genes.loci=merge(genes, loci.short, by= "GenomicLocus", all.x = T, all.y = F)


# check if the chr of genes and loci match

genes.loci$match.chr=genes.loci$chr.x %in% genes.loci$chr.y
table(genes.loci$match.chr) # 146/146 TRUE


# check if the loci and the genes overlap

genes.loci$overlap=0
genes.loci[genes.loci$start.x < genes.loci$start.y & genes.loci$end.x > genes.loci$start.y,"overlap"]=1
genes.loci[genes.loci$start.x < genes.loci$end.y & genes.loci$end.x > genes.loci$start.y,"overlap"]=1
table(genes.loci$overlap) # 141/146 overlap


# for genes that do not overlap, check the distance from the locus

genes.loci$delta.start=genes.loci$end.x- genes.loci$start.y
genes.loci$delta.end=genes.loci$start.x- genes.loci$end.y
genes.loci$out.range=TRUE
genes.loci[genes.loci$overlap==1, "out.range"]=FALSE
table(genes.loci$out.range)
genes.loci[abs(genes.loci$delta.start) <= 10000, "out.range"]=FALSE
table(genes.loci$out.range) # (141 + 3)/146 overlap or are upstream of < 10000 bp
genes.loci[abs(genes.loci$delta.end) <= 10000, "out.range"]=FALSE
table(genes.loci$out.range) # (141 + 3 + 2)/146 overlap or are upstream or downstream of < 10000 bp


# classify genes

genes.loci$match.type=0
genes.loci[genes.loci$overlap==1, "match.type"]="overlap"
genes.loci[genes.loci$match.type==0 & abs(genes.loci$delta.start) <=10000, "match.type"]="upstream_10kb"
genes.loci[genes.loci$match.type==0 & abs(genes.loci$delta.end) <=10000, "match.type"]="downstream_10kb"
table(genes.loci$match.type) # final classification

