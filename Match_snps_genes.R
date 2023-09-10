# load the data

genes=read.table("bsc05476@mn1.bsc.es:/gpfs/projects/bsc05/eduardmas/FUMA/FUMA_jobs/FUMA_job135670/genes.txt", header = T)
head(genes)
snps=read.table("bsc05476@mn1.bsc.es:/gpfs/projects/bsc05/eduardmas/FUMA/FUMA_jobs/FUMA_job135670/snps.txt", header = T)
head(snps)

# make a shorter version of the snps dataframe

snps.short=snps[,c("GenomicLocus","chr","pos")]
head(snps.short)

# merge snps and genes using the genomic locus

genes.snps=merge(genes, snps.short, by= c("GenomicLocus","chr"), all=TRUE)

# check if the chr of genes and snps match

genes.snps$match.chr=genes.snps$chr.x %in% genes.snps$chr.y
table(genes.snps$match.chr)

# check if the snps and the genes overlap

genes.snps$overlap=0 
genes.snps$overlap=ifelse([genes.snps$start < genes.snps$pos & genes.snps$end > genes.snps$pos],1,0)
table(genes.snps$overlap)

# for genes that do not overlap, check the distance from the locus

genes.snps$delta.start=0
genes.snps$delta.end=0
genes.snps$delta.start=genes.snps$start - genes.snps$pos
genes.snps$delta.end=genes.snps$end - genes.snps$pos
genes.snps$delta=0
genes.snps$delta=ifelse([genes.snps$delta.start < 10000 | genes.snps$delta.end < 10000],1,0)
genes.snps$list=unique(genes.snps$ensg) # to have a unique row for each gene
genes.snps$list[is.na(genes.snps$list)] <- 0 # eliminate NA rows and change them for 0
genes.mapping <- data.frame("genes"=genes.snps$list) # create a data frame with unique genes
rownames(genes.mapping) <- c("g"=genes.snps$ensg) # change row names to ensg ids
for (g in genes.mapping){
	genes.mapping$overlap <- genes.snps$overlap
	genes.mapping$delta <- genes.snps$delta
}
genes.mapping$mapped.genes=ifelse([genes.mapping$overlap >= 1 | genes.mapping$delta >= 1],1,0)
print(sum(genes.mapping$mapped.genes))







