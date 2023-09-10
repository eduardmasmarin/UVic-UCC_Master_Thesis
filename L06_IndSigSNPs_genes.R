# Load data: IndSigSNPs and genes

IndSigSNPs <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/IndSigSNPs.txt",header=T)
head(IndSigSNPs)
genes <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/genes.txt",header=T)
head(genes)

IndSigSNPs_short <- IndSigSNPs[,c("No","GenomicLocus","rsID","chr","pos","p")]
head(IndSigSNPs_short)
tail(IndSigSNPs_short)

genes_short <- genes[,c("ensg","symbol","chr","start","end","minGwasP","GenomicLocus")]
head(genes_short)
tail(genes_short)

# Check that chromosomes of IndSigSNPs are coincident in IndSigSNPs and genes

IndSigSNPs_genes <- merge(IndSigSNPs_short,genes_short,by="GenomicLocus",all.x=T,all.y=F)
head(IndSigSNPs_genes)
tail(IndSigSNPs_genes)

IndSigSNPs_genes$chr_match <- IndSigSNPs_genes$chr.x == IndSigSNPs_genes$chr.y
table(IndSigSNPs_genes$chr_match)

# Check whether the IndSigSNP is contained within mapped genes

IndSigSNPs_genes$IndSigSNP_match <- IndSigSNPs_genes$pos >= IndSigSNPs_genes$start & IndSigSNPs_genes$pos <= IndSigSNPs_genes$end
table(IndSigSNPs_genes$IndSigSNP_match)

# Check whether the IndSigSNP is contained within 10 kb from either start or end of mapped genes

IndSigSNPs_genes$IndSigSNP_10kb_match <- IndSigSNPs_genes$pos >= (IndSigSNPs_genes$start-10000) & IndSigSNPs_genes$pos <= (IndSigSNPs_genes$end+10000)
table(IndSigSNPs_genes$IndSigSNP_10kb_match)
#### Genes are not mapped by positional mapping considering IndSigSNPs.

# Separate the IndSigSNPs which are contained within either genes or a 10 kb distance from either gene start or end

IndSigSNP_mapped_genes <- IndSigSNPs_genes[which(IndSigSNPs_genes$IndSigSNP_match==TRUE),]
head(IndSigSNP_mapped_genes)
tail(IndSigSNP_mapped_genes)
IndSigSNP_mapped_genes

IndSigSNP_mapped_genes_10kb <- IndSigSNPs_genes[which(IndSigSNPs_genes$IndSigSNP_10kb_match==TRUE),]
head(IndSigSNP_mapped_genes_10kb)
tail(IndSigSNP_mapped_genes_10kb)
IndSigSNP_mapped_genes_10kb

