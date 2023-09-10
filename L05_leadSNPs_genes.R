# Load data: leadSNPs and genes

leadSNPs <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/leadSNPs.txt",header=T)
head(leadSNPs)
genes <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/genes.txt",header=T)
head(genes)

leadSNPs_short <- leadSNPs[,c("No","GenomicLocus","rsID","chr","pos","p","nIndSigSNPs")]
head(leadSNPs_short)
tail(leadSNPs_short)

genes_short <- genes[,c("ensg","symbol","chr","start","end","minGwasP","GenomicLocus")]
head(genes_short)
tail(genes_short)

# Check that chromosomes of leadSNPs are coincident in leadSNPs and genes

leadSNPs_genes <- merge(leadSNPs_short,genes_short,by="GenomicLocus",all.x=T,all.y=F)
head(leadSNPs_genes)
tail(leadSNPs_genes)

leadSNPs_genes$chr_match <- leadSNPs_genes$chr.x == leadSNPs_genes$chr.y
table(leadSNPs_genes$chr_match)

# Check whether the leadSNP is contained within mapped genes

leadSNPs_genes$leadSNP_match <- leadSNPs_genes$pos >= leadSNPs_genes$start & leadSNPs_genes$pos <= leadSNPs_genes$end
table(leadSNPs_genes$leadSNP_match)

# Check whether the leadSNP is contained within 10 kb from either start or end of mapped genes

leadSNPs_genes$leadSNP_10kb_match <- leadSNPs_genes$pos >= (leadSNPs_genes$start-10000) & leadSNPs_genes$pos <= (leadSNPs_genes$end+10000)
table(leadSNPs_genes$leadSNP_10kb_match)
#### Genes are not mapped by positional mapping considering leadSNPs.


# Separate the leadSNPs which are contained within either genes or a 10 kb distance from either gene start or end

leadSNP_mapped_genes <- leadSNPs_genes[which(leadSNPs_genes$leadSNP_match==TRUE),]
head(leadSNP_mapped_genes)
tail(leadSNP_mapped_genes)
leadSNP_mapped_genes

leadSNP_mapped_genes_10kb <- leadSNPs_genes[which(leadSNPs_genes$leadSNP_10kb_match==TRUE),]
head(leadSNP_mapped_genes_10kb)
tail(leadSNP_mapped_genes_10kb)
leadSNP_mapped_genes_10kb