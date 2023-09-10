# Load data: snps and genes

snps <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/snps.txt",header=T)
head(snps)
genes <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/genes.txt",header=T)
head(genes)

snps_short <- snps[,c("rsID","chr","pos","gwasP","or","r2","GenomicLocus","nearestGene","dist")]
head(snps_short)
tail(snps_short)

genes_short <- genes[,c("ensg","symbol","chr","start","end","minGwasP","GenomicLocus")]
head(genes_short)
tail(genes_short)

# Check that chromosomes of SNPs are coincident in snps and genes

snps_genes <- merge(snps_short,genes_short,by="GenomicLocus",all.x=T,all.y=F)
head(snps_genes)
tail(snps_genes)

snps_genes$chr_match <- snps_genes$chr.x == snps_genes$chr.y
table(snps_genes$chr_match)

# Check whether the SNP is contained within mapped genes

snps_genes$SNP_match <- snps_genes$pos >= snps_genes$start & snps_genes$pos <= snps_genes$end
table(snps_genes$SNP_match)
 
# Check whether the SNP is contained within 10 kb from either start or end of mapped genes

snps_genes$SNP_10kb_match <- snps_genes$pos >= (snps_genes$start-10000) & snps_genes$pos <= (snps_genes$end+10000)
table(snps_genes$SNP_10kb_match)

# Check whether either the start or the end of the gene is within 10kb distance from the SNP (20 kb + the SNP)

snps_genes$gene_10kb_match <- snps_genes$start <= (snps_genes$pos+10000) | snps_genes$end >= (snps_genes$pos-10000)
table(snps_genes$gene_10kb_match)
### Genes are mapped from total candidate SNPs… because genes are mapped from Genomic Locus and the Genomic Locus is built from candidate SNPs extreme positions.


# Separate the SNPs which are contained within either genes or a 10 kb distance from either gene start or end

SNP_mapped_genes <- snps_genes[which(snps_genes$SNP_match==TRUE),]
head(SNP_mapped_genes)
tail(SNP_mapped_genes)
SNP_mapped_genes

SNP_mapped_genes_10kb <- snps_genes[which(snps_genes$SNP_10kb_match==TRUE),]
head(SNP_mapped_genes_10kb)
tail(SNP_mapped_genes_10kb)
SNP_mapped_genes_10kb
