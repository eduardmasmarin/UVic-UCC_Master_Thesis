# Load data: GenomicRiskLoci and genes

GenomicRiskLoci <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/GenomicRiskLoci.txt",header=T)
head(GenomicRiskLoci)
genes <- read.table("/Users/eduardmas/Downloads/FUMA_job161043/genes.txt",header=T)
head(genes)
genes

GenomicRiskLoci_short <- GenomicRiskLoci[,c("GenomicLocus","rsID","chr","pos","start","end","nSNPs")]
head(GenomicRiskLoci_short)
tail(GenomicRiskLoci_short)

genes_short <- genes[,c("ensg","symbol","chr","start","end","GenomicLocus")]
head(genes_short)
tail(genes_short)

# Check that chromosomes of GenomicLocus are coincident in GenomicRiskLoci and genes

GenomicRiskLoci_genes <- merge(GenomicRiskLoci_short,genes_short,by="GenomicLocus",all.x=T,all.y=F)
head(GenomicRiskLoci_genes)
tail(GenomicRiskLoci_genes)

GenomicRiskLoci_genes$chr_match <- GenomicRiskLoci_genes$chr.x == GenomicRiskLoci_genes$chr.y
table(GenomicRiskLoci_genes$chr_match)

# Check overlapping between genes and GenomicRiskLoci

GenomicRiskLoci_genes$overlap_start <- GenomicRiskLoci_genes$start.y >= GenomicRiskLoci_genes$start.x & GenomicRiskLoci_genes$start.y <= GenomicRiskLoci_genes$end.x 
GenomicRiskLoci_genes$overlap_end <-  GenomicRiskLoci_genes$end.y >= GenomicRiskLoci_genes$start.x & GenomicRiskLoci_genes$end.y <= GenomicRiskLoci_genes$end.y
table(GenomicRiskLoci_genes$overlap_start)
table(GenomicRiskLoci_genes$overlap_end)

GenomicRiskLoci_genes$overlap <- (GenomicRiskLoci_genes$start.y >= GenomicRiskLoci_genes$start.x & GenomicRiskLoci_genes$start.y <= GenomicRiskLoci_genes$end.x) | (GenomicRiskLoci_genes$end.y >= GenomicRiskLoci_genes$start.x & GenomicRiskLoci_genes$end.y <= GenomicRiskLoci_genes$end.y)
table(GenomicRiskLoci_genes$overlap)

# For genes which do not overlap, determine their distance from GenomicLocus

GenomicRiskLoci_genes$distance_start <- GenomicRiskLoci_genes$start.x - GenomicRiskLoci_genes$end.y
GenomicRiskLoci_genes$distance_end <- GenomicRiskLoci_genes$start.y - GenomicRiskLoci_genes$end.x
head(GenomicRiskLoci_genes)

genes_not_overlap <- GenomicRiskLoci_genes[which(GenomicRiskLoci_genes$overlap==FALSE),]
head(genes_not_overlap)
tail(genes_not_overlap)
genes_not_overlap

# For genes which do not overlap, check if they are within 10 kb distance from either the start or the end of the GenomicLocus

genes_within_10kb <- genes_not_overlap$distance_start <= 10000 | genes_not_overlap$distance_end <= 10000
table(genes_within_10kb)