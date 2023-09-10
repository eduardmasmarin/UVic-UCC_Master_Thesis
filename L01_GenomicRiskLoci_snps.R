# Load data: GenomicRiskLoci and snps

GenomicRiskLoci <- read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135659/GenomicRiskLoci.txt",header=T)
head(GenomicRiskLoci)
snps <- read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135659/snps.txt",header=T)
head(snps)

GenomicRiskLoci_short <- GenomicRiskLoci[,c("GenomicLocus","rsID","chr","pos","start","end","nSNPs")]
head(GenomicRiskLoci_short)
tail(GenomicRiskLoci_short)

snps_short <- snps[,c("rsID","chr","pos","GenomicLocus")]
head(snps_short)
tail(snps_short)

# Check that chromosomes of GenomicLocus are coincident in GenomicRiskLoci and snps

GenomicRiskLoci_snps <- merge(GenomicRiskLoci_short,snps_short,by="GenomicLocus",all.x=T,all.y=F)
head(GenomicRiskLoci_snps)
tail(GenomicRiskLoci_snps)

GenomicRiskLoci_snps$chr_match <- GenomicRiskLoci_snps$chr.x == GenomicRiskLoci_snps$chr.y
table(GenomicRiskLoci_snps$chr_match)

# Check that all snps assigned to a GenomicLocus are contained within the GenomicLocus

GenomicRiskLoci_snps$snp_in_GenomicLocus <- GenomicRiskLoci_snps$pos.y >= GenomicRiskLoci_snps$start & GenomicRiskLoci_snps$pos.y <= GenomicRiskLoci_snps$end 
table(GenomicRiskLoci_snps$snp_in_GenomicLocus)

# Among all snps of a GenomicLocus, determine how many of them coincide with either the start or the end of the GenomicLocus

GenomicRiskLoci_snps$either_start_or_end <- GenomicRiskLoci_snps$pos.y == GenomicRiskLoci_snps$start | GenomicRiskLoci_snps$pos.y == GenomicRiskLoci_snps$end
table(GenomicRiskLoci_snps$either_start_or_end)

GenomicRiskLoci_snps$snp_start <- GenomicRiskLoci_snps$pos.y == GenomicRiskLoci_snps$start
table(GenomicRiskLoci_snps$snp_start)

GenomicRiskLoci_snps$snp_end <- GenomicRiskLoci_snps$pos.y == GenomicRiskLoci_snps$end
table(GenomicRiskLoci_snps$snp_end)

head(GenomicRiskLoci_snps)
tail(GenomicRiskLoci_snps)

snp_start <- GenomicRiskLoci_snps[which(GenomicRiskLoci_snps$snp_start==TRUE),]
head(snp_start)
tail(snp_start)
snp_start

snp_end <- GenomicRiskLoci_snps[which(GenomicRiskLoci_snps$snp_end==TRUE),]
head(snp_end)
tail(snp_end)
snp_end