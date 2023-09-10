# Load data: GenomicRiskLoci and IndSigSNPs

GenomicRiskLoci <- read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135659/GenomicRiskLoci.txt",header=T)
head(GenomicRiskLoci)
IndSigSNPs <- read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135659/IndSigSNPs.txt",header=T)
head(IndSigSNPs)

GenomicRiskLoci_short <- GenomicRiskLoci[,c("GenomicLocus","rsID","chr","pos","p","start","end","nSNPs")]
head(GenomicRiskLoci_short)
tail(GenomicRiskLoci_short)

IndSigSNPs_short <- IndSigSNPs[,c("No","GenomicLocus","rsID","chr","pos","p")]
head(IndSigSNPs_short)
tail(IndSigSNPs_short)

# Check that chromosomes of GenomicLocus are coincident in GenomicRiskLoci and IndSigSNPs

GenomicRiskLoci_IndSigSNPs <- merge(GenomicRiskLoci_short,IndSigSNPs_short,by="GenomicLocus",all.x=T,all.y=F)
head(GenomicRiskLoci_IndSigSNPs)
tail(GenomicRiskLoci_IndSigSNPs)

GenomicRiskLoci_IndSigSNPs$chr_match <- GenomicRiskLoci_IndSigSNPs$chr.x == GenomicRiskLoci_IndSigSNPs$chr.y
table(GenomicRiskLoci_IndSigSNPs$chr_match)

# Check which IndSigSNPs assigned to a GenomicLocus are contained within the GenomicLocus

GenomicRiskLoci_IndSigSNPs$IndSigSNP_in_GenomicLocus <- GenomicRiskLoci_IndSigSNPs$pos.y >= GenomicRiskLoci_IndSigSNPs$start & GenomicRiskLoci_IndSigSNPs$pos.y <= GenomicRiskLoci_IndSigSNPs$end 
table(GenomicRiskLoci_IndSigSNPs$IndSigSNP_in_GenomicLocus)

# Among each IndSigSNPs of every GenomicLocus, determine how many of them coincide with either the start or the end of the GenomicLocus

GenomicRiskLoci_IndSigSNPs$either_start_or_end <- GenomicRiskLoci_IndSigSNPs$pos.y == GenomicRiskLoci_IndSigSNPs$start | GenomicRiskLoci_IndSigSNPs$pos.y == GenomicRiskLoci_IndSigSNPs$end
table(GenomicRiskLoci_IndSigSNPs$either_start_or_end)

GenomicRiskLoci_IndSigSNPs$IndSigSNP_start <- GenomicRiskLoci_IndSigSNPs$pos.y == GenomicRiskLoci_IndSigSNPs$start
table(GenomicRiskLoci_IndSigSNPs$IndSigSNP_start)

GenomicRiskLoci_IndSigSNPs$IndSigSNP_end <- GenomicRiskLoci_IndSigSNPs$pos.y == GenomicRiskLoci_IndSigSNPs$end
table(GenomicRiskLoci_IndSigSNPs$IndSigSNP_end)

head(GenomicRiskLoci_IndSigSNPs)
tail(GenomicRiskLoci_IndSigSNPs)

IndSigSNP_start <- GenomicRiskLoci_IndSigSNPs[which(GenomicRiskLoci_IndSigSNPs$IndSigSNP_start==TRUE),]
head(IndSigSNP_start)
tail(IndSigSNP_start)
IndSigSNP_start

IndSigSNP_end <- GenomicRiskLoci_IndSigSNPs[which(GenomicRiskLoci_IndSigSNPs$IndSigSNP_end==TRUE),]
head(IndSigSNP_end)
tail(IndSigSNP_end)
IndSigSNP_end

# Determine which IndSigSNPs coincide with the GenomicRiskLoci in both P-value and position

GenomicRiskLoci_IndSigSNPs$same_position <- GenomicRiskLoci_IndSigSNPs$pos.x == GenomicRiskLoci_IndSigSNPs$pos.y
table(GenomicRiskLoci_IndSigSNPs$same_position)

GenomicRiskLoci_IndSigSNPs$same_pvalue <- GenomicRiskLoci_IndSigSNPs$p.x == GenomicRiskLoci_IndSigSNPs$p.y
table(GenomicRiskLoci_IndSigSNPs$same_pvalue)

GenomicRiskLoci_IndSigSNPs$same_position_pvalue <- (GenomicRiskLoci_IndSigSNPs$pos.x == GenomicRiskLoci_IndSigSNPs$pos.y) & (GenomicRiskLoci_IndSigSNPs$p.x == GenomicRiskLoci_IndSigSNPs$p.y)
table(GenomicRiskLoci_IndSigSNPs$same_position_pvalue)

# Check what happens with the remaining IndSigSNPs that are not leadSNPs

IndSigSNPs_different_position_and_pvalue <- GenomicRiskLoci_IndSigSNPs[which(!((GenomicRiskLoci_IndSigSNPs$pos.x == GenomicRiskLoci_IndSigSNPs$pos.y) & (GenomicRiskLoci_IndSigSNPs$p.x == GenomicRiskLoci_IndSigSNPs$p.y))),]
head(IndSigSNPs_different_position_and_pvalue)
tail(IndSigSNPs_different_position_and_pvalue)
IndSigSNPs_different_position_and_pvalue
IndSigSNPs_different_position_and_pvalue$rsID.y
IndSigSNPs_different_position_and_pvalue$rsID.x

IndSigSNPs_different_position_and_pvalue$rsID_match <- IndSigSNPs_different_position_and_pvalue$rsID.x == IndSigSNPs_different_position_and_pvalue$rsID.y
table(IndSigSNPs_different_position_and_pvalue$rsID_match)

IndSigSNPs_not_GenomicLocus <- IndSigSNPs_different_position_and_pvalue[,c("GenomicLocus","rsID.x","pos.x","p.x","rsID.y","pos.y","p.y")]
IndSigSNPs_not_GenomicLocus

IndSigSNPs_not_GenomicLocus$distance <- abs(IndSigSNPs_not_GenomicLocus$pos.y - IndSigSNPs_not_GenomicLocus$pos.x)

