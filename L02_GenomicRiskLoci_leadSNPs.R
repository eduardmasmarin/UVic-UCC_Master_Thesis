# Load data: GenomicRiskLoci and leadSNPs

GenomicRiskLoci <- read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135659/GenomicRiskLoci.txt",header=T)
head(GenomicRiskLoci)
leadSNPs <- read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135659/leadSNPs.txt",header=T)
head(leadSNPs)

GenomicRiskLoci_short <- GenomicRiskLoci[,c("GenomicLocus","rsID","chr","pos","p","start","end","nSNPs")]
head(GenomicRiskLoci_short)
tail(GenomicRiskLoci_short)

leadSNPs_short <- leadSNPs[,c("No","GenomicLocus","rsID","chr","pos","p","nIndSigSNPs")]
head(leadSNPs_short)
tail(leadSNPs_short)

# Check that chromosomes of GenomicLocus are coincident in GenomicRiskLoci and leadSNPs

GenomicRiskLoci_leadSNPs <- merge(GenomicRiskLoci_short,leadSNPs_short,by="GenomicLocus",all.x=T,all.y=F)
head(GenomicRiskLoci_leadSNPs)
tail(GenomicRiskLoci_leadSNPs)

GenomicRiskLoci_leadSNPs$chr_match <- GenomicRiskLoci_leadSNPs$chr.x == GenomicRiskLoci_leadSNPs$chr.y
table(GenomicRiskLoci_leadSNPs$chr_match)

# Check which leadSNPs assigned to a GenomicLocus are contained within the GenomicLocus

GenomicRiskLoci_leadSNPs$leadSNP_in_GenomicLocus <- GenomicRiskLoci_leadSNPs$pos.y >= GenomicRiskLoci_leadSNPs$start & GenomicRiskLoci_leadSNPs$pos.y <= GenomicRiskLoci_leadSNPs$end 
table(GenomicRiskLoci_leadSNPs$leadSNP_in_GenomicLocus)

# Among each leadSNPs of every GenomicLocus, determine how many of them coincide with either the start or the end of the GenomicLocus

GenomicRiskLoci_leadSNPs$either_start_or_end <- GenomicRiskLoci_leadSNPs$pos.y == GenomicRiskLoci_leadSNPs$start | GenomicRiskLoci_leadSNPs$pos.y == GenomicRiskLoci_leadSNPs$end
table(GenomicRiskLoci_leadSNPs$either_start_or_end)

GenomicRiskLoci_leadSNPs$leadSNP_start <- GenomicRiskLoci_leadSNPs$pos.y == GenomicRiskLoci_leadSNPs$start
table(GenomicRiskLoci_leadSNPs$leadSNP_start)

GenomicRiskLoci_leadSNPs$leadSNP_end <- GenomicRiskLoci_leadSNPs$pos.y == GenomicRiskLoci_leadSNPs$end
table(GenomicRiskLoci_leadSNPs$leadSNP_end)

head(GenomicRiskLoci_leadSNPs)
tail(GenomicRiskLoci_leadSNPs)

leadSNP_start <- GenomicRiskLoci_leadSNPs[which(GenomicRiskLoci_leadSNPs$leadSNP_start==TRUE),]
head(leadSNP_start)
tail(leadSNP_start)
leadSNP_start

leadSNP_end <- GenomicRiskLoci_leadSNPs[which(GenomicRiskLoci_leadSNPs$leadSNP_end==TRUE),]
head(leadSNP_end)
tail(leadSNP_end)
leadSNP_end

# Determine which leadSNPs coincide with the GenomicRiskLoci in both P-value and position

GenomicRiskLoci_leadSNPs$same_position <- GenomicRiskLoci_leadSNPs$pos.x == GenomicRiskLoci_leadSNPs$pos.y
table(GenomicRiskLoci_leadSNPs$same_position)

GenomicRiskLoci_leadSNPs$same_pvalue <- GenomicRiskLoci_leadSNPs$p.x == GenomicRiskLoci_leadSNPs$p.y
table(GenomicRiskLoci_leadSNPs$same_pvalue)

GenomicRiskLoci_leadSNPs$same_position_pvalue <- (GenomicRiskLoci_leadSNPs$pos.x == GenomicRiskLoci_leadSNPs$pos.y) & (GenomicRiskLoci_leadSNPs$p.x == GenomicRiskLoci_leadSNPs$p.y)
table(GenomicRiskLoci_leadSNPs$same_position_pvalue)

# Check what happens with the remaining leadSNPs (55-21=34) which do not constitute a GenomicLocus

leadSNPs_different_position_and_pvalue <- GenomicRiskLoci_leadSNPs[which(!((GenomicRiskLoci_leadSNPs$pos.x == GenomicRiskLoci_leadSNPs$pos.y) & (GenomicRiskLoci_leadSNPs$p.x == GenomicRiskLoci_leadSNPs$p.y))),]
head(leadSNPs_different_position_and_pvalue)
tail(leadSNPs_different_position_and_pvalue)
leadSNPs_different_position_and_pvalue
leadSNPs_different_position_and_pvalue$rsID.y
leadSNPs_different_position_and_pvalue$rsID.x

leadSNPs_different_position_and_pvalue$rsID_match <- leadSNPs_different_position_and_pvalue$rsID.x == leadSNPs_different_position_and_pvalue$rsID.y
table(leadSNPs_different_position_and_pvalue$rsID_match)

leadSNPs_not_GenomicLocus <- leadSNPs_different_position_and_pvalue[,c("GenomicLocus","rsID.x","pos.x","p.x","rsID.y","pos.y","p.y")]
leadSNPs_not_GenomicLocus

leadSNPs_not_GenomicLocus$distance <- abs(leadSNPs_not_GenomicLocus$pos.y - leadSNPs_not_GenomicLocus$pos.x)
