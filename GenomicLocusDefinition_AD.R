atopic.dermatitis_sumstats <- read.table("/Users/eduardmas/Library/Mobile Documents/com~apple~CloudDocs/Documents/BSC/FUMA/FUMA_GWAS_summary_statistics/AD_interview_meta_fixed_F-M_gc_EU_PanUKBB_chr1-23_raw_ALL.FUMA.txt", header=T)
head(atopic.dermatitis_sumstats)

# Select the variants that are beyond the borders of the genomic loci and have P-value < 0.05

# Genomic Locus #01

leadvar_GL1 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position==152285861),]
leadvar_GL1

before_start_GL1 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position<150050001 & atopic.dermatitis_sumstats$position>(150050001-250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
tail(before_start_GL1)

after_end_GL1 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position>154428283 & asthma_sumstats$position<(154428283+250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
head(after_end_GL1)


# Genomic Locus #12

leadvar_GL12 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position==6094697),]
leadvar_GL12

before_start_GL12 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position<6032661 & atopic.dermatitis_sumstats$position>(6032661-250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
tail(before_start_GL12)

after_end_GL12 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position>6196839 & asthma_sumstats$position<(6196839+250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
head(after_end_GL12)


# Genomic Locus #16

leadvar_GL16 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position==76293527),]
leadvar_GL16

before_start_GL16 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position<75276801 & atopic.dermatitis_sumstats$position>(75276801-250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
tail(before_start_GL16)

after_end_GL16 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position>77072589 & asthma_sumstats$position<(77072589+250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
head(after_end_GL16)


# Genomic Locus #19

leadvar_GL19 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position==11223454),]
leadvar_GL19

before_start_GL19 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position<11031741 & atopic.dermatitis_sumstats$position>(11031741-250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
tail(before_start_GL19)

after_end_GL19 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position>11285678 & asthma_sumstats$position<(11285678+250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
head(after_end_GL19)


# Genomic Locus #21

leadvar_GL21 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position==62306855),]
leadvar_GL21

before_start_GL21 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position<62240322 & atopic.dermatitis_sumstats$position>(62240322-250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
tail(before_start_GL21)

after_end_GL21 <- atopic.dermatitis_sumstats[which(atopic.dermatitis_sumstats$position>62401720 & asthma_sumstats$position<(62401720+250000) & atopic.dermatitis_sumstats$all_p.value<0.05),]
head(after_end_GL21)