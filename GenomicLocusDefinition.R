asthma_sumstats <- read.table("/Users/eduardmasmarin/Desktop/Input.txt", header=T)
head(asthma_sumstats)

# Select the variants that are beyond the borders of the genomic loci and have P-value < 0.05

leadvar_GL1 <- asthma_sumstats[which(asthma_sumstats$position==13023741),]
leadvar_GL1

before_start_GL1 <- asthma_sumstats[which(asthma_sumstats$position<12997469 & asthma_sumstats$position>(12997469-250000) & asthma_sumstats$all_p.value<0.05),]
tail(before_start_GL1)

after_end_GL1 <- asthma_sumstats[which(asthma_sumstats$position>13024079 & asthma_sumstats$position<(13024079+250000) & asthma_sumstats$all_p.value<0.05),]
head(after_end_GL1)

before_start_GL2 <- asthma_sumstats[which(asthma_sumstats$position<49077072 & asthma_sumstats$position>(49077072-250000) & asthma_sumstats$all_p.value<0.05),]
tail(before_start_GL2)

after_end_GL2 <- asthma_sumstats[which(asthma_sumstats$position>49401542 & asthma_sumstats$position<(49401542+250000) & asthma_sumstats$all_p.value<0.05),]
head(after_end_GL2)

lead_vars <- asthma_sumstats[which(asthma_sumstats$all_p.value<1e-06),]
tail(lead_vars)
min(lead_vars$all_p.value)

lead_var_GL3 <- asthma_sumstats[which(asthma_sumstats$all_p.value==1.41e-10),]
lead_var_GL3   # 153,744,507  rs9887660
which(grepl(1.41e-10, asthma_sumstats$all_p.value))

before_start_GL3 <- asthma_sumstats[which(asthma_sumstats$position<153744507 & asthma_sumstats$position>(153744507-250000) & asthma_sumstats$all_p.value<0.05),]
tail(before_start_GL3)

after_end_GL3 <- asthma_sumstats[which(asthma_sumstats$position>153744507 & asthma_sumstats$position<(153744507+250000) & asthma_sumstats$all_p.value<0.05),]
head(after_end_GL3)

TMSB4X_rs850637 <- asthma_sumstats[which(asthma_sumstats$rsid=="rs850637"),]
TMSB4X_rs850637$all_p.value

PPP1R3F_rs5953283 <- asthma_sumstats[which(asthma_sumstats$rsid=="rs5953283"),]
PPP1R3F_rs5953283$all_p.value









