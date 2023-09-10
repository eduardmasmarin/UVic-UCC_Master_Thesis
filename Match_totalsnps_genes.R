library("dplyr")
library("data.table")

# load the data

genes=read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135670/genes.txt", header = T)
head(genes)
snps=read.table("/Users/eduardmasmarin/Documents/BSC/BSC_COMPUTATIONAL_GENOMICS/BSC_COMPLEX_GENETIC_DISEASES/FUMA/FUMA_jobs/FUMA_job135670/snps.txt", header = T)
head(snps)

# make a shorter version of the snps dataframe

genes_short=genes[,c("GenomicLocus","chr","start","end")]
snps_short=snps[,c("GenomicLocus","chr","pos")]
head(genes_short)
head(snps_short)

# merge snps and genes using the genomic locus

genes_short$delta_start = genes_short$start - 10000
genes_short$delta_end = genes_short$end + 10000
genes_short_1 = genes_short %>% filter(genes_short$GenomicLocus == 1)
head(genes_short_1)
tail(genes_short_1)
snps_short_1 = snps_short %>% filter(snps_short$GenomicLocus == 1)
head(snps_short_1)
tail(snps_short_1)



for (i in 1:length(snps_short_1$pos)){
  if (genes_short$GenomicLocus == snps_short_1$GenomicLocus){
  (genes_short$pos_snps = snps_short_1$pos[i]) & (genes_short$GenomicLocus == genes_short$GenomicLocus_snps & genes_short$pos_snps > genes_short$delta_start & genes_short$pos_snps < genes_short$delta_end)
  }
}
head(genes_short$GenomicLocus)

warnings()



# if (genes_short$GenomicLocus == snps_short$GenomicLocus)
# for (genes_short$pos_snps in snps_short_1$pos){
# genes_short$pos_snps = 150050001
# genes_short$pos_snps = snps_short_1$pos

# head(genes_short)
# genes_short$GenomicLocus == genes_short$GenomicLocus_snps & genes_short$pos_snps > genes_short$delta_start & genes_short$pos_snps < genes_short$delta_end 
# head(genes_short$GenomicLocus)


# variant_pass = snps_short_1$pos > genes_short_1$delta_start & snps_short_1$pos < genes_short_1$delta_end


