# load the data (adjust the path if necessary!)

genes=read.table("bsc05476@mn1.bsc.es:/gpfs/projects/bsc05/eduardmas/FUMA/FUMA_jobs/FUMA_job136397/genes_136397_PM.txt", header = T)
head(genes)
variant=read.table("bsc05476@mn1.bsc.es:/gpfs/projects/bsc05/eduardmas/FUMA/FUMA_jobs/FUMA_job136397/IndSigSNPs.txt", header = T) # ¡¡¡ "p" column must be previously renamed as "minGwasP" !!!
head(variant)

# make a shorter version of the variant dataframe

variant.short=variant[,c(,"chr","pos","minGwasP")]
head(variant.short)


# merge variant and genes using the P-value

genes.variant=merge(genes, variant.short, by= "minGwasP", all.x = T, all.y = F)


# check if the chr of genes and variant match

genes.variant$match.chr=genes.variant$chr.x %in% genes.variant$chr.y
table(genes.variant$match.chr) # 146/146 TRUE


# check if the variant and the genes overlap

genes.variant$overlap=0
genes.variant[genes.variant$start.x < genes.variant$pos.y & genes.variant$end.x > genes.variant$pos.y,"overlap"]=1
table(genes.variant$overlap) # 141/146 overlap


# for genes that do not overlap, check the distance from the locus

genes.variant$delta.start=genes.variant$start.x- genes.variant$pos.y
genes.variant$delta.end=genes.variant$pos.y- genes.variant$end.x
genes.variant$out.range=TRUE
genes.variant[genes.variant$overlap==1, "out.range"]=FALSE
table(genes.variant$out.range)
genes.variant[genes.variant$delta.start <= 10000, "out.range"]=FALSE
table(genes.variant$out.range) # (141 + 3)/146 overlap or are upstream of < 10000 bp
genes.variant[genes.variant$delta.end <= 10000, "out.range"]=FALSE
table(genes.variant$out.range) # (141 + 3 + 2)/146 overlap or are upstream or downstream of < 10000 bp


# classify genes

genes.variant$match.type=0
genes.variant[genes.variant$overlap==1, "match.type"]="overlap"
genes.variant[genes.variant$match.type==0 & genes.variant$delta.start <=10000, "match.type"]="upstream_10kb"
genes.variant[genes.variant$match.type==0 & genes.variant$delta.end <=10000, "match.type"]="downstream_10kb"
table(genes.variant$match.type) # final classification