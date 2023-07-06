library(DESeq2)

## Load data
ps_CR <- readRDS(file = "ps_CR_clean.rds")
ps_genus <- readRDS(file = "ps_genus.rds")
CR_filt_asvs <- readRDS(file = "CR_filt_asvs.RDS")
CR_filt_genera <- readRDS(file = "CR_filt_gen.RDS")

## ASVs (raw)
ds_asv <- phyloseq_to_deseq2(ps_CR, ~ ShadeSun*Herbicide)
ds_asv <- DESeq(ds_asv, quiet = T)
res_ds_asv <- results(ds_asv)
res_ds_asv <- res_ds_asv[order(res_ds_asv$pvalue),]
sum(res_ds_asv$padj < 0.1, na.rm=TRUE) 
res_ds_asv_shade <- results(ds_asv, contrast=c("ShadeSun", "Sun", "Shade"))
res_ds_asv_shade <- res_ds_asv_shade[order(res_ds_asv_shade$pvalue),]
sum(res_ds_asv_shade$padj < 0.1, na.rm=TRUE) 
res_ds_asv_Herbicide <- results(ds_asv, contrast=c("Herbicide", "Conventional", "Reduced"))
res_ds_asv_Herbicide <- res_ds_asv_Herbicide[order(res_ds_asv_Herbicide$pvalue),]
sum(res_ds_asv_Herbicide$padj < 0.1, na.rm=TRUE) 

## ASVs (filtered)
ds_asv_filt <- phyloseq_to_deseq2(CR_filt_asvs, ~ ShadeSun*Herbicide)
ds_asv_filt <- DESeq(ds_asv_filt, quiet = T)
res_ds_asv_filt <- results(ds_asv_filt)
res_ds_asv_filt <- res_ds_asv_filt[order(res_ds_asv_filt$pvalue),]
sum(res_ds_asv_filt$padj < 0.1, na.rm=TRUE) 
res_ds_asv_filt_shade <- results(ds_asv_filt, contrast=c("ShadeSun", "Sun", "Shade"))
res_ds_asv_filt_shade <- res_ds_asv_filt_shade[order(res_ds_asv_filt_shade$pvalue),]
sum(res_ds_asv_filt_shade$padj < 0.1, na.rm=TRUE) 
res_ds_asv_filt_Herbicide <- results(ds_asv_filt, contrast=c("Herbicide", "Conventional", "Reduced"))
res_ds_asv_filt_Herbicide <- res_ds_asv_filt_Herbicide[order(res_ds_asv_filt_Herbicide$pvalue),]
sum(res_ds_asv_filt_Herbicide$padj < 0.1, na.rm=TRUE) 

## Genera (raw)
ds_genus <- phyloseq_to_deseq2(ps_genus, ~ ShadeSun*Herbicide)
ds_genus <- DESeq(ds_genus, quiet = T)
res_ds_genus <- results(ds_genus)
res_ds_genus <- res_ds_genus[order(res_ds_genus$pvalue),]
sum(res_ds_genus$padj < 0.1, na.rm=TRUE) 
res_ds_genus_shade <- results(ds_genus, contrast=c("ShadeSun", "Sun", "Shade"))
res_ds_genus_shade <- res_ds_genus_shade[order(res_ds_genus_shade$pvalue),]
sum(res_ds_genus_shade$padj < 0.1, na.rm=TRUE) 
row.names(res_ds_genus_shade)[1:9]
res_ds_genus_Herbicide <- results(ds_genus, contrast=c("Herbicide", "Conventional", "Reduced"))
res_ds_genus_Herbicide <- res_ds_genus_Herbicide[order(res_ds_genus_Herbicide$pvalue),]
sum(res_ds_genus_Herbicide$padj < 0.1, na.rm=TRUE) 

## Genera (filtered)
ds_genus_filt <- phyloseq_to_deseq2(CR_filt_genera, ~ ShadeSun*Herbicide)
ds_genus_filt <- DESeq(ds_genus_filt, quiet = T)
res_ds_genus_filt <- results(ds_genus_filt)
res_ds_genus_filt <- res_ds_genus_filt[order(res_ds_genus_filt$pvalue),]
sum(res_ds_genus_filt$padj < 0.1, na.rm=TRUE) # 5
res_ds_genus_filt_shade <- results(ds_genus_filt, contrast=c("ShadeSun", "Sun", "Shade"))
res_ds_genus_filt_shade <- res_ds_genus_filt_shade[order(res_ds_genus_filt_shade$pvalue),]
sum(res_ds_genus_filt_shade$padj < 0.1, na.rm=TRUE) # 6
row.names(res_ds_genus_filt_shade)[1:6]
res_ds_genus_filt_Herbicide <- results(ds_genus_filt, contrast=c("Herbicide", "Conventional", "Reduced"))
res_ds_genus_filt_Herbicide <- res_ds_genus_filt_Herbicide[order(res_ds_genus_filt_Herbicide$pvalue),]
sum(res_ds_genus_filt_Herbicide$padj < 0.1, na.rm=TRUE)

save.image("DESeq2_diff_abund.RData")
