library(phyloseq)
library(selbal)

ps_CR <- readRDS(file = "ps_CR_clean.rds")
ps_genus <- readRDS(file = "ps_genus.rds")
CR_filt_asvs <- readRDS(file = "CR_filt_asvs.RDS")
CR_filt_genera <- readRDS(file = "CR_filt_gen.RDS")

otu_table(ps_CR) <- otu_table(ps_CR)[,order(colSums(otu_table(ps_CR)), decreasing = T)]
taxa_names(ps_CR) <- 1:6873

#x=microbiome data matrix
y <- factor(sample_data(ps_CR)$ShadeSun) #y=response variable vector
z <- factor(sample_data(ps_CR)$Herbicide) #covar= cofounding variables dataframe 

#ASVs (raw and filtered are the same, due to removal of ASVs with >80% zeros)
x <- prune_taxa(colSums(otu_table(ps_CR)==0)<(24*.8), ps_CR) #424 taxa
bal.shade.asvs <- selbal.cv(x = otu_table(x), y = y, covar = data.frame(z))
bal.herb.asvs <- selbal.cv(x = otu_table(x), y = z, covar = data.frame(y))

saveRDS(bal.herb.asvs, "bal.herb.asvs.all.RDS")
saveRDS(bal.shade.asvs, "bal.shade.asvs.all.RDS")

#Genera (raw)
x <- prune_taxa(colSums(otu_table(ps_genus)==0)<(24*.8), ps_genus)
bal.shade.gen <- selbal.cv(x  = otu_table(x), y = y, covar= data.frame(z))
bal.herb.gen <- selbal.cv(x = otu_table(x), y = z, covar = data.frame(y))

#Genera (filtered)
bal.shade.gen.filt <- selbal.cv(x = otu_table(CR_filt_genera), y = y,covar = data.frame(z), n.iter=20) # default: n.iter=10
bal.herb.gen.filt <- selbal.cv(x = otu_table(CR_filt_genera),y = z, covar = data.frame(y), n.iter=20)
