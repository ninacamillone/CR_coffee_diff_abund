#### set up ####
ps_CR <- readRDS(file = "ps_CR_clean.rds")

library(phyloseq)

#### Rarefy ####
## rarefying ASVs
ps_CR_rare <- rarefy_even_depth(ps_CR, rngseed = 123)
sample_sums(ps_CR_rare) # check library size
rarePlot <- rarecurve(as.data.frame(otu_table(ps_CR)), step = 100, cex=0.5)
abline(v=2989)
abline(v=9612)
#one sample maxes out species much more quickly than the others, so we'll remove it
ps_CR_rare <- rarefy_even_depth(prune_samples(sample_sums(ps_CR)>2989, ps_CR),rngseed = 123)
rowSums(otu_table(ps_CR_rare)) # check final library size

## rarefied genera
# group rarefied ASVs by genus
ps_genus_rare <- tax_glom(ps_CR_rare, taxrank = "Genus") # 313 genera
taxa_names(ps_genus_rare) <- tax_table(ps_genus_rare)[,6]
saveRDS(ps_genus_rare, file = "ps_genus_rare.RDS")

#### Estimate absolute abundance ####

#ASVs
ps_CR_abs <- transform_sample_counts(ps_CR, function(x) x/sum(x)) # convert to relative abundance
otu_table(ps_CR_abs) <- otu_table(ps_CR_abs)*sample_data(ps_CR)$qPCR #multiply by qPCR counts
sample_sums(ps_CR_abs)==sample_data(ps_CR)$qPCR # check all true
saveRDS(ps_CR_abs, file="ps_CR_abs.rds")

#genera
genus_abs <- tax_glom(ps_CR_abs, taxrank = "Genus")
taxa_names(genus_abs) <- as.data.frame(tax_table(genus_abs))$Genus
saveRDS(genus_abs, file="ps_genus_abs.rds")
