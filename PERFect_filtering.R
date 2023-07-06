library(phyloseq)
library(PERFect)

ps_CR <- readRDS(file = "ps_CR_clean.rds")
ps_genus <- readRDS(file = "ps_genus.rds")
ps_CR_abs <- readRDS(file = "ps_CR_abs")
ps_genus_abs <- readRDS(file = "ps_genus_abs.rds")
ps_CR_rare <- readRDS(file = "ps_CR_rare.RDS")
ps_genus_rare <- readRDS(file = "ps_genus_rare.RDS")

## ASVs: simultaneous method
res_sim <- PERFect_sim(otu_table(ps_CR))
CR_filt_asvs <- phyloseq(otu_table(res_sim$filtX),sample_data(ps_CR))

## Genera: permutational method
CR_filt_genera <- ps_genus
# simultaneous method (preliminary)
res_sim_genera <- PERFect_sim(otu_table(CR_filt_genera))
pvals  <- c(1, round(res_sim_genera$pvals, 2))
pvals.sort  <- sort.int(pvals, decreasing = TRUE, index.return = TRUE)
otu_table(CR_filt_genera) <- otu_table(CR_filt_genera)[, pvals.sort$ix]
# permutational method
res_perm_genera <- PERFect_perm(X = otu_table(CR_filt_genera), Order = "pvals", 
                          pvals_sim = res_sim_genera, algorithm = "full")
CR_filt_genera <- phyloseq(otu_table(res_perm_genera$filtX),sample_data(ps_CR))

## absolute ASVs 
ASV_abs_filt <- ps_CR_abs
res_sim_abs <- PERFect_sim(otu_table(ASV_abs_filt))
ASV_abs_filt  <- phyloseq(otu_table(ASV_abs_filt),sample_data(ps_CR))

## absolute genera 
genus_abs_filt <- ps_genus_abs
# simultaneous method
res_sim_genera_abs <- PERFect_sim(otu_table(genus_abs_filt))
pvals  <- c(1, round(res_sim_genera_abs$pvals, 2))
pvals.sort  <- sort.int(pvals, decreasing = TRUE, index.return = TRUE)
otu_table(genus_abs_filt) <- otu_table(genus_abs_filt)[, pvals.sort$ix]
# permutational method
res_perm_genera_abs <- PERFect_perm(X = otu_table(genus_abs_filt), Order = "pvals", 
                                     pvals_sim = res_sim_genera_abs, algorithm = "full")
genus_abs_filt <- res_perm2_genera_abs$filtX
genus_abs_filt  <- phyloseq(otu_table(genus_abs_filt),sample_data(ps_CR))

## rarefied ASVs 
ASV_rare_filt <- ps_CR_rare 
res_sim_rare <- PERFect_sim(otu_table(ASV_rare_filt))
ASV_rare_filt <- res_sim_rare$filtX
ASV_rare_filt  <- phyloseq(otu_table(ASV_rare_filt),sample_data(ps_CR))

## rarefied genera 
genus_rare_filt <- ps_genus_rare
# simultaneous method
res_sim_genera_rare <- PERFect_sim(otu_table(genus_rare_filt))
dim(res_sim_genera_rare$filtX) # 22
pvals  <- c(1, round(res_sim_genera_rare$pvals, 2))
pvals.sort  <- sort.int(pvals, decreasing = TRUE, index.return = TRUE)
otu_table(genus_rare_filt) <- otu_table(genus_rare_filt)[, pvals.sort$ix]
# permutational method
res_perm_genera_rare <- PERFect_perm(X = otu_table(genus_rare_filt), Order = "pvals", 
                                    pvals_sim = res_sim_genera_rare, algorithm = "full")
genus_rare_filt <- res_perm_genera_rare$filtX
genus_rare_filt  <- phyloseq(otu_table(genus_rare_filt),sample_data(ps_CR))

# Save outputs
saveRDS(CR_filt_genera, file = "CR_filt_gen.RDS")
saveRDS(CR_filt_asvs, file = "CR_filt_asvs.RDS")
saveRDS(ASV_abs_filt, file = "ASV_abs_filt.RDS")
saveRDS(genus_abs_filt, file = "genus_abs_filt.RDS")
saveRDS(ASV_abs_filt, file = "ASV_rare_filt.RDS")
saveRDS(genus_abs_filt, file = "genus_rare_filt.RDS")
