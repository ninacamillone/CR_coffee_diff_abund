library(phyloseq)
library(PERFect)

ps_CR <- readRDS(file = "ps_CR_clean.rds")
ps_genus <- readRDS(file = "ps_genus.rds")

# ASVs: simultaneous method
res_sim <- PERFect_sim(otu_table(ps_CR))
CR_filt_asvs <- phyloseq(otu_table(res_sim$filtX),sample_data(ps_CR))

# Genera: permutational method
CR_filt_genera <- ps_genus
# simultaneous method (preliminary)
res_sim_genera <- PERFect_sim(otu_table(CR_filt_genera))
pvals  <- c(1, round(res_sim_genera$pvals, 2))
pvals.sort  <- sort.int(pvals, decreasing = TRUE, index.return = TRUE)
otu_table(CR_filt_genera) <- otu_table(CR_filt_genera)[, pvals.sort$ix]

# permutational method
res_perm_genera <- PERFect_perm(X = otu_table(CR_filt_genera), Order = "pvals", 
                          pvals_sim = res_sim_genera, algorithm = "full")
resSummary(X = otu_table(CR_filt_genera), filtX = res_perm_genera$filtX)

CR_filt_genera <- phyloseq(otu_table(res_perm_genera$filtX),sample_data(ps_CR))

# Save outputs
saveRDS(CR_filt_genera, file = "CR_filt_gen.RDS")
saveRDS(CR_filt_asvs, file = "CR_filt_asvs.RDS")
