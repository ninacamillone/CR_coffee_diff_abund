library(phyloseq)
library(ANCOMBC)

ps_CR <- readRDS(file = "ps_CR_clean.rds")
ps_genus <- readRDS(file = "ps_genus.rds")
CR_filt_asvs <- readRDS(file = "CR_filt_asvs.RDS")
CR_filt_genera <- readRDS(file = "CR_filt_gen.RDS")

#ASVs (raw)
set.seed(123)
out_asv <- ancombc2(data = ps_CR, assay_name = NULL,
                    fix_formula = "ShadeSun + Herbicide", p_adj_method = "holm",
                    prv_cut = 0.10,
                    group = NULL, struc_zero = FALSE, verbose = TRUE,
                    global = FALSE, pairwise = FALSE, 
                    dunnet = FALSE, trend = FALSE,
                    iter_control = list(tol = 1e-5, max_iter = 20, 
                                        verbose = FALSE),
                    em_control = list(tol = 1e-5, max_iter = 100),
                    lme_control = NULL)

#ASVs (filtered)
set.seed(123)
out_asv_filt <- ancombc2(data = CR_filt_asvs, assay_name = NULL,
                         fix_formula = "ShadeSun + Herbicide", p_adj_method = "holm",
                         prv_cut = 0.10,
                         group = NULL, struc_zero = FALSE, verbose = TRUE,
                         global = FALSE, pairwise = FALSE, 
                         dunnet = FALSE, trend = FALSE,
                         iter_control = list(tol = 1e-5, max_iter = 20, 
                                             verbose = FALSE),
                         em_control = list(tol = 1e-5, max_iter = 100),
                         lme_control = NULL)

set.seed(123)
out_genera_filt <- ancombc2(data = CR_filt_genera, assay_name = NULL, tax_level = NULL,
                            fix_formula = "ShadeSun + Herbicide", p_adj_method = "holm",
                            prv_cut = 0.10,group = NULL, struc_zero = FALSE, verbose = TRUE,
                            global = FALSE, pairwise = FALSE, 
                            dunnet = FALSE, trend = FALSE,
                            iter_control = list(tol = 1e-5, max_iter = 20,verbose = FALSE),
                            em_control = list(tol = 1e-5, max_iter = 100),
                            lme_control = NULL)
res_genera_filt <- out_genera_filt$res
res_genera_filt <- res_genera_filt[res_genera_filt$diff_ShadeSunSun|res_genera_filt$diff_HerbicideReduced,]
