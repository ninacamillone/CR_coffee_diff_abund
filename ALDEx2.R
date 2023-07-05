library(ALDEx2)

ps_CR <- readRDS(file = "ps_CR_clean.rds")
ps_genus <- readRDS(file = "ps_genus.rds")
CR_filt_asvs <- readRDS(file = "CR_filt_asvs.RDS")
CR_filt_genera <- readRDS(file = "CR_filt_gen.RDS")

#create model
mm <- data.frame(sample_data(ps_CR)$ShadeSun,sample_data(ps_CR)$Herbicide)
colnames(mm) <- c("Shade", "Herbicide")
mm <- model.matrix(~ Shade*Herbicide, mm)

#ASVs (raw)
x.asv <- aldex.clr(t(otu_table(ps_CR)), mm)
x.asv.glm.test <- aldex.glm(x.asv, mm)

#ASVs (filtered)
x.asv.filt <- aldex.clr(t(otu_table(CR_filt_asvs)), mm)
x.asv.filt.glm.test <- aldex.glm(x.asv.filt, mm)

#Genera (raw)
x.gen <- aldex.clr(t(otu_table(ps_genus)), mm)
x.gen.glm.test <- aldex.glm(x.gen, mm)

#Genera (filtered)
x.gen.filt <- aldex.clr(t(otu_table(CR_filt_genera)), mm)
x.gen.filt.glm.test <- aldex.glm(x.gen.filt, mm)
