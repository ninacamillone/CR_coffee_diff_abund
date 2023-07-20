#### load data if not already in environment ####
# rarefied data 
ps_CR_rare # ASVS
ASV_rare_filt # filtered ASVs
ps_genus_rare # genera
genus_rare_filt # filtered genera
# absolute abundance data (filtered and unfiltered)
ps_CR_abs # ASVS
ASV_abs_filt # filtered ASVs
ps_genus_abs # genera
genus_abs_filt # filtered genera

#### Wilcoxon on rarefied data ####

## ASVs (unfiltered)
wt.asv.shade <- apply(otu_table(ps_CR_rare), 2, 
                      function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.asv.shade <- data.frame(wt.asv.shade, p.adjust(wt.asv.shade, method = "holm"))
colnames(wt.asv.shade)<-c("p","p.holm")
sum(wt.asv.shade$p.holm<0.1) # check how many are significant

wt.asv.herb <- apply(otu_table(ps_CR_rare), 2, 
                     function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.asv.herb <- data.frame(wt.asv.herb, p.adjust(wt.asv.herb, method = "holm"))
colnames(wt.asv.herb)<-c("p","p.holm")
sum(wt.asv.herb$p.holm<0.1) # check how many are significant

## ASVs (filtered)
wt.asv.shade.filt <- apply(otu_table(ASV_rare_filt), 2, 
                           function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.asv.shade.filt <- data.frame(wt.asv.shade.filt, p.adjust(wt.asv.shade.filt, method = "holm"))
colnames(wt.asv.shade.filt)<-c("p","p.holm")
sum(wt.asv.shade.filt$p.holm<0.1) # check how many are significant

wt.asv.herb.filt <- apply(otu_table(ASV_rare_filt), 2, 
                          function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.asv.herb.filt <- data.frame(wt.asv.herb.filt, p.adjust(wt.asv.herb.filt, method = "holm"))
colnames(wt.asv.herb.filt)<-c("p","p.holm")
sum(wt.asv.herb.filt$p.holm<0.1) # check how many are significant

## Genera (unfiltered)
wt.gen.shade <- apply(otu_table(ps_genus_rare), 2, 
                      function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.gen.shade <- data.frame(wt.gen.shade, p.adjust(wt.gen.shade, method = "holm"))
colnames(wt.gen.shade)<-c("p","p.holm")
sum(wt.gen.shade$p.holm<0.1) # check how many are significant

wt.gen.herb <- apply(otu_table(ps_genus_rare), 2, 
                     function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.gen.herb <- data.frame(wt.gen.herb, p.adjust(wt.gen.herb, method = "holm"))
colnames(wt.gen.herb)<-c("p","p.holm")
sum(wt.gen.herb$p.holm<0.1) # check how many are significant

## Genera (filtered)
wt.gen.shade.filt <- apply(otu_table(genus_rare_filt), 2, 
                function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.gen.shade.filt <- data.frame(wt.gen.shade.filt, p.adjust(wt.gen.shade.filt, method = "holm"))
colnames(wt.gen.shade.filt)<-c("p","p.holm")
sum(wt.gen.shade.filt$p.holm<0.05) # check how many are significant

wt.gen.herb.filt <- apply(otu_table(genus_rare_filt), 2, 
                      function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.gen.herb.filt <- data.frame(wt.gen.herb.filt, p.adjust(wt.gen.herb.filt, method = "holm"))
colnames(wt.gen.herb.filt)<-c("p","p.holm")
sum(wt.gen.herb.filt$p.holm<0.05) # check how many are significant

#### Wilcoxon on absolute data ####

## ASVs (unfiltered)
wt.asv.abs.shade <- apply(otu_table(ps_CR_abs), 2, 
                      function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.asv.abs.shade <- data.frame(wt.asv.shade, p.adjust(wt.asv.shade, method = "holm"))
colnames(wt.asv.abs.shade)<-c("p","p.holm")
sum(wt.asv.abs.shade$p.holm<0.1) # check how many are significant

wt.asv.herb <- apply(otu_table(ps_CR_abs), 2, 
                     function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.asv.herb <- data.frame(wt.asv.herb, p.adjust(wt.asv.herb, method = "holm"))
colnames(wt.asv.herb)<-c("p","p.holm")
sum(wt.asv.herb$p.holm<0.1) # check how many are significant

## ASVs (filtered)
wt.asv.shade.abs.filt <- apply(otu_table(ASV_abs_filt), 2, 
                           function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.asv.shade.abs.filt <- data.frame(wt.asv.shade.filt, p.adjust(wt.asv.shade.filt, method = "holm"))
colnames(wt.asv.shade.abs.filt)<-c("p","p.holm")
sum(wt.asv.shade.abs.filt$p.holm<0.1) # check how many are significant

wt.asv.herb.abs.filt <- apply(otu_table(ASV_abs_filt), 2, 
                          function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.asv.herb.abs.filt <- data.frame(wt.asv.herb.filt, p.adjust(wt.asv.herb.filt, method = "holm"))
colnames(wt.asv.herb.abs.filt)<-c("p","p.holm")
sum(wt.asv.herb.abs.filt$p.holm<0.1) # check how many are significant

## Genera (unfiltered)
wt.gen.abs.shade <- apply(otu_table(ps_genus_abs), 2, 
                      function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.gen.abs.shade <- data.frame(wt.gen.abs.shade, p.adjust(wt.gen.abs.shade, method = "holm"))
colnames(wt.gen.abs.shade)<-c("p","p.holm")
sum(wt.gen.abs.shade$p.holm<0.1) # check how many are significant

wt.gen.abs.herb <- apply(otu_table(ps_genus_abs), 2, 
                     function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.gen.abs.herb <- data.frame(wt.gen.abs.herb, p.adjust(wt.gen.abs.herb, method = "holm"))
colnames(wt.gen.abs.herb)<-c("p","p.holm")
sum(wt.gen.abs.herb$p.holm<0.1) # check how many are significant

## Genera (filtered)
wt.gen.shade.abs.filt <- apply(otu_table(genus_abs_filt), 2, 
                           function(x) wilcox.test(x~sample_data(ps_CR)$ShadeSun)$p.value)
wt.gen.shade.abs.filt <- data.frame(wt.gen.shade.abs.filt, p.adjust(wt.gen.shade.abs.filt, method = "holm"))
colnames(wt.gen.shade.abs.filt)<-c("p","p.holm")
sum(wt.gen.shade.abs.filt$p.holm<0.1) # check how many are significant

wt.gen.herb.abs.filt <- apply(otu_table(genus_abs_filt), 2, 
                          function(x) wilcox.test(x~sample_data(ps_CR)$Herbicide)$p.value)
wt.gen.herb.abs.filt <- data.frame(wt.gen.herb.abs.filt, p.adjust(wt.gen.herb.abs.filt, method = "holm"))
colnames(wt.gen.herb.abs.filt)<-c("p","p.holm")
sum(wt.gen.herb.abs.filt$p.holm<0.1) # check how many are significant
