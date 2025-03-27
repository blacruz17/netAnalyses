library(curatedMetagenomicData)
library(tidyverse)
library(phyloseq)
library(NetCoMi)
library(mia)
library(genefilter)

rm(list = ls())

# GET PHYLOSEQ OBJECTS #########################################################
#### KarlssonFH_2013 ###########################################################
physeq.karlsson <- sampleMetadata %>%
  filter(study_name == "KarlssonFH_2013")  %>%
  returnSamples("relative_abundance", counts = FALSE)
physeq.karlsson <- makePhyloseqFromTreeSE(physeq.karlsson,
                                     assay.type = "relative_abundance")

(physeq.k.ctl <- subset_samples(physeq.karlsson, disease == "healthy"))
(physeq.k.igt <- subset_samples(physeq.karlsson, disease == "IGT"))
(physeq.k.t2d <- subset_samples(physeq.karlsson, disease == "T2D"))

rm(physeq_karlsson)

# Filtrado 0.01% en 10% de las muestras
(physeq.k.ctl.f <- filter_taxa(physeq.k.ctl,
                               flist =  filterfun(kOverA(k = 4, A = 0.01)),
                               prune = TRUE)) 
saveRDS(physeq.k.ctl.f, "physeq_k_ctl_f.RDS")

(physeq.k.igt.f <- filter_taxa(physeq.k.igt,
                                 flist =  filterfun(kOverA(k = 5, A = 0.01)),
                                 prune = TRUE)) 
saveRDS(physeq.k.igt.f, "physeq_k_igt_f.RDS")

(physeq.k.t2d.f <- filter_taxa(physeq.k.t2d,
                               flist =  filterfun(kOverA(k = 5, A = 0.01)),
                               prune = TRUE)) 
saveRDS(physeq.k.t2d.f, "physeq_k_t2d_f.RDS")

#### QinJ_2012: T2D ############################################################
physeq.qinj <- sampleMetadata %>%
  filter(study_name == "QinJ_2012")  %>%
  returnSamples("relative_abundance", counts = FALSE)
physeq.qinj <- makePhyloseqFromTreeSE(physeq.qinj,
                                     assay.type = "relative_abundance")

(physeq.qj.ctl <- subset_samples(physeq.qinj, study_condition == "control"))
(physeq.qj.t2d <- subset_samples(physeq.qinj, study_condition == "T2D"))

rm(physeq.qinj)

# Filtrado 0.01% en 10% de las muestras
(physeq.qj.ctl.f <- filter_taxa(physeq.qj.ctl,
                                flist =  filterfun(kOverA(k = 17, A = 0.01)),
                                prune = TRUE)) 
saveRDS(physeq.qj.ctl.f, "physeq_qj_ctl_f.RDS")

(physeq.qj.t2d.f <- filter_taxa(physeq.qj.t2d,
                                 flist =  filterfun(kOverA(k = 17, A = 0.01)),
                                 prune = TRUE)) 
saveRDS(physeq.qj.t2d.f, "physeq_qj_t2d_f.RDS")
