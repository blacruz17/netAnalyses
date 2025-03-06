library(curatedMetagenomicData)
library(tidyverse)
library(phyloseq)
library(NetCoMi)
library(mia)
library(genefilter)

rm(list = ls())

# HELPER EXPORT FUNCTION #######################################################

exportNet <- function(net, physeq, props, filename){
  # EDGES
  # Create edge object from the edge list exported by netConstruct()
  edges <- net$edgelist1 %>%
    select(v1, v2, asso) %>%
    rename(Source = v1,
           Target = v2,
           Weight = asso) %>%
    mutate(Type = "Undirected",
           Sign = if_else(Weight < 0, "Negative", "Positive"))
  
  # NODES
  # get taxonomic info:
  taxa <- data.frame(physeq@tax_table) %>%
    rownames_to_column(var = "Label")
  
  # get mean abundances:
  mean_ab <- data.frame(
    "Abundance" = apply(data.frame(physeq@otu_table), 1, mean)) %>%
    rownames_to_column(var = "Label")
  
  # get clusters:
  clusters <- data.frame("Cluster" =  props$clustering$clust1) %>%
    rownames_to_column(var = "Label")
  
  # get hubs:
  hubs <- data.frame(isHub = rep(1, length(props$hubs$hubs1)),
                     "Label" = props$hubs$hubs1)
  
  # join all of them together:
  metadata <- taxa %>%
    left_join(mean_ab, by = "Label") %>%
    left_join(clusters, by = "Label") %>%
    left_join(hubs, by = "Label") %>%
    mutate(isHub = if_else(is.na(isHub), 0, 1))
  
  
  # WRITE CSV FILES:
  write_csv(edges,
            file = paste0(filename, "_edges.csv"),
            quote = "needed")
  write_csv(metadata,
            paste0(filename, "_metadata.csv"),
            quote = "needed")
}

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