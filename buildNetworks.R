library(tidyverse)
library(NetCoMi)

rm(list = ls())

# Data Import #########
# karlsson
physeq.k.igt.f <- readRDS("physeq_k_igt_f.RDS")
physeq.k.t2d.f <- readRDS("physeq_k_t2d_f.RDS")
physeq.k.ctl.f <- readRDS("physeq_k_ctl_f.RDS")

# qinj
physeq.qj.t2d.f <- readRDS("physeq_qj_t2d_f.RDS")
physeq.qj.ctl.f <- readRDS("physeq_qj_ctl_f.RDS")

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

# CONSTRUIR REDES SPIECEASI ####################################################

library(doParallel)
num_cores <- 4 # cambiar segun donde se ejecute
cl <- makeCluster(num_cores)
registerDoParallel(cl)

getNetwork <- function(physeq.obj, filename) {
  mynet <- netConstruct(data = physeq.obj,
                               dataType = "counts",
                               measure = "spieceasi",
                               sparsMethod = "none",
                               nboot = 1000,
                               cores = 20L)
  myprops <- netAnalyze(mynet,
                               centrLCC = TRUE,
                               clustMethod = "cluster_fast_greedy",
                               hubPar = c("degree", "betweenness"),
                               hubQuant = 0.9,
                               weightDeg = FALSE, normDeg = FALSE)
  saveRDS(mynet, paste0("net_", filename, ".rds"))
  saveRDS(myprops, paste0("props_", filename, ".rds"))
  exportNet(mynet, physeq.obj, myprops, filename)
}

# Karlsson
getNetwork(physeq.k.ctl.f, "karlsson_control")
getNetwork(physeq.k.igt.f, "karlsson_igt")
getNetwork(physeq.k.t2d.f, "karlsson_t2d")

# QinJ
getNetwork(physeq.qj.ctl.f, "qinj_control")
getNetwork(physeq.qj.t2d.f, "qinj_t2d")

stopCluster(cl)
