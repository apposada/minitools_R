# Network Analysis scripts

# Data loading

#Transcription factors
class_color_df <- colors_all_tfclasses
colnames(class_color_df) <- c("TFclass","color")

TFclasses <- class_color_df$TFclass
TFclasses_colors <- class_color_df$color
names(TFclasses_colors) <- TFclasses

spp1_TFs <- spp1_tfs
colnames(spp1_TFs)[2] <- "TFclass"

spp1_TFs <- merge(
  spp1_TFs,
  class_color_df,
  by.x = 2,
  by.y = 1,
  all.x = T
)[,c(2,1,3)]

# Effector genes
spp1_TFEG <- data.frame(
  id = spp1_all_gene_names,
  TFEG = ifelse(spp1_all_gene_names %in% spp1_TFs$id, "TF", "EG")
)

# Functional categories (COG)
spp1_funcat <- read.table(
  "/path/to/.emapper.annotations.COG.tsv",
  col.names = c("id","funcat")
)

# Trans-developmental / housekeeping
spp1_TDHK <- data.frame(
  id = spp1_all_gene_names,
  TDHK =  ifelse(spp1_all_gene_names %in% list_spp1_transdev_exclusive_nohk, "td", ifelse(spp1_all_gene_names %in% list_spp1_housekeep_exclusive_notd, "hk", "none"))
)

spp1_attributes_list <- list(spp1_TFs,spp1_funcat,spp1_TDHK,spp1_TFEG)

#' ---------------------------------------------------------
#' SAMPLE 1
spp1_sample1_nw <- LoadNetworkData("/path/to/ananse/sample1.network")

spp1_sample1_nw2 <- FilterNetwork(spp1_sample1_nw,q=0.95)
spp1_sample1_graph <- GenerateNetwork(spp1_sample1_nw,q = 0.95)
spp1_sample1_parsenetwork <- ParseNetwork(spp1_sample1_graph, spp1_attributes_list)
spp1_sample1_graph2 <- spp1_sample1_parsenetwork[[1]]
spp1_sample1_df_attr <- spp1_sample1_parsenetwork[[2]]
spp1_sample1_stats <- NetworkStats(spp1_sample1_nw2, spp1_sample1_graph2, spp1_sample1_df_attr, N = 10)
NetworkPlots(
  stats = spp1_sample1_stats,
  nw = spp1_sample1_nw,
  tfcol = TFclasses_colors,
  pdf = F
)

#' ---------------------------------------------------------
#' SAMPLE 2
spp1_sample2_nw <- LoadNetworkData("/path/to/ananse/sample1.network")

spp1_sample2_nw2 <- FilterNetwork(spp1_sample2_nw,q=0.95)
spp1_sample2_graph <- GenerateNetwork(spp1_sample2_nw,q = 0.95)
spp1_sample2_parsenetwork <- ParseNetwork(spp1_sample2_graph, spp1_attributes_list)
spp1_sample2_graph2 <- spp1_sample2_parsenetwork[[1]]
spp1_sample2_df_attr <- spp1_sample2_parsenetwork[[2]]
spp1_sample2_stats <- NetworkStats(
  spp1_sample2_nw2, spp1_sample2_graph2, spp1_sample2_df_attr, N = 10
   )
NetworkPlots(
  stats = spp1_sample2_stats,
  nw = spp1_sample2_nw,
  tfcol = TFclasses_colors,
  pdf = F
)

#' ---------------------------------------------------------
#' NETWORK COMPARISONS

sample1_vs_sample2 <- compareNetworks(
  nw_a = spp1_sample1_nw2,
  nw_b = spp1_sample2_nw2,
  graph_a = spp1_sample1_graph2,
  graph_b = spp1_LG_graph2,
  influence = F,
  top = 0.9,
  id2go = spp1_geneID2GO,
  gene_universe = spp1_all_gene_names
)

ananse_sample1_comparedto_sample2 <- 
  read.table("/path/to/ananse/influence.txt",header=T)
ananse_sample1_comparedto_sample2 <- influ_table(ananse_sample1_comparedto_sample2,spp1_TFs)
