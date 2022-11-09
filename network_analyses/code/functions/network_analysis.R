# Network Analysis in Ptychodera

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
  "/home/ska/aperpos/Def_spp1/outputs/eggnog/20221014_ptyFlav3_CYi_longest.cds.fa.transdecoder.pep_eggnog_meNOG_all.emapper.annotations.COG.tsv",
  col.names = c("id","funcat")
)

# Trans-developmental / housekeeping
spp1_TDHK <- data.frame(
  id = spp1_all_gene_names,
  TDHK =  ifelse(spp1_all_gene_names %in% list_spp1_transdev_exclusive_nohk, "td", ifelse(spp1_all_gene_names %in% list_spp1_housekeep_exclusive_notd, "hk", "none"))
)

spp1_attributes_list <- list(spp1_TFs,spp1_funcat,spp1_TDHK,spp1_TFEG)

#' ---------------------------------------------------------
#' EARLY BLASTULA
spp1_EB_nw <- LoadNetworkData("/home/ska/aperpos/Def_spp1/outputs/ananse/20210901_EB.network")
spp1_EB_nw$tf <- gsub("TCONS", "TCONS_", spp1_EB_nw$tf)
spp1_EB_nw$tg <- gsub("TCONS", "TCONS_", spp1_EB_nw$tg)
spp1_EB_nw2 <- FilterNetwork(spp1_EB_nw,q=0.95)
spp1_EB_graph <- GenerateNetwork(spp1_EB_nw,q = 0.95)
spp1_EB_parsenetwork <- ParseNetwork(spp1_EB_graph, spp1_attributes_list)
spp1_EB_graph2 <- spp1_EB_parsenetwork[[1]]
spp1_EB_df_attr <- spp1_EB_parsenetwork[[2]]
spp1_EB_stats <- NetworkStats(spp1_EB_nw2, spp1_EB_graph2, spp1_EB_df_attr, N = 10)
NetworkPlots(
  stats = spp1_EB_stats,
  nw = spp1_EB_nw,
  tfcol = TFclasses_colors,
  pdf = F
)

#' ---------------------------------------------------------
#' EARLY GASTRULA
spp1_EG_nw <- LoadNetworkData("/home/ska/aperpos/Def_spp1/outputs/ananse/20220907_EG.network")
spp1_EG_nw$tf <- gsub("TCONS", "TCONS_", spp1_EG_nw$tf)
spp1_EG_nw$tg <- gsub("TCONS", "TCONS_", spp1_EG_nw$tg)
spp1_EG_nw2 <- FilterNetwork(spp1_EG_nw,q=0.95)
spp1_EG_graph <- GenerateNetwork(spp1_EG_nw,q = 0.95)
spp1_EG_parsenetwork <- ParseNetwork(spp1_EG_graph, spp1_attributes_list)
spp1_EG_graph2 <- spp1_EG_parsenetwork[[1]]
spp1_EG_df_attr <- spp1_EG_parsenetwork[[2]]
spp1_EG_stats <- NetworkStats(spp1_EG_nw2, spp1_EG_graph2, spp1_EG_df_attr, N = 10)
NetworkPlots(
  stats = spp1_EG_stats,
  nw = spp1_EG_nw,
  tfcol = TFclasses_colors,
  pdf = F
)

#' ---------------------------------------------------------
#' MID GASTRULA
spp1_MG_nw <- LoadNetworkData("/home/ska/aperpos/Def_spp1/outputs/ananse/20220919_MG.network")
spp1_MG_nw$tf <- gsub("TCONS", "TCONS_", spp1_MG_nw$tf)
spp1_MG_nw$tg <- gsub("TCONS", "TCONS_", spp1_MG_nw$tg)
spp1_MG_nw2 <- FilterNetwork(spp1_MG_nw,q=0.95)
spp1_MG_graph <- GenerateNetwork(spp1_MG_nw,q = 0.95)
spp1_MG_parsenetwork <- ParseNetwork(spp1_MG_graph, spp1_attributes_list)
spp1_MG_graph2 <- spp1_MG_parsenetwork[[1]]
spp1_MG_df_attr <- spp1_MG_parsenetwork[[2]]
spp1_MG_stats <- NetworkStats(spp1_MG_nw2, spp1_MG_graph2, spp1_MG_df_attr, N = 10)
NetworkPlots(
  stats = spp1_MG_stats,
  nw = spp1_MG_nw,
  tfcol = TFclasses_colors,
  pdf = F
)

#' ---------------------------------------------------------
#' LATE GASTRULA
spp1_LG_nw <- LoadNetworkData("/home/ska/aperpos/Def_spp1/outputs/ananse/20210908_LG.network")
spp1_LG_nw$tf <- gsub("TCONS", "TCONS_", spp1_LG_nw$tf)
spp1_LG_nw$tg <- gsub("TCONS", "TCONS_", spp1_LG_nw$tg)
spp1_LG_nw2 <- FilterNetwork(spp1_LG_nw,q=0.95)
spp1_LG_graph <- GenerateNetwork(spp1_LG_nw,q = 0.95)
spp1_LG_parsenetwork <- ParseNetwork(spp1_LG_graph, spp1_attributes_list)
spp1_LG_graph2 <- spp1_LG_parsenetwork[[1]]
spp1_LG_df_attr <- spp1_LG_parsenetwork[[2]]
spp1_LG_stats <- NetworkStats(spp1_LG_nw2, spp1_LG_graph2, spp1_LG_df_attr, N = 10)
NetworkPlots(
  stats = spp1_LG_stats,
  nw = spp1_LG_nw,
  tfcol = TFclasses_colors,
  pdf = F
)


#' ---------------------------------------------------------
#' NETWORK COMPARISONS

EB_vs_EG <- compareNetworks(
  nw_a = spp1_EB_nw2,
  nw_b = spp1_EG_nw2,
  graph_a = spp1_EB_graph2,
  graph_b = spp1_LG_graph2,
  influence = F,
  top = 0.9,
  id2go = spp1_geneID2GO,
  gene_universe = spp1_all_gene_names
)

EG_vs_MG <- compareNetworks(
  nw_a = spp1_EG_nw2,
  nw_b = spp1_MG_nw2,
  graph_a = spp1_EB_graph2,
  graph_b = spp1_LG_graph2,
  influence = F,
  top = 0.9,
  id2go = spp1_geneID2GO,
  gene_universe = spp1_all_gene_names
)


MG_vs_LG <- compareNetworks(
  nw_a = spp1_MG_nw2,
  nw_b = spp1_LG_nw2,
  graph_a = spp1_EB_graph2,
  graph_b = spp1_LG_graph2,
  influence = F,
  top = 0.9,
  id2go = spp1_geneID2GO,
  gene_universe = spp1_all_gene_names
)

EB_vs_LG <- compareNetworks(
  nw_a = spp1_EB_nw2,
  nw_b = spp1_LG_nw2,
  graph_a = spp1_EB_graph2,
  graph_b = spp1_LG_graph2,
  influence = F,
  top = 0.9,
  id2go = spp1_geneID2GO,
  gene_universe = spp1_all_gene_names
)

ananse_EB_comparedto_EG <- 
  read.table("/home/ska/aperpos/Def_spp1/outputs/ananse/20220927_EB_EG_influence.txt",header=T)
ananse_EB_comparedto_EG$factor <- sub("TCONS","TCONS_",ananse_EB_comparedto_EG$factor)
ananse_EB_comparedto_EG <- influ_table(ananse_EB_comparedto_EG,spp1_TFs)

ananse_EG_to_MG <- 
  read.table("/home/ska/aperpos/Def_spp1/outputs/ananse/20220927_MG_EG_influence.txt",header=T)
ananse_EG_to_MG$factor <- sub("TCONS","TCONS_",ananse_EG_to_MG$factor)
ananse_EG_to_MG <- influ_table(ananse_EG_to_MG,spp1_TFs)

ananse_MG_to_LG <- 
  read.table("/home/ska/aperpos/Def_spp1/outputs/ananse/20220927_LG_MG_influence.txt",header=T)
ananse_MG_to_LG$factor <- sub("TCONS","TCONS_",ananse_MG_to_LG$factor)
ananse_MG_to_LG <- influ_table(ananse_MG_to_LG,spp1_TFs)

ananse_EB_comparedto_LG <- 
  read.table("/home/ska/aperpos/Def_spp1/outputs/ananse/20220927_EB_LG_influence.txt",header=T)
ananse_EB_comparedto_LG$factor <- sub("TCONS","TCONS_",ananse_EB_comparedto_LG$factor)
ananse_EB_comparedto_LG <- influ_table(ananse_EB_comparedto_LG,spp1_TFs)

ananse_EB_to_LG <- 
  read.table("/home/ska/aperpos/Def_spp1/outputs/ananse/20210909_ANANSE_EG_EB",header=T)
ananse_EB_to_LG$factor <- sub("TCONS","TCONS_",ananse_EB_to_LG$factor)
ananse_EB_to_LG <- influ_table(ananse_EB_to_LG,spp1_TFs)
