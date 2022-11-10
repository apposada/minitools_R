#' Functions ANANSE
require(igraph)

#' Load network
LoadNetworkData <- function(x) {
  nw <- fread(
    text = gsub(
      "—", # replace dash line from ananse rows
      "\t", # with a tab, thus "tf" "\t" "tg"
      readLines(x) #taken from: 
    ),
    skip = 1 # skip the colnames
  )
  colnames (nw) <- c("tf","tg","prob")
  # output of function: a table/dataframe in the structure of igraph
  return(nw)
}

#' Filter network
FilterNetwork <- function(nw,q = 0.85) {
  q_filt <- nw$prob > quantile(nw$prob, q)
  ngenes <- length(
    unique(
      c(nw$tf,nw$tg)
    )
  )
  nw_filt <- nw[q_filt,]
  ngenes_filt <- length(
    unique(
      c(nw_filt$tf,nw_filt$tg)
    )
  )
  print(
    paste0(
      "Original network is ",
      ngenes,
      " genes. ",
      "Filtered network is ",
      ngenes_filt,
      " genes."
    )
  )
  # output of function: a filtered network in the shape of a data frame
  return(nw_filt)
}

#' Generate network
GenerateUndirectedNetwork <- function(nw, q = 0.85) {
  colnames(nw) <- c("gene1","gene2","value")
  q_filt <- nw$value > quantile(nw$value, q)
  ngenes <- length(
    unique(
      c(nw$gene1,nw$gene2)
    )
  )
  nw_G <- graph.data.frame( # the initial igraph object
    d = nw[q_filt,],
    directed = F
  )
  E(nw_G)$width <- nw[q_filt,]$value
  # output of function: a filtered network in the shape of an igraph object
  return(nw_G)
}

GenerateDirectedNetwork <- function(nw, q = 0.85) {
  colnames(nw) <- c("tf","target","prob")
  q_filt <- nw$prob > quantile(nw$prob, q)
  ngenes <- length(
    unique(
      c(nw$tf,nw$target)
    )
  )
  nw_G <- graph.data.frame( # the initial igraph object
    d = nw[q_filt,],
    directed = T
  )
  E(nw_G)$width <- nw[q_filt,]$prob
  # output of function: a filtered network in the shape of an igraph object
  return(nw_G)
}

GenerateNetwork <- function(nw, q = 0.85, directed = FALSE){
  if(directed == FALSE ){
    nw_G <- GenerateUndirectedNetwork(nw = nw, q = q)
  } else {
    nw_G <- GenerateDirectedNetwork(nw = nw, q = q)
  }
  return(nw_G)
}



#' Parse the network
ParseNetwork <- function(
  graph,
  list_attr = list(
    tflist, efflist, classANDcolor, td_hk, gfams, funcat, markerlist, ...
    )
  ) {
  # 1.1 asigna lo q tenga que asignar de datos para igraph (TF class, color..)
  df_attr <- data.frame(
    id = V(graph)$name,
    index = seq(
      1:length(V(graph))
    )
  )

  # merge data
  for (i in list_attr){
    df_attr <- merge(
      df_attr,
      i,
      by.x = 1,
      by.y = 1,
      all.x = T
    ) 
  }

  # cleanup data
  df_attr[is.na(df_attr)] <- ""
  df_attr <- df_attr[order(df_attr$index),] #IMPORTANT
  df_attr <- df_attr[!duplicated(df_attr$id),]
  rownames(df_attr) <- NULL

  # slap together with the network
  for (i in colnames(df_attr)[-1]){
    graph <- set_vertex_attr(graph, paste0(i), value = df_attr[[i]])
  }
  #' output of the function: a graph with additional 
  #' attributes to the different nodes and edges
  output <- list(
    graph,
    df_attr
  )
  return(output)
}

#' Network stats and metrics
NetworkStats <- function(
  nw, graph, att, N = 10, C = 1,
  geneid2go, gene_universe
  ){

  print("Basic Stats") # add option of directed or undirected
  numgenes_network <- length(V(graph))
  # top N connected TFs measured as # of stemming edges
  connectionxtf <- c(table(nw$tf))
  numtfs <- length(connectionxtf)
  meanconnectionxtf <- mean(connectionxtf)
  # topEmmitters <- names(rev(sort(table(nw$tf)))[1:N])

  # top N connected TGs measured as # of incoming edges
  connectionxtg <- c(table(nw$tg))
  numtgs <- length(connectionxtg)
  meanconnectionxtg <- mean(connectionxtg)
  # topReceivers <- names(rev(sort(table(nw$tg)))[1:N])

  # Num outgoing vs incoming interactions
  nw_InOut <- computeInVsOut(nw)
  topEmitters <- nw_InOut$id[rev(order(nw_InOut$emit))][1:N]
  topReceivers <- nw_InOut$id[rev(order(nw_InOut$recv))][1:N]

  # Num TFs that self-regulate in the network
  selfreg <- length(which(nw$tf == nw$tg))

  # Divide strongest and weakest components
  print("Components")
  graph_components <- data.frame(
    id=names(components(graph,mode=c("strong"))$membership),
    member=components(graph,mode=c("strong"))$membership
  )

  # new subgraph with only genes from largest CC
  mainCC <- names(
    rev(sort(table(
      graph_components$member
    )))[C])
  print(mainCC)

  graph_strong <- induced_subgraph(
    graph,
    graph_components$id[graph_components$member == mainCC],
    impl = "auto"
  )

  ccsizes <- table(graph_components$member)[1:10]
  mainccsize <- max(table(graph_components$member))

  mainCC_genes <- graph_components$id[graph_components$member == mainCC]

  print("Centrality and category-based metrics")
  #' number TF genes
  numTFs <- length(V(graph_strong)[V(graph_strong)$TFEG == "TF"])
  #' number effector genes
  numEGs <- length(V(graph_strong)[V(graph_strong)$TFEG == "EG"])
  #' ratio TF/effectors in the network
  ratioTFEG <- numTFs/numEGs

  print("Centrality")
  #' centrality of all TFs --> top central TFs
  print("Centrality of TFs")
  tf_central <- getcentral(graph, att, "TFclass")
  centralTFs <-  rev(sort(unlist(tf_central)))[1:N]

  #' numero y centrality de transdevs en las networks,
  print("Centrality of Trans-Devs")
  td_central <- getcentral(graph, att, "TDHK")
  centralTDs <- rev(sort(unlist(td_central$td)))[1:N]

  #' number incoming edges per funcat family
  print("Edges of Functional Categories")
  funcat_edges <- merge(
    nw_InOut,
    att[,c("id","funcat")],
    by.x = 1,
    by.y = 1,
    all.x = T
  )
  #' centrality per funcat family
  print("Centrality of Functional Categories")
  funcat_central <- getcentral(graph,att,"funcat")

  #' GOs
  print("Getting Gene Ontology")
  GO_outputs <- getGOs( # re-define how these categories are defined
    list(
      tfs = V(graph_strong)$name[
        V(graph_strong)$name %in% att$id[att$TFclass != ""]
        ],
      active_tfs = unique(nw$tf),
      egs = V(graph_strong)$name[
        !(V(graph_strong)$name %in% att$id[att$TFclass != ""])
        ],
      target_genes = unique(nw$tg),
      whole_graph = V(graph)$name
      ),
    gene_universe = gene_universe,
    gene2GO = geneid2go,
    alg = "elim"
  )

  print("Generating output")
  report <- list(
    numgenes_network = numgenes_network,
    connection_per_TF = connectionxtf,
    num_active_TFs = numtfs,
    mean_connection_per_TF = meanconnectionxtf,
    connection_per_TG = connectionxtg,
    num_TGs = numtgs,
    mean_connection_per_TG = meanconnectionxtg,
    top_Emitters = topEmitters,
    In_Out_per_Gene = nw_InOut,
    top_Receivers = topReceivers,
    num_selfregulated_TFs = selfreg,
    genes_per_component = graph_components,
    Connected_component_sizes = ccsizes,
    main_Connected_Component = mainCC,
    main_Connected_Component_graph = graph_strong,
    main_Connected_Component_size = mainccsize,
    main_Connected_Component_genes = mainCC_genes,
    num_total_TFs = numTFs,
    num_Effector_Genes = numEGs,
    ratio_TF_per_EG = ratioTFEG,
    Centrality_per_TFclass = tf_central,
    most_central_TFs = centralTFs,
    Centrality_per_TransDev_HK = td_central,
    most_central_transdevs = centralTDs,
    edges_per_funcat_family = funcat_edges,
    Centrality_per_funcat_family = funcat_central,
    GO_tables = GO_outputs$GOtable
  )
  return(report)
  print("Done.")
}

#' Generate tables with the stat reports
digestStats <- function(stats){
  #
  #
}

#' Plot the basic metrics of the networks
NetworkPlots <- function(stats, nw, tfcol, layout = TRUE, colpal = default_colpal, pdf = TRUE, pdfname = NULL, ...){
  if (layout == TRUE) {
    layout(
      matrix(
        c(1,2,3,4,5,6,7,8,9,10,10,11,12,12,13),
        nrow = 5,
        byrow = TRUE
      )
    )
  }

  density_network(nw)
  plot_numgenes_network(stats)
  plot_typegenes_network(stats)
  plot_connections(stats)
  # plot_connections_per_gene(stats)
  plot_genes_ratioconnections(stats)
  plot_ratio_connections(stats)
  plot_selfreg(stats)
  plot_size_CC(stats)
  plot_centrality_TDHK(stats)
  plot_centrality_TFclass(stats,  f = tfcol)
  plot_top_central_tfs(stats)
  plot_behavior_per_category2(stats)
  plot_top_central_transdev(stats)
  # plot_GOs(stats$GO_plots)
  par(mfrow=c(1,1))
}


#' Network Comparison
compareNetworks <- function(nw_a,nw_b,graph_a,graph_b,
                            influence = NULL, col_a = "green",
                            col_b = "purple", geneset_interest = NULL, 
                            top = 0.9, tfs, gene_universe, id2go){
  
  require(VennDiagram)

  #common_genes
  print("Defining exclusive and common targets")
  a_top_tgt <- names(table(nw_a$tg)[table(nw_a$tg) >= quantile(table(nw_a$tg),top)])
  b_top_tgt <- names(table(nw_b$tg)[table(nw_b$tg) >= quantile(table(nw_b$tg),top)])
  
  a_top_tgt_no_b <- a_top_tgt[!(a_top_tgt %in% b_top_tgt)]
  b_top_tgt_no_a <- b_top_tgt[!(b_top_tgt %in% a_top_tgt)]
  
  list_targets <- list(
    targets_exclusive_a = a_top_tgt_no_b,
    targets_exclusive_b = b_top_tgt_no_a,
    targets_common_ab = a_top_tgt[a_top_tgt %in% b_top_tgt]
  )
  
  numtgts_exclusive_a <- length(a_top_tgt) - length(which(a_top_tgt %in% b_top_tgt))
  numtgts_exclusive_b <- length(b_top_tgt) - length(which(a_top_tgt %in% b_top_tgt))
  numtgts_common <- length(which(a_top_tgt %in% b_top_tgt))
  
  # plots
  
  boxplot(
    list(
      c(table(nw_a$tg)),
      c(table(nw_b$tg))
    ),
    col=c(
      'green',
      'purple'
    ),
    main="connections per target gene",
    xaxt="n"
  )
  axis(1,at=c(1,2),c("a","b"))
  
  venn <- venn.diagram(
    x = list(a_top_tgt, b_top_tgt),
    category.names = c("A" , "B"),
    filename = NULL,
    cex = 1,
    cat.cex = 1,
    lwd = 0.5,
    fill=c("green", 'purple'),
    
  )
  plot(0,type="n", 
       xlab = "",
       ylab = "",
       main = "Targets shared between networks",
       axes = F
  )
  grid::grid.draw(venn)
  
  # Genes of interest
  
  if ( is.null(geneset_interest) ) { 
    print("Number of genes of interest in targets")
    genes_interest_innetwork = c(
      "A" = length(which(a_top_tgt_no_b %in% geneset_interest)),
      "B" = length(which(b_top_tgt_no_a %in% geneset_interest))
    )
    barplot(
      height=c(
        genes_interest_innetwork[1]/length(a_top_tgt_no_b)*100,
        genes_interest_innetwork[2]/length(b_top_tgt_no_a)*100
      ),
      col = c(
        col_a,
        col_b
      ),
      ylim=c(0,40),
      main="Genes of interest\nin network",
      names=c("a","b"),
      ylab="gene percent",
      las=1
    )
  } else {
    genes_interest_innetwork = "No genes of interest provided."
  }
  
  # GOs
  print("getting GO terms")
  GOs_targets <- getGOs(
    list_targets,
    gene_universe = gene_universe,
    gene2GO = id2go,
    alg = "elim"
  )
  
  # centrality
  print("gene centrality across networks")
  common_genes_in_nw <-
    V(graph_a)$name [
      V(graph_a)$name %in% V(graph_b)$name
    ]
  
  centrality_per_network <- data.frame(
    id = common_genes_in_nw,
    centrality_a = closeness(
      graph_a,
      vids = common_genes_in_nw,
      mode = c("all"),
      weights = NULL,
      normalized = F
    ),
    centrality_b = closeness(
      graph_b,
      vids = common_genes_in_nw,
      mode = c("all"),
      weights = NULL,
      normalized = F
    )
  )
  
  plot(
    centrality_per_network$centrality_a,
    centrality_per_network$centrality_b
  )
  
  # Influence
  if ( influence == TRUE ){
    print("Influence")
    influence_results <- influ_table(influence,tfs)
  } else {
    influence = "no Influence provided."
  }
  
  res = list(
    top_targets_a = a_top_tgt,
    top_targets_b = b_top_tgt,
    targets_exclusive_a = a_top_tgt_no_b,
    targets_exclusive_b = b_top_tgt_no_a,
    connects_per_tgt_gene_a = table(nw_a$tg),
    connects_per_tgt_gene_b = table(nw_b$tg),
    numtgts_exclusive_a = numtgts_exclusive_a,
    numtgts_exclusive_b = numtgts_exclusive_b,
    numtgts_common = numtgts_common,
    genes_interest_innetwork,
    GOs_targets = GOs_targets,
    centrality_per_network = centrality_per_network,
    influence_results = influence_results
  )
  
  return(res)
}

# 3. funcion para plotear la centrality de los TFs over time (fuzz/broadlines depicting quantile or median or whatever)


#' 1.10 concordancia motivos atac homer y motivos de los tfs de la network (los + centrales, o los q tienen mas targets, etc.)
#' Calculate number (barplot) and percentage (stacked, pie) of genes of each module in the network, grab color from module names
#' Calculate number of markers of each cell type in the whole network (using a table of gene-tissue/cell/whatever)
#' Calculate number of markers of each cell type in the EFFECTOR genes of the network (using a table of gene-tissue/cell/whatever)
#' Generate a Heatmap with rows = TFs, Cols = TGs Done i think?? it doesnt see well
