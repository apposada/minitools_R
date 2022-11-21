# Setup
options(scipen=999)
require(vroom)
require(reshape2)
require(ComplexHeatmap)
require(circlize)
require(viridis)
require(colorspace)
require(RColorBrewer)
require(dplyr)
require(stringr)
source('./code/functions/get_best_scores.R')
source('./code/functions/tfclasses_proportions.R')
source('./code/functions/barplotcluster.R')
source('./code/functions/correg_tf_modules.R')

TF_whole_analysis <- function(
  x_cpm, x_tfs, x_wg_kME, topClasses = 30,
  tf_colors = colors(), sample_colors = colors(),
  user_seed = 4343, filter_top_tfs = FALSE,
  clu_method = "ward.D2", cor_method = "spearman",
  barplot_clu_method = "ward.D2"
  min_tfs_class = 4, min_rho = 0.5, 
  ){
    
    # Get info from data provided
    ngenes = nrow(x_cpm)
    samples = colnames(x_cpm)
    nsamples = ncol(x_cpm)
    tfclasses = sort(unique(x_tfs$class))

    # top TF classes
    x_tfs_scores <- get_best_scores(x_tfs[,1:2], x_tfs_cpm[,1:nsamples])
    x_tfs_scores_top <- x_tfs_scores$class[c(1:topClasses)] 
    x_tfs_scores_top <- sort(x_tfs_scores_top) # sort alphabetically
    
    # Parse TFs from cpm table
    x_tfs_cpm <- x_cpm[
        rownames(x_cpm) %in% x_tfs$id,
        1:nsamples
    ]
    colnames(x_tfs_cpm) <- colnames(x_cpm)[1:nsamples]
    
    ntfs <- nrow(x_tfs_cpm)
    
    x_tfs_cpm <- as.data.frame(x_tfs_cpm)
    
    # CV , standard deviation / mean
    x_tfs_cpm$cv <- apply(
        x_tfs_cpm[,1:nsamples],
        1,
        function(x){
        sd(x)/mean(x)
        }
    )
    
    x_tfs_cv <- data.frame(
        id = rownames(x_tfs_cpm),
        cv = x_tfs_cpm$cv
    )
    
    # sort by CV
    x_tfs_cpm <- x_tfs_cpm[rev(order(x_tfs_cpm$cv)),1:nsamples] #remove cv column

    # FILTER
    if ( filter_top_tfs != FALSE ){
        x_tfs_cpm_top <- x_tfs_cpm[1:filter_top_tfs,]
    } else {
        x_tfs_cpm_top <- x_tfs_cpm
    }
    x_tfs_cpm_top <- x_tfs_cpm_top[order(rownames(x_tfs_cpm_top)),]
    
    # The Heatmap
    
    # Annotations of dendrograms
    x_ha_clu <- sample_colors[1:nsamples]
    names(x_ha_clu) <- samples
    
    clu_ha = HeatmapAnnotation( # Color annotation for samples
        cluster = names(x_ha_clu),
        col=list(cluster=x_ha_clu)
        )
    
    # x_ha = rowAnnotation(class = x_top_tfs_id_to_name[,2],col=list(class=x_ha_col)) # class of transcription factors
    #top_annotation=tfs_ha
    
    # Matrices of data
    x_tfs_cor <- cor(
        scale(
        t(
            log(x_tfs_cpm_top[,1:nsamples]+1)
            )
        ), 
        method = cor_method
        )
    
    x_tfs_zsco_expr <-  t(
        scale(
        t(
            x_tfs_cpm_top[,1:nsamples]
        )
        )
    )
    
    # Color palettes
    x_hm_cor_col <- colorRamp2(
        c( # breaks:
        seq(-0.4,0.2,len=10), # clipped
        seq(0.3,0.5,len=10) # clipped
        ),
        colorRampPalette( #color:
        rev(brewer.pal(7,"RdYlBu"))
        )(20)
    )


    x_hm_expr_col <- colorRamp2(c(1:6,6.5),rev(sequential_hcl(7,"YlGnBu")))
    
    x_hm_cor <- Heatmap(
        x_tfs_cor,
        name = cor_method,
        clustering_method_columns = clu_method,
        clustering_method_rows = clu_method,
        col=x_tfs_heat_col,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5)
    )

    x_hm_expr <- Heatmap(
        x_tfs_zsco_expr+2,
        name="norm\nexpr",
        col=x_hm_expr_col,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        cluster_rows=T,
        # clustering_distance_columns = "kendall",
        clustering_method_columns = "single",
        clustering_method_rows = "ward.D2",
        cluster_columns=F,
        #left_annotation=x_ha,
        top_annotation=clu_ha,
    )
    
    # Proportions of TF classes
    #data
    x_tfs_cpm_scoring_top <- 
        merge(
        x_tfs_cpm[,1:nsamples],
        x_tfs,
        by.x = 0,
        by.y = 1,
        )
    
    rownames(x_tfs_cpm_scoring_top) <- x_tfs_cpm_scoring_top$Row.names
    x_tfs_cpm_scoring_top$Row.names <-  NULL
    x_tfs_cpm_scoring_top <- x_tfs_cpm_scoring_top[x_tfs_cpm_scoring_top$class %in% x_tfs_scores_top,]
    # here a call to tfclasses_proportions.R
    x_tfs_proportions <- tfclasses_proportions(x_tfs_cpm_scoring_top, nsamples, ...)

    set.seed(user_seed)
    x_tfclass_cols = sample(c(sample(colors(),15),rainbow(15)))
    # call here ggplot2 with the same types of barplots as ptychodera
    
    barplot(x_tfs_proportions$relative_exp_pergene_perclass, col = x_tfclass_cols, las = 2) #we should investigate this
    
    barplot_cluster(t(x_tfs_proportions$relative_exp_pergene_perclass), clu_method = barplot_clu_method) # provide custom colors, return tree as well

    numgenes_hm <- Heatmap(
        name = "numTFs",
        x_tfs_proportions$numgenes,
        col=colorRamp2( c(0,1,2,5,10,20), cividis(6)),
        clustering_method_columns = "ward.D2",
    )
    
    rel_exp_pergene_perclass_hm_cividis <- Heatmap(
        name = "% exp/TFgene/class",
        x_tfs_proportions$relative_exp_pergene_perclass,
        col=colorRamp2(c(0,0.02,0.1,0.15,0.2,0.3), cividis(6)),
        clustering_method_columns = "ward.D2",
    )
    
    rel_exp_pergene_perclass_hm <- Heatmap( # add celltype annotation of cluster colors by Jordi
        name = "# exp/TFgene/class"
        log(t(x_tfs_proportions$expression_pergene_perclass)+1),
        clustering_method = "ward.D2",
        # col = colorRampPalette(rev(brewer.pal(7,"RdYlBu")))(20)
        col = colorRamp2(
            c(0,1,2,4,5,6),
            rev(brewer.pal(7,"RdYlBu")),
    )

    exp_pergene_perclass_hm <- Heatmap(
        name = "# exp/TFgene/class",
        log(t(x_tfs_proportions$expression_pergene_perclass)+1),
        col = colorRamp2(
            c(0,1,2,4,5,6),
            c("#f9ecce","#ffe8b0","#F7BA3C","#F5A100","#C32200","#A10706")
            ),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2"
    )

    ## Corregulation TFs ~ Modules

    if ( x_wg_kME != FALSE){
        x_correg_tf_modules <- correg_tf_modules(x_tfs_cpm = x_tfs_cpm,x_tfs = x_tfs,x_wg_kME = x_wg_kME, min_tfs_class = min_tfs_class, min_rho = min_rho)
    } 
  
    # list of tabular quantitative data
    tables = list(
        tfs_cpm = x_tfs_cpm,
        tfs_cpm_top = x_tfs_cpm_top,
        proportions = x_tfs_proportions
    )
    # list of different types of data, smaller
    sets = list(
        tfs_CV = x_tfs_cv
    )
    
    # list of plots and grid objects for plotting
    plots = list()



}
