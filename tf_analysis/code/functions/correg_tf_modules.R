#' Correg TF classes and network modules
require(ComplexHeatmap)
require(circlize)
require(gplots)

correg_tf_modules <- function(x_tfs_cpm, x_tfs, x_wg_kME, min_tfs_class = 4, min_rho = 0.5){
    # Data tidy
    x_modules <- sort(gsub("ME", "", colnames(x_wg_kME)))

    # keep classes with more than 4 genes per class
    x_class_correg <- names(
        table(x_tfs$class)[
            table(x_tfs$class) > min_tfs_class
            ]
    )

    x_tfs_correg <- x_tfs[
        x_tfs$class %in% x_class_correg,
    ]

    x_tfs_modules_cor <- x_wg_kME[
        rownames(x_wg_kME) %in% x_tfs$id,
    ]

    colnames(x_tfs_modules_cor) <- x_modules

    x_tfs_modules <- matrix(
        0,
        nrow = length(x_class_correg),
        ncol = length(x_modules)
    )

    rownames(x_tfs_modules) <- x_class_correg

    colnames(x_tfs_modules) <- x_modules

    x_tfs_modules_abs <- x_tfs_modules

    for (i in x_modules) {
        h <- which(x_modules == i)

        for (j in x_class_correg) {
        n <- which(x_class_correg == j)
        all_tfs_class_j <- x_tfs$id[x_tfs$class == j]

        cors_jclass_in_imodule <-
            x_tfs_modules_cor[
                rownames(x_tfs_modules_cor) %in% all_tfs_class_j,
                colnames(x_tfs_modules_cor) == i
                ] # raise an eyebrow

        good_tfs_jclass_in_imodule <-
            length(
                cors_jclass_in_imodule[cors_jclass_in_imodule > min_rho]
                ) 

        x_tfs_modules[n, h] <-
            good_tfs_jclass_in_imodule / length(cors_jclass_in_imodule)
        }
        x_tfs_modules_abs[n, h] <-
            good_tfs_jclass_in_imodule
    }

    # Heatmap parameters
    viridis_pastel <- c("#ffee61","#96e88e","#5dc9ac","#4da2ba","#6b6eab","#552761")
    col_heat_cor_tfs <- colorRamp2(
        seq(0, min_rho, len=20),
        colorRampPalette(rev(viridis_pastel))(20)
    )

    modules_ha_colors <- gplots::col2hex(x_modules)
    names(modules_ha_colors) <- x_modules

    modules_ha <- HeatmapAnnotation(
        modules = colnames(x_tfs_modules),
        col = list(modules = modules_ha_colors), annotation_name_side = "left"
        )

    # Heatmap
    correg_hm <- Heatmap(
        name = paste0("% cor(Sp) > ", min_rho),
        x_tfs_modules,
        col = col_heat_cor_tfs,
        heatmap_legend_param = list(
        title = expression(
            paste0("% ", rho, " > ", min_rho)
            ), # check paste0 does not break expression()
        title_position = "leftcenter-rot",
        at = c(0,0.2,0.4,0.6), # must adapt to min_rho parameter
        labels = c("0","20","40","60+") # must adapt to min_rho parameter
        ),
        row_names_side = "right",
        row_names_gp = gpar(cex = 0.7),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        cluster_rows = TRUE,
        column_title = "Coexpression of TF classes and\ngene modules",
        column_title_gp = gpar(cex = 0.8,fontface = "bold"),
        # left_annotation=tfs_hb # TF annotation by color
        top_annotation = modules_ha # WGCNA modules
    )

    res <- list(
        correg_tfs_modules_fraction = x_tfs_modules,
        correg_tfs_modules_total = x_tfs_modules_abs,
        correg_hm = correg_hm
    )

}

