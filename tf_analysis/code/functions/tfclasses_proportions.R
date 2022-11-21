tfclasses_proportions <- function(x_tfs_cpm_scoring_top, nsamples, ...) {

    # STEP 1 : Aggregate expressed genes per class

    #create presence absence matrix of genes expressed
    x_tfs_cpm_ifelse <-  data.frame(
        ifelse(x_tfs_cpm_scoring_top[,1:nsamples] > 0, 1, 0),
        class = x_tfs_cpm_scoring_top$class
        )

    x_tfs_ngenes <-
    aggregate(
        x_tfs_cpm_ifelse[, 1:nsamples], # count amt of class "$class" expressed
        by = list(class = x_tfs_cpm_ifelse$class),
        FUN = sum
    )

    rownames(x_tfs_ngenes) <- 
    x_tfs_ngenes$class

    x_tfs_ngenes <- x_tfs_ngenes[
    order(x_tfs_ngenes$class), -1 # order, important for later divisions
    ]

    colnames(x_tfs_ngenes) <- colnames(x_tfs_cpm_scoring_top)
    x_tfs_ngenes[is.na(x_tfs_ngenes)] <- 0

    relative_ngenes <- 
        apply(x_tfs_ngenes, 2, function(x) {
        x / sum(x)
        })

    ### STEP 2 : aggregate tf cpms by tf class

    x_tfs_numcpms <-
        aggregate(
            x_tfs_cpm_scoring_top[, 1:nsamples], # aggregate cpms
            by = list(class = x_tfs_cpm_scoring_top$class),
            FUN = sum
            )

    rownames(x_tfs_numcpms) <- 
        x_tfs_numcpms$class

    x_tfs_numcpms <- x_tfs_numcpms[
        order(x_tfs_numcpms$class), # order, important for later divisions
    ]
    relative_numcpms <- 
        apply(x_tfs_numcpms, 2, function(x) { x / sum(x)})


    # STEP3: aggregate on a level of expression by gene, family basis

    x_tfs_expngenes <- 
        x_tfs_numcpms / x_tfs_ngenes
    x_tfs_expngenes <- apply(
        x_tfs_expngenes, 2, function(x) {
            ifelse(is.nan(x), 0, x)
            }
        )

    relative_exp_pergene_perclass <- 
        apply(x_tfs_expngenes, 2, function(x) {
        x / sum(x)
        })

    res <- list(
        numgenes = x_tfs_ngenes,
        relative_numgenes = relative_ngenes,
        numcpms = x_tfs_numcpms,
        relative_numcpms = relative_numcpms,
        expression_pergene_perclass = x_tfs_expngenes,
        relative_exp_pergene_perclass = relative_exp_pergene_perclass
    )

    return(res)
}


