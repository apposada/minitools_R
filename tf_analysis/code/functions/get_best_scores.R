get_best_scores = function (x,y) { # function accepting (i) a tf list of gene <--> tfclass ; (ii) an expression table with names in rownames and columns of values
    exp <- y
    for (i in 1:nrow(x)){
      gene <- x$id[i]
      genexp <- unlist(c(exp[grep(gene,rownames(exp)),]))
      genemean <- mean (genexp,na.rm=T)
      genesum <- sum(genexp)
      coefv <- sd(genexp,na.rm=T)/genemean
      x$mean[i] <- genemean
      x$counts[i] <- genesum
      x$cv[i] <- coefv
    }
    x$score <- (log(x$mean) + log(x$counts)) ^ x$cv # can be tweaked..
    x$score[is.nan(x$score)] <- 0
    x <- x[complete.cases(x),]
    tfs_scores <-
      aggregate(x$score,
                by = list(Category = x$class),
                FUN = sum)
    tfs_scores$score <- log(tfs_scores$x)
    tfs_scores <- tfs_scores[complete.cases(tfs_scores),]
    tfs_scores <- tfs_scores[rev(order(tfs_scores$score)),]
    return(tfs_scores)
  }
  