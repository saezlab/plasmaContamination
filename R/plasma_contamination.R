#' @title plasma_contamination
#'
#' @description This function performs an enrichment analysis to test whether a given
#' proteomic sample is contaminated with markers from erythrocyte, platelets
#' or coagulation. The signatures are taken from:
#' https://www.sciencedirect.com/science/article/pii/S2405471216300722
#'
#' @param df as data frame with row names as gene symbols, and then each column
#' should be a sample or the t-values of a differential analysis
#' @import decoupleR
#' @import ggplot2
#' @export
#' @return a data frame with normalized enrichment (mean) scores for the
#' contamination signatures. The higher the score, the higher the potential
#' contamination
#' @examples
#' library(plasmaContamination)
#' data(test_data)
#'
#' conta_scores <- plasma_contamination(test_data)
plasma_contamination <- function(df, n_perm = 1000, plot = T, verbose = T, signatures = NULL, manual_colors = NULL) {
  if(is.null(signatures))
  {
    data("reg_sets")
    signatures <- reg_sets
  }
  if(dim(df)[2] == 1)
  {
    if(verbose)
    {
      print("Computing scores for 1 sample/contrast...")
    }
    df <- df[complete.cases(df),,drop = F]
    df$decouplerIsGreat <- df[,1]
    means <- decoupleR::decouple(mat = as.matrix(df), network = signatures, .source = "set", .target = "feature",
                                 statistics = c("wmean"),
                                 args = list(
                                   mean = list(.mor = "mor", .likelihood = "likelihood", times=n_perm)
                                 )
    )
    means <- as.data.frame(means)
    means <- means[grepl("norm_wmean",means$statistic),]
    means <- means[!grepl("^decouplerIsGreat$",means$condition),]
    means <- means[,-c(1,2,6)]
    names(means)[1] <- "set"

    score_df <- means
  } else
  {
    score_list <- list()
    for(i in 1:length(df[1,]))
    {
      if(verbose)
      {
        print(paste0("Computing scores for sample ", i, " out of ",length(df[1,]),"..."))
      }
      sub_df <- df[,i,drop = F]
      sub_df <- sub_df[complete.cases(sub_df),,drop = F]
      sub_df$decouplerIsGreat <- sub_df[,1]
      means <- decoupleR::decouple(mat = as.matrix(sub_df), network = signatures, .source = "set", .target = "feature",
                                   statistics = c("wmean"),
                                   args = list(
                                     mean = list(.mor = "mor", .likelihood = "likelihood", times=n_perm)
                                   )
      )
      means <- as.data.frame(means)
      means <- means[grepl("norm_wmean",means$statistic),]
      means <- means[!grepl("^decouplerIsGreat$",means$condition),]
      means <- means[,-c(1,2,6)]
      names(means)[1] <- "set"
      score_list[[i]] <- means
    }
    score_df <- as.data.frame(do.call(rbind,score_list))
  }
  if(plot)
  {
    if(is.null(manual_colors))
    {
      plot(ggplot(score_df, aes(y = score, x = condition, group = set, fill = set)) +
             geom_bar(position = "dodge", stat = "identity", color="black") +
             theme_minimal() +
             theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
             theme(plot.margin = unit(c(1,1,1,2), "cm")))
    } else
    {
      plot(ggplot(score_df, aes(y = score, x = condition, group = set, fill = set)) +
             geom_bar(position = "dodge", stat = "identity", color="black") +
             theme_minimal() +
             theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
             theme(plot.margin = unit(c(1,1,1,2), "cm")) +
             scale_fill_manual(values=manual_colors))

    }


  }
  return(score_df)
}
