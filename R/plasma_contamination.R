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
plasma_contamination <- function(df, n_perm = 1000, plot = T) {
  if(dim(df)[2] == 1)
  {
    df$decouplerIsGreat <- df[,1]
  }
  means <- decoupleR::decouple(mat = as.matrix(df), network = reg_sets, .source = "set", .target = "gene",
                              statistics = c("mean"),
                              args = list(
                                mean = list(.mor = "mor", .likelihood = "likelihood", times=n_perm)
                              )
  )
  means <- as.data.frame(means)
  means <- means[grepl("normalized_mean",means$statistic),]
  means <- means[!grepl("^decouplerIsGreat$",means$condition),]
  means <- means[,-c(1,2,6)]
  names(means)[1] <- "set"

  contamination_scores <- means

  if(plot)
  {
    plot(ggplot(contamination_scores, aes(y = score, x = condition, group = set, fill = set)) +
           geom_bar(position = "dodge", stat = "identity", color="black") +
           theme_minimal())
  }
  return(contamination_scores)
}
