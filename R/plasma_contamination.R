#' Ipsum...
#'
#' @param df ipsum...
#' @import decoupleR
#' @export
#' @return ipsum...
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
