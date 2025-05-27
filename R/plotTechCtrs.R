#' qPCR Plot Technical Controls
#'
#' This function imports assay info, sample info and results file from QuantStudio5 384-well instrument.
#'
#' @param df Dataframe imported from importFiles.
#' @param techCtrs Character string that defines the names of your technical controls, e.g. 'water', 'RT-', etc. Default string is within 'R/defaultTechCtrs.R'
#' @param outputPath Path to save files, plots and the report
#' @return Plots the technical controls (e.g., water, RT-) and saves plots and table in outuput path.
#' @export
#'
plotTechCtrs <- function(df, techCtrs = defaultTechCtrs, outputPath){
# Filter out your technical controls
dfCtrl <- filter(df, .data$Group %in% tolower(techCtrs))
write_csv(dfCtrl, file = paste0(normalizePath(outputPath), "/tech_ctr_table.csv"))

if (all(is.na(dfCtrl$Cq))){
  print("All Cq values are NA. Skipping tech control plots...")
  return(list(p1 = NULL, dfCtrl = dfCtrl))
} else {
  # Plot Cts from controls
  p1 <- dfCtrl %>%
    ggplot(aes(x=.data$Name, y=.data$Cq)) +
    geom_boxplot() +
    geom_sina(color="#AA4499") +
    facet_wrap(facets = ~.data$Gene) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(p1, filename = paste0(normalizePath(outputPath), "/tech_ctrl_plot.png"),
         width = 6,
         height = 6)

  return(list(p1 = p1, dfCtrl = dfCtrl))
  }

}
