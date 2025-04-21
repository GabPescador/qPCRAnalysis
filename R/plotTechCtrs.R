#' qPCR Plot Technical Controls
#'
#' This function imports assay info, sample info and results file from QuantStudio5 384-well instrument.
#'
#' @param df Dataframe imported from importFiles.
#' @param techCtrs Path to Assay info, Sample info and Results file from QuantStudio 5 384-well instrument
#' @param outputPath Path to save files, plots and the report
#' @return Plots the technical controls (e.g., water, RT-) and saves plots and table in outuput path.
#' @export
#'
plotTechCtrs <- function(df, techCtrs, outputPath){
# Filter out your technical controls
dfCtrl <- filter(df, .data$Name %in% tolower(techCtrs))
write_csv(dfCtrl, file = paste0(normalizePath(outputPath), "/tech_ctr_table.csv"))

# Plot Cts from controls
p1 <- dfCtrl %>%
  ggplot(aes(x=.data$Name, y=.data$Cq)) +
  geom_boxplot() +
  geom_sina(color="#AA4499") +
  facet_wrap(facets = ~.data$Gene) +
  theme_bw()

ggsave(p1, filename = paste0(normalizePath(outputPath), "/tech_ctrl_plot.png"),
       width = 6,
       height = 6)

return(list(p1 = p1, dfCtrl = dfCtrl))
}
