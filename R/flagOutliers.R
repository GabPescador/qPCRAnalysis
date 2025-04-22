#' qPCR Flag Technical Outliers
#'
#' This function flags technical outliers based on the Ct standard deviation from sdCutOff. As deviations in technical
#' replicates typically arise from pipetting error, we can exclude them from skewing the data points.
#'
#' @param df Dataframe imported from importFiles.
#' @param techCtrs Character string for technical controls (e.g., water, RT-). Default is within defaultTechCtrs.R
#' @param outputPath Path to save files, plots and the report
#' @param sdCutOff Defines the sd cutoff of technical replicates to be accounted in the analysis. Default is 0.2
#' @param assayInfo Path to assay info csv file
#' @param sampleInfo Path to sample info csv file
#' @return Flags possible outliers based on sdCutOff, filters them out, recalculates Ct means and plots before and after Ct values.
#' @export
#'
flagOutliers <- function(df,
                         techCtrs = defaultTechCtrs,
                         outputPath,
                         sdCutOff = 0.2,
                         assayInfo,
                         sampleInfo){
# Filter out your experimental datapoints
dfExp <- filter(df, !.data$Group %in% tolower(techCtrs))

# Calculating square distances of the mean
dfExp$dMean <- (dfExp$Cq-dfExp$`Cq Mean`)^2

# Making a column that gives rank order of the dMean column for each target/sample combination
dfExp <- dfExp %>%
  group_by(.data$Sample, .data$Target) %>%
  mutate(MeanRank = rank(.data$dMean))

# Making an outlier column
dfExp$Outlier <- ifelse((dfExp$MeanRank > 2 & dfExp$`Cq SD` >= sdCutOff),
                        "Yes",
                        "No")

# Saving table of outlier technical replicates
write_csv(dfExp,
          file = paste0(normalizePath(outputPath), "/flagged_results_table.csv"))

# Plotting outliers
p2 <- dfExp %>%
  ggplot(aes(x=.data$Sample, y=.data$Cq, color=.data$Outlier)) +
  geom_point() +
  facet_wrap(facets = ~.data$Gene) +
  scale_color_manual(values = c("#DDCC77", "#AA4499")) +
  coord_cartesian(ylim = c(20, 30))

ggsave(p2, filename = paste0(normalizePath(outputPath), "/outlier_plot.png"),
       width = 6,
       height = 6)

print("Recalculating means without technical outliers...")
# Recalculating means without technical outliers
means_df <- filter(dfExp, .data$Outlier == "No")

means_df <- means_df %>%
  group_by(.data$Target, .data$Sample) %>%
  summarise(`Cq Mean` = mean(.data$Cq), `Cq SD` = sd(.data$Cq))

if (is.character(assayInfo) && file.exists(assayInfo)) {
  assay_info <- readr::read_csv(assayInfo) # If it's a path to a csv file use read_csv
} else if (is.data.frame(assayInfo)) {
  assay_info <- assayInfo # If it is a loaded dataframe just import as a new object
} else {
  stop("`assayInfo` must be a path to a CSV file or a data frame.")
}

if (is.character(sampleInfo) && file.exists(sampleInfo)) {
  sample_info <- readr::read_csv(sampleInfo)
} else if (is.data.frame(sampleInfo)) {
  sample_info <- sampleInfo
} else {
  stop("`sampleInfo` must be a path to a CSV file or a data frame.")
}

means_df2 <- merge(means_df, assay_info, by="Target")
means_df2 <- merge(means_df2, sample_info, by="Sample")

# Saving the new table without outliers
write_csv(means_df2,
          file = paste0(normalizePath(outputPath), "/cleaned_results_table.csv"))

print("Plotting Cts without outliers...")
# Plot cts without outliers
p3 <- means_df2 %>%
  ggplot(aes(x=.data$Group, y=.data$`Cq Mean`, color=.data$Group)) +
  geom_boxplot() +
  geom_sina() +
  facet_wrap(facets = ~Gene) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave(p3, filename = paste0(normalizePath(outputPath), "/clean_ct_plot.png"),
       width = 6,
       height = 6)

return(list(p2 = p2, dfExp = dfExp,
            p3 = p3, means_df2 = means_df2))

}
