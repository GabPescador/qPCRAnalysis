#' qPCR delta delta Ct calculation
#'
#' This function calculates the delta delta Ct values based on the specified housekeeping gene.
#'
#' @param df Dataframe imported from importFiles.
#' @param hsk Defines a character string with the name of your housekeeping gene used to normalize the qPCR data. Default is 'cdk2ap2'
#' @param ctr Defines a character string with the name of your control group used to normalize delta-delta-Cts. Default is 'pcas13d'
#' @param outputPath Path to save files, plots and the report
#' @return Outputs a dataframe with delta delta Ct values based on each target.
#' @export
#'
deltaDeltaCtCalc <- function(df, hsk = "cdk2ap2", ctr = "pcas13d", outputPath){
# Getting means of housekeeper
hsk_means <- filter(df, .data$Gene %in% tolower(hsk)) %>%
  select(.data$`Cq Mean`, .data$Name) %>%
  mutate(hsk_cq = .data$`Cq Mean`) %>%
  select(.data$Name, .data$hsk_cq)

means_df3 <- merge(df, hsk_means, by="Name")

# Calculating delta CT
means_df3$dCT <- means_df3$`Cq Mean` - means_df3$hsk_cq

# Calculating delta delta CT
# First get means of delta CT
dCT_means <- means_df3 %>%
  group_by(.data$Group, .data$Gene) %>%
  mutate(dCT_mean = mean(.data$dCT))
# Then, select the means of control group and put it back in the means table
means_df4 <- dCT_means %>%
  filter(.data$Group %in% tolower(ctr)) %>%
  select(.data$Gene, .data$Group, .data$dCT_mean) %>%
  merge(x = means_df3, by="Gene")
# Calculate ddCT
means_df4$ddCT <- means_df4$dCT - means_df4$dCT_mean
# Calculate Fold Change
means_df4$FoldChange <- 2^(-means_df4$ddCT)

# Cleaning table
means_df4 <- means_df4 %>%
  select(!.data$Group.y) %>%
  rename("Group" = "Group.x")

# Save final table
write_csv(means_df4,
          file = paste0(normalizePath(outputPath), "/FoldChange_table.csv"))

return(means_df4)
}
