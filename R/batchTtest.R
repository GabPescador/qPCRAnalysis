#' qPCR Batch T-test
#'
#' This function performs T-test between your ctr group and all other
#'
#' @param df Dataframe imported from importFiles.
#' @param ctr Defines a character string with the name of your control group used to perform T-test in comparison to. Default is 'pcas13d'
#' @param hsk Defines a character string with the housekeeping gene used in the experiment. Default is "cdk2ap2"
#' @param outputPath Path to save files, plots and the report
#' @return Dataframe with T-test results from all pairwise comparisons.
#' @export
#'
batchTtest <- function(df, ctr = "pcas13d", hsk = "cdk2ap2", outputPath){

  # Getting rid of hsk gene since it will always have 0 variance
  df1 <- filter(df, !.data$Gene %in% hsk)

results <- list()
results_df_final <- data.frame()
control_group <- ctr
groups <- unique(filter(df1, !.data$Group %in% ctr)$Group)
for (gene in unique(df1$Gene)){
  gene_df <- filter(df1, .data$Gene == gene)

  for (group in groups) {
    control_data <- filter(gene_df, .data$Group == control_group)$FoldChange
    test_data <- filter(gene_df, .data$Group == group)$FoldChange

    if (length(control_data) > 1 && length(test_data) > 1) {
      test_result <- t.test(control_data, test_data)

      results[[paste(control_group, "vs", group)]] <- list(
        control = control_group,
        gene = gene,
        test = group,
        p.value = test_result$p.value,
        mean.control = mean(control_data),
        mean.test = mean(test_data),
        t.statistic = test_result$statistic)

      # Convert results to data.frame
      results_df <- do.call(rbind, lapply(names(results), function(x) {
        data.frame(Comparison = x,
                   Gene = results[[x]]$gene,
                   Control = results[[x]]$control,
                   Test = results[[x]]$test,
                   PValue = results[[x]]$p.value,
                   MeanControl = results[[x]]$mean.control,
                   MeanTest = results[[x]]$mean.test,
                   TStatistic = results[[x]]$t.statistic)
      }))

      results_df_final <- rbind(results_df_final, results_df)

    } else {
      print(paste0(gene, ": not enough data points for: ", control_group, " vs ", group))
    }
  }
}
# # Convert results to data.frame
# results_df <- do.call(rbind, lapply(names(results), function(x) {
#   data.frame(Comparison = x,
#              Gene = results[[x]]$gene,
#              Control = results[[x]]$control,
#              Test = results[[x]]$test,
#              PValue = results[[x]]$p.value,
#              MeanControl = results[[x]]$mean.control,
#              MeanTest = results[[x]]$mean.test,
#              TStatistic = results[[x]]$t.statistic)
# }))

# Again for some reason rows are getting duplicated
results_df_final <- unique(results_df_final)

# Save final table
write_csv(results_df_final,
          file = paste0(normalizePath(outputPath), "/Ttest_table.csv"))

return(results_df_final)

}
