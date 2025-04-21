#' qPCR Plot Fold Changes per Targeted Gene
#'
#' This function calculates mean and standard error of the mean, then plots the FoldChanges compared to your control group.
#'
#' @param df Dataframe imported from importFiles.
#' @param ctr Defines a character string with the name of your control group used to normalize delta-delta-Cts. Default is 'pcas13d'
#' @param outputPath Path to save files, plots and the report
#' @param allGroups TRUE or FALSE to determine if plots will include all conditions or only conditions that match the target name. Useful when you name your groups from knockdowns based on the target gene name. Default is TRUE.
#' @return Plots the technical controls (e.g., water, RT-) and saves plots and table in outuput path.
#' @export
#'
plotFoldChange <- function(df, ctr = "pcas13d", outputPath, allGroups = TRUE){
# Plots per gene
means_df5 <- df %>%
  group_by(.data$Group, .data$Gene, .data$Class_gene) %>%
  summarise(FC_Mean = mean(.data$FoldChange), FC_SD = std.error(.data$FoldChange))

# This will keep uninjected and defines ctr groups first in the plots
means_df5$Group <- factor(x = means_df5$Group,
                          levels = c(ctr, unique(filter(means_df5,
                                                        !.data$Group %in% c(ctr))$Group)))
if(allGroups == TRUE){
  temp2 <- list()
  for(i in unique(means_df5$Gene)) {
    cat("\nGene:", i, "\n")

      selection <- unique(as.character(means_df5$Group))
      selection <- unique(c(selection, tolower(ctr)))

    temp <-  means_df5 %>%
      filter(.data$Gene == i & .data$Group %in% selection) %>%
      ggplot(aes(x=.data$Group, y=.data$FC_Mean, fill=.data$Group)) +
      geom_bar(stat="identity", color="black") +
      geom_errorbar(aes(ymin=.data$FC_Mean-.data$FC_SD, ymax=.data$FC_Mean+.data$FC_SD), width=.2, color="black") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      scale_fill_manual(values = c(rep("black", sum((selection %in% c(tolower(ctr))))),
                                   rep("#AA4499", sum(!(selection %in% c(tolower(ctr))))))) +
      theme_classic() +
      theme(legend.position = "none",
            text = element_text(size = 10),
            axis.text = element_text(size = 10),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 12, hjust = 0.5)) +
      xlab("") +
      ylab("Fold Change") +
      ggtitle(i)

    print(temp)

    ggsave(plot = temp,
           filename = paste0(normalizePath(outputPath), "/", i, "_qPCR.png"),
           width = 4,
           height = 4)

    temp2 <- append(temp2, list(temp))
  }
} else {
temp2 <- list()
for(i in unique(means_df5$Gene)) {
  cat("\nGene:", i, "\n")

  # For some reason it does not put all groups from hsk genes... Will need to figure it out
  if(any(filter(means_df5, .data$Gene == i)$Class_gene == "exp")){
    selection <- unique(c(stringr::str_subset(means_df5$Group, i)))
    print("exp-only: subsetting Group")
  } else {
    selection <- unique(as.character(means_df5$Group))
    print("non-exp or mixed: keeping all Groups")
  }

  selection <- unique(c(selection, tolower(ctr)))

  temp <-  means_df5 %>%
    filter(.data$Gene == i & .data$Group %in% selection) %>%
    ggplot(aes(x=.data$Group, y=.data$FC_Mean, fill=.data$Group)) +
    geom_bar(stat="identity", color="black") +
    geom_errorbar(aes(ymin=.data$FC_Mean-.data$FC_SD, ymax=.data$FC_Mean+.data$FC_SD), width=.2, color="black") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = c(rep("black", sum((selection %in% c(tolower(ctr))))),
                                 rep("#AA4499", sum(!(selection %in% c("uninjected", tolower(ctr))))))) +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 12, hjust = 0.5)) +
    xlab("") +
    ylab("Fold Change") +
    ggtitle(i)

  print(temp)

  ggsave(plot = temp,
         filename = paste0(normalizePath(outputPath), "/", i, "_qPCR.png"),
         width = 4,
         height = 4)

  temp2 <- append(temp2, list(temp))

          }
      }
return(temp2)
}
