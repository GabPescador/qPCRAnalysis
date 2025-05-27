#' qPCR Analysis Report
#'
#' This function process your qPCR data and combines all results and processing into an HTML report.
#' All middle steps are also saved in the output folder path.
#'
#' @param hsk Defines a character string with the name of your housekeeping gene used to normalize the qPCR data. Default is 'cdk2ap2'
#' @param ctr Defines a character string with the name of your control group used to normalize delta-delta-Cts. Default is 'cas13d'
#' @param sdCutOff Defines the sd cutoff of technical replicates to be accounted in the analysis. Default is 0.2
#' @param inputPath Path to Assay info, Sample info and Results file from QuantStudio 5 384-well instrument
#' @param outputPath Path to save files, plots and the report
#' @param techCtrs Character string that defines the names of your technical controls, e.g. 'water', 'RT-', etc. Default string is within 'R/defaultTechCtrs.R'
#' @param allGroups TRUE or FALSE to determine if plots will include all conditions or only conditions that match the target name. Useful when you name your groups from knockdowns based on the target gene name. Default is TRUE.
#' @param skiprows Number of rows to skip from the Results file. Default is 25
#' @return Plots, tables and HTML report are saved in the output folder path.
#' @export
#'
qPCRreport <- function(hsk = "cdk2ap2",
                 ctr = "pcas13d",
                 sdCutOff = 0.2,
                 inputPath,
                 outputPath,
                 techCtrs = defaultTechCtrs,
                 allGroups = TRUE,
                 skiprows = 25){

  # Checking that output directory exists, and if not create it
  # Check if the directory exists
  if (!dir.exists(outputPath)) {
    # If it doesn't exist, create the directory
    dir.create(outputPath, recursive = TRUE)  # recursive = TRUE ensures that any missing parent directories are also created
    print(paste("Directory created:", outputPath))
  } else {
    print(paste("Directory already exists:", outputPath))
  }

  ################### Importing module ###################
  print("Importing files...")
  l0 <- importFiles(inputPath = inputPath,
                    skiprows = skiprows)

  ################### Control module ###################
  print("Plotting technical controls...")
  l1 <- plotTechCtrs(l0$df2, techCtrs = techCtrs, outputPath = outputPath)

  ################### Checking for outliers module ###################
  print("Plotting possible technical outliers...")
  l2 <- flagOutliers(l0$df2, techCtrs = techCtrs, outputPath = outputPath,
                     assayInfo = l0$assay_info, sampleInfo = l0$sample_info)

  ################### Calculating delta-deltaCt module ###################
  print("Calculating delta-deltaCt values...")
  means_df4 <- deltaDeltaCtCalc(l2$means_df2, hsk = hsk, ctr = ctr, outputPath = outputPath)

  table5 <- batchTtest(means_df4, ctr = ctr, hsk = hsk, outputPath = outputPath)

  print("Plotting fold change results per target...")
  temp2 <- plotFoldChange(means_df4, ctr = ctr, outputPath = outputPath, allGroups = allGroups)

  p4 <- cowplot::plot_grid(plotlist = temp2)
  print("Generating report...")

  plot1 <- l1$p1
  plot2 <- l2$p2
  plot3 <- l2$p3
  plot4 <- p4
  table1 <- l1$dfCtrl
  table2 <- l2$dfExp
  table3 <- l2$means_df2
  table4 <- means_df4

  template <- system.file("rmd", "html_report_template.Rmd", package = "qPCRAnalysis")

  rmarkdown::render(
    input = template,
    output_file = "report.html",
    output_dir = normalizePath(outputPath),
    params = list(plot1 = plot1,
                  plot2 = plot2,
                  plot3 = plot3,
                  plot4 = plot4,
                  table1 = table1,
                  table2 = table2,
                  table3 = table3,
                  table4 = table4,
                  table5 = table5),
    envir = new.env(parent = globalenv())
  )

  message("HTML report saved to: ", normalizePath(outputPath))
  print(paste0("HTML report saved to: ", normalizePath(outputPath)))
}
