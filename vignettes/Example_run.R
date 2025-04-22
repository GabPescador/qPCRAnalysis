devtools::load_all("./")

test <- importFiles("/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/")
l1 <- plotTechCtrs(test$df2, techCtrs = defaultTechCtrs, outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/")
l2 <- flagOutliers(test$df2, techCtrs = defaultTechCtrs,
                   outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/",
                   assayInfo = test$assay_info, sampleInfo = test$sample_info)
df4 <- deltaDeltaCtCalc(l2$means_df2,outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/")
hsk = "cdk2ap2"


# toy <- data.frame(Gene = c(rep("znf281b", 7)),
#            Group = c(rep("znf281b_g3", 3),
#                      rep("pcas13d", 3),
#                      "znf281b_g2"),
#            Class_gene = c(rep("exp", 7)),
#            FoldChange = c(runif(3, 0.1, 0.2),
#                           runif(3, 0.9, 1),
#                           0.2))
#
# test <- batchTtest(toy, outputPath = outputPath)

a <- plotFoldChange(df4, outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/")

outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/"
ctr = "pcas13d"

qPCRreport(hsk = "cdk2ap2",
     ctr = "pcas13d",
     sdCutOff = 0.2,
     inputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/",
     outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/",
     techCtrs = defaultTechCtrs)

inputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2024_02_13_qPCR_znf281b-g1-3_pCas13d/20230212_znf281b-g1-3_KDs_Copy/"
outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2024_02_13_qPCR_znf281b-g1-3_pCas13d/test/"

qPCRreport(hsk = "cdk2ap2",
           ctr = "pcas13d",
           sdCutOff = 0.2,
           inputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2024_02_13_qPCR_znf281b-g1-3_pCas13d/20230212_znf281b-g1-3_KDs_Copy/",
           outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2024_02_13_qPCR_znf281b-g1-3_pCas13d/test/",
           techCtrs = defaultTechCtrs)

l0 <- importFiles(inputPath = inputPath)
l1 <- plotTechCtrs(l0$df2, techCtrs = defaultTechCtrs, outputPath = outputPath)
# dfCtrl <- filter(l0$df2, .data$Name %in% tolower(defaultTechCtrs))

a <- read_csv(paste0(outputPath, "FoldChange_table.csv"))
a <- a %>%
  filter(str_detect(Gene, "znf281b|taf15") & !Group %in% defaultTechCtrs)

