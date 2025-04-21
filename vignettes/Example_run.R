devtools::load_all("./")

test <- importFiles("/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/")
l1 <- plotTechCtrs(test$df2, techCtrs = defaultTechCtrs, outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/")
l2 <- flagOutliers(test$df2, techCtrs = defaultTechCtrs,
                   outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/",
                   assayInfo = test$assay_info, sampleInfo = test$sample_info)
df4 <- deltaDeltaCtCalc(l2$means_df2,outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/")
a <- plotFoldChange(df4, outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/")

outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/"
ctr = "pcas13d"

qPCRreport(hsk = "cdk2ap2",
     ctr = c("pcas13d","uninjected"),
     sdCutOff = 0.2,
     inputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/",
     outputPath = "/n/sci/SCI-004255-ZFPROT/gd2417/Candidates_KDs/qPCR/2023_06_12_qPCR_Candidates_2/test/",
     techCtrs = defaultTechCtrs)
