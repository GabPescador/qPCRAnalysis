# qPCRAnalysis: A package to analyze RT-qPCR data from QuantStudio Software.
It is specifically designed for the Bazzini Lab at Stowers Institute for Medical Research.

## Instalation
You can install the released version from GitHub with:
```
# install.packages("devtools")
devtools::install_github("GabPescador/qPCRAnalysis")
```

## Example
This is a basic example which generates the full report based on your inputs.
```
library(qPCRAnalysis)
qPCRreport(
     hsk = "cdk2ap2", # Housekeeping gene
     ctr = "pcas13d", # Control sample for normalization
     sdCutOff = 0.2, # Cutoff for technical replicates variation
     inputPath = system.file("extdata", "", package = "qPCRAnalysis"),
     outputPath = "./",
     techCtrs = defaultTechCtrs # Technical control sample names
     )
```

One example report and all output files can be found in:
```
system.file("extdata/output", "", package = "qPCRAnalysis")
```

