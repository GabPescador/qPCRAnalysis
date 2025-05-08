# qPCRAnalysis: A package to analyze RT-qPCR data from QuantStudio Software.
It is specifically designed for the Bazzini Lab at Stowers Institute for Medical Research.

## Instalation
You can install the released version from GitHub with:
```
# install.packages("devtools")
devtools::install_github("GabPescador/qPCRAnalysis")
```

## Required files
You should have a folder that contains 3 .csv files:

1. Results.csv - the output from QuantStudio software with your Cq values
2. sample_info.csv - a table with your sample information
3. assay_info.csv - a table with your assays (oligo pairs) information

It is important that the name of the files contain the names in the examples, e.g., results file can be named as: "20250401_target_Results.csv" since the function will look for files that have "Results" in the file names to import.

### Results.csv example

*Note:* The import function will skip all the metadata generated from the equipment run (the "# File Name:", etc rows), if the number of rows changes the import will not work properly.

```
# File Name: Run_428277102.eds
# Comment:
# Operator: RemoteUser
# Barcode: C3C01IO4
# Instrument Type: QuantStudio™ 7 Pro System
# Block Type: 384-Well Block
# Instrument Name: QS7Pro-2778721020029
# Instrument Serial Number: 2778721020029
# Heated Cover Serial Number:
# Block serial number:
# Run Start Date/Time: 2024-02-12 09:52:38 PM UTC
# Run End Data/Time: 2024-02-12 11:30:22 PM UTC
# Run Duration: 97 minutes 43 seconds
# Sample Volume: 10.0
# Cover Temperature: 105.0
# Passive Reference: ROX
# PCR Stage/Step Number: Stage 2 Step 2
# Melt Stage Number: 3
# Quantification Cycle Method: CT
# Analysis Date/Time: 2024-02-13 06:45:43 PM UTC
# Software Name and Version: Design & Analysis Software v2.7.0
# Plugin Name and Version: Primary Analysis v1.8.0
# Reduce dye signal crosstalk by algorithm: No (Default)
# Exported By: guest
# Exported On: 2024-02-13 06:46:04 PM UTC
Well	Well Position	Omit	Sample	Target	Task	Reporter	Quencher	Amp Status	Amp Score	Curve Quality	Result Quality Issues	Cq	Cq Confidence	Cq Mean	Cq SD	Auto Threshold	Threshold	Auto Baseline	Baseline Start	Baseline End	Tm1	Tm2	Tm3	Tm4
1	A1	FALSE	s01	1	UNKNOWN	SYBR	None	INCONCLUSIVE	1.171677	NA	NA	22.38276809	0.6045832	22.35764	0.2016855	FALSE	0.1730069	TRUE	3	17	81.25739	NA	NA	NA
2	A2	FALSE	s01	1	UNKNOWN	SYBR	None	AMP	1.193406	NA	NA	22.14456235	0.6381591	22.35764	0.2016855	FALSE	0.1730069	TRUE	3	17	80.99167	NA	NA	NA
3	A3	FALSE	s01	1	UNKNOWN	SYBR	None	AMP	1.164003	NA	NA	22.54557773	0.5889061	22.35764	0.2016855	FALSE	0.1730069	TRUE	3	18	81.25739	NA	NA	NA
4	A4	FALSE	s17	1	UNKNOWN	SYBR	None	NO_AMP	0.000000	NA	NA	UNDETERMINED	0.0000000	NA	NA	FALSE	0.1730069	TRUE	3	39	88.03301	83.78163	79.39740	73.95033
5	A5	FALSE	s17	1	UNKNOWN	SYBR	None	NO_AMP	0.000000	NA	NA	UNDETERMINED	0.0000000	NA	NA	FALSE	0.1730069	TRUE	3	39	76.85406	80.43716	83.75482	71.54580
6	A6	FALSE	s17	1	UNKNOWN	SYBR	None	AMP	1.073138	NA	NA	37.24113127	0.9399698	37.24113	NA	FALSE	0.1730069	TRUE	3	32	75.26159	87.07249	89.85934	NA
```

### sample_info.csv example

Sample should be the same number as the Sample column from Results.csv, but without the "s0".

```
Sample	Name	Group
1	uninjected_rep1	uninjected
2	uninjected_rep2	uninjected
3	uninjected_rep3	uninjected
4	pCas13d_rep1	pCas13d
5	pCas13d_rep2	pCas13d
6	pCas13d_rep3	pCas13d
```

### assay_info.csv example

Target should be the same number as the Target column from "Results.csv".

```
Target	Oligos	Gene	Class_gene
1	946/947	cdk2ap2	hsk
2	5209/5210	znf_3prime	exp
3	944/945	taf15	exp
4	4615/4616	znf_cds	exp
```

## Example

This is a basic example which generates the full report based on your input files.

- "hsk" should match the housekeeping "Gene" column in assay_info.csv
- "ctr" should match the control "Group" column in sample_info.csv

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
"/inst/extdata/output"
```

