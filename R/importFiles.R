#' qPCR Import Files
#'
#' This function imports assay info, sample info and results file from QuantStudio5 384-well instrument and
#' merges all tables into a single dataframe with all information.
#'
#' @param inputPath Path to Assay info, Sample info and Results file from QuantStudio 5 384-well instrument
#' @return Imports tables and merges them creating a final dataframe with all sample and target info.
#' @export
#'
importFiles <- function(inputPath){
# First imports results files and makes columns cleaner
raw_results_file <- list.files(path = inputPath,
                               pattern = "Results", full.names = TRUE)
df1 <- read_csv(raw_results_file,
                skip = 25, # This is based on the output from QuantStudio5
                col_select = c("Sample", "Target", "Cq", "Cq Mean", "Cq SD"))
df1$Cq <- as.numeric(df1$Cq)
df1$Sample <- as.numeric(str_split_fixed(df1$Sample, "s", n = Inf)[,2])

# Read sample info file and fixes possible errors with lower and upper cases
sample_info_file <- list.files(path = inputPath,
                               pattern = "sample_info", full.names = TRUE)
sample_info <- read_csv(sample_info_file)
sample_info$Name <- tolower(sample_info$Name)
sample_info$Group <- tolower(sample_info$Group)
colnames(sample_info) <- c("Sample", "Name", "Group")

# Read assay info file and fixes possible errors with lower and upper cases
assay_info_file <- list.files(path = inputPath,
                              pattern = "assay_info", full.names = TRUE)
assay_info <- read_csv(assay_info_file)
assay_info$Gene <- tolower(assay_info$Gene)
assay_info$Class_gene <- tolower(assay_info$Class_gene)
colnames(assay_info) <- c("Target", "Oligos", "Gene", "Class_gene")

# Finalize by merging all information into one table
df2 <- merge(df1, sample_info, by="Sample")
df2 <- merge(df2, assay_info, by="Target")

return(list(df2 = df2, sample_info = sample_info, assay_info = assay_info))
}
