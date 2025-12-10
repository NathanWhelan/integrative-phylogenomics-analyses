library(dplyr)

##Command line structure: Rscript chi-square-test_noGQS-limit.R <GLS-table> <GQS-table>



# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if arguments were provided and assign file paths
if (length(args) == 2) {
  gls_file <- args[1]
  gqs_file <- args[2]
  cat(paste("Using command line arguments:", gls_file, "and", gqs_file, "\n"))
} else {
  cat(paste("No command line arguments provided. \n"))
  cat("Usage example: Rscript chi-square-test_noGQS-limit.R <GLS_File> <GQS_File>\n")
  quit()
}

# Load GLS Data

cat("--- Reading GLS Data ---\n")
gls_data <- read.table(
  gls_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
# Clean 'gene_id' for GLS data -
gls_data <- gls_data %>%
  mutate(
    gene_id = trimws(gene_id), # Crucial: remove any leading/trailing whitespace
    gene_id = sub("\\..*", "", gene_id) # Removes everything from the first dot onwards
  )
# Display the initial structure of the GLS data
cat("GLS data loaded and 'gene_id' cleaned for merging. Dimensions:", nrow(gls_data), "rows,", ncol(gls_data), "columns.\n")


# Load GQS Data
# The GQS file needs the header from the GLS file. 

header_names <- c("gene_id", "tree_supported", "tr1_value", "tr2_value", "df_value")

cat("\n--- Reading GQS Data ---\n")
gqs_data <- read.table(
  gqs_file,
  header = FALSE, 
  sep = "\t",
  stringsAsFactors = FALSE
)

# Assign the header names
colnames(gqs_data) <- header_names


# The 'gene_id' in GQS data contains file extensions (e.g., .faa.mafft.clipkit.treefile).
gqs_data <- gqs_data %>%
  mutate(
    gene_id = sub("\\..*", "", gene_id) # Removes everything from the first dot onwards
  )


# Display the initial structure of the GQS data
cat("GQS data loaded and header applied. Dimensions:", nrow(gqs_data), "rows,", ncol(gqs_data), "columns.\n")


# Add Absolute Value Columns 

# For GLS data
gls_data <- gls_data %>%
  mutate(
    abs_df_GLS = abs(df_GLS)
  )



cat("\nAbsolute value column 'abs_df_GLS' added to GLS data.\n")

# For GQS data
gqs_data <- gqs_data %>%
  mutate(
    abs_df_GQS = abs(df_value)
  )

cat("Absolute value column 'abs_df_GQS' added to GQS data.\n")

# Merge Data for Combined Calculation ---
merged_data <- inner_join(
  gls_data, 
  gqs_data, 
  by = "gene_id", 
  suffix = c("_GLS", "_GQS") 
)

cat("\nData merged successfully based on 'gene_id'. Total rows:", nrow(merged_data), "\n")


# Calculate the Number of Loci Meeting All Criteria 


result_count_tr1 <- merged_data %>%
  filter(
    # Criteria 1: Tree support in BOTH files must be "tr1"
    tree_supported_GLS == "tr1" & tree_supported_GQS == "tr1",
    
    # Criteria 2: New GLS column (abs_df_GLS) is 2 or greater
    abs_df_GLS >= 2,
    
    # Criteria 3: New GQS column (abs_df_GQS) is 0.00000001 or greater
    abs_df_GQS >= 0.00000001
  ) %>%
  # Count the rows that meet all the conditions
  nrow()
result_count_tr2 <- merged_data %>%
  filter(
    # Criteria 1: Tree support in BOTH files must be "tr1"
    tree_supported_GLS == "tr2" & tree_supported_GQS == "tr2",
    
    # Criteria 2: New GLS column (abs_df_GLS) is 2 or greater
    abs_df_GLS >= 2,
    
    # Criteria 3: New GQS column (abs_df_GQS) is 0.00000001 or greater
    abs_df_GQS >= 0.0000001
  ) %>%
  # Count the rows that meet all the conditions
  nrow()



total <- result_count_tr1 + result_count_tr2

# Perform the binomial test

binom_result <- binom.test(result_count_tr1, total, p = 0.5, alternative = "two.sided")


# Output Results 
cat("\n======================================================\n")
cat("Summary of Counts:\n")
cat("  - Count (tr1 criteria):", result_count_tr1, "\n")
cat("  - Count (tr2 criteria):", result_count_tr2, "\n")
cat("  - Total Count:", total, "\n")

cat("\nBinomial Test Results (tr1 vs. tr2):\n")
cat("  - Null Hypothesis (p=0.5): The proportion of tr1 is equal to the proportion of tr2 in the combined subset.\n")



cat("\n======================================================\n")
cat("Final Calculation Result:\n")
cat("Number of rows where:\n")
cat("  - Tree support ('tree_supported') is 'tr1' in BOTH GLS and GQS files", gls_file, "\n")
cat("  - Absolute delta in GLS ('abs_df_GLS') is >= 2.0, AND\n")
cat("  - Absolute delta in GQS ('abs_df_GQS') is >= 0.0.00000001\n")
cat("  - P-Value:", binom_result$p.value, "\n")
cat("\nResult Count tr1:", result_count_tr1, "\n")
cat("======================================================\n")

