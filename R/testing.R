library(testthat)

source("R/01_multipleAccessionHandling.r")

sample <- read.csv("data/sample_data.csv")
data_multiacc  <-  read.csv("data/sample_lfq_data.csv")
data_single_acc <- read.csv("data/quant_data.csv")

# test for count accession
expect_equal(
  object = countAccession(data_multiacc, "Accession", delimiter = ";"),
  expected = 10
)

expect_equal(
  object = countAccession(data_single_acc, "Accession", delimiter = ";"),
  expected = 1
)

# test for cleaning
expect_equal(
  object = unname(cleanAccessionColumn(data_multiacc[2,], "Accession",delimiter = ";", split_char = "\\|")[1,1]),
  expected = "Q3ZC07,A0A3Q1M558,G8JKX4,P60712,Q5E9B5,A0A452DJE4,A0A3Q1NKP5"
)

#############
# adding #remove to be removed



colnames(quant_data)[1] <- "accession"

gene_data <- quant_data
colName <- "accession"
delimiter = ";"
split_char = "\\|"
information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence"
verbose = FALSE

#getGeneInformationforMultipleAccession
getGeneInformationforMultipleAccession(
  gene_data_separated[1:10,], #remove [1:10,]
  colName, 
  delimiter, 
  information, 
  verbose_ = verbose
)


# populateGeneNames
a <- populateGeneNames(quant_data[1:10,], 
                  "accession", 
                  delimiter = ";",
                  split_char = "\\|",
                  information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence", 
                  verbose = FALSE)

 # populateGeneNamesfromMultipleAccession
populateGeneNamesfromMultipleAccession(
  gene_data, 
  colName, 
  delimiter, 
  information, 
  verbose = verbose
)
# 
# Error in `pivot_longer()`:                                                                                                                                                                
#   ! Names must be unique.
# ✖ These names are duplicated:
#   * "accession" at locations 1 and 38.
# ℹ Use argument `names_repair` to specify repair strategy.
# Run `rlang::last_trace()` to see where the error occurred.
# Called from: signal_abort(cnd)
# Warning message:
#   Debug mode activated: adding variables `accession_ok`, `accession_pieces`, and `accession_remainder`. 
# Error during wrapup: unimplemented type (29) in 'eval'
# 
# Error: no more error handlers available (recursive errors?); invoking 'abort' restart
# Error during wrapup: INTEGER() can only be applied to a 'integer', not a 'unknown type #29'
# Error: no more error handlers available (recursive errors?); invoking 'abort' restart

# getGeneInformationforMultipleAccession



# Error in `filter()`:                                                                                                                                                                      
#   ℹ In argument: `protein_name == "merged"`.
# Caused by error:
#   ! object 'protein_name' not found

# getGeneUniProt
## protein_name is becomin Protein_name
# as.vector(gene_data_separated[[current_col]])


# Add this before the separate_wider_delim
print("Number of max accessions:")
print(max_accessions)

print("Sample of data to be separated:")
print(head(gene_data[[colName]]))

# Also check the actual number of commas in each entry
comma_counts <- stringr::str_count(gene_data[[colName]], ",")
print("Distribution of comma counts:")
print(table(comma_counts))

### understand and include the last part from the normal function that is deleted


