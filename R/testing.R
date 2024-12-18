library(testthat)

source("R/01_multipleAccessionHandling.r")

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



