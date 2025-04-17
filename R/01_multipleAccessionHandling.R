# Gene Data Processing Toolkit
# A comprehensive R package for cleaning, processing, and enriching gene identification data
# using UniProt API and advanced data manipulation techniques

# Required Libraries
library(tidyverse)
library(httr)
library(jsonlite)
library(progress)
library(stringr)

# 1. Gene ID Cleaning and Extraction Functions ---------------------------

# helper function
#' Clean and Extract Gene Accessions
#' 
#' Extracts clean accession IDs from complex string formats
#' 
#' @param input_string A string containing gene accessions
#' @param delimiter Delimiter separating different gene entries (default: ";")
#' @param split_char Character used to split each gene entry (default: "|")
#' @return Vector of clean accession IDs
#' @export
cleanAccessionExtractor <- function(input_string,
                                    delimiter, 
                                    split_char
) {
  # Handle NA or empty input
  if (is.na(input_string) || input_string == "") {
    return(NA_character_)
  }
  
  # Split the input string into individual gene entries
  gene_entries <- unlist(strsplit(input_string, delimiter))
  
  # Extract first part of each entry (before split_char)
  clean_accessions <- sapply(gene_entries, function(entry) {
    # Split the entry and take the first part
    parts <- unlist(strsplit(entry, split_char))
    return(parts[1])
  })
  
  # Remove any potential NA values and get unique accessions
  clean_accessions <- unique(clean_accessions[!is.na(clean_accessions)])
  
  # Return as comma-separated string or NA if no valid accessions
  if (length(clean_accessions) > 0) {
    return(paste(clean_accessions, collapse = ","))
  } else {
    return(NA_character_)
  }
}

# Wrapper function for use in data manipulation
#' Clean Accessions in a Data Frame Column
#' 
#' Applies accession cleaning to a specific column in a data frame
#' 
#' @param data Data frame containing gene data
#' @param column Name of the column with accession information
#' @param delimiter Delimiter separating different gene entries
#' @param split_char Character used to split each gene entry
#' @return Data frame with cleaned accession column
#' @export
cleanAccessionColumn <- function(data, 
                                 column, 
                                 delimiter = ";", 
                                 split_char = "\\|") {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  if (!column %in% names(data)) {
    stop("Specified column does not exist in the data frame")
  }
  
  # Apply cleaning function to the specified column
  data %>%
    mutate(!!column := sapply(.[[column]], 
                              cleanAccessionExtractor, 
                              delimiter = delimiter, 
                              split_char = split_char))
}

# 2. Accession Counting and Processing Functions --------------------------

#' Count Maximum Number of Gene IDs in a Column
#' 
#' Determines the maximum number of gene IDs in a single cell of a given column
#' 
#' @param gene_data Data frame containing gene data
#' @param colName Column name with gene IDs
#' @param delimiter Regular expression delimiter between gene IDs
#' @param endsWith Logical, whether delimiter is at the end of the string
#' @param verbose Logical, whether to display additional information
#' @return Maximum number of gene IDs in a cell
#' @export
#' Count Maximum Number of Delimited Entries in a Column
#' 
#' Determines the maximum number of entries in a column when split by a delimiter
#' 
#' @param gene_data Data frame containing gene data
#' @param colName Column name with entries to count
#' @param delimiter Regular expression delimiter between entries
#' @param endsWith Logical, whether to count delimiter at the end of the string
#' @return Maximum number of entries in a cell
#' @export
countAccession <- function(gene_data, colName, delimiter = "", endsWith = FALSE) {
  # Input validation
  if (!is.data.frame(gene_data)) {
    stop("Input must be a data frame.")
  }
  
  if (!colName %in% names(gene_data)) {
    stop("Specified column does not exist in the data frame.")
  }
  
  # Handle empty delimiter case
  if (nchar(delimiter) == 0) {
    return(1)  # If no delimiter, assume single entry
  }
  
  # Use vectorized approach
  counts <- sapply(gene_data[[colName]], function(x) {
    # If the string is NA or empty, return 0
    if (is.na(x) || nchar(x) == 0) {
      return(0)
    }
    
    # Count matches
    matches <- gregexpr(delimiter, x)[[1]]
    
    # Adjust count based on endsWith parameter
    if (matches[1] == -1) {
      # No delimiter found
      count <- 1
    } else {
      # Count delimiters and adjust
      count <- length(matches)
      
      if (!endsWith) {
        # Add 1 to count for number of entries
        count <- count + 1
      }
    }
    
    return(count)
  })
  
  # Return maximum count
  max(counts, na.rm = TRUE)
}

# 3. UniProt Data Retrieval and Processing Functions ---------------------

#' Retrieve Merged Gene Names (helper function)
#' 
#' Identifies gene IDs that have been merged or replaced in UniProt
#' 
#' @param accession_ Vector of gene accession numbers
#' @return Data frame with original and new accession numbers
#' @export
getMergedGeneNames <- function(accession_) {
  new_accession_list <- c()
  base_url <- "https://www.uniprot.org/uniprot/"
  
  new_accession_list <- sapply(accession_, function(acc) {
    query_url <- paste0(base_url, acc, ".json")
    response <- GET(query_url)
    
    if (status_code(response) != 200) {
      data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      return(ifelse(!is.null(data$path), sub(".*/", "", data$path), NA))
    }
    
    return(NA)
  })
  
  data.frame(accession = accession_, new_Accession = new_accession_list)
}

#' Fetch Individual UniProt Data (helper function)
#' 
#' Retrieves data for a single accession from UniProt
#' 
#' @param accession_ Single UniProt accession number
#' @param information Desired information fields
#' @param base_url Base UniProt API URL
#' @param verbose Logical, whether to show detailed information
#' @return Data frame with UniProt information or NULL
#' @export
fetchUniprotData <- function(accession_, 
                             information, 
                             base_url = "https://rest.uniprot.org/uniprotkb/search?query=accession:", 
                             verbose_ = FALSE) {
  # Construct and encode query URL
  query_url <- paste0(base_url, accession_, "&format=tsv&fields=", information)
  query_url <- URLencode(query_url)
  
  # Make API request with error handling
  tryCatch({
    response <- httr::GET(query_url)
    
    # Validate response
    if (response$status_code != 200) {
      if (verbose_) {
        message("Failed to fetch data for accession: ", accession_, 
                ". HTTP status: ", response$status_code)
      }
      return(NULL)
    }
    
    # Parse response
    protein_data <- read.csv(text = httr::content(response, "text", encoding = "UTF-8"), 
                             header = TRUE, 
                             sep = "\t", 
                             stringsAsFactors = FALSE)
    
    return(protein_data)
  }, 
  error = function(e) {
    if (verbose_) {
      message("Error processing accession ", accession_, ": ", e$message)
    }
    return(NULL)
  })
}

#' Retrieve Gene Information from UniProt (wrapper function)
#' 
#' Fetches detailed information about gene IDs from the UniProt database
#' 
#' @param accession_list Vector of UniProt accession numbers
#' @param information Comma-separated string of desired information fields
#' @param verbose Logical, whether to show progress and additional information
#' @return Data frame with retrieved gene information
#' @export
getGeneUniProt <- function(accession_list, 
                           information, 
                           verbose = FALSE) {
  message("Processing accessions, please wait...")
  
  # Initialize progress bar
  progress_bar <- progress_bar$new(total = length(accession_list))
  progress_bar$tick(0)
  
  # Retrieve and combine data for each accession
  combined_data <- lapply(accession_list, function(accession_) {
    protein_data <- fetchUniprotData(accession_, information, verbose_ = verbose)
    progress_bar$tick()
    
    if (!is.null(protein_data)) {
      protein_data[1, ]
    } else {
      NULL
    }
  }) %>% 
    bind_rows()
  
  return(combined_data)
}

# 4. Multiple Accession Processing Functions -----------------------------

#' Get Gene Information for Multiple Accessions (helper))
#' 
#' Retrieves and processes UniProt information for multiple accession numbers
#' 
#' @param gene_data_separated Data frame with separated gene data
#' @param colName Column name containing gene IDs
#' @param delimiter Delimiter separating multiple gene IDs
#' @param information UniProt information fields to retrieve
#' @param verbose Logical, whether to show detailed information
#' @return Processed UniProt information
#' @export
getGeneInformationforMultipleAccession <- function(gene_data_separated, 
                                                   colName, 
                                                   delimiter, 
                                                   information, 
                                                   verbose_ = FALSE) {
  # Initialize empty data frame for UniProt names
  uniprot_cols <- unlist(strsplit(information, ","))
  uniprot_names <- data.frame(matrix(nrow = 0, ncol = length(uniprot_cols)))
  colnames(uniprot_names) <- uniprot_cols
  
  # Count number of Accession columns
  no_accession <- sum(grepl("accession", colnames(gene_data_separated)))
  
  # Process accessions iteratively
  for (n in 1:(no_accession-1)) {
    current_col <- paste0(colName, ".", n)
    
    # Retrieve UniProt information
    current_uniprot_names <- getGeneUniProt(
      as.vector(gene_data_separated[[current_col]]),
      information
    )
    
    # Capitalize column names
    colnames(current_uniprot_names) <- uniprot_cols # change  str_to_title(uniprot_cols)
    rownames(current_uniprot_names) <- NULL
    
    # Combine results
    uniprot_names <- rbind(uniprot_names, current_uniprot_names)
    
    # Find entries with NA Gene_primary
    na_gene_primary <- which(is.na(uniprot_names$gene_primary))
    
    # if (length(na_gene_primary) == 0) break
    
    # Prepare next iteration's data
    gene_data_separated <- gene_data_separated[gene_data_separated[[current_col]] %in% uniprot_names$accession[na_gene_primary], ]
    gene_data_separated <- gene_data_separated[!is.na(gene_data_separated[[paste0(colName, ".", n+1)]]), ]
    
    if (nrow(gene_data_separated) == 0) break
    
    uniprot_names = uniprot_names %>%
      filter(get(colName) %in% setdiff(uniprot_names[[colName]], gene_data_separated[[current_col]]))
  }
  
  # Handle merged genes
  merged_accessions <- uniprot_names %>% 
    filter(protein_name == "merged")
  
  if (nrow(merged_accessions) > 0) {
    # Retrieve and process merged accessions
    merged_info <- getMergedGeneNames(merged_accessions$accession)
    merged_uniprot_info <- getGeneUniProt(merged_info$new_Accession, information) # debug the column new_Accession
    
    colnames(merged_uniprot_info) <- uniprot_cols # changed: str_to_title(uniprot_cols)
    
    # Update UniProt names
    uniprot_names <- uniprot_names %>% 
      filter(protein_name != "merged") %>%
      rbind(merged_uniprot_info)
  }
  
  return(uniprot_names)
}


#' Process Gene Information for Multiple Accessions (wrapper)
#' 
#' Handles gene data with multiple accession numbers in a single cell
#' 
#' @param gene_data Data frame with gene data
#' @param colName Column name containing gene IDs
#' @param delimiter Delimiter separating multiple gene IDs
#' @param information UniProt information fields to retrieve
#' @param verbose Logical, whether to show detailed information
#' @return Processed gene data with UniProt information
#' @export
populateGeneNamesfromMultipleAccession <- function(gene_data, 
                                                   colName, 
                                                   delimiter,
                                                   split_char,
                                                   information, 
                                                   verbose = FALSE) {
  # Validate inputs
  if (!is.data.frame(gene_data) || 
      !colName %in% names(gene_data) || 
      nchar(delimiter) < 1 ||
      nchar(split_char)<1) {
    stop("Invalid input parameters. Check data frame, column name, and delimiter.")
  }
  
  # data cleaning
  gene_data[[colName]]<- unname(cleanAccessionColumn(gene_data, 
                                              colName, 
                                              delimiter = delimiter, 
                                              split_char = split_char)[,1])
  
  # counting number of accessions
  max_accessions <- stringr::str_count(gene_data[[colName]], ",") %>% max() + 1 # dependent of clean data
  #print(max_accessions)
  
  gene_data_separated <- gene_data %>%
    separate_wider_delim(
      cols = !!colName, 
      delim = ",", #dependent on cleandata
      names = paste0(colName, ".", 1:max_accessions),
      too_few = "align_start" #,
      # too_many = "debug"
    )
  # print(gene_data_separated[, 1:12])
  
  # Retrieve UniProt information
  uniprot_info <- getGeneInformationforMultipleAccession(
    gene_data_separated, 
    colName, 
    delimiter, 
    information, 
    verbose_ = verbose
  )
  
  # Process and merge results
  processed_data <- gene_data_separated %>%
    pivot_longer(
      cols = starts_with(paste0(colName, ".")), 
      names_to = paste0(colName, "_no"),
      names_prefix = paste0(colName, "."),
      values_to = colName,
      values_drop_na = TRUE
    ) %>%
    select(-paste0(colName, "_no"))
  
  # Join UniProt information with original data
  uniprot_info <- uniprot_info %>%
    mutate(accession = as.character(accession))
  processed_data <- processed_data %>%
    mutate(accession = as.character(accession)) # was giving it as a logical 
  
  enriched_data <- uniprot_info %>%
    left_join(processed_data, by = "accession")
  
  # Remove deleted entries
  final_data <- enriched_data %>%
    filter(!(is.na(gene_primary) & protein_name == "deleted"))
  
  return(final_data)
}


# 5. Main Wrapper Function -----------------------------------------------

#' Populate Gene Names
#' 
#' Comprehensive function to process and enrich gene data with UniProt information
#' 
#' @param gene_data Data frame with gene data
#' @param colName Column name containing gene IDs
#' @param delimiter Delimiter separating multiple gene IDs
#' @param information UniProt information fields to retrieve
#' @param verbose Logical, whether to show detailed information
#' @return List containing processed quantitative data and UniProt information
#' @export
populateGeneNames <- function(gene_data, 
                              colName, 
                              delimiter = ";",
                              split_char = "\\|",
                              information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence", 
                              verbose = FALSE) {
  # Determine processing method based on number of accessions
  accession_count <- countAccession(gene_data, colName, delimiter)
  
  if (accession_count > 1) {
    # Multiple accessions
    processed_data <- populateGeneNamesfromMultipleAccession(
      gene_data = gene_data,      
      colName = colName,
      delimiter = delimiter,
      split_char = split_char,
      information = information,  
      verbose = verbose
    )
    
    quant_col_starts <- (length(strsplit(information, ",")[[1]])+1)
    quant_col_ends <- ncol(processed_data)
    
    quant_data <- processed_data[, c(1, quant_col_starts:quant_col_ends)]
    uniprot_data <- processed_data[, 1:length(strsplit(information, ",")[[1]])]
    
  } else {
    # Single accession
    uniprot_data <- getGeneUniProt(
      gene_data[[colName]], 
      information, 
      verbose = verbose
    )
    
    colnames(uniprot_data) <- str_to_title(strsplit(information, ",")[[1]])
    rownames(uniprot_data) <- NULL
    
    quant_data <- gene_data
  }
  
  return(list(quant_data, uniprot_data))
}

# # Export key functions for package use
# export(getCleanAccession)
# export(populateGeneNames)
# export(getGeneUniProt)
# export(getMergedGeneNames)
