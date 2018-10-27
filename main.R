# install necessary packages (if not installed already)
load.libraries <- c("tidyverse", "data.table", "webchem", 
                    "Rcpi", "ChemmineOB", "BBmisc", "seqinr")
install.lib <- load.libraries[!load.libraries %in% installed.packages()]
for(libs in install.lib) {
  if (libs %in% c("Biobase", "Rcpi")) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Rcpi", dependencies = c("Imports", "Enhances"))
  } else {
    install.packages(libs, dependences = TRUE)
  }
}
  
sapply(load.libraries, require, character = TRUE)


# read the data
dtc_data <- fread(file = "data/DTC_data.csv",
                  select = c("compound_id", "target_id", "standard_type",
                             "standard_relation", "standard_value"),
                  showProgress = TRUE,
                  na.strings = "")
fwrite(dtc_data, file = "data/DTC_data_selected.csv")


# set the standard_value as numeric (which will introduce NA values due to empty strings)
dtc_data$standard_value <- as.numeric(dtc_data$standard_value)


# make sure we have no NA values in specific variables
dtc_data <- dtc_data %>%
  filter(
    !is.na(compound_id) &
      !is.na(target_id) &
      !is.na(standard_value) &
      !is.na(standard_type)
  )


# we are only interested in PKD measurements (i.e. the target variable)
dtc_data <- dtc_data %>% 
  filter(standard_type == 'PKD') %>% 
  select(-standard_type)


# rename some variables
#dtc_data <- dtc_data %>% rename(inchi_key = standard_inchi_key, pkd = standard_value)
dtc_data <- dtc_data %>%
  rename(pkd = standard_value)


# deal with standard_relation: '<', '>'
dtc_data <- dtc_data %>%
  filter(standard_relation == '=') %>% 
  select(-standard_relation)


# there are three problematic records in dtc_data that have these values
# - "Q9JHJ5, P23979"                    -->  split to two records   (pkd value is the same for these records)
# - "Q15303, P21860, P04626, P00533"    -->  split to four records  (pkd value is the same for these records)
# - "P28223, P28335, P41595"            -->  split to three records (pkd value is the same for these records)
dtc_data <- dtc_data %>% 
  separate_rows(target_id, sep=",\\s+")


# unique compound ids
compound_IDs <- dtc_data %>% 
  select(compound_id) %>% 
  unique()
compound_IDs <- compound_IDs$compound_id
target_IDs <- dtc_data %>%
  select(target_id) %>% 
  unique()
target_IDs <- target_IDs$target_id


fwrite(dtc_data, file = "data/DTC_data_final.csv")


# get features!

# steps:
# A] For each compound in compound_IDs:
#   1. Rcpi::getMolFromChEMBL('<compound_id>')    # to get 'Mol' string
#   2. Write 'Mol' string to file: MOL/<compound_id>.mol
get_compounds <- function(compound_id) {
  if ("" != compound_id){
    mol <- getMolFromChEMBL(compound_id)
    if ("" != mol) {
      mol_id <- paste("data/MOL/", compound_id, ".mol", sep = "")
      sink(mol_id)
      cat(mol)
      sink()
    }
  }
  
}

walk(compound_IDs, ~get_compounds(.x))


#   3. Rcpi::convMolFormat('MOL/<compound_id>.mol', 'SMILES/<compound_id>.smiles', 'mol', 'smiles')
#           Alternative: Rcpi::convMolFormat('MOL/<compound_id>.mol', 'SMILES/<compound_id>.smiles', 'mol', 'smi')
#           # convert MOL format to SMILES format
#           Alternative to all previous steps: use webchem package to generate SMILES from InChiKeys found in dtc_data
#                webchem::cs_inchi_smiles() function  -->  then write SMILES into file 'SMILES/<compound_id>.smiles'
draw_smiles <- function(compound_id) {
  if ("" != compound_id) {
    mol_id <- paste("data/MOL/", compound_id, ".mol", sep = "")
    if (file.exists(mol_id)) {
      smile_id <-
        paste("data/SMILES/", compound_id, ".smiles", sep = "")
      convMolFormat(mol_id, smile_id, "mol", "smiles")
    }
  }
}

walk(compound_IDs, ~draw_smiles(.x))


# ------------------------------
# Problem: compound ID 'CHEMBL1823872' has no downloadable MOL file
# Solution: InchiKey from dtc_data, convert to InChi, then convert to SMILES
CHEMBL1823872_SMILES <- 
  webchem::cs_inchikey_inchi('NHXLMOGPVYXJNR-ATOGVRKGSA-N') %>% 
  webchem::cs_inchi_smiles()
CHEMBL1823872_SMILES_filename <- 
  "data/SMILES/CHEMBL1823872.smiles"
fileConn<-file(CHEMBL1823872_SMILES_filename)
writeLines(CHEMBL1823872_SMILES, fileConn)
close(fileConn)
# ------------------------------


#   3. Rcpi::extractDrugAIO(readMolFromSmi(filename, 'mol'), warn = FALSE)    # to get features for compound
#   4. Add to file, 'compoundFeatures.txt', where each row represents a unique compound
#   5. Features post-processing: remove constant-valued features & scale to range [0,1] using min-max normalization, etc.

extract_features <- function(compound_id) {
  if ("" != compound_id) {
    smile_id <-
      paste("data/SMILES/", compound_id, ".smiles", sep = "")
    if (file.exists(smile_id)) {
      return(extractDrugAIO(readMolFromSmi(smile_id), warn = FALSE))
    }
  }
}

mol_features <- map_df(compound_IDs, ~extract_features(.x))


# Features post-processing: remove constant-valued features 
mol_features <- mol_features[sapply(mol_features, function(x) length(unique(na.omit(x)))) > 1]
fwrite(mol_features, "data/compoundFeatures.csv", quote = FALSE, row.names = FALSE)


# Features post-processing:  scale to range [0,1] using min-max normalization, etc.
mol_features_normalized <- normalize(mol_features, method = "range")
fwrite(mol_features_normalized, "data/compoundFeaturesNormalized.csv",
       quote = FALSE, row.names = FALSE)


# ------------------------
# B] For targets:
#   1. go to: https://www.uniprot.org/uploadlists/
#   2. Provide protein identifiers (i.e. target_IDs) to text box in link, and click 'Submit' (default settings)
#   3. Click 'Download' with format "FASTA (canonical)"
#   4. Figure out a way to convert FASTA format to RAW protein sequence (get an R package that does this?)
target_sequences <- read.fasta("data/targets.fasta", seqtype = "AA")
sequences <- getSequence(target_sequences, as.string = TRUE)
#   5. Put all protein sequences (one in each row) in a file called 'targetSequences.txt'
write_sequence <- function(sequence) {
  sink("data/targetSequences.txt", append = TRUE)
  cat(sequence[[1]])
  cat("\n")
  sink()
}
walk(sequences, ~write_sequence(.x))
#   6. Submit 'targetSequences.txt' to http://137.132.97.65/cgi-bin/profeat2016/protein/profnew.cgi
#          (default settings)
#   7. Load 'targetFeatures.txt' and fix target IDs
#          take string after 1st '|' separator 
target_features <- fread("data/targetFeatures.out", sep = ",")

# ------------------------
# C] merge all the training data into a single data frame
#   1. have the following data frames ready
#       a) compound_id, target_id, pkd
#       b) compound_id, followed by the compound features
#       c) target_id, followed by the target features
#   2. INNER JOIN: a+b
#   3. INNER JOIN: (a+b) + c


# ------------------------
# Repeat the above for the test set (which is referred to as "Round 1 Submission Template" on the competition website). 
# Check this link for the test set:  https://www.synapse.org/#!Synapse:syn15667962/wiki/583675


# You may consider exploring the functions that Rcpi has for generating descriptors (the ones starting with the prefix, 'extract')

