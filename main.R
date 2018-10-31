## install necessary packages (if not installed already)
load.libraries <- c("tidyverse", "data.table", "webchem", 
                    "Rcpi", "ChemmineOB", "BBmisc", "seqinr", "rdetools")
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


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Create smaller file to deal with!    [TRAINING DATA]


## read the data
dtc_data <- fread(file = "data/DTC_data.csv",
                  select = c("compound_id", "target_id", "standard_type",
                             "standard_relation", "standard_value"),
                  showProgress = TRUE,
                  na.strings = "")


fwrite(dtc_data, file = "data/train/DTC_data_selected.csv")


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Filter out the data that we don't need!    [TRAINING DATA]


## clear environment
rm(list = ls())


## read the data
dtc_data <- read_csv("data/train/DTC_data_selected.csv")


## make sure we have no NA values in these variables
dtc_data <- 
    dtc_data %>%
    filter(!is.na(compound_id) &
           !is.na(target_id) &
           !is.na(standard_value) &
           !is.na(standard_type))


## we are only interested in PKD measurements (i.e. the target variable)
dtc_data <- 
    dtc_data %>% 
    filter(standard_type == 'PKD') %>% 
    select(-standard_type)


## rename some variables
dtc_data <- 
    dtc_data %>%
    rename(pkd = standard_value)


## deal with standard_relation: '<', '>'
dtc_data <- 
    dtc_data %>%
    filter(standard_relation == '=') %>% 
    select(-standard_relation)


## three problematic records in dtc_data have these target_id values
## - "Q9JHJ5, P23979"                  --> split to 2 records (same pkd value)
## - "Q15303, P21860, P04626, P00533"  --> split to 4 records (same pkd value)
## - "P28223, P28335, P41595"          --> split to 3 records (same pkd value)
dtc_data <- 
    dtc_data %>% 
    separate_rows(target_id, sep=",\\s+")


## problematic records: compound-target pairs appearing more than once with 
## different pkd values
dtc_data <- 
    dtc_data %>% 
    group_by(compound_id,target_id) %>% 
    summarise(pkd = mean(pkd))


## unique compound and target ids
compound_IDs <- unique(dtc_data$compound_id)
target_IDs <- unique(dtc_data$target_id)


fwrite(dtc_data, file = "data/train/DTC_data_final.csv")


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## generate compound structure files!    [TRAINING DATA]


## clear environment
rm(list = ls())


## read the data
dtc <- read_csv("data/train/DTC_data_final.csv")


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)

# steps:
# A] For each compound in compound_IDs:
#   1. Rcpi::getMolFromChEMBL('<compound_id>')    # to get 'Mol' string
#   2. Write 'Mol' string to file: MOL/<compound_id>.mol
#   3. Rcpi::convMolFormat('MOL/<compound_id>.mol', 'SMILES/<compound_id>.smiles', 'mol', 'smiles')
#           Alternative: Rcpi::convMolFormat('MOL/<compound_id>.mol', 'SMILES/<compound_id>.smiles', 'mol', 'smi')
#           # convert MOL format to SMILES format
#           Alternative to all previous steps: use webchem package to generate SMILES from InChiKeys found in dtc_data
#                webchem::cs_inchi_smiles() function  -->  then write SMILES into file 'SMILES/<compound_id>.smiles'
#   3. Rcpi::extractDrugAIO(readMolFromSmi(filename, 'mol'), warn = FALSE)    # to get features for compound
#   4. Add to file, 'compoundFeatures.txt', where each row represents a unique compound
#   5. Features post-processing: remove constant-valued features & scale to range [0,1] using min-max normalization

get_compounds <- function(compound_id) {
  if ("" != compound_id){
    mol <- getMolFromChEMBL(compound_id)
    if ("" != mol) {
      mol_id <- paste("data/train/MOL/", compound_id, ".mol", sep = "")
      sink(mol_id)
      cat(mol)
      sink()
    }
  }
}
walk(compound_IDs, ~get_compounds(.x))


draw_smiles <- function(compound_id) {
  if ("" != compound_id) {
    mol_id <- paste("data/train/MOL/", compound_id, ".mol", sep = "")
    if (file.exists(mol_id)) {
      smile_id <-
        paste("data/train/SMILES/", compound_id, ".smiles", sep = "")
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
  "data/train/SMILES/CHEMBL1823872.smiles"
fileConn<-file(CHEMBL1823872_SMILES_filename)
writeLines(CHEMBL1823872_SMILES, fileConn)
close(fileConn)
# ------------------------------


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get compound features!    [TRAINING DATA]


## clear environment
rm(list = ls())


## read the data
dtc <- read_csv("data/train/DTC_data_final.csv")


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


extract_features <- function(compound_id) {
  if ("" != compound_id) {
    smile_id <-
      paste("data/train/SMILES/", compound_id, ".smiles", sep = "")
    if (file.exists(smile_id)) {
      return(extractDrugAIO(readMolFromSmi(smile_id), warn = FALSE))
    }
  }
}
mol_features <- map_df(compound_IDs, ~extract_features(.x))


## Features post-processing: remove constant-valued features 
mol_features <- mol_features[sapply(mol_features, function(x) length(unique(na.omit(x)))) > 1]
fwrite(cbind(compound_IDs, mol_features), 
       "data/train/compoundFeatures.csv", 
       quote = FALSE, 
       row.names = FALSE)


## fill NA's with column mean (if any)
replaceWithMean <- function(x) {
    x <- replace_na(x, mean(x, na.rm = TRUE))
}
mol_features_normalized <- map_df(mol_features, replaceWithMean)


## Features post-processing: scale to range [0,1] using min-max normalization
mol_features_normalized <- normalize(mol_features_normalized, method = "range")
fwrite(cbind(compound_IDs, mol_features_normalized), 
       "data/train/compoundFeaturesNormalized.csv",
       quote = FALSE, 
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP TO BE DONE BY HAND!

## --------------------------
## print target IDs on console
# for (target_ID in target_IDs) {
#     cat(target_ID, '\n')
# }
##
## Copy output and paste in this link:
## https://www.uniprot.org/uploadlists/
##
## default settings --> submit --> 'Download' with format "FASTA (canonical)"
## --------------------------

## output of this step is:  "data/train/targets.fasta"


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


# ------------------------
# B] For targets:
#   1. go to: https://www.uniprot.org/uploadlists/
#   2. Provide protein identifiers (i.e. target_IDs) to text box in link, and click 'Submit' (default settings)
#   3. Click 'Download' with format "FASTA (canonical)"
#   4. Figure out a way to convert FASTA format to RAW protein sequence (get an R package that does this?)
#   5. Put all protein sequences (one in each row) in a file called 'targetSequences.txt'
#   6. Submit 'targetSequences.txt' to http://137.132.97.65/cgi-bin/profeat2016/protein/profnew.cgi
#          (default settings)
#   7. Load 'targetFeatures.txt' and fix target IDs
#          take string after 1st '|' separator 
# ------------------------


## get target sequences file!    [TRAINING DATA]


## clear environment
rm(list = ls())


## read the data
dtc <- read_csv("data/train/DTC_data_final.csv")


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


target_sequences <- read.fasta("data/train/targets.fasta", seqtype = "AA")
sequences <- getSequence(target_sequences, as.string = TRUE)
sequences <- paste(target_IDs, as.vector(unlist(sequences)), sep = ",")


write_sequence <- function(sequence, file_name) {
  sink(file_name, append = TRUE)
  splitParts <- strsplit(sequence, ',')
  cat('>', splitParts[[1]][1], '\n')
  cat(splitParts[[1]][2], '\n\n')
  # cat(sequence[[1]])
  # cat("\n")
  sink()
}
walk(sequences, ~write_sequence(.x, "data/train/targets_modifiedIDs.fasta"))


# ------------------------------
## to check equivalence: target_IDs == IDs from fasta file?
# fastaTargetIDs <-
#     strsplit(names(target_sequences), split = '|', fixed = TRUE) %>%
#     sapply(function(x){ return(x[2]) })
# all(target_IDs == fastaTargetIDs)
# ------------------------------

# ------------------------------
## Fallen code:
# walk(sequences, ~write_sequence(.x, "data/train/targetSequences.txt"))
# target_features <- read_csv("data/train/targetFeatures.out", sep = ",")
# target_features$Feature <- substr(target_features$Feature, start = 4, stop = 9)
# fwrite(target_features, "data/train/targetFeatures.out", row.names = FALSE, quote = FALSE)
# ------------------------------


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP TO BE DONE BY HAND!

## --------------------------
## submit the file, "data/train/targets_modifiedIDs.fasta", to the link below:
## http://137.132.97.65/cgi-bin/profeat2016/protein/profnew.cgi
## 
## use default settings, output csv file and download
## --------------------------

## output of this step is:  "data/train/targetFeatures.out"


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get target features!    [TRAINING DATA]


## clear environment
rm(list = ls())


## read the data
dtc <- read_csv("data/train/DTC_data_final.csv")


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


## normalize?
target_features <- read_csv("data/train/targetFeatures.out")
target_features <- rename(target_features, target_id = Feature)
## ---------------------
## Check equivalance of target IDs?
# all(target_IDs == target_features$target_id)
## ---------------------
target_features <- target_features %>% select(-target_id)


## Features post-processing: remove constant-valued features 
target_features <- target_features[sapply(target_features, function(x) length(unique(na.omit(x)))) > 1]
fwrite(cbind(target_IDs, target_features), 
       "data/train/targetFeatures.csv", 
       quote = FALSE, 
       row.names = FALSE)


## fill NA's with column mean (if any)
replaceWithMean <- function(x) {
    x <- replace_na(x, mean(x, na.rm = TRUE))
}
target_features_normalized <- map_df(target_features, replaceWithMean)


## Features post-processing: percentage-style variables to be divided by 100
## These are the variables starting with '[G1', '[G2' and '[G4'
target_features_normalized[,1:419]    <- 
    target_features_normalized[,1:419]    / 100
target_features_normalized[,690:1188] <- 
    target_features_normalized[,690:1188] / 100


## Features post-processing: scale to range [0,1] using min-max normalization
remainingFeatures <- 1:ncol(target_features_normalized)
remainingFeatures <- setdiff(remainingFeatures, c(1:419, 690:1188))
target_features_normalized[,remainingFeatures] <- 
    normalize(target_features_normalized[,remainingFeatures], method = "range")
fwrite(cbind(target_IDs, target_features_normalized), 
       "data/train/targetFeaturesNormalized.csv",
       quote = FALSE, 
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Form proper dataset    [TRAINING DATA]


## clear environment
rm(list = ls())


# ------------------------
# C] merge all the training data into a single data frame
#   1. have the following data frames ready
#       a) compound_id, target_id, pkd
#       b) compound_id, followed by the compound features
#       c) target_id, followed by the target features
#   2. INNER JOIN: a+b
#   3. INNER JOIN: (a+b) + c
# ------------------------

dtc <- read_csv("data/train/DTC_data_final.csv")
compound_features <- read_csv("data/train/compoundFeatures.csv")
compound_features <- rename(compound_features,  compound_id = compound_IDs)
target_features <- read_csv("data/train/targetFeatures.csv")
target_features <- rename(target_features, target_id = target_IDs)
full_dataset <- 
    dtc %>% 
    inner_join(compound_features, by = 'compound_id') %>% 
    inner_join(target_features, by = 'target_id')
fwrite(full_dataset, "data/train/train.csv")


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Form proper dataset    [TEST DATA]


## clear environment
rm(list = ls())

# ------------------------
# Repeat the above for the test set (which is referred to as "Round 1 Submission Template" on
# the competition website). 
# Check this link for the test set:  https://www.synapse.org/#!Synapse:syn15667962/wiki/583675
# ------------------------


## read the test smiles
compound_smiles <- read_csv("data/round_1_template.csv")
compound_smiles <- compound_smiles$Compound_SMILES


## write it into smiles format
for (i in 1:length(compound_smiles)) {
    if (i == 1) {
        cat(compound_smiles[i], '\n', 
            file = "data/test/compoundSmiles.smiles", 
            append = FALSE)
    } else {
        cat(compound_smiles[i], '\n', 
            file = "data/test/compoundSmiles.smiles", 
            append = TRUE)
    }
}


## extract features
compound_features_test <- 
    extractDrugAIO(molecules = readMolFromSmi("data/test/compoundSmiles.smiles", 'mol'), 
                   warn = FALSE)


## keep only features that appear in the training set as well
compound_features_train <- read_csv("data/train/compoundFeatures.csv")
compound_features_test <- compound_features_test[,colnames(compound_features_train)[2:166]]


## Features post-processing: remove constant-valued features 
# compound_features_test <- compound_features_test[sapply(compound_features_test, function(x) length(unique(na.omit(x)))) > 1]
fwrite(compound_features_test, 
       "data/test/compoundFeatures.csv", 
       quote = FALSE, 
       row.names = FALSE)
## NOTE: compound features are saved WITHOUT their compound IDs


# Features post-processing: scale to range [0,1] using min-max normalization
compound_features_test_normalized <- normalize(compound_features_test, method = "range")
fwrite(compound_features_test_normalized, 
       "data/test/compoundFeaturesNormalized.csv", 
       quote = FALSE, 
       row.names = FALSE)


# get test target features
target_uniprot <- read_csv("data/round_1_template.csv")
target_uniprot <- target_uniprot$UniProt_Id

# fwrite(target_uniprot, 
#        "data/test/targetUniprotIDs.txt", 
#        row.names = FALSE, 
#        col.names = FALSE, 
#        quote = FALSE)

## --------------------------------
## print target IDs on console
# for (target_ID in target_uniprot) {
#     cat(target_ID, '\n')
# }
## --------------------------------


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP
## Uniprot get sequences in FASTA format
## Input:   target UniProt IDs from console or from 'targetUniprotIDs.txt'
## Output:  "data/test/targets.fasta"


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get sequences and merge them with their corresponding UniProt IDs
target_sequences <- read.fasta("data/test/targets.fasta", seqtype = "AA")
sequences <- getSequence(target_sequences, as.string = TRUE)
sequences <- paste(target_uniprot, as.vector(unlist(sequences)), sep = ",")


## write fasta file containing UniProt IDs and their sequences
write_sequence <- function(sequence, file_name) {
  sink(file_name, append = TRUE)
  splitParts <- strsplit(sequence, ',')
  cat('>', splitParts[[1]][1], '\n')
  cat(splitParts[[1]][2], '\n\n')
  # cat(sequence[[1]])
  # cat("\n")
  sink()
}
walk(sequences, ~write_sequence(.x, "data/test/targets_modifiedIDs.fasta"))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP
## PROFEAT: get features using FASTA file
## Input:   "data/test/targets_modifiedIDs.fasta"
## Output:  "data/test/targetFeatures.out"


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


#  Load 'testTargetFeatures.txt' and fix target IDs
#          take string after 1st '|' separator 
target_features_test <- read_csv("data/test/targetFeatures.out")
target_features_test <- rename(target_features_test, target_id = Feature)
## ---------------------
## Check equivalance of target IDs?
# all(target_uniprot == target_features_test$target_id)
## ---------------------
target_features_test <- target_features_test %>% select(-target_id)


## keep only features that appear in the training set as well
target_features_train <- read_csv("data/train/targetFeatures.csv")
target_features_test <- target_features_test[,colnames(target_features_train)[2:1432]]


## Features post-processing: remove constant-valued features 
# target_features <- target_features[sapply(target_features, function(x) length(unique(na.omit(x)))) > 1]
fwrite(target_features_test, 
       "data/test/targetFeatures.csv", 
       quote = FALSE, 
       row.names = FALSE)
## NOTE: target features are saved WITHOUT their target IDs


## fill NA's with column mean (if any)
replaceWithMean <- function(x) {
    x <- replace_na(x, mean(x, na.rm = TRUE))
}
target_features_test_normalized <- map_df(target_features_test, replaceWithMean)


## Features post-processing: percentage-style variables to be divided by 100
## These are the variables starting with '[G1', '[G2' and '[G4'
target_features_test_normalized[,1:419]    <- 
    target_features_test_normalized[,1:419]    / 100
target_features_test_normalized[,690:1188] <- 
    target_features_test_normalized[,690:1188] / 100


# Features post-processing: scale to range [0,1] using min-max normalization
remainingFeatures <- 1:ncol(target_features_test_normalized)
remainingFeatures <- setdiff(remainingFeatures, c(1:419, 690:1188))
target_features_test_normalized[,remainingFeatures] <- 
    normalize(target_features_test_normalized[,remainingFeatures], method = "range")
fwrite(target_features_test_normalized, 
       "data/test/targetFeaturesNormalized.csv",
       quote = FALSE, 
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


##
## Afterthought:
## 1. Rcpi has functions for generating descriptors for proteins (the ones starting with the prefix, 'extract')
## 2. Normalization should be performed on the training and testing together!
## 3. functions of the data.table package (e.g. 'fread')  coerce the data into 
##    data structures that require knowledge of non-standard syntax for 
##    indexing and retrieving data from within 
## 4. PROFEAT web server seems to prefer FASTA files to files containing raw
##    protein sequences
## 

