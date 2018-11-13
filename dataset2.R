

## path of the training and test datasets to be created
path <- 'data/dataset2/'


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## install necessary packages (if not installed already)
load.libraries <- c('tidyverse', 'data.table', 'webchem', 
                    'Rcpi', 'ChemmineOB', 'BBmisc', 'seqinr', 'rdetools')
install.lib <- load.libraries[!load.libraries %in% installed.packages()]
for(libs in install.lib) {
  if (libs %in% c('Biobase', 'Rcpi')) {
    source('https://bioconductor.org/biocLite.R')
    biocLite('Rcpi', dependencies = c('Imports', 'Enhances'))
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
dtc_data <- fread(file = 'data/DTC_data.csv',
                  select = c('compound_id', 'target_id', 'standard_type',
                             'standard_relation', 'standard_value',
                             'standard_units'),
                  showProgress = TRUE,
                  na.strings = '')


fwrite(dtc_data, file = paste(path, 'train/DTC_data_selected.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Filter out the data that we don't need!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
dtc_data <- read_csv(paste(path, 'train/DTC_data_selected.csv', sep = ''))


## make sure we have no NA values in these variables
dtc_data <-
  dtc_data %>%
  filter(
    !is.na(compound_id) &
      !is.na(target_id) &
      !is.na(standard_value) &
      !is.na(standard_type) &
      toupper(standard_type) == 'KD' &
      toupper(standard_units) == 'NM' &
      standard_relation == '='
  ) %>%
  select(-standard_type, -standard_relation, -standard_units) %>% 
  rename(pkd = standard_value)


## three problematic records in dtc_data have these target_id values
## - 'Q9JHJ5, P23979'                  --> split to 2 records (same pkd value)
## - 'Q15303, P21860, P04626, P00533'  --> split to 4 records (same pkd value)
## - 'P28223, P28335, P41595'          --> split to 3 records (same pkd value)
dtc_data <- 
    dtc_data %>% 
    separate_rows(target_id, sep=',\\s+')


## problematic records: compound-target pairs appearing more than once with 
## different pkd values
dtc_data <- 
    dtc_data %>% 
    group_by(compound_id,target_id) %>% 
    summarise(pkd = mean(pkd))

dtc_data$pkd <- -log10(dtc_data$pkd * 10 ^ -9)

fwrite(dtc_data, file = paste(path, 'train/DTC_data_final.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


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


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## generate compound structure files!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


## get compound MOL files
get_compounds <- function(compound_id) {
  mol_id <- paste(path, 'train/MOL/', compound_id, '.mol', sep = '')
  if ('' != compound_id & !file.exists(mol_id)) {
    mol <- getMolFromChEMBL(compound_id)
    if ('' != mol) {
      sink(mol_id)
      cat(mol)
      sink()
    }
  }
}
walk(compound_IDs, ~ get_compounds(.x))


## get compound SMILES files
draw_smiles <- function(compound_id) {
  if ('' != compound_id) {
    mol_id <- paste(path, 'train/MOL/', compound_id, '.mol', sep = '')
    if (file.exists(mol_id)) {
      smile_id <-
        paste(path, 'train/SMILES/', compound_id, '.smiles', sep = '')
      convMolFormat(mol_id, smile_id, 'mol', 'smiles')
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
  paste(path, 'train/SMILES/CHEMBL1823872.smiles', sep = '')
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
rm(list = setdiff(ls(), 'path'))


## read the data
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


extract_features <- function(compound_id) {
  if ('' != compound_id) {
    smile_id <-
      paste(path, 'train/SMILES/', compound_id, '.smiles', sep = '')
    if (file.exists(smile_id)) {
      return(extractDrugAIO(readMolFromSmi(smile_id), warn = FALSE))
    }
  }
}
mol_features <- map_df(compound_IDs, ~extract_features(.x))


## Features post-processing: remove constant-valued features + write to file
mol_features <- 
    mol_features[sapply(mol_features, function(x) length(unique(na.omit(x))))>1]
fwrite(cbind(compound_IDs, mol_features), 
       paste(path, 'train/compoundFeatures.csv', sep = ''), 
       quote = FALSE, 
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------

# B] For targets:
#   1. go to: https://www.uniprot.org/uploadlists/
#   2. Provide protein identifiers (i.e. target_IDs) to text box in link, and click 'Submit' (default settings)
#   3. Click 'Download' with format 'FASTA (canonical)'
#   4. Figure out a way to convert FASTA format to RAW protein sequence (get an R package that does this?)
#   5. Put all protein sequences (one in each row) in a file called 'targetSequences.txt'
#   6. Submit 'targetSequences.txt' to http://137.132.97.65/cgi-bin/profeat2016/protein/profnew.cgi
#          (default settings)
#   7. Load 'targetFeatures.txt' and fix target IDs
#          take string after 1st '|' separator 


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
## default settings --> submit --> 'Download' with format 'FASTA (canonical)'
## --------------------------

## output of this step is:  paste(path, 'train/targets.fasta', sep = '')


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get target sequences file!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


## simplify protein IDs in FASTA file before submitting to PROFEAT
target_sequences <- read.fasta(paste(path, 'train/targets.fasta', sep = ''), seqtype = 'AA')
sequences <- getSequence(target_sequences, as.string = TRUE)
sequences <- paste(target_IDs, as.vector(unlist(sequences)), sep = ',')
write_sequence <- function(sequence, file_name) {
  sink(file_name, append = TRUE)
  splitParts <- strsplit(sequence, ',')
  cat('>', splitParts[[1]][1], '\n')
  cat(splitParts[[1]][2], '\n\n')
  # cat(sequence[[1]])
  # cat('\n')
  sink()
}
walk(sequences, ~write_sequence(.x, paste(path, 'train/targets_modifiedIDs.fasta', sep = '')))


# ------------------------------
## to check equivalence: target_IDs == IDs from fasta file?
# fastaTargetIDs <-
#     strsplit(names(target_sequences), split = '|', fixed = TRUE) %>%
#     sapply(function(x){ return(x[2]) })
# all(target_IDs == fastaTargetIDs)
# ------------------------------


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP TO BE DONE BY HAND!

## --------------------------
## submit the file, paste(path, 'train/targets_modifiedIDs.fasta', sep = ''), 
## to the link below:
## http://137.132.97.65/cgi-bin/profeat2016/protein/profnew.cgi
## 
## use default settings, output csv file and download
## --------------------------

## output of this step is:  paste(path, 'train/targetFeatures.out', sep = '')


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get target features!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


## Features post-processing: remove constant-valued features + write to file
target_features <- read_csv(paste(path, 'train/targetFeatures.out', sep = ''))
target_features <- rename(target_features, target_id = Feature)
## ---------------------
## Check equivalance of target IDs?
# all(target_IDs == target_features$target_id)
## ---------------------
target_features <- target_features %>% select(-target_id)
target_features <- 
    target_features[sapply(target_features, function(x) length(unique(na.omit(x))))>1]
fwrite(cbind(target_IDs, target_features), 
       paste(path, 'train/targetFeatures.csv', sep = ''),
       quote = FALSE,
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## repeat above steps for test data    [TEST DATA]

## ------------------------
## Repeat the above for the test set (which is referred to as 'Round 1 
## Submission Template' on the competition website). 
## Check this link for the test set:  
##   https://www.synapse.org/#!Synapse:syn15667962/wiki/583675
## ------------------------


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the test smiles
compound_smiles <- read_csv('data/round_1_template.csv')
compound_smiles <- compound_smiles$Compound_SMILES


## write it into smiles format
for (i in 1:length(compound_smiles)) {
    if (i == 1) {
        cat(compound_smiles[i], '\n', 
            file = paste(path, 'test/compoundSmiles.smiles', sep = ''), 
            append = FALSE)
    } else {
        cat(compound_smiles[i], '\n', 
            file = paste(path, 'test/compoundSmiles.smiles', sep = ''), 
            append = TRUE)
    }
}


## extract features
compound_features_test <- 
    extractDrugAIO(molecules = readMolFromSmi(paste(path, 'test/compoundSmiles.smiles', sep = ''), 'mol'), 
                   warn = FALSE)


## keep only features that appear in the training set as well
compound_features_train <- read_csv(paste(path, 'train/compoundFeatures.csv', sep = ''))
compound_features_test <- compound_features_test[,colnames(compound_features_train)[2:166]]
fwrite(compound_features_test, 
       paste(path, 'test/compoundFeatures.csv', sep = ''), 
       quote = FALSE, 
       row.names = FALSE)
## NOTE: compound features are saved WITHOUT their compound IDs


# get test target features
target_uniprot <- read_csv('data/round_1_template.csv')
target_uniprot <- target_uniprot$UniProt_Id

# fwrite(target_uniprot, 
#        paste(path, 'test/targetUniprotIDs.txt', sep = ''), 
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
## Output:  paste(path, 'test/targets.fasta', sep = '')


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get sequences and merge them with their corresponding UniProt IDs
target_sequences <- read.fasta(paste(path, 'test/targets.fasta', sep = ''), seqtype = 'AA')
sequences <- getSequence(target_sequences, as.string = TRUE)
sequences <- paste(target_uniprot, as.vector(unlist(sequences)), sep = ',')


## write fasta file containing UniProt IDs and their sequences
write_sequence <- function(sequence, file_name) {
  sink(file_name, append = TRUE)
  splitParts <- strsplit(sequence, ',')
  cat('>', splitParts[[1]][1], '\n')
  cat(splitParts[[1]][2], '\n\n')
  # cat(sequence[[1]])
  # cat('\n')
  sink()
}
walk(sequences, ~write_sequence(.x, paste(path, 'test/targets_modifiedIDs.fasta', sep = '')))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP
## PROFEAT: get features using FASTA file
## Input:   paste(path, 'test/targets_modifiedIDs.fasta', sep = '')
## Output:  paste(path, 'test/targetFeatures.out', sep = '')


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


#  Load 'testTargetFeatures.txt' and fix target IDs
#          take string after 1st '|' separator 
target_features_test <- read_csv(paste(path, 'test/targetFeatures.out', sep = ''))
target_features_test <- rename(target_features_test, target_id = Feature)
## ---------------------
## Check equivalance of target IDs?
# all(target_uniprot == target_features_test$target_id)
## ---------------------
target_features_test <- target_features_test %>% select(-target_id)


## keep only features that appear in the training set as well
target_features_train <- read_csv(paste(path, 'train/targetFeatures.csv', sep = ''))
target_features_test <- target_features_test[,colnames(target_features_train)[2:1432]]
fwrite(target_features_test, 
       paste(path, 'test/targetFeatures.csv', sep = ''), 
       quote = FALSE, 
       row.names = FALSE)
## NOTE: target features are saved WITHOUT their target IDs


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## fill NAs in compound/target features + normalization
source('normalizer.R')


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Form proper dataset    [TRAINING DATA]

# ------------------------
# C] merge all the training data into a single data frame
#   1. have the following data frames ready
#       a) compound_id, target_id, pkd
#       b) compound_id, followed by the compound features
#       c) target_id, followed by the target features
#   2. INNER JOIN: a+b
#   3. INNER JOIN: (a+b) + c
# ------------------------


## clear environment
rm(list = setdiff(ls(), 'path'))


## generate training set for training prediction models
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))
compound_features <- read_csv(paste(path, 'train/compoundFeaturesNormalized.csv', sep = ''))
compound_features <- rename(compound_features,  compound_id = compound_IDs)
target_features <- read_csv(paste(path, 'train/targetFeaturesNormalized.csv', sep = ''))
target_features <- rename(target_features, target_id = target_IDs)
full_dataset <- 
    dtc %>% 
    inner_join(compound_features, by = 'compound_id') %>% 
    inner_join(target_features, by = 'target_id')
fwrite(full_dataset, paste(path, 'train/train.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Form proper dataset    [TEST DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## generate test set
compound_features <- read_csv(paste(path, 'test/compoundFeaturesNormalized.csv', sep = ''))
target_features <- read_csv(paste(path, 'test/targetFeaturesNormalized.csv', sep = ''))
full_dataset <- cbind(compound_features %>% select(-compound_IDs),
                      target_features %>% select(-target_IDs))
fwrite(full_dataset, paste(path, 'test/test.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


##
## Afterthought:
## 1. Rcpi has functions for generating descriptors for proteins (the ones 
##    starting with the prefix, 'extract'). Maybe worth checking out?
## 2. Normalization should be performed on the training and test sets together!
## 3. NA imputation should be performed on the training and test sets together!
## 4. functions of the data.table package (e.g. 'fread')  coerce the data into 
##    data structures that require knowledge of non-standard syntax for 
##    indexing and retrieving data from within. Better to avoid if possible.
## 5. PROFEAT web server seems to prefer FASTA files to files containing raw
##    protein sequences
## 

