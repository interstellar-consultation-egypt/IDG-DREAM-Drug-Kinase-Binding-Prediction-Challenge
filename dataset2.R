

## path of the training and test datasets to be created
path <- 'data/dataset2/'


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## install necessary packages (if not installed already)
load.libraries <- c('tidyverse', 
                    'data.table', 
                    'webchem', 
                    'Rcpi', 
                    'ChemmineOB', 
                    'BBmisc', 
                    'seqinr', 
                    'rdetools')
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
                             'standard_inchi_key', 'standard_units'),
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
      standard_relation == '=' & 
      !is.na(standard_inchi_key)
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
    group_by(compound_id,target_id,standard_inchi_key) %>% 
    summarise(pkd = mean(pkd))

dtc_data$pkd <- -log10(dtc_data$pkd * 10 ^ -9)

fwrite(dtc_data, file = paste(path, 'train/DTC_data_final.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Get compound SMILES!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## stuff we'll need below
compound_IDs <- unique(dtc$compound_id)
target_IDs <- unique(dtc$target_id)


## get SMILES (this step takes 4+ hours!)
compound_IDs_plus_inchikeys <- 
    dtc %>% 
    select(compound_id,standard_inchi_key) %>% 
    unique() %>% 
    mutate(smiles = rep('', length(compound_id)))
for (i in 1:nrow(compound_IDs_plus_inchikeys)) {
    current_inchi_key <- compound_IDs_plus_inchikeys$standard_inchi_key[i]
    compound_IDs_plus_inchikeys$smiles[i] <-
        webchem::cs_inchikey_inchi(current_inchi_key) %>%
        webchem::cs_inchi_smiles()
}


## remove any compounds for which SMILES could not be generated
compound_IDs_plus_inchikeys <- 
    compound_IDs_plus_inchikeys %>% 
    filter(!is.na(smiles)) %>% 
    filter(smiles != '')


## remove any compounds with SMILES that are too short or too long
compound_IDs_plus_inchikeys <- 
    compound_IDs_plus_inchikeys %>% 
    filter(map_dbl(smiles, nchar) > 20) %>% 
    filter(map_dbl(smiles, nchar) < 125)


## remove compounds having identical SMILES 
## (if more than 1 compound have the same SMILES, remove them)
identical_smiles <- 
    compound_IDs_plus_inchikeys %>% 
    group_by(smiles) %>% 
    summarise(count = n()) %>% 
    filter(count > 1)
compound_IDs_plus_inchikeys <- 
    compound_IDs_plus_inchikeys %>% 
    filter(!(smiles %in% identical_smiles$smiles))


## write final set of smiles to file
allSmilesFilename <- paste(path, 'train/compoundSmiles.smi', sep = '')
write_lines(compound_IDs_plus_inchikeys$smiles, allSmilesFilename)


## for safekeeping
allSmilesWithCompFilename <- 
    paste(path, 'train/compoundSmilesWithCompIDs.csv', sep = '')
write_lines('compound_id,compound_smiles',
            path = allSmilesWithCompFilename)
write_lines(paste(compound_IDs_plus_inchikeys$compound_id,
                  compound_IDs_plus_inchikeys$smiles,
                  sep =','),
            path = allSmilesWithCompFilename,
            append = T)


## removed compounds (from above statement) to be also removed from dtc
dtc <- 
    dtc %>% 
    filter(compound_id %in% compound_IDs_plus_inchikeys$compound_id)
fwrite(dtc, file = paste(path, 'train/DTC_data_final.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP TO BE DONE BY HAND!

## Generate compound features!    [TRAINING DATA]

## Compound features are obtained via the python package: Mordred. It is 
## available from the following link:
##      https://github.com/mordred-descriptor/mordred
## 
## The command used to generate the features was:
##      python -m mordred compoundSmiles.smi -o compoundFeaturesMordred.csv
## 
## To generate the features, the mordred package needed to be installed first, 
## which was done as follows:
##      Install Rdkit: conda install -c conda-forge rdkit
##      Install MORDRED: conda install -c rdkit -c mordred-descriptor mordred
##          if doesn't work: pip install mordred
## 


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Ensure consistency!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## load data
colTypesSpecs <- paste('c', paste(rep('d',1613), collapse=''), sep='')
compound_features_train <- 
    read_csv(paste(path, 'train/compoundFeaturesMordred.csv', sep = ''), 
             col_types = colTypesSpecs)
compound_IDs_plus_smiles <- 
    read_csv(paste(path, 'train/compoundSmilesWithCompIDs.csv', sep = ''))
dtc <- 
    read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## remove two variables: 'Lipinski', 'GhoseFilter'
compound_features_train <- 
    compound_features_train %>% 
    select(-Lipinski, -GhoseFilter)


## get compound (ChEMBL) IDs that exist in the compound features file
compound_IDs_plus_smiles <- 
    compound_IDs_plus_smiles %>% 
    filter(compound_smiles %in% compound_features_train$name)


## filter with compound IDs in 'dtc'
dtc <- 
    dtc %>% 
    filter (compound_id %in% compound_IDs_plus_smiles$compound_id)


## prepare compound feature file with finalized formatting
## (there is still some post-processing to be done later)
compound_features_train <- 
    compound_features_train %>% 
    rename(compound_smiles = name)
compound_IDs_plus_smiles %>% 
    inner_join(compound_features_train, by = 'compound_smiles') %>% 
    select(-compound_smiles) %>% 
    fwrite(file = paste(path, 'train/compoundFeatures.csv', sep = ''))


## update 'DTC_data_final.csv'
fwrite(dtc, file = paste(path, 'train/DTC_data_final.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Features post-processing!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## remove constant-valued features + write to file
colTypesSpecs <- paste('c', paste(rep('d',1611), collapse=''), sep='')
compound_features_train <- 
    read_csv(paste(path, 'train/compoundFeatures.csv', sep = ''), 
             col_types = colTypesSpecs)
compound_features_train <- 
    compound_features_train[sapply(compound_features_train, 
                                   function(x) length(unique(na.omit(x))))>1]
fwrite(compound_features_train, 
       paste(path, 'train/compoundFeatures.csv', sep = ''), 
       quote = FALSE, 
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Get Target IDs!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## load data
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))


## will need this for next manual step!
write_lines(dtc$target_id %>% unique() %>% as.vector(), 
            path = paste(path, 'train/targetIDs.txt', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP TO BE DONE BY HAND!

## Get Target sequences!    [TRAINING DATA]

## Copy target IDs from 'train/targetIDs.txt' and paste them in this link:
## https://www.uniprot.org/uploadlists/
##
## default settings --> submit --> 'Download' with format 'FASTA (canonical)'

## output of this step is:  'train/targets.fasta'


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## get target sequences file!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## stuff we'll need below
target_IDs <- read_lines(paste(path, 'train/targetIDs.txt', sep = ''))


## parse sequences and their Uniprot IDs
target_sequences <- read.fasta(paste(path, 'train/targets.fasta', sep = ''), seqtype = 'AA')
sequences <- 
    getSequence(target_sequences, as.string = TRUE) %>% 
    unlist %>% 
    as.vector
getTargetId <- function(x) {
    x <- strsplit(x, '|', fixed = TRUE)[[1]][2]
    return(x)
}
seqHeaders <- 
    getName(target_sequences) %>% 
    sapply(getTargetId) %>% 
    as.vector

## consistency check
sequences <- sequences[seqHeaders %in% target_IDs]
seqHeaders <- seqHeaders[seqHeaders %in% target_IDs]


## the amino acid X can cause problems in feature generation later
Xsequences <- grepl('X', sequences)
sequences <- sequences[!Xsequences]
seqHeaders <- seqHeaders[!Xsequences]
sequences <- paste(seqHeaders, sequences, sep = ',')


## simplify protein IDs in FASTA file before submitting to PROFEAT
write_sequence <- function(sequence, file_name) {
  sink(file_name, append = TRUE)
  splitParts <- strsplit(sequence, ',')
  cat('>', splitParts[[1]][1], '\n')
  cat(splitParts[[1]][2], '\n\n')
  # cat(sequence[[1]])
  # cat('\n')
  sink()
}
walk(sequences, 
     ~write_sequence(.x, 
                     paste(path, 'train/targets_modifiedIDs.fasta', sep = '')))


## save protein Uniprot IDs along with their sequences
write_lines(sequences, 
            paste(path, 'train/targets_IDs_sequences.csv', sep = ''))


## filter with target IDs in 'DTC_data_final.csv'
dtc <- read_csv(paste(path, 'train/DTC_data_final.csv', sep = ''))
dtc <- 
    dtc %>% 
    filter (target_id %in% seqHeaders)


## update 'DTC_data_final.csv'
fwrite(dtc, file = paste(path, 'train/DTC_data_final.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Get target features!    [TRAINING DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
seq_data <- read_csv(paste(path, 'train/targets_IDs_sequences.csv', sep = ''),
                      col_names = F)
sequences <- seq_data[,2] %>% unlist %>% as.vector
target_IDs <- seq_data[,1] %>% unlist %>% as.vector


## generate features (this step takes less than an hour for 1260 proteins)
first_time <- T
for (s in sequences) {
    ## get features of current sequence
    sFeatures <- 
        extractProtAAC(s) %>% 
        append(extractProtDC(s)) %>% 
        append(extractProtGeary(s)) %>% 
        append(extractProtCTDC(s)) %>% 
        append(extractProtCTDT(s)) %>% 
        append(extractProtCTDD(s)) %>% 
        append(extractProtQSO(s)) %>% 
        append(extractProtAPAAC(s)) %>% 
        append(extractProtSOCN(s)) %>% 
        append(extractProtCTriad(s))
    
    if (first_time) {
        target_features <- sFeatures
        first_time <- F
    } else {
        target_features <- rbind(target_features, sFeatures)
    }
}


## Features post-processing: remove constant-valued features + write to file
target_features <- as_data_frame(target_features)
target_features <- 
    target_features[sapply(target_features, 
                           function(x) length(unique(na.omit(x))))>1]
fwrite(cbind(target_IDs, target_features), 
       paste(path, 'train/targetFeatures.csv', sep = ''),
       quote = FALSE,
       row.names = FALSE)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Get compound SMILES!    [TEST DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the test smiles
compound_smiles <- read_csv('data/round_1_template.csv')
compound_smiles <- compound_smiles$Compound_SMILES


## write final set of smiles to file
allSmilesFilename <- paste(path, 'test/compoundSmiles.smi', sep = '')
write_lines(compound_smiles, allSmilesFilename)


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## MANUAL STEP TO BE DONE BY HAND!

## Generate compound features!    [TRAINING DATA]

## Compound features are obtained via the python package: Mordred. It is 
## available from the following link:
##      https://github.com/mordred-descriptor/mordred
## 
## The command used to generate the features was:
##      python -m mordred compoundSmiles.smi -o compoundFeaturesMordred.csv
## 
## To generate the features, the mordred package needed to be installed first, 
## which was done as follows:
##      Install Rdkit: conda install -c conda-forge rdkit
##      Install MORDRED: conda install -c rdkit -c mordred-descriptor mordred
##          if doesn't work: pip install mordred
## 


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Compound features post-processing!    [TEST DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## load data
colTypesSpecs <- paste('c', paste(rep('d',1426), collapse=''), sep='')
compound_features_train <- 
    read_csv(paste(path, 'train/compoundFeatures.csv', sep = ''), 
             col_types = colTypesSpecs)
colTypesSpecs <- paste('c', paste(rep('d',1613), collapse=''), sep='')
compound_features_test <- 
    read_csv(paste(path, 'test/compoundFeaturesMordred.csv', sep = ''), 
             col_types = colTypesSpecs) %>% 
    rename(compound_smiles = name)


## remove two variables: 'Lipinski', 'GhoseFilter'
compound_features_test <- 
    compound_features_test %>% 
    select(-Lipinski, -GhoseFilter)


## keep only features that appear in the training set as well
numFeatures <- ncol(compound_features_train)
compound_features_test <- 
    compound_features_test %>% 
    select(-compound_smiles)
compound_features_test <- 
    compound_features_test[,colnames(compound_features_train)[2:numFeatures]]
fwrite(compound_features_test, 
       paste(path, 'test/compoundFeatures.csv', sep = ''), 
       quote = FALSE, 
       row.names = FALSE)
## NOTE: compound features are saved WITHOUT their compound IDs


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Get Targets IDs!    [TEST DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


# get test target features
target_uniprot <- read_csv('data/round_1_template.csv')
target_uniprot <- target_uniprot$UniProt_Id %>% unique
write_lines(target_uniprot, 
            path = paste(path, 'test/targetIDs.txt', sep = ''))


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


## get target sequences file!    [TEST DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## stuff we'll need below
target_IDs <- read_lines(paste(path, 'test/targetIDs.txt', sep = ''))


## parse sequences and their Uniprot IDs
target_sequences <- read.fasta(paste(path, 'test/targets.fasta', sep = ''), 
                               seqtype = 'AA')
sequences <- 
    getSequence(target_sequences, as.string = TRUE) %>% 
    unlist %>% 
    as.vector
getTargetId <- function(x) {
    x <- strsplit(x, '|', fixed = TRUE)[[1]][2]
    return(x)
}
seqHeaders <- 
    getName(target_sequences) %>% 
    sapply(getTargetId) %>% 
    as.vector

## consistency check
sequences <- sequences[seqHeaders %in% target_IDs]
seqHeaders <- seqHeaders[seqHeaders %in% target_IDs]
sequences <- paste(seqHeaders, sequences, sep = ',')


## simplify protein IDs in FASTA file before submitting to PROFEAT
write_sequence <- function(sequence, file_name) {
  sink(file_name, append = TRUE)
  splitParts <- strsplit(sequence, ',')
  cat('>', splitParts[[1]][1], '\n')
  cat(splitParts[[1]][2], '\n\n')
  # cat(sequence[[1]])
  # cat('\n')
  sink()
}
walk(sequences, 
     ~write_sequence(.x, 
                     paste(path, 'test/targets_modifiedIDs.fasta', sep = '')))


## save protein Uniprot IDs along with their sequences
write_lines(sequences, 
            paste(path, 'test/targets_IDs_sequences.csv', sep = ''))


##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------


## Get target features!    [TEST DATA]


## clear environment
rm(list = setdiff(ls(), 'path'))


## read the data
seq_data <- read_csv(paste(path, 'test/targets_IDs_sequences.csv', sep = ''),
                      col_names = F)
sequences <- seq_data[,2] %>% unlist %>% as.vector
target_IDs <- seq_data[,1] %>% unlist %>% as.vector


## generate features (this step takes less than an hour for 1260 proteins)
first_time <- T
for (s in sequences) {
    ## get features of current sequence
    sFeatures <- 
        extractProtAAC(s) %>% 
        append(extractProtDC(s)) %>% 
        append(extractProtGeary(s)) %>% 
        append(extractProtCTDC(s)) %>% 
        append(extractProtCTDT(s)) %>% 
        append(extractProtCTDD(s)) %>% 
        append(extractProtQSO(s)) %>% 
        append(extractProtAPAAC(s)) %>% 
        append(extractProtSOCN(s)) %>% 
        append(extractProtCTriad(s))
    
    if (first_time) {
        target_features <- sFeatures
        first_time <- F
    } else {
        target_features <- rbind(target_features, sFeatures)
    }
}


## write to file
target_features <- as_data_frame(target_features)
fwrite(target_features,    # cbind(target_IDs, target_features),
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
compound_features <- 
    read_csv(paste(path, 'train/compoundFeaturesNormalized.csv', sep = ''))
target_features <- 
    read_csv(paste(path, 'train/targetFeaturesNormalized.csv', sep = ''),
             na = character())
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
compound_features <- 
    read_csv(paste(path, 'test/compoundFeaturesNormalized.csv', sep = ''))
target_features <- 
    read_csv(paste(path, 'test/targetFeaturesNormalized.csv', sep = ''),
             na = character())
target_features$target_IDs <- 
    read_lines(paste(path, 'test/targetIDs.txt', sep = ''))
target_IDs <- 
    read_csv('data/round_1_template.csv') %>% 
    select(UniProt_Id) %>% 
    rename(target_IDs = UniProt_Id)
target_features <- inner_join(target_IDs,target_features, by = 'target_IDs')
full_dataset <- cbind(compound_features %>% select(-compound_id),
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

