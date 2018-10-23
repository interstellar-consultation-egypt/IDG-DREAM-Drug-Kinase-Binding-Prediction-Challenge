# install necessary packages (if not installed already)
load.libraries <- c("tidyverse", "data.table", "webchem", "Rcpi")
install.lib <- load.libraries[!load.libraries %in% installed.packages()]
for(libs in install.lib) {
  if (libs %in% c("Biobase", "Rcpi")) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Rcpi")
  } else {
    install.packages(libs, dependences = TRUE)
  }
}
  
sapply(load.libraries, require, character = TRUE)


# read the data
dtc_data <- fread(file = "data/DTC_data.csv",
                  select = c("compound_id", "target_id", "standard_type",
                             "standard_relation", "standard_value"),
                  showProgress = TRUE)


# keep only specific columns
dtc_data <- dtc_data %>%
  select(compound_id, # compound
         target_id, # target
         standard_type,  # measurement type
         standard_relation,  # '=', '>', '<', etc.
         standard_value)       # the value


# make sure we have no NA values in specific variables
dtc_data <- dtc_data %>% filter(!is.na(compound_id) & !is.na(target_id) & !is.na(standard_value) & !is.na(standard_type))


# we are only interested in PKD measurements (i.e. the target variable)
dtc_data <- dtc_data %>% filter(standard_type == 'PKD') %>% select(-standard_type)


# rename some variables
#dtc_data <- dtc_data %>% rename(inchi_key = standard_inchi_key, pkd = standard_value)
dtc_data <- dtc_data %>% rename(pkd = standard_value)


# deal with standard_relation: '<', '>'
dtc_data <- dtc_data %>% filter(standard_relation == '=') %>% select(-standard_relation)
# dtc_data <- dtc_data %>% mutate(standard_value = ifelse(standard_relation == '<', standard_value - 0.5, standard_value))
# dtc_data <- dtc_data %>% mutate(standard_value = ifelse(standard_relation == '>', standard_value + 0.5, standard_value))


# unique compound ids
compound_IDs <- dtc_data %>% select(compound_id) %>% unique()
compound_IDs <- compound_IDs$compound_id
target_IDs <- dtc_data %>% select(target_id) %>% unique()
target_IDs <- target_IDs$target_id


# dealing with records containing multiple targets in target_id
target_IDs <- c(target_IDs, 'Q9JHJ5', 'P23979', 'Q15303', 'P21860', 'P04626', 'P00533', 'P28223', 'P28335', 'P41595')
target_IDs <- target_IDs[!grepl(',',target_IDs)]    # remove multi-target records
target_IDs <- unique(target_IDs)
# TO DO:
# there are three problematic records in dtc_data that have these values
# - "Q9JHJ5, P23979"                    -->  split to two records   (pkd value is the same for these records)
# - "Q15303, P21860, P04626, P00533"    -->  split to four records  (pkd value is the same for these records)
# - "P28223, P28335, P41595"            -->  split to three records (pkd value is the same for these records)


# get features!

# steps:
# A] For each compound in compound_IDs:
#   1. Rcpi::getMolFromChEMBL('<compound_id>')    # to get 'Mol' string
#   2. Write 'Mol' string to file: MOL/<compound_id>.mol
#   3. Rcpi::convMolFormat('MOL/<compound_id>.mol', 'SMILES/<compound_id>.smiles', 'mol', 'smiles')
#           Alternative: Rcpi::convMolFormat('MOL/<compound_id>.mol', 'SMILES/<compound_id>.smiles', 'mol', 'smi')
#           # convert MOL format to SMILES format
#           Alternative to steps all previous steps: use webchem package to generate SMILES from InChiKeys found in dtc_data
#                webchem::cs_inchi_smiles() function  -->  then write SMILES into file 'SMILES/<compound_id>.smiles'
#   3. Rcpi::extractDrugAIO(readMolFromSmi(filename, 'mol'), warn = FALSE)    # to get features for compound
#   4. Add to file, 'compoundFeatures.txt', where each row represents a unique compound
#   5. Features post-processing: remove constant-valued features & scale to range [0,1] using min-max normalization, etc.
#
#
# B] For targets:
#   1. go to: https://www.uniprot.org/uploadlists/
#   2. Provide protein identifiers (i.e. target_IDs) to text box in link, and click 'Submit' (default settings)
#   3. Click 'Download' with format "FASTA (canonical)"
#   4. Figure out a way to convert FASTA format to RAW protein sequence (get an R package that does this?)
#   5. Put all protein sequences (one in each row) in a file called 'targetSequences.txt'
#   6. Submit 'targetSequences.txt' to http://137.132.97.65/cgi-bin/profeat2016/protein/profnew.cgi
#   7. You will probably need my help at this point  -->  call me!
#
#
# Repeat the above for the test set (which is referred to as "Round 1 Submission Template" on the competition website). 
# Check this link for the test set:  https://www.synapse.org/#!Synapse:syn15667962/wiki/583675


# You may consider exploring the functions that Rcpi has for generating descriptors (the ones starting with the prefix, 'extract')


