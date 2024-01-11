# Loading packages
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table","dplyr","ggplot2","scales")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

##############################################################
# Set variables for local testing
# setwd("/home/au543206/GenomeDK/Trf_models/workflow/01_distribution_data/05_Taxon_match/") # local
# setwd("/home/owrisberg/Trf_models/workflow/01_distribution_data/05_Taxon_match") # srun
# occurences_input <- "../04_common_format/gbif_common_format.rds"
# input_file_taxonomy <- "../04_common_format/wcvp_names_apg_aligned.rds" 
# renaming_file <- "gbif_taxon_matched.rds"
# output_file <- "gbif_renamed.rds"
# raw_occurences <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/02_data_parsing/gbif_parsed.rds") # srun 
##############################################################

# # Command Line arguments for script
args <- commandArgs(trailingOnly = TRUE)
occurences_input <- as.character(args[1]) # here you define the name of your occurences file
input_file_taxonomy <- as.character(args[2]) # here you define the name of your input file for taxonomy (i.e WCVP)
renaming_file <- as.character(args[3]) # Here you define the name of the file from the taxonomy matcher.
output_file <- as.character(args[4]) # Here you define the name of the output file
raw_occurences <- as.character(args[5]) # Here you define the name of the raw occurences file

###############################################################
# Load the taxonomic data from wcvp
taxonomy <- readRDS(input_file_taxonomy)

# Load the renaming file
renaming <- readRDS(renaming_file)

# Load the occurences file
occurences <- readRDS(occurences_input)

# Appending the name of the species which is pointed to by the accepted_plant_name_id column in the renaming file from wcvp.
renaming$accepted_plant_name_id_name <- taxonomy$taxon_name[match(renaming$accepted_plant_name_id, taxonomy$plant_name_id)]

# Examining if there are any NA's in the renaming file or the occurenced file.
#occurrences_na <- sum(is.na(occurences$accepted_plant_name_id))
renaming_na <- sum(is.na(renaming$accepted_plant_name_id))
renaming_not_na <- sum(!is.na(renaming$accepted_plant_name_id))

# Printing the amount of rows with NA's in the occurrences file
# if (occurrences_na > 0 ){
#   cat("There are ", occurrences_na, " NA's in the occurrences file \n")
# } else {
#   cat("There are no NA's in the occurrences file \n")
# }

# Printing the amount of rows with NA's in the renaming file
if (renaming_na > 0) {
  cat("There are ", renaming_na, " NA's in the renaming file \n")
} else {
  cat("There are no NA's in the renaming file \n")
}

# Printing the amount of rows with no NA's in the renaming file
if (renaming_not_na > 0) {
  cat("There are ", renaming_not_na, "Not NA's in the renaming file \n")
} else {
  cat("There are only NA's in the renaming file \n")
}

# # Printing the rows with NA's in the renaming file
# print(renaming[is.na(renaming$accepted_plant_name_id),])

# Removin all occurrences where the accepted_plant_name_id is NA
renaming <- renaming[!is.na(renaming$accepted_plant_name_id),]
renaming_na <- sum(is.na(renaming$accepted_plant_name_id))

# Printing the amount of rows with NA's in the renaming file
if (renaming_na > 0) {
  cat("There are ", renaming_na, " NA's in the renaming file \n")
} else {
  cat("There are no NA's in the renaming file \n")
}

# Using the match function to rename the names in the occurrences file
# Can you make the below file print out the where there are NA's in the subsscribted assignment?
# for (i in unique(renaming$accepted_plant_name_id)){
#   #cat("Were currently looking for a match for" ,renaming$accepted_plant_name_id_name[i])
#   if(i %in% occurences$accepted_plant_name_id){
# 	cat( i," is found in occurrences \n")
#   }else{
# 	cat( i," is not found in occurrences \n")
#   }
  #occurences$taxon_name[match(renaming$accepted_plant_name_id[i], occurences$accepted_plant_name_id)] <- renaming$accepted_plant_name_id_name[i]
#}

# This is the code which crashes with the following error message:
#Error in occurences$taxon_name[match(renaming$accepted_plant_name_id,  : 
#  NAs are not allowed in subscripted assignments
# Is it because these two accepted_plant_name_id is not the same ID?

# Here I am renaming the names in the occurrences file using the renaming file.
# This is probably most easily done using the match function, but I somehow have to take into account that there are some names in the occurrences file which are not in the renaming file.
print("Head of occurrences")
head(occurences)

# Removing names from occurences which are not in renaming
length(which(occurences$occurrence_species_name %in% renaming$occurrence_species_name)) # names in renaming

length(which(!occurences$occurrence_species_name %in% renaming$occurrence_species_name)) # names not in renaming

occurences_no_na <- occurences[occurences$occurrence_species_name %in% renaming$occurrence_species_name,]

# Now I should be able to assign new names using match
occurences_no_na$new_name <- renaming$accepted_plant_name_id_name[match(occurences_no_na$occurrence_species_name, renaming$occurrence_species_name)]

dim(occurences_no_na[which(occurences_no_na$new_name != occurences_no_na$occurrence_species_name),]) # We correctly changed 19604 names.


# Doing some checks to see how much is recovered
if(length(unique(occurences_no_na$new_name)) == length(unique(renaming$accepted_plant_name_id_name))){
  cat("The number of unique names in the occurrences file is the same as the number of unique names in the renaming file \n")
}else{
  cat("The number of unique names in the occurrences file is not the same as the number of unique names in the renaming file \n")
}

# Some simple stats
cat("There are ", length(unique(occurences$taxon_name)), " unique names in the occurrences file \n")


# I now have a file with all the unique plant name ID's from the occurrences downloaded from GBIF and the corresponding correct name from WCVP.
# Now I need to use this file to rename all species names in the file with the raw occurrence data. 


################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

cat("\n\n")

cat("Appending the names from the renaming file to the raw_occurences file \n")
cat("Unfurtunately there are some names in the raw_occurences file which are not found in the renaming file \n")
cat("The reason for this is that these names are not accepted in the GBIF taxonomy and the taxonomy look up on GBIF did not return any valid names for me to match to the WCVP taxonomy \n")
cat("These names can be found in the file bad_names.rds \n")

## Removing occurrences from raw_occurences that do not have a species name found in the occurrences_no_na file
raw_occurences_no_na <- raw_occurences[raw_occurences$species %in% occurences_no_na$occurrence_species_name,]

# Adding the wcvp_taxon_name to the raw_occurences_no_na file
raw_occurences_no_na$wcvp_taxon_name <- occurences_no_na$new_name[match(raw_occurences_no_na$species, occurences_no_na$occurrence_species_name)]

# Checking if the number of unique names in the raw_occurences_no_na file is the same as the number of unique names in the renaming file
cat("The number of unique names in the raw_occurences_no_na file is ", length(unique(raw_occurences_no_na$wcvp_taxon_name)), "\n")
cat("The number of unique names in the renaming file is ", length(unique(renaming$accepted_plant_name_id_name)), "\n")

#Saving the raw_occurences_no_na file
saveRDS(raw_occurences_no_na, "gbif_renamed.rds")


################################################################################################################################
################################################################################################################################
################################################################################################################################
# Load the raw occurrences file
# raw_occurences <- readRDS("/home/owrisberg/Trf_models/workflow/01_distribution_data/02_data_parsing/gbif_parsed.rds") # srun 
# head(raw_occurences)
# dim(raw_occurences)
# length(unique(raw_occurences$species)) #367815 unique species names in the raw_occurences file

# # Find which species names are not found in the occurrences_no_na file
# length(which(!unique(raw_occurences$species) %in% occurences_no_na$occurrence_species_name))

# # Can I somehow print some of these names?
# bad_dtf <- raw_occurences[which(!raw_occurences$species %in% occurences_no_na$occurrence_species_name),]
# bad_names_unique <- unique(raw_occurences$species[which(!raw_occurences$species %in% occurences_no_na$occurrence_species_name)]) #This should be the names which are not found in the occurrences_no_na file
# bad_names_dups <- raw_occurences$species[which(!raw_occurences$species %in% occurences_no_na$occurrence_species_name)]
# counts_names <- table(bad_names_dups)
# counts_names <- as.data.frame(counts_names)

# # The total number of occurrences of which have a species name which is not found in the occurrences_no_na file
# sum(counts_names$Freq) #16.553.921 occurrences of species names which are not found in the occurrences_no_na file
# # 7.755.451 of which have no name at all
# # 8.798.470 of which have a name which is not found in the occurrences_no_na file
# # The number of species not found in the occurrences_no_na file
# length(counts_names$bad_names) # 45.539 unqiue species names which are not found in the occurrences_no_na file 


# # Can I order the data frame by the number of occurrences?
# counts_names <- counts_names[order(-counts_names$Freq),]
# counts_names[1:100,] 

# # So the test is now: If i remove the bad names from the raw_occurences file, then I should be able to rename all the species names in the raw_occurences file using the match function.
# raw_occurences_no_bad <- subset(raw_occurences, !species %in% bad_names_unique)
# dim(raw_occurences_no_bad)

# # appending the wcvp_taxon_name to the raw_occurences_no_bad file
# raw_occurences_no_bad$wcvp_taxon_name <- occurences_no_na$new_name[match(raw_occurences_no_bad$species, occurences_no_na$occurrence_species_name)]


