# Loading packages
library(data.table)
library(dplyr)
library(ggplot2)

# Command Line arguments
args <- commandArgs(trailingOnly = TRUE)
occurences <- as.character(args[1]) # here you define the name of your occurences file
input_file_taxonomy <- as.character(args[2]) # here you define the name of your input file for taxonomy (i.e WCVP)
renaming_file <- as.character(args[3]) # Here you define the name of the file from the taxonomy matcher.
output_file <- as.character(args[3]) # Here you define the name of the output file


# Load the taxonomic data from wcvp
taxonomy <- readRDS(input_file_taxonomy)

# Load the renaming file
renaming <- readRDS(renaming_file)

# Load the occurences file
occurences <- readRDS(occurences)

# Appending the name of the species which is pointed to by the accepted_plant_name_id column in the renaming file from wcvp.
renaming$accepted_plant_name_id_name <- taxonomy$species[match(renaming$accepted_plant_name_id, taxonomy$plant_name_id)]


# Examining if there are any NA's in the renaming file or the occurenced file.
occurrences_na <- sum(is.na(occurences$accepted_plant_name_id))
renaming_na <- sum(is.na(renaming$accepted_plant_name_id))
renaming_not_na <- sum(!is.na(renaming$accepted_plant_name_id))

# Printing the amount of rows with NA's in the occurrences file
if (occurrences_na > 0 ){
  cat("There are ", occurrences_na, " NA's in the occurrences file \n")
} else {
  cat("There are no NA's in the occurrences file \n")
}

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
occurences$taxon_name[match(renaming$tip, occurences$tip)] <- renaming$accepted_plant_name_id_name

# Doing some checks to see how much is recovered
if(length(unique(occurences$taxon_name)) == length(unique(renaming$accepted_plant_name_id_name))){
  cat("The number of unique names in the occurrences file is the same as the number of unique names in the renaming file \n")
}else{
  cat("The number of unique names in the occurrences file is not the same as the number of unique names in the renaming file \n")
}

# Some simple stats
cat("There are ", length(unique(occurences$accepted_plant_name_id_name)), " unique names in the occurrences file \n")
cat("There are ", length(occurences$accepted_plant_name_id_name), " occurences in the occurrences file \n")

# Can i print how many species have 1,2,3,4,5,6 .etc occurences? with the species having more than X occurences not being included in the count of X occurences?
# This would probably be easiest with a sort of for loop and the duplicated function.
# But how would I know what the max number of occurences is? I would have to find the max number of occurences and then use that as the max number in the for loop.
# In order to find the max number of occurences I will use the table function.

#Counting the species
species_counts <- table((occurences$accepted_plant_name_id_name))

# finding the species with the highest count.
most_common_species <- names(species_counts[which.max(species_counts)])

# Making a bar plot showing the number of species X the number of occurences using ggplot2
# I will use the species_counts object to make the plot.
occurence_plot <- ggplot(data = species_counts, aes(x = species_counts)) + geom_bar() + labs(x = "Number of occurences", y = "Number of species")

# Saving the plot as a pdf
ggsave("occurence_plot.pdf", plot = occurence_plot, width = 10, height = 10, units = "cm")



# Now I want to somehow link the accepted_plant_name_id and accepted_plant_id_name up to the occurrences so that my occurrences file has the accepted plant name in it.
# This needs to be done by using the renaming file as a look up table. as it links the species names in the occurrences with the accepted plant names.

# Here is a list of the columns which comes from the common_format_occurrences file
# Tip, taxon_name, family, author, split_length, genus_hybrid, species_hybrid, species, genus, taxon_rank, infra_name, comment, usable, id.

# Here is a list of the columns which comes from the wcvp file during the taxon_matcher step
# Accepted_plant_name_id, id, family.apg

