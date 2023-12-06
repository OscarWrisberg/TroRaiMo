# builds a common format based on the Smith & Brown ALLMB tree tip labels from GBIF to feed into the taxonomy matcher # nolint: line_length_linter.

#This R script builds a common format based on the Smith & Brown ALLMB tree tip labels from GBIF to feed into the taxonomy matcher.
#It takes in three command line arguments: the name of the input file containing occurrence data, the name of the input file containing taxonomic
# lookup data from GBIF, and the name of the output file to be created. The script loads the taxonomic lookup data, selects the required columns,
# and removes lines that do not include all the required columns. It then creates a dataframe from the selected columns.
# The script also loads the occurrence data and tests if all the taxon IDs from the taxonomy lookup are found in the occurrences.
# Finally, it joins the occurrences with the taxonomy data gathered and outputs the resulting dataframe to a file with the specified name.

#######################################################################################
#                              Setting up testing                                     #
#######################################################################################

# # setting the wd
# setwd("/home/au543206/GenomeDK/Biome_estimation/workflow/05_create_common_format") # Mounted drive
# setwd("/home/owrisberg/Biome_estimation/workflow/04_gbif_lookup/") # Drive when using SRUN


# input_file_occurrences <- "/home/owrisberg/Biome_estimation/workflow/03_datasplit/gbif_parsed.rds" # here you define the name of your input file
# input_file_taxonomy <-  "gbif_parsed_taxon_data.rds"# here you define the name of your input file
# output_file <-  "gbif_common_format.rds"# Here you define the name of the output file


###############################################################################
#                                 Set up                                      #
###############################################################################

rm(list = setdiff(ls(), lsf.str()))
library(phytools)
library(data.table)
library(dplyr)

# Command Line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file_occurrences <- as.character(args[1]) # here you define the name of your input file
input_file_taxonomy <- as.character(args[2]) # here you define the name of your input file
output_file <- as.character(args[3]) # Here you define the name of the output file


cat("loading taxonomy file \n")

# Load the taxonomic look up data from GBIF
taxonomy <- readRDS(input_file_taxonomy)

##################################################################################
#                                    GBIF                                        #
##################################################################################

# transform list to dataframe
## extract data from list
taxonomy_dat <- sapply(taxonomy, "[[", "data")

## removing lines that do not include all required columns
cat("Removing lines which do not include all the required columns \n")
cols <- c("key", # gbifID
          "scientificName", # scientificName
          "authorship", # authorship
          "taxonomicStatus", #taxonomicStatus
          "rank", # taxonRank ?
          "family", # family
          "genus", # genus
          "species") # species

# This is a function which checks if all the columns are in the data frame
cat("Calling the function \n")
func1 <- function(x) {
  all(cols %in% names(x))
}
remo <- which(unlist(lapply(taxonomy_dat, func1)) == FALSE)
non_null_res <- taxonomy_dat[-remo]
rm(taxonomy_dat)


# Selecting columns
cat("Creating a dataframe from select columns of the taxonomy data \n")
taxonid_tax <- sapply(non_null_res, "[[", "key")
scientificname_tax <- sapply(non_null_res, "[[", "scientificName")
authorship_tax <- sapply(non_null_res, "[[", "authorship")
taxonomicstatus_tax <- sapply(non_null_res, "[[", "taxonomicStatus")
rank_tax <- sapply(non_null_res, "[[", "rank")
family_tax <- sapply(non_null_res, "[[", "family")
genus_tax <- sapply(non_null_res, "[[", "genus")
species_tax <- sapply(non_null_res, "[[", "species")


# Adding columns together to create a dataframe
taxonomy_df <- data.frame(taxonid_tax,
  scientificname_tax,
  authorship_tax,
  taxonomicstatus_tax,
  rank_tax,
  family_tax,
  genus_tax,
  species_tax
)

rm(non_null_res)


##################################################################################
#                                 Occurrences                                    #
##################################################################################

cat("Loading the occurrences data \n")
occurences <- readRDS(input_file_occurrences)

# This tests if all the taxon id's from the taxonomy look up is found in the occurrences.
cat(paste0("Are all the taxonID's in taxonomy_df found in the acceptedNameUsageID in occurrence data: ", all(taxonomy_df$taxonid_tax %in% occurences$acceptedNameUsageID)), "\n")

# This test if all the acceptedNameUsageID's are in the taxon look up database
cat(paste0("Are all the acceptedNameUsageID's in occurrence data found in the taxon_df: ", all(occurences$acceptedNameUsageID %in% taxonomy_df$taxonid_tax)), "\n")

# Printing a random subset of the occurences data which are not found in the taxon_df
if (all(occurences$acceptedNameUsageID %in% taxonomy_df$taxonid_tax) == FALSE) {
  cat("These are examples of the acceptedNameUsage's where the acceptedNameUsageID are not in the taxon_df \n")
  cat("For more infotion, see the file called acceptedNameUsageID_not_in_taxon_df.txt \n")
  write.table(occurences[!occurences$acceptedNameUsageID %in% taxonomy_df$taxonid_tax, ], "acceptedNameUsageID_not_in_taxon_df.txt", sep = "\t", row.names = FALSE)
  rand_ind <- sample(nrow(occurences[!occurences$acceptedNameUsageID %in% taxonomy_df$taxonid_tax]), size = 5)
  cat("Random sample \n")
  print(occurences[!occurences$acceptedNameUsageID %in% taxonomy_df$taxonid_tax, ][rand_ind, c("family", "genus", "species", "scientificName")])
  cat("\n \n")
  rm(rand_ind)
}


# Here were joining the occurrences with the taxonomy data gathered
# This means were joing the taxonomy data from GBIF with the occurrences data from GBIF.
# Maybe I should remove the all.x = true? This would mean that I remove all the rows which are not in the taxonomy data. # , all.x=TRUE
occurrences_taxonomy <- merge(occurences, distinct(taxonomy_df), by.x = "acceptedNameUsageID", by.y = "taxonid_tax")

# This is a sample from occurence_taxonomy
cat("This is a sample from occurence_taxonomy \n")
rand_ind <- sample(nrow(occurrences_taxonomy), size = 5)
print(occurrences_taxonomy[rand_ind, ])
cat("\n \n")
rm(rand_ind)


#################################################################################
#                              set up common format                             #
#################################################################################

cat("Setting up the framework for the common format \n")
split_length <- unlist(lapply(strsplit(as.character(occurrences_taxonomy$species_tax), split = " "), length))
cat("Number of rows in split_length is ", length(split_length), " \n")

#Checking the number of rows for each of the variables that I want to include in my dataframe
cat("Checking the number of rows for each of the variables that I want to include in my dataframe \n")
list_of_vars <- c("species", "scientificname_tax", "family_tax", "authorship_tax", "genus_tax", "species_tax", "rank_tax")
for (i in list_of_vars){
  cat(paste0("The number of rows in ", i, " is ", length(occurrences_taxonomy[[i]]), "\n"))
}

# Selecting columns from the occurrences_taxonomy data frame
input <- data.frame(occurrence_species_name = occurrences_taxonomy$species,
  species_tax = occurrences_taxonomy$species_tax,
  taxon_name = occurrences_taxonomy$scientificname_tax,
  family = occurrences_taxonomy$family_tax,
  author = occurrences_taxonomy$authorship_tax,
  split_length = split_length,
  genus_hybrid = NA,
  species_hybrid = NA,
  species = sub(".* ", "", occurrences_taxonomy$species_tax),
  genus = occurrences_taxonomy$genus_tax,
  taxon_rank = tolower(occurrences_taxonomy$rank_tax),
  infra_name = NA,
  comment = NA,
  usable = NA
)

## Cleaning
cat("Doing some cleaning \n")

# remove sp. species
cat("Removing sp. species \n")
input <- input[!grepl("sp\\.", input$occurrence_species_name), ]

# remove aggregates
cat("Removing aggregates \n")
input <- input[!grepl("aggr$", input$occurrence_species_name), ]

# order the dataframe by split length
cat("Ordering the dataframe by split length \n")
input <- input[order(input$split_length), ]

# Removing duplicate rows
cat("Removing duplicate rows \n")
input <- input[!duplicated(input$occurrence_species_name), ]

cat("Assigning an id column \n")
input$id <- c(1:nrow(input)) # THIS IS WHERE THE ID COLUMN IS FROM!!!! and it is just the row number. # nolint

# create split list
# which is a list of the tip names split by space and ordered by ID
cat("Creating a split list \n")
split_list <- strsplit(as.character(input$species_tax), split = " ")
cat("Naming the split list based on the id column \n")
names(split_list) <- input$id

cat("This is the unique values in split_length in the input dataframe \n")
print(unique(input$split_length))

###############################################################
#                  Vectorize + loop mix                       #
###############################################################
cat(" Vectorize and loop mix \n")

## split length == 0
cat("Split length == 0 \n")
ind <- which(input$split_length == 0)
cat("The length of indexes with split length == 0 is ", length(ind), "\n")
input$usable[input$split_length == 0] <- "no"

## split length == 1
cat("Split length == 1 \n")
ind <- which(input$split_length == 1)
cat("The length of indexes with split length == 1 is ", length(ind), "\n")
input$usable[input$split_length == 1] <- "no"

## split length == 2
cat("Split length == 2 \n")
ind <- which(input$split_length == 2)
cat("The length of indexes with split length == 2 is ", length(ind), "\n")
ind_hyb <- which(grepl("^×", input$genus))
input$genus[ind_hyb] <- gsub("^×", "", sapply(split_list[ind_hyb], "[[", 1))
input$genus_hybrid[ind_hyb] <- "x"
input$usable[ind_hyb] <- "no"

## split length == 3
cat("Split length == 3 \n")
ind <- which(input$split_length == 3)
cat("The length of indexes with split length == 3 is ", length(ind), "\n")
if (length(ind >= 1)) { # nolint
  for (i in 1:length(ind)){ # nolint: seq_linter.
    print("Printing split list")
    print(split_list[[ind[i]]])
    if (split_list[[ind[i]]][2] == "sp.") {
      input$usable[ind[i]] <- "no"
    }
    if (split_list[[ind[i]]][2] == "x") { # does not occur but you never know
      if (grepl("[A-Z]", split_list[[ind[i]]][3])) {
        input$genus_hybrid[ind[i]] <- "x"
        input$usable[ind[i]] <- "no"
      }else {
        input$species_hybrid[ind[i]] <- "x"
      }
    }
    if (split_list[[ind[i]]][1] != "x" && split_list[[ind[i]]][2] != "x" && !grepl("[A-Z]", split_list[[ind[i]]][3])) {
      input$infra_name[ind[i]] <- split_list[[ind[i]]][[3]]
    }
    if (split_list[[ind[i]]][2] != "x" && grepl("[A-Z]", split_list[[ind[i]]][3])) {
    }
    if (grepl("sp\\.", split_list[[ind[i]]][[2]])) {
      input$usable[ind[i]] <- "no"
    }
  }
}

cat("How many indexes have usable == no \n")
nrow(which(input$usable == "no"))

# remove all unusable taxa
#cat("Removing all unusable taxa \n")
cat("nrow before removing unusable taxa ", nrow(input), "\n")
if (any(!is.na(input$usable))) {
  cat("There are some unusable taxa \n")
  input <- input[which(input$usable != "no"), ] # why does this line of code remove all the data?
} else {
  cat("There are no unusable taxa \n")
}

cat("nrow after removing unusable taxa ", nrow(input), "\n")

# Remove the usable column
cat("Removing the usable column \n")

cat("nrow before removing the usable column ", nrow(input), "\n")

# Attemping to remove the usable column from the dataframe without removing all the entries in the dataframe.
input <- input[,-which(names(input) == "usable")]

cat("nrow after removing the usable column ", nrow(input), "\n")

# Removing the duplicate columns
cat("Removing the duplicate columns \n")
cat("nrow before removing the duplicate columns ", nrow(input), "\n")
input <- distinct(input)
cat("nrow after removing the duplicate columns ", nrow(input), "\n")

# Saving the file
saveRDS(input, output_file)

# How do I remove a column from a dataframe in R
# You can use the following syntax to remove a column from a data frame:
# mydataframe <- mydataframe[,-columnnumber]