#!/usr/bin/env Rscript
# read Trachaeophyta gbif data. 120gb file

library(data.table)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1]) # here you define the name of your input file
chunk_size <-as.numeric(args[2]) # here you have to define the size of each chunk
#chunk_nr <- as.numeric(args[3]) # here you define the number of the chunk
total_lines <- as.numeric(args[3]) # here you define the total number of lines in the file 
output_file <- as.character(args[4]) # Here you define the name of the output file


#wd <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(wd)

# data prep ------------------------------------------------------------------
# ----------------------------------------------------------------------------
# -- the filtered Trachaeophyta download has 19,044,905 occurrences. Handle wisely! --
# ----------------------------------------------------------------------------

# There are 19,044,905 lines in the file with occurrences and habitat description.
# This number can be divided by 11 to achieve a whole number
# process the data chunk by chunk since it is huge.

# These are the variables which will be kept after the parsing.
#keep <- c("gbifID", "establishmentMeans", "habitat", "eventRemarks", "decimalLatitude", "decimalLongitude", "scientificName", "family", "genus", "taxonRank", "datasetKey", "species")


# define start and end based on the number of chunks and the chunk size
# end_of_chunk = chunk_nr*chunk_size
# start_of_chunk = end_of_chunk-chunk_size

# cat(paste0("This is the start of this chunk ", start_of_chunk, "\n\n"))
# cat(paste0("This is the end of this chunk ", end_of_chunk, "\n\n"))


## get chunk ----
#if(chunk_nr == 1){ # first chunk gets col names
    # dat <- read_tsv(input_file,  
    #                 n_max = chunk_size, skip=0, show_col_types = FALSE, skip_empty_rows = FALSE) # Remember to change the keep_num variable n_ma
dat <- fread(input_file, header=TRUE, sep="\t",verbose=TRUE, quote="")
    
    
# }else{ # add col names to other chunks for easy processing
#     # dat <- read_tsv(input_file,  
#     #                 n_max = chunk_size, skip=start_of_chunk, show_col_types = FALSE, skip_empty_rows = FALSE)
#     dat <- fread(input_file, header=FALSE, # Remember to change the keep_num variable
#                  sep="\t",verbose=TRUE, quote="", nrows = chunk_size, skip=start_of_chunk)
    
#names(dat) = keep
#	}

cat(paste0("fread is done reading in the file at: ", Sys.time(),"\n"))


cat(paste0("Commence saving file as RDS at: ", Sys.time(),"\n"))

saveRDS(dat,output_file)

cat(paste0("Finished saving file as RDS at: ", Sys.time(),"\n"))

rm(dat)
