#!/usr/bin/env Rscript
# read Trachaeophyta gbif data. 120gb file

library(data.table)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1]) # here you define the name of your input file
output_file <- as.character(args[2]) # Here you define the name of the output file

dat <- fread(input_file, header=TRUE, sep="\t",verbose=TRUE, quote="")
    

cat(paste0("fread is done reading in the file at: ", Sys.time(),"\n"))


cat(paste0("Commence saving file as RDS at: ", Sys.time(),"\n"))

saveRDS(dat,output_file)

cat(paste0("Finished saving file as RDS at: ", Sys.time(),"\n"))

rm(dat)
