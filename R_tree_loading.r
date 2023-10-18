#!/usr/bin/env Rscript

# Loading packages
library(data.table)
library(readr)
library(ape)

# Loading the command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1]) # here you define the name of your input file
output_file_tips <- as.character(args[2]) # here you define the name of your output file

# Loading the Smith and Brown 2018 GBMB data
tree <- read.tree(input_file)

# Creating an object which is the tip names in this tree
tips <- tree$tip.label

#! save the tips object as a csv file
write.table(tips[,2], output_file_tips, sep = "\t", row.names = FALSE, col.names = FALSE)


