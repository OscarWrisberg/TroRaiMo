#wd <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(wd)

#############################################################################
######################## Command line arguments #############################
#############################################################################

args <- commandArgs(trailingOnly = TRUE)
input_file <- as.character(args[1]) # here you define the name of your input file
output_file <- as.character(args[2]) # Here you define the name of the output file

library(rgbif)
library(data.table)

print(paste0("The input file is ", input_file))


# Prep unmatched tip labels for taxonomy matching ------------------------
# Load Smith & Brown tip labels not covered by NCBI (-->GBIF)
sb <- readRDS(input_file)

cat(paste0("This is the number of unique names  ",length(unique(sb$acceptedNameUsageID)), "\n"))

total_ids <- length(unique(sb$acceptedNameUsageID))


#GBIF download to complete tax info --------------------------------------
#THIS TAKES A WHILE!
#Carefull, If i change this from acceptedNameUsageID to something else, then I have to look up lots more ID's
# But by using the acceptedNameUsageID I kinda remove the option for Wcvp to have a different opinion on the correct name for something that isent currently accepted.

ids <- unique(sb$acceptedNameUsageID)

res <- list()
for(i in 1:length(ids)){
  tryCatch({q
    res[[i]] <- name_usage(ids[i], data = "all")
    if(!i%%1000) cat(i,"\r")

    # failsafe
    if(i %in% c(100,50000, 100000, 150000, 200000,250000,300000,350000,4000000)){saveRDS(res, "midway_gbif.rds")}
    if(i==length(ids)){saveRDS(res, output_file)}
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}




