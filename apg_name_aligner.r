
# Command Line arguments
args <- commandArgs(trailingOnly = TRUE)
apg_in <- as.character(args[1]) # here you define the name of your input file
wcp_in <- as.character(args[2]) # here you define the name of your input file
output_file <- as.character(args[3]) # Here you define the name of the output file


# Loading packages
library(data.table)
library(dplyr)


# Loading name files
apg <- fread(apg_in)
wcp <- fread(wcp_in)


#Families in wcvp
sort(unique(wcp$family)) # Length 459

#Families in wcvp which are not in the file apg synonym families
unique(wcp$family[!wcp$family %in% apg$Syn_Fam])

# Incertae_sedis (of unknown placement)
wcp <- wcp[which(wcp$family != "Incertae_sedis"),]

# Isoetaceae - ( probably a problem with ö)
wcp[which(wcp$family == "Isoetaceae"),]$family <- "Isoëtaceae"

# Osmundaceae - a fern family so it is probably also right

# Gigaspermaceae = moss family
wcp <- wcp[which(wcp$family != "Gigaspermaceae"),]

# Ripogonaceae - Could also be named Rhipogonaceae
wcp[which(wcp$family == "Ripogonaceae"),]$family <- "Rhipogonaceae"

# Tiganophytaceae - includes monotypic genus, but should be alright


# Fixing and removal of names / families
wcp <- wcp[!wcp$family=="Pseudotubulare",]

unique(wcp$family[!wcp$family %in% apg$Syn_Fam])

# This should update family column of our WCVP data to be right according to APGIV
wcp$family.apg <- apg$Acc_Fam[match(wcp$family, apg$Syn_Fam)]

saveRDS(wcp, output_file)


