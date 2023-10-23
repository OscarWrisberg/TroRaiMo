
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
#sort(unique(wcp$family)) # Length 459

#Families in wcvp which are not in the file apg synonym families
#unique(wcp$family[!wcp$family %in% apg$Syn_Fam])

print(apg)

# I need to remove the ""'s around the family name in the wcvp file
# I think this is the problem.
wcp$family <- gsub('"', '', wcp$family)

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

cat("Families in wcvp which are not in the file apg synonym families \n")
unique(wcp$family[!wcp$family %in% apg$Syn_Fam])

# This should update family column of our WCVP data to be right according to APGIV
wcp$family.apg <- apg$Acc_Fam[match(wcp$family, apg$Syn_Fam)]

#Why is there some NA's in the family.apg column?
# There should'nt be any NA's in the family.apg column.
no_na_family.apg <- sum(!is.na(wcp$family.apg))
if ( no_na_family.apg == nrow(wcp)) {
  cat("There are no NA's in the family.apg column \n")
} else {
  cat("There are ", nrow(wcp) - no_na_family.apg, " NA's in the family.apg column, out of a total of ", nrow(wcp), " rows \n")
}

# Printing the rows where the family column is not NA but the family.apg column is NA.
row_fam_no_fam.apg <- wcp[!is.na(wcp$family) & is.na(wcp$family.apg),]
cat("Total rows where the family column is not NA but the family.apg column is NA ", nrow(row_fam_no_fam.apg), "\n")
cat("The rows where the family column is not NA but the family.apg column is NA \n")
print(row_fam_no_fam.apg)


saveRDS(wcp, output_file)


