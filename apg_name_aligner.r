# Loading packages
# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table","dplyr","tidyr")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Files for local testing
# setwd("/home/au543206/GenomeDK/Trf_models/data/")
# apg_in <- "/home/au543206/GenomeDK/Trf_models/TroRaiMo/apgweb_parsed.csv"
# wcp_in <- "wcvp_names.csv"
# output_file <- "wcvp_names_apg_aligned.rds"


# Command Line arguments
args <- commandArgs(trailingOnly = TRUE)
apg_in <- as.character(args[1]) # here you define the name of your input file
wcp_in <- as.character(args[2]) # here you define the name of your input file
output_file <- as.character(args[3]) # Here you define the name of the output file


# Loading name files
apg <- fread(apg_in)
wcp <- fread(wcp_in)


#Families in wcvp
#sort(unique(wcp$family)) # Length 459

#Families in wcvp which are not in the file apg synonym families
#unique(wcp$family[!wcp$family %in% apg$Syn_Fam])

# This might solve why the family.apg column has NA's in it.
apg <- as.data.frame(apg)
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

cat("Families in apg$syn which are not in the wcvp file \n")
unique(apg$Syn_Fam[!apg$Syn_Fam %in% wcp$family])


# This should update family column of our WCVP data to be right according to APGIV
wcp$family.apg <- apg$Acc_Fam[match(wcp$family, apg$Syn_Fam)]

# Trying with a for loop instead of a match function.
# for(family in unique(wcp$family)){
#   if(family %in% apg$Syn_Fam){
#     wcp[which(wcp$family == family),]$family.apg <- apg[which(apg$Syn_Fam == family),]$Acc_Fam
#     cat("Family ", family, " has been updated to ", apg[which(apg$Syn_Fam == family),]$Acc_Fam, "\n")
#   }else{
#     cat("Family ", family, " is not in the apg file \n")  }
# }

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


