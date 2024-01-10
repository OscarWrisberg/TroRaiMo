
# This sctipt is going to loop through a list of orders in a tree which are non-monophyletic
# It is going to find the most recent common ancestor (MRCA) of the tips in the order.
# Because the orders are non-monophyletic the decendants of the MRCA will be both species in the order and species outside the order.

# This script will then find out how many tips are descendants of the MRCA which are not supposed to be in the order.

# Depending on the results. This script will then consider pruning the rogue tips which are in the order but not supposed to be in the order.
# or it will remove tips which are supposed to be in the order but are placing elsewhere in the tree.


#########################################################################################################################
################################################## Loading packages  ####################################################
#########################################################################################################################

# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("data.table", "ape", "phytools", "geiger", "castor", "MonoPhy")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


###########################################################################################################################
##################  Testing the code by runnin it on GDK through VScode and its built in terminal  ########################
###########################################################################################################################

# setwd("/home/au543206/GenomeDK/Trf_models/data") # Set working directory when local
# wcvp <- readRDS("../workflow/02_adding_orders/wcvp_names_apg_aligned.rds")  # Read the WCVP names file into a data frame
# tree <- read.tree("GBMB_pruned.tre") # Read the GBMB pruned tree
# output_path <- "../workflow/02_adding_orders/pruning/"
# apg <- fread("../TroRaiMo/apgweb_parsed.csv")
# tips_families <- fread(paste0(output_path,"tips_families.txt"))
# non_monophyletic_orders <- fread("../workflow/02_adding_orders/pruning/non_mono_order.txt", header = FALSE, sep = "\t")
# output_file <- "Orders_which_could_not_be_solved.txt"

###########################################################################################################################
############################# Getting command line file names for workflow ################################################
###########################################################################################################################

# # Setting the wd for the script
setwd("/home/owrisberg/Trf_models/data") # Set working directory when remote

# Command line arguments
input_file_tree <- commandArgs(trailingOnly = TRUE)[1]
input_file_wcvp <- commandArgs(trailingOnly = TRUE)[2]
output_file <- commandArgs(trailingOnly = TRUE)[3]
output_path <- commandArgs(trailingOnly = TRUE)[4]
apg <- commandArgs(trailingOnly = TRUE)[5]

# Read the WCVP names file into a data frame
cat("Opening ", input_file_wcvp, "\n")  
wcvp <- readRDS(input_file_wcvp)

# Read the GBMB tree
cat("Opening ", input_file_tree, "\n")
tree <- read.tree(input_file_tree)

# Loading the apgweb_parsed.csv file
cat("Loading the apgweb_parsed.csv file \n")
apg <- fread(apg)

 # Loading tips families so I dont have to wait so fucking long..
tips_families <- fread(paste0(output_path,"tips_families.txt"),)


# Loading the list of non-monophyletic orders
cat("Loading the list of non-monophyletic orders \n")
non_monophyletic_orders <- fread("../workflow/02_adding_orders/pruning/non_mono_order.txt", header = FALSE, sep = "\t")


###################################################################################
###########################  Basic clean up  ######################################
###################################################################################

# Make tree bifurcating
tree <- multi2di(tree)

# Remove "_"
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

# Dropping the x order
non_monophyletic_orders <- non_monophyletic_orders[which(non_monophyletic_orders$V1 != "x"),]

# Removing some orders of Ferns which I am not interested in and which are causing problems
fern_orders <- c("Gleicheniales", "Hymenophyllales", "Lycopodiales", "Ophioglossales", "Schizaeales", "Marsileales", "Polypodiales-eupolypod_I", "eupolypod_II", "Polypodiales")
non_monophyletic_orders <- non_monophyletic_orders[which(!non_monophyletic_orders$V1 %in% fern_orders),]

# Find unique families
unique_families <- unique(tips_families$families)
unique_families <- as.character(unique_families)
cat("This is the number of unique_families \n")
cat(length(unique_families), "\n")

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))
non_mono_family <- character(0)

# Dropping tips which are outlier taxa based on MonoPhy
outlier_taxa <- c("Androya decaryi","Pteleocarpa lamponga", "Tectaria heracleifolia", "Streptopus parviflorus","Streptopus amplexifolius","Streptopus koreanus",
"Streptopus lanceolatus","Streptopus ovalis","Streptopus obtusatus","Streptopus parasimplex")
tree <- drop.tip(tree, tip = outlier_taxa)

# Fixing a wrong family in tips_families
tips_families[which(tips_families$name == "Saussurea japonica"),"families"] <- "Asteraceae"

####################################################################################
####################  Finding the order for each family  ###########################
####################################################################################

# Find order function
find_order <- function(fams, apg) {
  fam_list <- character(0)
  orders <- character(0)

  for (family in unique(fams)) {
    order <- apg[which(apg$Syn_Fam == family), "Clade"] # Finding the order of that family in APG file
    order <- as.character(order[1]) # Selecting the order of the family
    cat("family ", family, "Order", order, "\n") # Printing the family and order
    fam_list <- c(fam_list, family)
    orders <- c(orders, order)
  }

  df_orders <- data.frame(family = fam_list, order = orders)
  return(df_orders)
}

# Running the function
cat("Finding the order for each of the families \n")
family_orders <- find_order(unique_families, apg)
length(family_orders$order)

# Merging the tips_families and family_orders data frames
cat("Merging the tips_families and family_orders data frames \n")
tips_family_orders <- merge(tips_families, family_orders, by.x = "families", by.y = "family")
unique(tips_family_orders$order)

# saving the tips_families_orders
cat("Saving the tips_families_orders \n")
write.table(tips_family_orders, paste0(output_path, "tips_family_orders.txt"), sep = "\t", row.names = FALSE)


##########################################################################################################################
############################  Defining function for finding the largest clades############################################
##########################################################################################################################

# Function for finding the largest monophyletic clade of the tips in the order
find_largest_clade <- function(tips_in_order, tree) {
  biggest_subtree <- list()
  rogue_species <- list()
  length_biggest_subtree <- 0

  # Looping through all the tips in the order
  for (i in seq_along(tips_in_order)) {

    # Progress bar
    if (!i %% 100) cat("Percentage done", format(round((i / length(tips_in_order)) * 100, 2), nsmall = 2), " at ", format(Sys.time(), '%H:%M:%S'), "\n")
    tip <- tips_in_order[i]

    # Find the node of the tip because getParent only works on nodes
    node <- which(tree$tip.label == tip)[1]

    # Get the parent node
    pnode <- getParent(tree, node)

    # Check if the parent node is the root
    if (pnode != 0) {
      subtree <- get_subtree_at_node(tree, pnode - Ntip(tree))$subtree

      # Check if there are any tips in the subtree & if they are all found in the order
      if (length(subtree$tip.label) > 0 && all(subtree$tip.label %in% tips_in_order)) {

        # While all the tip labels are in the order continue "diving" into the tree until you come to a species which is not in the order
        while (all(subtree$tip.label %in% tips_in_order)) {
          last_tree <- subtree
          last_pnode <- pnode

          pnode <- getParent(tree, pnode)
          subtree <- get_subtree_at_node(tree, pnode - Ntip(tree))$subtree
        }

        # Select the last tree which contained only tips in the order

        # Checking if the tips in the subtree are all in the order
        if (all(last_tree$tip.label %in% tips_in_order)) {
          # cat("All the tips in the subtree are in the order \n")
        } else {
          # cat("Not all the tips in the subtree are in the order \n")
          break
        }

        # When you find a node where not all the tips are in the order.
        # If the clade is bigger than the biggest clade found so far then save it as the biggest clade
        # And update the rogue species in order to see what species are breaking the order

        # Calculate the proportion of tips in the order
        proportion_in_order <- length(last_tree$tip.label[which(last_tree$tip.label %in% tips_in_order)]) / length(last_tree$tip.label)

        # Check if the proportion is at least 90%
        if (proportion_in_order >= 0.9) {
          if (length(last_tree$tip.label) > length(biggest_subtree$tip.label)) {
            biggest_subtree <- last_tree
            length_biggest_subtree <- length(biggest_subtree$tip.label)
            rogue_species <- subtree$tip.label[which(!subtree$tip.label %in% tips_in_order)]
          }
        }
      }
    }
  }

  cat("Length of biggest subtree is: ", length_biggest_subtree, "\n")

  # Check if biggest_subtree is non-empty before returning
  if (length(biggest_subtree) > 0) {
    return(list(biggest_subtree = biggest_subtree, rogue_species = rogue_species, subtree = subtree))
  } else {
    return(NULL)
  }
}


#################################################################################################################################################
############################ Looping through all non monophyletic orders to see if we can make them monophyletic#################################
#################################################################################################################################################

	# Thought process behind this loop

	# First we find the Most Recent Common Ancestor (MRCA) of all the tips in the order which we are investigating.
	# Then because these orders are non-monophyletic the MRCA will have descendants which are both in the order and outside the order.
	# We find the number of tips which are outside the order (Rogue tips)

	# The proportion of Rogue tips in the tree created from the MRCA will define how we continue.
	# I have decided that the proportion of Rogue tips I am willing to tolerate is 10 %

	# if the number of Rogue tips is smaller than 10% of the number of tips in the order
	# Then it means we have a few tips which are placed in the order clade but which are not supposed to be in the order clade.
	# Then we can just justifiably prune the rogue tips from the order clade in order to get a monophyletic tree of the order.

	# if the number of rogue tips is larger than 10% of the number of tips in the order
	# The problem is likely that we have a tip which SHOULD be placed in the order clade but is not.
	# and the pruning of this species would probably result in a monophyletic order.
	# In order to deal with these species which should be in the order but which are not
	# We can loop through the tree to find the largest monophyletic clade which is in the order.
	# If this largest clade then contains at least 90 % of the tips in the order we can save it as the monophyletic tree of the order.

	# Orders where the largest clade does not contain at least 90 % of the tips in the order will be saved to a file
	# Potentially I can somehow later join some orders in order to get a monophyletic clade of two or more orders.

# Create a dataframe to store the orders which I could not solve
no_solvable_tips_family <- data.frame(order = character(0), max_clade_found = numeric(0), total_tips_in_order = numeric(0),max_clade_coverage = numeric(0), rogue_tips = character(0))

# Looping through the non-monophyletic orders
for (i in seq_along(non_monophyletic_orders[[1]])) {

	# Finding the order
	order <- non_monophyletic_orders[[1]][i]

	# Determining if all the tips are present and have an order
	#tips_in_order_complete <- tips_family_orders$name[complete.cases(tips_family_orders$order) & tips_family_orders$order == order]

	# Finding the tips in the order
	cat("Starting on ",order, " \n")
	tips_in_order_1 <- tips_family_orders$name[which(tips_family_orders$order == order)] # Selecting the tips which are in the selected order
	cat("Number of tips in the order: ", length(tips_in_order_1), "   ") 

	# Finding all the tips in the tree which are in the order
	tips_in_order <- tips_in_order_1[which(tips_in_order_1 %in% tree$tip.label)]

	#Counting the number of tips in the order which are in the tree
	cat("Number of tips in the order which are in the tree: ", length(tips_in_order), "  \n")
	
	# Finding the MRCA of the tips in the order
	MRCA <- getMRCA(tree, tips_in_order)

	# Finding the tips which are descendants of the MRCA
	descendants <- tips(tree, MRCA)
	cat("The number of descendants from the MRCA of all tips in ",order," is ", length(descendants), " ")

	# finding the tips which are descendants of the MRCA and are not in the order
	rogue_tips <- descendants[which(!descendants %in% tips_in_order)]
	cat(" and the number of rogue tips is:", length(rogue_tips), "\n ")


	cat("Checking if the number of rogue tips is smaller than 10 % of the number of tips in the order \n")
	cat("and if the number of rogue tips is larger than 0 \n")
	cat("length(rogue_tips) : ", length(rogue_tips), "\n")
	cat("length(tips_in_order) : ", length(tips_in_order), "\n")
	cat("length(rogue_tips) <= 0.1 * length(tips_in_order) : " ,length(rogue_tips) <= 0.1 * length(tips_in_order), " \n" )
	cat("length(rogue_tips) > 0 : ", length(rogue_tips) > 0, "\n")

	# If the number of rogue tips is smaller than 10 % of the number of tips in the order
	if ( length(rogue_tips) <= (0.1 * length(tips_in_order)) & length(rogue_tips) > 0){
		cat("The number of rogue tips is smaller than 10 % of the number of tips in the order \n")
		cat("Pruning the rogue tips and extracting the order \n")
		tree <- drop.tip(tree, tip = rogue_tips)

		# Check if the tips in each order form a monophyletic clade
		if( is.monophyletic(tree, tips_in_order) == TRUE){
				order_tree <- drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% tips_in_order])
				 # Save the pruned tree to a file
  				write.tree(order_tree, paste0(output_path,"orders/", "pruned_tree_order", order, "_GBMB.tre"))
				cat("Done with ", order, "\n")
		}

	} else if (length(rogue_tips) == 0) { # Trees which have already had the rogue tips pruned
		cat("The number of rogue tips is 0 \n")
		cat("Extracting the order \n")
		order_tree <- drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% tips_in_order])
		write.tree(order_tree, paste0(output_path,"orders/", "pruned_tree_order", order, "_GBMB.tre"))
		next
	} else {
		# Here I will loop through the tips in the order and find the largest monophyletic clade which is in the order.

		result <- find_largest_clade(tips_in_order, tree)

		largest_clade <- result$biggest_subtree
		rogue_species_breaking <- result$rogue_species
		subtree <- result$subtree

		# If the largest clade contains atleast 90 % of the species in the order. I will save it
		if ( length(largest_clade$tip.label) >= 0.9 * length(tips_in_order) ) {
			cat("The largest monophyletic clade in the order contains atleat 90 % of the tips in the order \n")
			cat("Saving the largest clade \n")
  			write.tree(largest_clade, paste0(output_path,"orders/", "pruned_tree_order", order, "_GBMB.tre"))
			cat("Done with ", order, "\n")
			
		# If it dosent I will save the order to a file and save the number of tips in the largest clade and the total number of tips in the order
		# Later down the line I can then see if I can fix some of these orders by 
		} else {
			cat("Problem order is: ", order, "\n")

			if(length(rogue_tips) <= 64486) { # total number of tips in the tree
				cat("Assesing the monophyly of ", order, "\n")
			rogue_sub_tree <- ape::extract.clade(tree, MRCA)
			rogue_tips_orders <- tips_family_orders[which(tips_family_orders$name %in% rogue_sub_tree$tip.label)] # Selecting the tips which are in the MRCA tree
			rogue_tips_orders <- rogue_tips_orders[, c("name", "order")]
			monophy <- AssessMonophyly(rogue_sub_tree, rogue_tips_orders)

			print(monophy)
			
			height_per_species <- 0.15

			#Calculate the total height based on the number of species
			if ( length(rogue_sub_tree$tip.label) < 10) {
				total_height <- 50
			 } else {
				 total_height <- round(length(subtree$tip.label) * height_per_species)
			}
			cat("Total height is: ", total_height, "\n")


			# Open PDF file for plotting
			pdf(paste0(output_path, "rogue_monophy_", order, "_GBMB.pdf"), height = total_height)
			
			# Plot using PlotMonophyly
			PlotMonophyly(monophy, subtree, plot.type='monophyly', ladderize=TRUE, cex=0.5)
			
			# Close the PDF file
			dev.off()
			} else {
				cat("The number of rogue tips is too large to asses monophyly \n")
			}

			new_row <- data.frame(order = order,
				max_clade_found = length(largest_clade$tip.label),
				total_tips_in_order = length(tips_in_order),
				max_clade_coverage = length(largest_clade$tip.label)/length(tips_in_order),
				rogue_tips = paste(rogue_species_breaking, collapse = ","))

    		no_solvable_tips_family <- rbind(no_solvable_tips_family, new_row)

			cat("Done with ", order, "\n")
			
			}
		#write.tree(subtree, paste0(output_path, "Rogue_MRCA_tree_", order, "_GBMB.txt")) # biggest 
	}
}




cat("Printing output_file to: ",paste0(output_path,output_file), "\n")
# Writing out the no_solvable_tips_family data frame to a file
write.table(no_solvable_tips_family, paste0(output_path,output_file), sep = "\t", row.names = FALSE) #nolint



