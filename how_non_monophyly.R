
# This sctipt is going to loop through a list of orders in a tree which are non-monophyletic
# It is going to find the most recent common ancestor (MRCA) of the tips in the order.
# Because the orders are non-monophyletic the decendants of the MRCA will be both species in the order and species outside the order.

# This script will then find out how many tips are descendants of the MRCA which are not supposed to be in the order.

# Depending on the results. This script will then consider pruning the rogue tips which are in the order but not supposed to be in the order.
# or it will remove tips which are supposed to be in the order but are placing elsewhere in the tree.

# Loading packages
library(data.table)
library(ape)
library(phytools)
library(geiger)
library(castor)

###########################################################################################################################
#Testing the code by runnin it on GDK through VScode and its built in terminal
setwd("/home/owrisberg/Trf_models/data") # Set working directory when local
setwd("/home/owrisberg/Trf_models/data") # Set working directory when remove
wcvp <- readRDS("../workflow/02_adding_orders/wcvp_names_apg_aligned.rds")  # Read the WCVP names file into a data frame
tree <- read.tree("GBMB.tre") # Read the GBMB tree
path_out <- "../workflow/02_adding_orders/pruning/"

###########################################################################################################################

# Loading the tree
# cat("Loading the tree \n")
# tree <- read.tree(tree)

# Removing the _ and " from the tip labels
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub('"', '', tree$tip.label)  # nolint

# Loading the list of non-monophyletic orders
cat("Loading the list of non-monophyletic orders \n")
non_monophyletic_orders <- fread("../workflow/02_adding_orders/pruning/non_mono_order.txt", header = FALSE, sep = "\t")
# Dropping the x order
non_monophyletic_orders <- non_monophyletic_orders[which(non_monophyletic_orders$V1 != "x"),]
 
 # Loading tips families so I dont have to wait so fucking long..
tips_families <- fread("tips_families.txt")

# Find unique families
unique_families <- unique(tips_families$families)
unique_families <- as.character(unique_families)
cat("This is the number of unique_families \n")
cat(length(unique_families), "\n")

# Create a data frame to store the number of tips in each family
df_number_tips <- data.frame(family = character(0), number_tips = numeric(0))
non_mono_family <- character(0)

# Loading the apgweb_parsed.csv file
cat("Loading the apgweb_parsed.csv file \n")
apg <- fread("../TroRaiMo/apgweb_parsed.csv")

# Find order function
cat("Creating the find_order function \n")
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
cat("Running the function \n")
family_orders <- find_order(unique_families, apg)

# Merging the tips_families and family_orders data frames
cat("Merging the tips_families and family_orders data frames \n")
tips_family_orders <- merge(tips_families, family_orders, by.x = "families", by.y = "family")


#####################################################################
# Function for finding the largest monophyletic clade in an order that is not monophyletic
find_largest_clade <- function(tips_in_order, tree){
	biggest_subtree <- character(0)
	for(i in seq_along(tips_in_order)){
	tip <- tips_in_order[i]
	#print(tip)

    cat(which(tips_in_order == tip), "\r")

    node <- which(tree$tip.label==tip)

    pnode <- getParent(tree, node)

    subtree <- get_subtree_at_node(tree, pnode-Ntip(tree))$subtree

    # condition: any of the tips in this subtree not a species in this order
    if(any(!subtree$tip.label %in% tips_in_order)){
		# We need to save the subtree if it is bigger than the current biggest subtree
		if(length(subtree$tip.label) > length(biggest_subtree$tip.label)){
			biggest_subtree <- subtree
		}
    }

    # condition: all of the tips in this subtree an island endemic?
    if(all(subtree$tip.label %in% tips_in_order)){
      
      while(all(subtree$tip.label %in% tips_in_order)){
        #print("down the rabbit hole")
        last_tree <- subtree
        last_pnode <- pnode

        pnode <- getParent(tree, pnode)
        subtree <- get_subtree_at_node(tree, pnode-Ntip(tree))$subtree 
      }
      
	  biggest_subtree <- subtree
		# Here I should save the last tree if it is bigger than the current biggest subtree
    	}
	}
return(biggest_subtree)
}
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


rogue_tips_family <- data.frame(order = character(0), rogue_tips = numeric(0))

# Looping through the non-monophyletic orders
for (i in seq_along(non_monophyletic_orders[[1]])) {
	# Finding the order
	order <- non_monophyletic_orders[[1]][i]

	# Determining if all the tips are present and have an order
	#tips_in_order_complete <- tips_family_orders$name[complete.cases(tips_family_orders$order) & tips_family_orders$order == order]

	# Finding the tips in the order
	print(order)
	tips_in_order <- tips_family_orders$name[which(tips_family_orders$order == order)]
	cat("Number of tips in the order: ", length(tips_in_order), "   ")

	# Finding the MRCA of the tips in the order
	if(order == "Gentianales"){
		cat("Skipping Gentianales \n")
		next
	}else {
	MRCA <- getMRCA(tree, tips_in_order)
	}
	#print(MRCA)

	# Finding the tips which are descendants of the MRCA
	descendants <- tips(tree, MRCA)


	# finding the tips which are descendants of the MRCA and are not in the order
	rogue_tips <- descendants[which(descendants %in% tips_in_order == FALSE)]

	cat("Number of rogue tips is:", length(rogue_tips), "\n ")

	# if the number of rogue tips is smaller than some proportion of the number of tips in the order
	# Then we can just justifiably prune the rogue tips in order to get monophyletic orders

	# if the number of rogue tips is larger than some proportion of the number of tips in the order
	# The problem is likely that we have a tip which SHOULD be placed in the order but is not.
	# and the pruning of this species would probably result in a monophyletic order.

	if ( length(rogue_tips) <= 0.1 * length(tips_in_order) & length(rogue_tips) > 0){
		cat("The number of rogue tips is smaller than 10 % of the number of tips in the order \n")
		cat("Pruning the rogue tips and extracting the order \n")
		tree <- drop.tip(tree, tip = rogue_tips)

		if( is.monophyletic(tree, tips_in_order) == TRUE){
				order_tree <- drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% tips_in_order])
				 # Save the pruned tree to a file
  				write.tree(order_tree, paste0(path_out, "twice_pruned_tree_family_", order, "_GBMB.txt"))
				rogue_tips_family <- rbind(rogue_tips_family, data.frame(order = order, rogue_tips = c(rogue_tips)))
		}

	} else if (length(rogue_tips) == 0) { # Trees which have already had the rogue tips pruned
		next
	} else {
		# Here I will loop through the tips in the order and find the largest monophyletic clade which is in the order.
		#largest_clade <- find_largest_clade(tips_in_order, tree)
		#cat("The largest monophyletic clade in the", order, "is: ", length(largest_clade$tip.label), "\n")

		subtree <- get_subtree_at_node(tree, MRCA-Ntip(tree))$subtree
		# Saving the tree with all the descendants of the MRCA
		# Can i somehow 
		write.tree(subtree, paste0(path_out, "MRCA_tree_", order, "_GBMB.txt"))
		}
	}

tree
plot(tree)
