# Load Packages
library(ape)
library(phytools)
library(data.table)

# Loading the command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_in <- as.character(args[1]) # here you define the name of your input file
wcvp_file <- as.character(args[2]) # here you define the name of your output file
output_folder <- as.character(args[3]) # here you define the name of your output file

# Set the working directory
setwd(path_in)

# Load the wcvp file
wcvp <- fread(wcvp_file)

# I want to find out how many accepted species there are in each family and genus according to the wcvp file
# Create a subset of wcvp where taxon_status is "accepted" and taxon_rank is "species"
accepted_species <- subset(wcvp, taxon_status == "Accepted" & taxon_rank == "Species")

# Save the subset phylogenies to a folder.
folder <- output_folder

# List of orders
orders <- c("Zingiberales", "Laurales","Poales","Ranunculales","Rosales","Sapindales","Saxifragales","Myrtales","Malvales","Malpighiales","Lamiales","Gentianales","Fabales",
			"Ericales", "Apiales","Asterales","Asparagales","Caryophyllales","Arecales","Brassicales")

# List of families to join or run on their own
Zingiberales_trees <- list("Zingiberaceae", list("Marantaceae", "Cannaceae"), "Costaceae", list("Heliconiaceae", "Lowiaceae", "Strelitziaceae"))
Laurales_trees <- list("Lauraceae","Monimiaceae")
Poales_trees <- list("Poaceae","Cyperaceae","Bromeliaceae", "Restionaceae", list("Xyridaceae", "Eriocaulaceae"),"Juncaceae","Typhaceae")
Ranunculales_trees <- list("Menispermaceae","Berberidaceae","Ranunculaceae","Papaveraceae")
Rosales_trees <- list("Rosaceae","Urticaceae",list("Rhamnaceae","Barbeyaceae","Dirachmaceae","Elaeagnaceae"),"Moraceae","Ulmaceae","Cannabaceae")
Sapindales_trees <-list(list("Anacardiaceae","Burseraceae","Kirkiaceae"),"Sapindaceae","Rutaceae","Meliaceae","Simaroubaceae")
Saxifragales_trees <-list(list("Crassulaceae","Aphanopetalaceae","Halograceae","Penthoraceae","Tetracarpaeaceae"),list("Saxifragaceae","Iteaceae","Grossulariaceae"),list("Cercidiphyllaceae","Hamamelidaceae","Daphniphyllaceae","Altingiaceae","Paeoniaceae"))
Myrtales_trees <-list("Melastomataceae","Myrtaceae",list("Lythraceae","Onagraceae"),list("Alzateaceae","Crypteroniaceae","Penaeaceae"),"Combretaceae")
Malvales_trees <-list("Malvaceae","Thymelaeaceae",list("Dipterocarpaceae","Bixaceae","Cistaceae","Sarcoleanaceae","Muntingiaceae","Sphaerosepalaceae"))
Malpighiales_trees <- list(list("Salicaceae","Lacistemataceae"),"Euphorbiaceae",list("Chrysobalanaceae","Malpighiaceae","Caryocaraceae","Balanopaceae","Elatinaceae","Centroplacaceae","Dichapetalaceae","Putranjivaceae","Euphroniaceae","Lophopyxidaceae","Trigoniaceae"),list("Phyllanthaceae","Picodendraceae","Linaceae","Ixonanthaceae"),list("Ochnaceae","Clusiaceae","Erythroxylaceae","Podostemaceae","Bonnetiaceae","Rhizophoraceae","Calophyllaceae","Hypericaceae","Ctenolophonaceae","Irvingiaceae","Pandaceae"),"Passifloraceae",list("Violaceae","Goupiaceae"))
Lamiales_trees <- list(list("Verbenaceae","Schlegeliaceae","Lentibulariaceae","Thomandersiaceae"),"Lamiaceae",list("Acanthaceae","Martyniaceae","Pedaliaceae"),list("Gesneriaceae","Calceolariaceae"),"Bignoniaceae",list("Orobanchaceae","Phrymaceae","Mazaceae","Paulowniaceae"),"Scrophulariaceae","Plantaginaceae")
Gentiales_trees <- list("Rubiaceae","Apocynaceae",list("Loganiaceae","Gelsemiaceae"),"Gentianaceae")
Fabales_trees <- list("Fabaceae",list("Polygalaceae","Surianaceae"))
Ericales_trees <- list("Sapotaceae",list("Polemoniaceae","Lecythidaceae","Fouquieriaceae"),list("Ericaceae","Clethraceae","Cyrillaceae"),list("Pentaphylacaceae","Sladeniaceae"),"Primulaceae",list("Styracaceae","Diapensiaceae","Symplocaceae"),"Theaceae","Ebenaceae",list("Balsaminaceae","Marcgraviaceae","Tetrameristaceae"))
Apiales_trees <- list("Apiaceae","Araliaceae","Pittosporaceae")
Asterales_trees <- list("Asteraceae","Calyceraceae",list("Campanulaceae","Rousseaceae"),"Goodeniaceae","Menyanthaceae")
Asparagales_trees <- list("Asphodelaceae","Orchidaceae","Amaryllidaceae","Iridaceae","Asparagaceae")
Caryophyllales_trees <- list(list("Cactaceae","Molluginaceae","Didiereaceae","Anacompserotaceae","Basellaceae","Montiaceae","Halophytaceae","Portulacaceae","Talinaceae"),list("Plumbaginaceae","Polygonaceae","Frankeniaceae","Tamaricaceae"),list("Caryophyllaceae","Achatocarpaceae","Amaranthaceae"),list("Aizoaceae","Phytolaccaceae","Barbeuiaceae","Lophiocarpaceae","Gisekiaceae","Sarcobataceae"),list("Droseraceae","Ancistrocladaceae","Drosophyllaceae","Nepenthaceae","Dioncophyllaceae"))
Arecales_trees <- list("Arecaceae")
Brassicales_trees <- list("Brassicaceae","Resedaceae","Capparaceae","Cleomaceae")

# Combine all these lists into one list
Tree_list <- list(Zingiberales_trees, Laurales_trees, Poales_trees, Ranunculales_trees, Rosales_trees, Sapindales_trees, Saxifragales_trees, Myrtales_trees, Malvales_trees,
				  Malpighiales_trees, Lamiales_trees, Gentiales_trees, Fabales_trees, Ericales_trees, Apiales_trees, Asterales_trees, Asparagales_trees, Caryophyllales_trees,
				  Arecales_trees, Brassicales_trees)

# Loop through all the lists in Tree_list and make sub phylogenies for each of them. 
for (k in seq_along(Tree_list)) {
	family_trees <- Tree_list[[k]]

	family_phylos <- list()

	# Load the tree
	phy_tree <- read.tree(paste0("pruned_tree_order_",orders[k],"_GBMB.tre"))

	# Change _ to " " in the tip labels
	phy_tree$tip.label <- gsub("_", " ", phy_tree$tip.label)
	length(phy_tree$tip.label)

	# Convert phy_tree$tip.label to a data frame with column name "taxon_name"
	phy_df <- data.frame(taxon_name = phy_tree$tip.label)

	# Merge with wcvp
	matched <- merge(phy_df, wcvp, by = "taxon_name", all.x = TRUE)
	matched <- subset(matched, select = c("taxon_name", "family", "genus"))

	# Remove duplicate rows
	matched <- unique(matched)

	# Create a subset of the phylogeny for each family
	for (i in seq_along(family_trees)) {
		if (is.list(family_trees[[i]])) {
			family_subset <- matched[matched$family %in% family_trees[[i]], ]

			for (j in seq_along(family_trees[[i]])) {
				cat("The number of species in ", family_trees[[i]][[j]], " is ", length(family_subset[which(family_subset$family == family_trees[[i]][[j]]),1]), "\n")
			}

		} else {
			family_subset <- matched[matched$family == family_trees[[i]], ]
			cat("The number of species in ", family_trees[[i]], " is ", length(family_subset[which(family_subset$family == family_trees[[i]]),1]), "\n")
		}
		tips_family <- as.character(family_subset$taxon_name)

		# Prune the tree so it only contains the tips in the family
		family_phylo <- keep.tip(phy_tree, tip = tips_family)

		# Check if the number of tips in family_phylo is the same as the number of species in family_subset
		if (length(family_phylo$tip.label) != nrow(family_subset)) {
			cat("Error: The number of tips in family_phylo is not the same as the number of species in family_subset.\n")
			cat("Number of tips in family_phylo: ", length(family_phylo$tip.label), "\n")
			cat("Number of species in family_subset: ", nrow(family_subset), "\n")
			cat("Rows in family_subset not in family_phylo:\n")
			missing_rows <- family_subset[!(family_subset$taxon_name %in% family_phylo$tip.label), ]
			print(missing_rows)
			break
		}

		# Append the pruned phylogeny to the list
		family_phylos[[i]] <- family_phylo
	}

	# Save the subset phylogenies to a folder.
	for (l in seq_along(family_phylos)) {
		if (is.list(family_trees[[l]])) {
			family_names <- unlist(family_trees[[l]])
			tree_name <- paste0("family_phylo_", paste(family_names, collapse = "_"), ".tre")
		} else {
			tree_name <- paste0("family_phylo_", family_trees[[l]], ".tre")
		}
		
		write.tree(family_phylos[[l]], file = paste0(folder, tree_name))
	}
}

