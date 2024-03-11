# Setting Cran mirror
chooseCRANmirror(ind = 30)

#Packages
packages <- c("ape","phytools","data.table")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Loading the command line arguments
args <- commandArgs(trailingOnly = TRUE)
path_in <- as.character(args[1]) # here you define the name of your input file
output_folder <- as.character(args[2]) # here you define the name of your output file

# Set the working directory
setwd(path_in)

# Creating subphylogenies of some of the really big families.
# Asteraceae 
sp_pairs_for_mrca_asteraceae <- list(list("Duhaldea rubricaulis","Centaurea amblensis"),
								list("Abrotanella forsteroides","Grindelia oxylepis"),
								list("Crepis setosa","Platycarphella carlinoides"),
								list("Centaurea olympica","Macledium zeyheri"),
								list("Cyclolepis genistoides","Chaptalia tomentosa"),
								list("Dasyphyllum brevispinum","Doniophyton weddellii"),
								list("Pertya robusta","Ainsliaea pingbianensis")
								) # were losing 1 species

sp_pairs_for_mrca_orchidaceae <- list(list("Codonorchis lessonii","Orthoceras strictum"),
								   list("Notylia buchtienii","Thecopus maingayi"),
								   list("Phalaenopsis marriottiana","Bromheadia finlaysoniana"),
								   list("Masdevallia pinocchio", "Wullschlaegelia aphylla"),
								   list("Bulbophyllum clandestinum","Liparis rheedei"),
								   list("Eriodes barbata","Ancistrochilus rothschildianus"),
								   list("Callostylis rigida","Phreatia tahitensis"),
								   list("Coelogyne mayeriana","Arundina graminifolia"),
								   list("Sobralia ecuadorana","Nervilia cumberlegii"),
								   list("Epipactis helleborine","Neottia smallii"),
								   list("Phragmipedium warszewiczianum","Selenipedium aequinoctiale")
								   ) # Were losing 15 species

sp_pairs_for_mrca_fabaceae <- list(list("Trifolium macraei","Lennea viridiflora"),
								list("Psoralea asarina","Phylloxylon perrieri"),
								list("Crotalaria lebrunii","Dalbergia sericea"),
								list("Acacia crassa","Ceratonia siliqua"),
								list("Bauhinia grevei","Duparquetia orchidacea"),
								list("Bikinia aciculifera","Barnebydendron riedelii")
								) # were losing 72 species

sp_pairs_for_mrca_caryophyllaceae <- list(list("Schiedea jacobii","Telephium imperati"),
										list("Aphanisma blitoides","Phaulothamnus spinescens")
										) # were losing 0 species 


sp_pairs_for_mrca_lamiaceae <- list(list("Clinopodium betulifolium","Salvia rosmarinus"),
								list("Eriope foetida","Perillula reptans"),
								list("Teucrium odontites","Congea tomentosa")
								) # were losing 0 species

sp_pairs_for_mrca_rubiaceae <- list(list("Galium hirtum","Luculia grandifolia"),
								 list("Tarenna seemanniana","Sipaneopsis rupicola")
								 ) # were losing 0 species

sp_pairs_for_mrca_poaceae <- list(list("Cynodon transvaalensis","Sartidia jucunda"),
								list("Festuca rupicola", "Streptogyna americana")) # were losing 10 species

Sp_pairs_for_mrca_apiaceae <- list(list("Rhodosciadium pringlei","Arcuatopterus thalictrioideus"),
								list("Ferula dissecta","Conioselinum chinense"),
								list("Lilaeopsis carolinensis","Perideridia kelloggii"),
								list("Bupleurum gracilipes","Pleurospermopsis sikkimensis"),
								list("Eryngium smithii","Steganotaenia araliacea"),
								list("Azorella valentini","Klotzschia glaziovii")
								) # Were losing 149 species.

sp_pairs_for_mrca_euphorbiaceae <- list(list("Neoscortechinia kingii","Argythamnia candicans"),
									 list("Cladogelonium madagascariense","Suregada boiviniana"),
									 list("Nealchornea yapurensis","Maprounea africana"))


sp_pairs_for_mrca_ericaceae <- list(list("Erica interrupta","Cassiope selaginoides"),
								list("Disterigma codonanthum", "Prionotes cerinthoides"),
								list("Pyrola calliantha","Monotropsis odorata")
								) # were losing 27 tips


sp_pairs_for_mrca_apocynaceae <- list(list("Ditassa taxifolia","Periploca laevigata"),
								   list("Mandevilla sellowii","Rhabdadenia biflora"),
								   list("Tabernaemontana contorta","Dyera costulata")
								   ) # Were losing 0 species

# sp_pairs_for_mrca_melastomataceae <- c(c("Macrocentrum cristatum","Miconia roseopetala"),
# 										c("Rhynchanthera bracteata","Salpinga margaritacea")) were losing 60 sp


# Create a list of the sp_pairs
sp_pairs <- list(sp_pairs_for_mrca_asteraceae, sp_pairs_for_mrca_orchidaceae, sp_pairs_for_mrca_fabaceae, sp_pairs_for_mrca_caryophyllaceae, sp_pairs_for_mrca_lamiaceae,
				 sp_pairs_for_mrca_rubiaceae, sp_pairs_for_mrca_poaceae, Sp_pairs_for_mrca_apiaceae, sp_pairs_for_mrca_euphorbiaceae, sp_pairs_for_mrca_ericaceae,
				 sp_pairs_for_mrca_apocynaceae)


# Asteraceae
fam_list <- list("Asteraceae","Orchidaceae","Fabaceae","Caryophyllaceae","Lamiaceae","Rubiaceae","Poaceae","Apiaceae","Euphorbiaceae","Ericaceae","Apocynaceae")

# Now loop through all the sp_pairs and make sub phylogenies for each of them.
for (k in seq_along(sp_pairs)){
  list <- sp_pairs[[k]]
  file_name <- paste0("family_phylo_",fam_list[k],".tre")
  if (file.exists(file_name)) {
    print(file_name)
    phy_tree <- read.tree(file_name)
    phy_tree$tip.label <- gsub("_", " ", phy_tree$tip.label)
    print(phy_tree)
    for (i in seq_along(list)) {
      cat("Is Sp 1: ",list[[i]][[1]]," in the phylogeny?:",list[[i]][[1]] %in% phy_tree$tip.label,"\n")
      cat("Is Sp 2:",list[[i]][[2]]," in the phylogeny?: ",list[[i]][[2]] %in% phy_tree$tip.label,"\n")
      
      MRCA <- getMRCA(phy_tree, c(list[[i]][[1]],list[[i]][[2]]))

	  if (is.null(MRCA)) {
		cat("The MRCA is null\n")
	  } else {
		cat("The MRCA is: ", MRCA, "\n")
		sub_phylo <- extract.clade(phy_tree, MRCA)
		tree_name <- paste0("sub_phylo_",fam_list[k],"_",i,".tre")
		write.tree(sub_phylo, file = paste0(output_folder, tree_name))
	  }
    }
  } else {
    cat("File not found:", file_name, "\n")
  }
}