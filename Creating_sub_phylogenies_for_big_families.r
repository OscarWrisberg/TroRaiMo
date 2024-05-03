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
sp_pairs_for_mrca_asteraceae <- list(list("Doniophyton weddellii","Dasyphyllum brevispinum"),
								list("Cyclolepis genistoides","Chaptalia tomentosa"),
								list("Macledium zeyheri","Centaurea olympica"),
								list("Platycarphella carlinoides","Lychnophora granmogolensis"),
								list("Gundelia tournefortii","Crepis rhaetica"),
								list("Abrotanella forsteroides","Senecio laetevirens"),
								list("Doronicum stenoglossum","Artemisia eriopoda"),
								list("Cyathocline purpurea","Grindelia fraxinipratensis"),
								list("Duhaldea rubricaulis","Geigeria nianganensis"),
								list("Callilepis salicifolia","Helichrysum fulvescens"),
								list("Ericentrodea decomposita","Bidens amplectens"),
								list("Sabazia sarmentosa","Aldama flava"),
								list("Polymnia laevigata","Liatris provincialis")
								) # were losing 1 species

sp_pairs_for_mrca_orchidaceae <- list(list("Codonorchis lessonii","Pterygodium magnum"),
								   list("Chloraea cylindrostachya","Orthoceras strictum"),
								   list("Notylia buchtienii","Psychopsis sanderae"),
								   list("Maxillaria huancabambae","Cryptarrhena guatemalensis"),
								   list("Eulophia flava","Cyrtopodium andersonii"),
								   list("Phalaenopsis marriottiana","Bromheadia finlaysoniana"),
								   list("Masdevallia pinocchio", "Neocogniauxia hexaptera"),
								   list("Encyclia candollei", "Arpophyllum giganteum"),
								   list("Bulbophyllum clandestinum","Dendrobium mariae"),
								   list("Liparis_anopheles","Liparis rheedei"),
								   list("Eriodes barbata","Ancistrochilus rothschildianus"),
								   list("Callostylis rigida","Phreatia tahitensis"),
								   list("Coelogyne mayeriana","Arundina graminifolia"),
								   list("Sobralia ecuadorana","Nervilia cumberlegii"),
								   list("Epipactis helleborine","Neottia smallii"),
								   list("Phragmipedium warszewiczianum","Selenipedium aequinoctiale")
								   ) # Were losing 15 species

sp_pairs_for_mrca_fabaceae <- list(list("Trifolium macraei","Parochetus communis"),
								list("Astragalus lilacinus","Chesneya ferganensis"),
								list("Lotus dumetorum", "Lennea viridiflora"),
								list("Psoralea asarina","Platycyamus regnellii"),
								list("Derris glabra","Disynstemon_paullinioides"),
								list("Indigofera humbertiana","Phylloxylon perrieri"),
								list("Crotalaria lebrunii","Dalbergia sericea"),
								list("Acacia crassa","Senegalia catechu"),
								list("Inga acrocephala","Mariosousa heterophylla"),
								list("Hoffmannseggia repens","Hererolandia pearsonii"),
								list("Senna skinneri","Pterogyne nitens"), 
								list("Bauhinia grevei","Duparquetia orchidacea"),
								list("Bikinia aciculifera","Barnebydendron riedelii")
								) # were losing ca 200 species

sp_pairs_for_mrca_caryophyllaceae <- list(list("Schiedea jacobii","Telephium imperati"),
										list("Aphanisma blitoides","Phaulothamnus spinescens")
										) # were losing 0 species 


sp_pairs_for_mrca_lamiaceae <- list(list("Clinopodium betulifolium","Lepechinia ganderi"),
								list("Salvia purpurea","Salvia rosmarinus"),
								list("Eriope foetida","Hanceola exserta"),
								list("Syncolostemon macranthus","Isodon serra"),
								list("Teucrium odontites","Caryopteris incana"),
								list("Lamium multifidum","Scutellaria amoena"),
								list("Prostanthera incisa","Callicarpa longifolia")
								) # were losing 42 species

sp_pairs_for_mrca_rubiaceae <- list(list("Galium hirtum","Luculia grandifolia"),
								 list("Tarenna seemanniana","Hekistocarpa minutiflora"),
								 list("Uncaria hirsuta","Sipaneopsis rupicola")
								 ) # were losing 0 species

sp_pairs_for_mrca_poaceae <- list(list("Cynodon transvaalensis","Monachather paradoxus"),
								list("Trichanthecium rivale", "Bromuniola gossweileri"),
								list("Aristida ternipes", "Sartidia jucunda"),
								list("Festuca rupicola", "Bromus benekenii"),
								list("Austrostipa hemipogon", "Ampelodesmos mauritanicus"),
								list("Dendrocalamus sinicus", "Shibataea chinensis"),
								list("Oryza glaberrima", "Streptogyna americana"),
								) # were losing 60 species

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
								) # were losing 27 species


sp_pairs_for_mrca_apocynaceae <- list(list("Ditassa taxifolia","Eustegia minuta"),
									list("Huernia aspera","Rhyssolobium dumosum"),
								   list("Mandevilla sellowii","Rhabdadenia biflora"),
								   list("Tabernaemontana contorta","Dyera costulata")
								   ) # Were losing 36 species

sp_pairs_for_mrca_brassicaceae <- list(list("Streptanthus insignis","Ricotia lunaria"),
										list("Alyssum gallaecicum","Asperuginoides axillaris"),
										list("Rhammatophyllum gaudanense","Clausia robusta"),
										list("Erysimum bastetanum","Yinshania henryi")
										) # were losing 31 species

sp_pairs_for_mrca_amaryllidaceae <- list(list("Narcissus gaditanus","Agapanthus campanulatus"),
										 list("Allium cretaceum","Tulbaghia violacea")
										 ) # Were losing 0 species

sp_pairs_for_mrca_asparagaceae <- list(list("Scilla messeniaca", "Oziroe biflora"),
									   list("Anemarrhena asphodeloides","Hosta ventricosa"),
									   list("Lomandra hystrix","Asparagus pseudoscaber")
										) # were losing 17 species

sp_pairs_for_mrca_bromeliaceae <- list(list("Aechmea arenaria","Bakerantha caerulea"),
										list("Tillandsia lechneri", "Glomeropitcairnia erectiflora")
										) # were losing 10 species

sp_pairs_for_mrca_cyperaceae <- list(list("Carex paeninsulae","Blysmus rufus"),
										list("Cyperus luzulae", "Bulbostylis barbata"),
										list("Tetraria secans", "Carpha glomerata")
										) # were losing 36 species

sp_pairs_for_mrca_iridaceae <- list(list("Iris fosteriana", "Diplarrena latifolia"),
									list("Gladiolus hirsutus", "Patersonia sericea")
									) # were losing 1 species

sp_pairs_for_mrca_lauraceae <- list(list("Caryodaphnopsis laotica", "Persea barbujana"),
									list("Beilschmiedia kweichowensis", "Hypodaphnis zenkeri")
									) # were losing 0 species

sp_pairs_for_mrca_melastomataceae <- list(list("Macrocentrum cristatum", "Miconia roseopetala"),
										list("Rhynchanthera bracteata","Salpinga_margaritacea")
										)	# were losing 59 species

sp_pairs_for_mrca_myrtaceae <- list(list("Austromyrtus dulcis","Syzygium bullatum"),
									list("Arillastrum gummiferum","Stockwellia quadrifida"),
									list("Lindsayomyrtus racemoides", "Neofabricia sericisepala")
									) # were losing 34 species

sp_pairs_for_mrca_plantaginaceae <- list(list("Veronica oltensis","Lafuentea rotundifolia"),
										list("Penstemon mensarum","Russelia verticillata"),
										list("Gratiola_floridana","Ourisia_microphylla")
										) # were losing 9 species

sp_pairs_for_mrca_ranunculaceae <- list(list("Asteropyrum cavaleriei","Hamadryas magellanica"),
										list("Gymnaconitum gymnandrum","Nigella orientalis"),
										list("Enemion raddeanum","Thalictrum ichangense")
										) # were losing 14 species

sp_pairs_for_mrca_rosaceae <- list(list("Potentilla argaea","Potaninia mongolica"),
									list("Cliffortia juniperina","Hagenia abyssinica"),
									list("Rosa rubiginosa","Rosa bracteata"),
									list("Rubus perrobustus","Rubus biflorus"),
									list("Eriolobus trilobatus","Lyonothamnus floribundus")
									) # were losing 40 species

sp_pairs_for_mrca_acanthaceae_martyniaceae_pedaliaceae <- list(list("Ruellia beyrichiana","Neuracanthus_umbraticus"),
																list("Pachystachys killipii","Spathacanthus_parviflorus"),
																list("Sesamum indicum","Holubia saccata")
																) # were losing 17 species

sp_pairs_for_mrca_anacardiaceae_burseraceae_kirkiaceae <- list(list("Campnosperma gummiferum","Attilaea abalak"),
																list("Beiselia mexicana", "Protium varians")
																) # were losing 2 species

sp_pairs_for_mrca_caryophyllaceae_achatocarpaceae_amaranthaceae <- list(list("Schiedea jacobii","Telephium imperati"),
																list("Aphanisma blitoides", "Phaulothamnus spinescens"),
																list("Pentabrachion reticulatum","Croizatia brevipetiolata")
																) # were losing 0 species

sp_pairs_for_mrca_gesneriaceae_calceolariaceae <- list(list("Primulina pengii","Epithema tenue"),
														list("Nematanthus brasiliensis", "Napeanthus apodemus")
														) # were losing 31 species

sp_pairs_for_mrca_phyllanthaceae_picodendraceae_linaceae_ixonanthaceae <- list(list("Cyrillopsis paraensis","Roucheria monsalveae"),
																				list("Notoleptopus decaisnei","Glochidion eucleoides")
																				) # were losing 42 species


# Create a list of the sp_pairs
sp_pairs <- list(sp_pairs_for_mrca_asteraceae, sp_pairs_for_mrca_orchidaceae, sp_pairs_for_mrca_fabaceae, sp_pairs_for_mrca_caryophyllaceae, sp_pairs_for_mrca_lamiaceae,
				 sp_pairs_for_mrca_rubiaceae, sp_pairs_for_mrca_poaceae, Sp_pairs_for_mrca_apiaceae, sp_pairs_for_mrca_euphorbiaceae, sp_pairs_for_mrca_ericaceae,
				 sp_pairs_for_mrca_apocynaceae,sp_pairs_for_mrca_brassicaceae, sp_pairs_for_mrca_amaryllidaceae, sp_pairs_for_mrca_asparagaceae, sp_pairs_for_mrca_bromeliaceae,
				 sp_pairs_for_mrca_cyperaceae, sp_pairs_for_mrca_iridaceae, sp_pairs_for_mrca_lauraceae, sp_pairs_for_mrca_melastomataceae, sp_pairs_for_mrca_myrtaceae,
				 sp_pairs_for_mrca_plantaginaceae, sp_pairs_for_mrca_ranunculaceae, sp_pairs_for_mrca_rosaceae, sp_pairs_for_mrca_acanthaceae_martyniaceae_pedaliaceae,
				 sp_pairs_for_mrca_anacardiaceae_burseraceae_kirkiaceae, sp_pairs_for_mrca_caryophyllaceae_achatocarpaceae_amaranthaceae, sp_pairs_for_mrca_gesneriaceae_calceolariaceae,
				 sp_pairs_for_mrca_phyllanthaceae_picodendraceae_linaceae_ixonanthaceae)


# Family list
fam_list <- list("Asteraceae","Orchidaceae","Fabaceae","Caryophyllaceae","Lamiaceae","Rubiaceae","Poaceae","Apiaceae","Euphorbiaceae","Ericaceae","Apocynaceae",
 				"Brassicaceae", "Amaryllidaceae", "Asparagaceae", "Bromeliaceae", "Cyperaceae", "Iridaceae", "Lauraceae", "Melastomataceae", "Myrtaceae", "Plantaginaceae",
				"Ranunculaceae", "Rosaceae", "Acanthaceae_Martyniaceae_Pedaliaceae", "Anacardiaceae_Burseraceae_Kirkiaceae", "Caryophyllaceae_Achatocarpaceae_Amaranthaceae",
				"Gesneriaceae_Calceolariaceae", "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae")

# Now loop through all the sp_pairs and make sub phylogenies for each of them.
for (k in seq_along(sp_pairs)){ # loop through the list of sp_pairs
  list <- sp_pairs[[k]] # Select the list
  file_name <- paste0("family_phylo_",fam_list[k],".tre") # Define the file name
  if (file.exists(file_name)) { # Check if the file exists
    print(file_name)
    phy_tree <- read.tree(file_name) # Reading the tree
    phy_tree$tip.label <- gsub("_", " ", phy_tree$tip.label) # subbing the underscores for spaces

    for (i in seq_along(list)) { # looping through the lists within the list
      cat("Is Sp 1: ",list[[i]][[1]]," in the phylogeny?:",list[[i]][[1]] %in% phy_tree$tip.label,"\n") # checking if sp one is in the phylogeny
      cat("Is Sp 2:",list[[i]][[2]]," in the phylogeny?: ",list[[i]][[2]] %in% phy_tree$tip.label,"\n") # checking if sp two is in the phylogeny
      
      MRCA <- getMRCA(phy_tree, c(list[[i]][[1]],list[[i]][[2]])) # Finding the MRCA

	  if (is.null(MRCA)) { # Checking if there is an MRCA
		cat("The MRCA is null\n") 
	  } else {
		cat("The MRCA is: ", MRCA, "\n") # Printing the MRCA
		sub_phylo <- extract.clade(phy_tree, MRCA) # Extracting the clade
		tree_name <- paste0("sub_phylo_",fam_list[k],"_",i,".tre") # Defining the tree name
		write.tree(sub_phylo, file = paste0(output_folder, tree_name)) # Writing the subtree to a file.
	  }
    }
  } else {
    cat("File not found:", file_name, "\n")
  }
}