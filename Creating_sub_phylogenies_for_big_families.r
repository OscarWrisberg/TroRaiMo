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
sp_pairs_for_mrca_asteraceae <- list(list("Doniophyton weddellii","Dasyphyllum brevispinum"), #1
								list("Cyclolepis genistoides","Chaptalia tomentosa"), # 2
								list("Cardopatium corymbosum","Echinops ritrodes"), # 3
								list("Tugarinovia mongolica","Carlina sicula"), #  3
								list("Berardia lanuginosa","Saussurea przewalskii"), # 3
								list("Schischkinia_albispina","Centaurea olympica"), # 3 
								list("Platycarphella carlinoides","Berkheya carlinopsis"), # 4
								list("Cacosmia rugosa","Munnozia maronii"), # 4
								list("Moquinia racemosa","Lychnophora granmogolensis"), # 4
								list("Geropogon hybridus","Tragopogon samaritani"), # 5
								list("Erythroseris somalensis","Pilosella fusca"), # 5
								list("Faberia_thibetica","Crepis rhaetica"), # 5 
								list("Abrotanella forsteroides","Abrotanella spathulata"), # 6
								list("Gymnodiscus capillaris","Euryops brevilobus"), # 6
								list("Dolichoglottis lyallii", "Brachyglottis southlandica"), # 6
								list("Capelio caledonica","Tephroseris palustris"), # 6
								list("Hasteola robusta","Senecio laetevirens"), # 6 
								list("Doronicum stenoglossum","Artemisia eriopoda"), # 7
								list("Pachystegia insignis","Olearia_lirata"), # 8
								list("Mairia purpurata","Nardophyllum chiliotrichoides"), # 8
								list("Pteronia glomerata","Felicia cymbalarioides"), # 8
								list("Madagaster mandrarensis","Psiadia pollicina"), # 8
								list("Dichrocephala benthamii","Bellis pappulosa"), # 8
								list("Ceratogyne obionoides","Brachyscome staceae"), # 8
								list("Guynesomia scoparia", "Inulopsis scaposa"), # 8
								list("Erigeron primulifolius","Erigeron caucasicus"), # 8
								list("Tracyina rostrata", "Ericameria fasciculata"), # 8
								list("Gundlachia corymbosa","Solidago ptarmicoides"), # 8
								list("Townsendia formosa","Grindelia fraxinipratensis"), # 8 
								list("Duhaldea rubricaulis","Geigeria nianganensis"), # 9
								list("Callilepis salicifolia","Helichrysum fulvescens"), # 10
								list("Ericentrodea decomposita","Bidens amplectens"), # 11
								list("Sabazia sarmentosa","Aldama flava"), # 12
								list("Polymnia laevigata","Espeletia episcopalis"), #13
								list("Coulterella capitata","Pectis multiceps"), # 14
								list("Neurolaena lobata","Gaillardia amblyodon"),# 15
								list("Orochaenactis thysanocarpha","Hymenothrix loomisii"), # 16
								list("Guardiola platyphylla","Deinandra frutescens"), #17
								list("Hofmeisteria urenifolia","Stevia tomentosa"), # 18
								list("Ageratina glechonophylla","Liatris provincialis") # 19
								) # 

sp_pairs_for_mrca_orchidaceae <- list(list("Codonorchis lessonii","Pterygodium magnum"), # 1
								   list("Chloraea cylindrostachya","Orthoceras strictum"), # 2
								   list("Notylia buchtienii","Psychopsis sanderae"), # 3
								   list("Maxillaria huancabambae","Eriopsis rutidobulbon"), # 4 X
								   list("Stanhopea ecornuta","Braemia vittata"), # 4
								   list("Pescatoria coronaria","Cryptarrhena guatemalensis"), # 4
								   list("Eulophia flava","Cyrtopodium andersonii"), # 5
								   list("Phalaenopsis marriottiana","Bromheadia finlaysoniana"), # 6
								   list("Masdevallia pinocchio", "Lepanthes bradei"), # 7 
								   list("Pabstiella leucopyramis","Stelis montis-mortensis"), # 7
								   list("Scaphosepalum ursinum","Andinia vestigipetala"), # 7
								   list("Acianthera minima","Acianthera atroglossa"), # 7
								   list("Barbosella orbicularis","Restrepiella ovatipetala"), # 7
								   list("Epidendrum oxyglossum","Laelia eyermaniana"), # 8
								   list("Guarianthe skinneri","Meiracyllium gemma"), # 8
								   list("Encyclia bractescens","Oestlundia cyanocolumna"), # 8
								   list("Scaphyglottis longicaulis", "Nidema boothii"), # 8 X
								   list("Bulbophyllum clandestinum","Dendrobium mariae"), # 9
								   list("Liparis anopheles","Liparis rheedei"), # 10
								   list("Eriodes barbata","Ancistrochilus rothschildianus"), # 11
								   list("Callostylis rigida","Phreatia tahitensis"), # 12
								   list("Coelogyne mayeriana","Arundina graminifolia"), # 13
								   list("Sobralia ecuadorana","Nervilia cumberlegii"), # 14
								   list("Epipactis helleborine","Neottia smallii"), # 15
								   list("Phragmipedium warszewiczianum","Selenipedium aequinoctiale") # 16
								   ) # Were losing 15 species

sp_pairs_for_mrca_fabaceae <- list(list("Trifolium macraei","Parochetus communis"), # 1
								list("Astragalus lilacinus","Erophaca baetica"), # 2
								list("Onobrychis hajastana","Alhagi maurorum"), # 2
								list("Caragana stenophylla","Chesneya ferganensis"), # 2
								list("Lotus dumetorum", "Lennea viridiflora"), # 3
								list("Psoralea asarina","Toxicopueraria peduncularis"), # 4 
								list("Desmodium paniculatum","Campylotropis macrocarpa"), # 4
								list("Phaseolus zimapanensis","Otoptera burchellii"), # 4
								list("Derris glabra","Disynstemon paullinioides"), # 5
								list("Indigofera humbertiana","Phylloxylon perrieri"), 	# 6
								list("Crotalaria lebrunii","Amphiodon effusus"), # 7
								list("Dalbergia subcymosa","Salweenia wardii"), # 7
								list("Luetzelburgia purpurea","Vataireopsis surinamensis"), # 7
								list("Ateleia standleyana","Dalbergia sericea"), # 7 
								list("Acacia crassa","Senegalia catechu"), # 8
								list("Inga acrocephala","Mariosousa heterophylla"), # 9
								list("Hoffmannseggia repens","Hererolandia pearsonii"), # 10
								list("Senna skinneri","Pterogyne nitens"),  # 11
								list("Bauhinia grevei","Duparquetia orchidacea"), # 12
								list("Bikinia aciculifera","Barnebydendron riedelii") # 13
								) # 

sp_pairs_for_mrca_caryophyllaceae <- list(list("Schiedea jacobii","Telephium imperati"),
										list("Aphanisma blitoides","Phaulothamnus spinescens")
										) # were losing 0 species 


sp_pairs_for_mrca_lamiaceae <- list(list("Clinopodium betulifolium","Satureja spicigera"),
								list("Thymus mandschuricus","Micromeria croatica"),
								list("Nepeta meyeri","Cedronella canariensis"),
								list("Lepechinia heteromorpha","Lepechinia ganderi"),
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

sp_pairs_for_mrca_poaceae <- list(list("Cynodon transvaalensis","Eragrostiella brachyphylla"), # 1 X
								list("Muhlenbergia huegelii","Cleistogenes songorica"), # 1
								list("Sporobolus pumilus","Zoysia matrella"), # 1
								list("Eragrostis multicaulis","Viguierella madagascariensis"), # 1
								list("Rytidosperma penicillatum","Merxmuellera drakensbergensis"), # 1
								list("Isachne globosa","Monachather paradoxus"), # 1
								list("Trichanthecium rivale","Digitaria albescens"), # 2
								list("Andropogon bicornis","Reynaudia filiformis"), # 2
								list("Loudetia lanata","Bromuniola gossweileri"), # 2 
								list("Aristida ternipes", "Sartidia jucunda"), # 3
								list("Festuca rupicola", "Colpodium biebersteinianum"), # 4 X
								list("Agrostis gracililaxa","Tricholemma jahandiezii"), # 4
								list("Hordeum parodii","Bromus benekenii"), # 4
								list("Austrostipa hemipogon", "Ampelodesmos mauritanicus"), # 5
								list("Dendrocalamus sinicus","Buergersiochloa bambusoides"), # 6
								list("Fargesia yuanjiangensis", "Shibataea chinensis"), # 6 X
								list("Oryza glaberrima", "Streptogyna americana")
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
									 list("Nealchornea yapurensis","Maprounea africana")
									 )


sp_pairs_for_mrca_ericaceae <- list(list("Erica interrupta","Cassiope selaginoides"),
								list("Disterigma codonanthum", "Prionotes cerinthoides"),
								list("Pyrola calliantha","Monotropsis odorata")
								) # were losing 27 species


sp_pairs_for_mrca_apocynaceae <- list(list("Ditassa taxifolia","Eustegia minuta"),
									list("Huernia aspera","Heterostemma_acuminatum"),
									list("Hoya_blashernaezii","Rhyssolobium dumosum"),
								   list("Mandevilla sellowii","Rhabdadenia biflora"),
								   list("Tabernaemontana contorta","Dyera costulata")
								   ) # Were losing 36 species

sp_pairs_for_mrca_brassicaceae <- list(list("Streptanthus insignis","Ricotia lunaria"),
										list("Alyssum gallaecicum","Asperuginoides axillaris"),
										list("Rhammatophyllum gaudanense","Clausia robusta"),
										list("Erysimum bastetanum","Yinshania henryi")
										) # were losing 31 species

sp_pairs_for_mrca_amaryllidaceae <- list(list("Narcissus gaditanus","Ammocharis coranica"),
										 list("Allium cretaceum","Tulbaghia violacea")
										 ) # Were losing 3 species

sp_pairs_for_mrca_asparagaceae <- list(list("Scilla messeniaca", "Oziroe biflora"),
									   list("Anemarrhena asphodeloides","Hosta ventricosa"),
									   list("Lomandra hystrix","Asparagus pseudoscaber")
										) # were losing 17 species

sp_pairs_for_mrca_bromeliaceae <- list(
										list("Aechmea arenaria","Greigia sphacelata"),
										list("Puya compacta","Puya gilmartiniae"),
										list("Dyckia maracasensis","Pitcairnia riparia"),
										list("Tillandsia lechneri","Gregbrownia fulgens"),
										list("Vriesea bituminosa", "Vriesea pleiosticha")
										) 

sp_pairs_for_mrca_cyperaceae <- list(
										list("Carex paeninsulae","Carex ludwigii"),
										list("Scirpus flaccidifolius","Erioscirpus comosus"),
										list("Cyperus luzulae","Cyperus steudneri"),
										list("Isolepis prolifera","Scirpoides mexicana"),
										list("Schoenoplectiella wallichii","Fuirena abnormalis"),
										list("Eleocharis palustris","Arthrostylis aphylla"),
										list("Rhynchospora colorata","Rhynchospora aristata"),
										list("Tetraria secans", "Carpha glomerata")
										) 

sp_pairs_for_mrca_iridaceae <- list(list("Iris fosteriana","Iris unguicularis"),
									list("Moraea demissa","Dietes iridioides"),
									list("Sisyrinchium_fuscatum","Trimezia northiana"),
									list("Gladiolus hirsutus", "Patersonia sericea")
									) # were losing 2 species

sp_pairs_for_mrca_lauraceae <- list(list("Ocotea spixiana", "Persea barbujana"),
									list("Mezilaurus glabriantha","Williamodendron cinnamomeum"),
									list("Beilschmiedia kweichowensis", "Hypodaphnis zenkeri")
									) # were losing 0 species

sp_pairs_for_mrca_melastomataceae <- list(list("Macrocentrum cristatum", "Graffenrieda rotundifolia"),
										list("Eriocnema fulva","Miconia roseopetala"),
										list("Rhynchanthera bracteata","Pachyloma huberioides"),
										list("Acanthella sprucei","Macairea thyrsiflora"),
										list("Tristemma mauritianum","Pterogastra minor"),
										list("Triolena obliqua","Cambessedesia hilariana"),
										list("Sonerila cantonensis","Blakea purpusii"),
										list("Bertolonia mosenii","Salpinga margaritacea")
										)	# were losing 59 species

sp_pairs_for_mrca_myrtaceae <- list(list("Austromyrtus dulcis","Gossia grayi"),
									list("Myrcia insigniflora","Psidium grandifolium"),
									list("Backhousia leptopetala","Syzygium bullatum"),
									list("Arillastrum gummiferum","Stockwellia quadrifida"),
									list("Lindsayomyrtus racemoides", "Neofabricia sericisepala")
									) # were losing 34 species

sp_pairs_for_mrca_plantaginaceae <- list(list("Veronica oltensis","Lafuentea rotundifolia"),
										list("Penstemon mensarum","Russelia verticillata"),
										list("Gratiola floridana","Ourisia microphylla")
										) # were losing 9 species

sp_pairs_for_mrca_ranunculaceae <- list(list("Asteropyrum cavaleriei","Hamadryas magellanica"),
										list("Gymnaconitum gymnandrum","Nigella orientalis"),
										list("Enemion raddeanum","Thalictrum ichangense")
										) # were losing 14 species

sp_pairs_for_mrca_rosaceae <- list(list("Potentilla argaea","Potentilla hebiichigo"),
									list("Alchemilla crinita","Potaninia mongolica"),
									list("Cliffortia juniperina","Hagenia abyssinica"),
									list("Rosa rubiginosa","Rosa bracteata"),
									list("Rubus perrobustus","Rubus biflorus"),
									list("Eriolobus trilobatus","Lyonothamnus floribundus")
									) # were losing 40 species

sp_pairs_for_mrca_acanthaceae_martyniaceae_pedaliaceae <- list(
																list("Ruellia beyrichiana", "Physacanthus_nematosiphon"),
																list("Aphelandra_fasciculata","Kudoacanthus_albonervosus"),
																list("Barleria_lupulina","Neuracanthus umbraticus"),
																list("Pachystachys killipii","Spathacanthus parviflorus"),
																list("Sesamum indicum","Holubia saccata")
																) # were losing 17 species

sp_pairs_for_mrca_anacardiaceae_burseraceae_kirkiaceae <- list(list("Campnosperma gummiferum","Attilaea abalak"),
																list("Beiselia mexicana", "Protium varians")
																) # were losing 2 species

sp_pairs_for_mrca_caryophyllaceae_achatocarpaceae_amaranthaceae <- list(
																list("Schiedea jacobii","Drypis spinosa"),
																list("Thylacospermum caespitosum","Rabelera holostea"),
																list("Acanthophyllum pachystegium","Psammosilene tunicoides"),
																list("Heliosperma tommasinii","Agrostemma githago"),
																list("Paronychia canadensis","Drymaria cordata"),
																list("Aphanisma blitoides","Corispermum squarrosum"),
																list("Halimocnemis villosa","Allenrolfea occidentalis"),
																list("Polycnemum majus","Charpentiera ovata")
																) 

sp_pairs_for_mrca_gesneriaceae_calceolariaceae <- list(list("Primulina pengii","Microchirita tubulosa"),
														list("Streptocarpus erubescens","Streptocarpus strigosus"),
														list("Paraboea paniculata","Tribounia venosa"),
														list("Rhynchotechum parviflorum","Tetraphylloides roseus"),
														list("Monophyllaea glauca","Epithema tenue"),
														list("Nematanthus brasiliensis", "Napeanthus apodemus")
														) # were losing 35 species

sp_pairs_for_mrca_phyllanthaceae_picodendraceae_linaceae_ixonanthaceae <- list(list("Cyrillopsis paraensis","Roucheria monsalveae"),
																				list("Glochidion emarginatum","Glochidion cordatum"), # 2
																				list("Phyllanthus purpusii","Phyllanthus vakinankaratrae"), # 2
																				list("Breynia mollis","Synostemon bacciformis"), # 2 
																				list("Pentabrachion reticulatum","Croizatia brevipetiolata"),
																				list("Spondianthus preussii","Bischofia javanica")
																				) # were losing 42 species

####
sp_pairs_for_mrca_arecaceae <- list(
	list("Eugeissona tristis", "Raphia farinifera"),
	list("Iriartea deltoidea", "Wendlandiella gracilis"),
	list("Ceroxylon parvifrons","Pseudophoenix lediniana"),
	list("Chuniophoenix nana","Coccothrinax argentea")
) # were losing 1 species

sp_pairs_for_ericaceae_clethraceae_cyrillaceae <- list(
	list("Erica interrupta", "Daboecia azorica"),
	list("Rhododendron thomsonii", "Ceratiola ericoides"),
	list("Disterigma codonanthum","Prionotes cerinthoides")
) # were losing 114 species

sp_pairs_for_ochnaceae_clusiaceae_erythroxylaceae_podostemaceae_bonnetiaceae_rhizophoraceae_calophyllaceae_hypericaceae_ctenolophonaceae_irvingiaceae_pandaceae <- list(
	list("Medusagyne oppositifolia", "Testulea gabonensis"),
	list("Irvingia gabonensis","Carallia brachiata"),
	list("Terniopsis brevis","Hypericum fraseri"),
	list("Garcinia verrucosa","Dystovomita paniculata"),
	list("Bonnetia stricta","Endodesmia_calophylloides")
)

sp_pairs_for_orobanchaceae_phrymaceae_mazaceae_paulowniaceae <- list(
	list("Pedicularis lutescens","Leptorhabdos parviflora"),
	list("Castilleja lutescens","Lamourouxia rhinanthifolia"),
	list("Euphrasia tricuspidata","Pterygiella bartschioides"),
	list("Harveya bolusii","Orobanche boninsimae"),
	list("Erythranthe floribunda","Wightia speciosissima")
)

sp_pairs_for_papaveraceae <- list(
	list("Pteridophyllum racemosum","Dicentra uniflora"),
	list("Eschscholzia palmeri","Arctomecon humilis")
)

sp_pairs_for_plumbaginaceae_polygonaceae_frankeniaceae_tamaricaceae <- list(
	list("Podopterus mexicanus","Afrobrunnichia erecta"),
	list("Pteropyrum aucheri","Fagopyrum tataricum"),
	list("Psylliostachys suworowii","Aegialitis annulata"),
	list("Frankenia corymbosa","Tamarix usneoides")
)

sp_pairs_for_scrophulariaceae <- list(
	list("Selago aspera","Aptosimum indivisum"),
) # losing 1 species.

sp_pairs_for_zingiberaceae <- list(
	list("Epiamomum roseisquamosum","Tamijia flagellaris"),
	list("Hemiorchis rhodorrhachis","Globba substrigosa")
)

# Create a list of the sp_pairs
sp_pairs <- list(sp_pairs_for_mrca_asteraceae, sp_pairs_for_mrca_orchidaceae, sp_pairs_for_mrca_fabaceae, sp_pairs_for_mrca_caryophyllaceae, sp_pairs_for_mrca_lamiaceae,
				 sp_pairs_for_mrca_rubiaceae, sp_pairs_for_mrca_poaceae, Sp_pairs_for_mrca_apiaceae, sp_pairs_for_mrca_euphorbiaceae, sp_pairs_for_mrca_ericaceae,
				 sp_pairs_for_mrca_apocynaceae,sp_pairs_for_mrca_brassicaceae, sp_pairs_for_mrca_amaryllidaceae, sp_pairs_for_mrca_asparagaceae, sp_pairs_for_mrca_bromeliaceae,
				 sp_pairs_for_mrca_cyperaceae, sp_pairs_for_mrca_iridaceae, sp_pairs_for_mrca_lauraceae, sp_pairs_for_mrca_melastomataceae, sp_pairs_for_mrca_myrtaceae,
				 sp_pairs_for_mrca_plantaginaceae, sp_pairs_for_mrca_ranunculaceae, sp_pairs_for_mrca_rosaceae, sp_pairs_for_mrca_acanthaceae_martyniaceae_pedaliaceae,
				 sp_pairs_for_mrca_anacardiaceae_burseraceae_kirkiaceae, sp_pairs_for_mrca_caryophyllaceae_achatocarpaceae_amaranthaceae, sp_pairs_for_mrca_gesneriaceae_calceolariaceae,
				 sp_pairs_for_mrca_phyllanthaceae_picodendraceae_linaceae_ixonanthaceae, sp_pairs_for_mrca_arecaceae, sp_pairs_for_ericaceae_clethraceae_cyrillaceae,
				 sp_pairs_for_ochnaceae_clusiaceae_erythroxylaceae_podostemaceae_bonnetiaceae_rhizophoraceae_calophyllaceae_hypericaceae_ctenolophonaceae_irvingiaceae_pandaceae,
				 sp_pairs_for_orobanchaceae_phrymaceae_mazaceae_paulowniaceae, sp_pairs_for_papaveraceae, sp_pairs_for_plumbaginaceae_polygonaceae_frankeniaceae_tamaricaceae,
				 sp_pairs_for_scrophulariaceae, sp_pairs_for_zingiberaceae)


# Family list
fam_list <- list("Asteraceae","Orchidaceae","Fabaceae","Caryophyllaceae","Lamiaceae","Rubiaceae","Poaceae","Apiaceae","Euphorbiaceae","Ericaceae","Apocynaceae",
 				"Brassicaceae", "Amaryllidaceae", "Asparagaceae", "Bromeliaceae", "Cyperaceae", "Iridaceae", "Lauraceae", "Melastomataceae", "Myrtaceae", "Plantaginaceae",
				"Ranunculaceae", "Rosaceae", "Acanthaceae_Martyniaceae_Pedaliaceae", "Anacardiaceae_Burseraceae_Kirkiaceae", "Caryophyllaceae_Achatocarpaceae_Amaranthaceae",
				"Gesneriaceae_Calceolariaceae", "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae", "Arecaceae", "Ericaceae_Clethraceae_Cyrillaceae",
				"Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae",
				"Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae", "Papaveraceae", "Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae", "Scrophulariaceae", "Zingiberaceae")

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