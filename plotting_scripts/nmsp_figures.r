# Empirical figures for marginal ancestral probabilities.
# 
#
# Ignacio Quintero
#  t(-_-t)
# 
# 30 11 20
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-

source('./node_probabilities.r')
source('plotting_scripts/node_probabilities.r')
data_dir <- "/home/au543206/GenomeDK/Trf_models/workflow/04_results/Esse_output/Esse_done_results/"

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# prepare data
sf = list.files(data_dir, 
    full.names = TRUE, recursive = TRUE, pattern = '*.csv') # Find the file paths of all files with the '.log' extension
	

# sf = grep('h2', sf, value = TRUE)

# This for loop iterates over the elements in the 'sf' vector.
# Each element represents a file path.
for (i in 1:length(sf)) {

  sfi = sf[i]  # Assign the current file path to 'sfi'

  spl = strsplit(sfi, '/')[[1]]  # Split the file path by '/' and store the result in 'spl'

  cl  = spl[7]  # Extract the 7th element from 'spl' and assign it to 'cl'
  cli = strsplit(cl,'')[[1]][1]  # Split 'cl' into individual characters and assign the first character to 'cli'
  k   = spl[8]  # Extract the 8th element from 'spl' and assign it to 'k'
  h   = spl[9]  # Extract the 9th element from 'spl' and assign it to 'h'
  c   = spl[10]  # Extract the 10th element from 'spl' and assign it to 'c'
  v   = spl[11]  # Extract the 11th element from 'spl' and assign it to 'v'
  tn  = strsplit(spl[12], '_|\\.')[[1]][2]  # Split the 12th element from 'spl' by '_' or '.' and assign the second element to 'tn'

  # Check if the file path contains '_GO'
  # If it does, assign TRUE to 'igo', otherwise assign FALSE
  if (length(grep('_GO', sfi)) > 0) igo = TRUE else igo = FALSE
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Loop through each element in the 'sf' vector
for (i in 1:length(sf)) {

  sfi = sf[i]  # Get the current element

  spl = strsplit(sfi, '/')[[1]]  # Split the element by '/'

  cl  = spl[7]  # Get the 7th element of 'spl'
  cli = strsplit(cl,'')[[1]][1]  # Get the first character of 'cl'
  k   = spl[8]  # Get the 8th element of 'spl'
  h   = spl[9]  # Get the 9th element of 'spl'
  c   = spl[10]  # Get the 10th element of 'spl'
  v   = spl[11]  # Get the 11th element of 'spl'
  tn  = strsplit(spl[12], '_|\\.')[[1]][2]  # Split the 12th element of 'spl' by '_' or '.'

  # Check if the string '_GO' is present in 'sfi'
  if (length(grep('_GO', sfi)) > 0) igo = TRUE else igo = FALSE

  # Create the path to the parameters file
  rfi = paste0('~/data/esse/empirical_comp_raw/', 
    cl,'/', k,'/',h,'/',c,'/',v,'/r_',tn,'.log')

  # Assign the appropriate value to 'vc' based on the value of 'v'
  if (v == 'rb0') vc = 'no_covariates'
  if (v == 'rl')  vc = 'linear'
  if (v == 'r0')  vc = 'log_area'
  if (v == 'r1')  vc = 'log_deriv_area'
  if (v == 'rt')  vc = 'temp'
  if (v == 'rdt') vc = 'deriv_temp'

  # Create the writing directory path
  wts = paste0('~/data/esse/empirical_vetted/',
    cl,'/', k,'/',h,'/',c,'/',vc,'/')

  # Create the path to the tree file based on the values of 'igo' and 'cli'
  if (igo) {
    trf = paste0('~/data/esse/trees/',cli,'tree_GO_',tn,'.tre')
  } else {
    if (cli == 'r') {
      trf = paste0('~/data/esse/trees/',cli,'tree0_',tn,'.tre')
    } else {
      trf = paste0('~/data/esse/trees/',cli,'tree_',tn,'.tre')
    }
  }
}

  # check if not processed already
  if (!file.exists(paste0(wts,'d_',tn,'.rda'))) {

    d = prepare_nmsp(trf, sfi, rfi)

    save(d, file = paste0(wts,'d_',tn,'.rda'))
  }

  cat(i, 'out of', length(sf), '\n')
}



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# plot


cs0 = c("#3b6a00D0",
        "#032e42D0",
        "#627e9aD0")

cs0f = c("#3b6a00",
         "#032e42",
         "#627e9a")

cs2 = c("#3b6a00D0",
        "#6a994eD0",
        "#032e42D0",
        "#0077B6D0",
        "#627e9aD0",
        "#7f9990D0")

cs2f = c("#3b6a00",
         "#6a994e",
         "#032e42",
         "#0077B6",
         "#627e9a",
         "#7f9990")


lf = list.files('~/data/esse/empirical_vetted', 
    full.names = TRUE, recursive = TRUE, pattern = 'd.*.rda$')

for (i in 1:length(lf)) {

  lfi = lf[i]

  spl = strsplit(lfi, '/')[[1]]

  cl  = spl[7]
  cli = strsplit(cl,'')[[1]][1]
  k   = spl[8]
  h   = spl[9]
  c   = spl[10]
  v   = spl[11]
  tn  = strsplit(spl[12], '_|\\.')[[1]][2]

  wts = paste0('~/data/esse/empirical_vetted/',
    cl,'/', k,'/',h,'/',c,'/',v,'/')

  load(paste0(wts,'d_',tn,'.rda'))

  # which hidden state colors
  if (length(grep('1', h)) > 0) {
    cols  = cs0
    colsf = cs0f
  }
  if (length(grep('2', h)) > 0) {
    cols  = cs2
    colsf = cs2f
  }

  ns = .subset2(d, 'ns')

  if (!file.exists(paste0(wts,'lttstates_',tn,'.pdf'))) {
    pdf(file = paste0(wts,'lttstates_',tn,'.pdf'), width = 14, height = 4)
      par(mfrow = c(1, 3))
      plot_nmsp_ltt(d, colstates = colsf, lwd = 2)
      plot_nmsp_prop(d, colstates = colsf)
      if (ns == 3) plot_densities_ns3(d, colstates = colsf)
      if (ns == 6) plot_densities_ns6(d, colstates = colsf)
    dev.off()
  }

  if (!file.exists(paste0(wts,'tree_',tn,'.png'))) {
    png(file = paste0(wts,'tree_',tn,'.png'), width = 3000, height = 3000)
      plot_nmsp(d, wcircle = .0125, colstates = colsf, lrect = 0.04)
    dev.off()
  }

  if (ns > 3 && !file.exists(paste0(wts,'tree_hs_',tn,'.png'))) {
    png(file = paste0(wts,'tree_hs_',tn,'.png'), width = 3000, height = 3000)
      plot_nmsp_hs(d, wcircle = .0125, lrect = 0.04,
        colstates = c('#011f4b', '#6497b1'))
    dev.off()
  }

  cat(i, 'out of', length(lf),'\n')
}



