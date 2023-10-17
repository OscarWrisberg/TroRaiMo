'''
------------------------------------------------------------------------------------------------------------------------
This is gonna be my workflow for the Biome estimation of all vascular plant species using the occurence records of
Herbarium specimens.

The Idea is that if you download the repository then you can run the entire analysis using GWF on a cluster using slurm.
Ideally it should run the entire analysis as a pipeline where the output of the frist function is the input to the next
one and result in all the final data output.

Download the repository from Github into a folder for the project.
When you run the gwf run command, the pipeline should by itself download the necessary data from the internet, and perform
all the required steps in order. 

------------------------------------------------------------------------------------------------------------------------



------------------------------------------------------------------------------------------------------------------------
Author: Oscar Wrisberg
Date: 13/Sep/2023
------------------------------------------------------------------------------------------------------------------------
'''

from gwf import Workflow, AnonymousTarget
from inspect import getsourcefile
import os

# Starting workflow
gwf = Workflow()


####################################################################################################
####################################################################################################
###############################---- Defining functions for workflow ----############################
####################################################################################################
####################################################################################################


##############################################################
##########---- Function for downloading data ----#############
##############################################################

# Webpage for Smith and Brown combined trees
# https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip

# Webpage for trees for individual orders
# https://www-personal.umich.edu/~eebsmith/big_seed_plant_datasets/trees/

def download_data(path_out,
                  smb_doi,
                  kew_doi,
                  output_gbif,
                  output_kew):
    """This function should download all the necessary files for the project"""
    inputs = []
    outputs = [path_out+output_kew,path_out+output_gbif]
    options = {
        'cores': 5,
        'memory': '40g',
        'account':"Biome_estimation",
        'walltime': "04:00:00"
    }

    spec = '''

    # Webpage for Smith and Brown combined trees
    # https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip

    # Webpage for trees for individual orders
    # https://www-personal.umich.edu/~eebsmith/big_seed_plant_datasets/trees/   

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    #Going to data folder
    cd {path_out}

    # Writing expected output
    echo "Expected output is: \n \n {path_out}{output_gbif} \n \n and: \n \n {path_out}{output_kew} \n \n "


    # Checking if file has already been downloaded. If not download it at # https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip

    # First im extracting the name of the downloaded file.
    filename_gbif=$(echo "https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip" | awk -F "/" '{{print $NF}}')
    echo "Were trying to download $filename_gbif"

    #Then i am checking if file exists
    if [ -f $filename_SmB ]; then
        echo "$filename_SmB is already downloaded"
    else
        echo "starting download of Smith & Brown Angiosperm Phylogeny data at: "
        date
        wget {smb_doi}
        echo "Finished downloading Smith & Brown Angiosperm Phylogeny data at: "
        date
    fi

    # Unzipping GBIF data
    if [ -f {path_out}{output_gbif} ]; then
        echo "Files from 
        Smith & Brown Angiosperm Phylogeny have already been unzipped"
    else
        echo "\n  starting to unzip at Smith & Brown Angiosperm Phylogeny data at: "
        date
        unzip -o {path_out}$filename_gbif
        echo " Finished unzipping Smith & Brown Angiosperm Phylogeny data at :"
        date
    fi


    '''.format(path_out = path_out, smb_doi = smb_doi, kew_doi = kew_doi, output_gbif=output_gbif, output_kew=output_kew)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
##########---- Function for downloading data ----#############
##############################################################

def (path_out):
    """This function should download all the necessary files for the project"""
    inputs = []
    outputs = [path_out]
    options = {
        'cores': 5,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }

    spec = '''



    '''.format(path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


########################################################################################################################
####################################---- Running pipeline ----##########################################################
########################################################################################################################

# Setting up some global variables with relative paths 
script_dir = os.path.dirname(getsourcefile(Rm_cols)) + "/"
data_dir = os.path.normpath(os.path.join(script_dir, "../data/")) +"/"
workflow_dir = os.path.normpath(os.path.join(script_dir, "../workflow/")) +"/"



gwf.target_from_template(name = "Download_Data",
                          template=download_data(
                            path_out = data_dir,
                            smb_doi = ,
                            kew_doi = ,
                            output_gbif ="occurrence.txt",
                            output_kew = "wcvp_names.csv"
                          ))

