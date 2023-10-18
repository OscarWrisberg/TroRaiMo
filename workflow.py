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

# Before Running this workflow make sure to download GWF from here: 
# And activate the slurm backend using the following command: gwf config set backend slurm

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
                  output_smb,
                  output_kew):
    """This function should download all the necessary files for the project"""
    inputs = []
    outputs = [path_out+output_kew,path_out+output_smb]
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
    echo "Expected output is: \n \n {path_out}{output_smb} \n \n and: \n \n {path_out}{output_kew} \n \n "


    # Checking if file has already been downloaded. If not download it at # https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip

    # First im extracting the name of the downloaded file.
    filename_smb=$(echo "https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip" | awk -F "/" '{{print $NF}}')
    echo "Were trying to download $filename_smb"

    #Then i am checking if file exists
    if [ -f $filename_smb ]; then
        echo "$filename_smb is already downloaded"
    else
        echo "starting download of Smith & Brown Angiosperm Phylogeny data at: "
        date
        wget {smb_doi}
        echo "Finished downloading Smith & Brown Angiosperm Phylogeny data at: "
        date
    fi

    # Unzipping SmB data
    if [ -f {path_out}{output_smb} ]; then
        echo "Files from 
        Smith & Brown Angiosperm Phylogeny have already been unzipped"
    else
        echo "\n  starting to unzip at Smith & Brown Angiosperm Phylogeny data at: "
        date
        unzip -o {path_out}$filename_smb
        echo " Finished unzipping Smith & Brown Angiosperm Phylogeny data at :"
        date
        cd $filename_smb
        mv * ../.
    fi

    # Checking if file has already been downloaded. if not download it at  # "http://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip"

    # First im extracting the name of the downloaded file.
    filename_kew=$(echo "http://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip" | awk -F "/" '{{print $NF}}')
    echo "Were trying to download $filename_kew"
    

    #Then i am checking if the file exists    
    if [ -f $filename_kew ]; then
        echo " $filename_kew is already downloaded \n"
    else
        echo "starting download of Kew data at: "
        date
        wget {kew_doi}
        echo " Finished downloading Kew data at :"
        date
    fi

    # Unzipping Kew data
    if [ -f {path_out}{output_kew}]; then
        echo "Files from Kew have already been unzipped"
    else
        echo "\n  starting to unzip at Kew data at: "
        date
        unzip -o {path_out}$filename_kew
        echo " Finished unzipping Kew data at :"
        date
    fi


    '''.format(path_out = path_out, smb_doi = smb_doi, kew_doi = kew_doi, output_smb=output_smb, output_kew=output_kew)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
#############---- Loading the SmB tree tips ----##############
##############################################################
def Load_tree(input_file, output_file, path_in,path_out, script_dir):
    """Here I load the SmB tree and save it as an rds file for further analysis."""
    inputs = [path_in+input_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 1,
        'memory': '5g',
        'account':"Biome_estimation",
        'walltime': "00:10:00"
    }

    spec = '''

    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    cd {path_in}

    echo Starting the R script
    
    date

    # Using the R script to load the tree into R and save the tips as a csv file.
    Rscript --vanilla {script_dir}R_tree_loading.r {input_file} {output_file}

    echo Ended the R script

    date
    
    mv {output_file} {path_out}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in,script_dir = script_dir, path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
#############---- Loading the SmB tree tips ----##############
##############################################################
# def 


#######################################################################
#######---- Function for Running ESSE on entire SmB tree ----##########
#######################################################################

def Esse(path_in, path_out, script_dir, ):
    """Function for running ESSE on the Smith and Brown phylogeny"""
    inputs = [path_in]
    outputs = [path_out]
    options = {
        'cores': 5,
        'memory': '25g',
        'account':"Trf_models",
        'walltime': "03:00:00"
    }

    spec = '''

    # Going to input folder
    cd {path_in}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate julia_env

    echo " Starting to run Julia script \n "
    date

    julia {script_dir}Esse_test.jl

    echo " Finished running Julia script \n "
    date

    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


########################################################################################################################
####################################---- Running pipeline ----##########################################################
########################################################################################################################

# Setting up some global variables with relative paths 
script_dir = os.path.dirname(getsourcefile(download_data)) + "/"
data_dir = os.path.normpath(os.path.join(script_dir, "../data/")) +"/"
workflow_dir = os.path.normpath(os.path.join(script_dir, "../workflow/")) +"/"



gwf.target_from_template(name = "Download_Data",
                          template=download_data(
                            path_out = data_dir,
                            smb_doi = "https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip",
                            kew_doi = "http://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip",
                            output_smb ="GBMB.tre",
                            output_kew = "wcvp_names.csv"
                          ))

gwf.target_from_template(name = "Load_tree",
                          template=Load_tree(
                            input_file = "GBMB.tre",
                            output_file = "GBMB_tips.rds",
                            path_in = data_dir,
                            path_out = data_dir,
                            script_dir = script_dir
                          ))

gwf.target_from_template(name = "Esse",
                          template=Esse(
                            path_in = data_dir,
                            path_out = "/home/owrisberg/Trf_models/Esse_test",
                            script_dir = script_dir
                          ))
