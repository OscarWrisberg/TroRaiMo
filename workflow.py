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
        'account':"Trf_models",
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
        'cores': 10,
        'memory': '15g',
        'account':"Trf_models",
        'walltime': "00:30:00"
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
#########---- Adding random states to tips file ----##########
##############################################################
def sim_state_data(input_file, output_file, path_in,path_out, script_dir, nr_states):
    """This function should be used to simulate state data for the tips of the SmB tree
    by writing 0 or 1 for each state defined by nr_states and each tip in the input_file"""
    inputs = [path_in+input_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '20g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env

    # Going to input folder
    cd {path_in}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo Starting the Adding states script
    date

    # Loading the input file which is the file containing the tips of the SmB tree
    python3 {script_dir}adding_states.py {input_file} {nr_states} {output_file}

    echo Ended the Adding states script
    date

    mv {output_file} {path_out}

    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, nr_states = nr_states, input_file = input_file, output_file = output_file)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
#############---- Creating covariate files ----##############
##############################################################
def sim_covariate_data(input_file, output_file, path_in,path_out, script_dir, nr_states):
    """This function should be used to simulate the covariate data table through time for the states in the """
    inputs = [path_in+input_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '20g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate python3_env

    # Going to input folder
    cd {path_in}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo Starting the Adding states script
    date

    # Loading the input file which is the file containing the tips of the SmB tree
    python3 {script_dir}adding_states.py {input_file} {nr_states} {output_file}

    echo Ended the Adding states script
    date


    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, nr_states = nr_states, input_file = input_file, output_file = output_file)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#######################################################################
#######---- Function for Running ESSE on entire SmB tree ----##########
#######################################################################
"""
From the supplement material of the The build-up of the present-day tropical diversity of tetrapods: Line 215-220

_______________________________________________________________________________________________________________________
As a practical estimate, the inference procedure of the G+E+H ESSE model on mammals
(4001 species tree with 6 states and 20 free parameters) takes about 2 weeks of computation
time to estimate a MCMC chain with sufficient Effective Sample Sizes (ESS). Since there are
|K| - 1 states for |K| regions (see above), adding a region or adding more hidden states
doubles the number of states, with an exponential increase in the number of coupled equations
220 and parameters and thus of computational time.
_______________________________________________________________________________________________________________________

SOOOoooooo if 4001 species with 6 states (2 hidden, in 3 different regions (Tropical, Extratropical, Both)) takes 2 weeks, 
and a phylogeny with 50 species, 3 different regions and 2 hidden states takes 1 hour to run.
Then a phylogeny with 80000 species, 3 different regions and 2 hidden states takes 1600 hours to run. which is 66 days.....

This means I have to cut my tree into smaller pieces and run them separately.
cutting the tree into smaller 
I think I should aim at getting trees with around 2000 species in them in order to get a reasonable computation time.

I could do this by cutting the tree into smaller pieces based on either
1. The order of the species (i.e taxonomy of the species)
    I could potentially use the World checklist of vascular plants to infer the order of the species in the tree and then prune the tree based on this

2. The geographic location of the species (i.e. the region they are found in) (Tropical rainforest, , Both)
    I.e. by finding subtrees in the phylogeny with a proportion of species in Tropical rainforest over a certain threshold.

3. Time cut off ( i.e. cut the tree into smaller pieces by taking all the branches which )
    What time cut offs would make sense to use in regards to Tropical rainforest?
    And how would I find out if this results in some subtrees with more than 2000 species in them?

I think I need to try out all three and compare the results between them.

"""
##############################################################
##############---- Preparing WCVP-names file ----################
##############################################################
def apg_name_align(apg,wcp, output_file, path_in, script_dir, path_out):
    """Here I want to create a common format the wcvp and the previous file and update some family names to APGIV."""
    inputs = [script_dir+apg, path_in+wcp]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo This is the data dir \n
    echo {script_dir}

    # Starting Conda environment
    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    #Navigating to folder
    cd {path_in}


    Rscript --vanilla {script_dir}apg_name_aligner.r {script_dir}{apg} {wcp} {output_file}

    echo Done with Rscript

    mv {output_file} {path_out}

    '''.format(apg=apg,wcp = wcp, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
#############---- Approach 1 finding the orders ----##############
##############################################################
def finding_orders(input_file, output_file, path_in,path_out, script_dir, wcvp_file):
    """This function should be used to simulate the covariate data table through time for the states in the """
    inputs = [path_in+input_file, wcvp_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo Starting the Adding orders script
    date

    # Loading the input file which is the file containing the tips of the SmB tree
    Rscript --vanilla {script_dir}adding_order.R {input_file} {output_file} {wcvp_file} {path_out}


    echo Ended the Adding orders script
    date


    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, input_file = input_file, output_file = output_file, wcvp_file = wcvp_file)


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
                            output_file = "GBMB_tips.txt",
                            path_in = data_dir,
                            path_out = data_dir,
                            script_dir = script_dir
                          ))

gwf.target_from_template(name = "APG_preb",
                             template = apg_name_align(
                                 wcp = "wcvp_names.csv",
                                 apg ="apgweb_parsed.csv",
                                 output_file = "wcvp_names_apg_aligned.rds",
                                 path_in = data_dir,
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"02_adding_orders/"
                             ))
    
gwf.target_from_template(name = "slicing_families",
                        template=finding_orders(
                            input_file = "GBMB.tre",
                            output_file = "GBMB_orders.txt",
                            path_in = data_dir,
                            path_out = workflow_dir+"02_adding_orders/",
                            script_dir = script_dir,
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds"
                            ))




# gwf.target_from_template(name = "Esse",
#                           template=Esse(
#                             path_in = data_dir,
#                             path_out = "/home/owrisberg/Trf_models/Esse_test",
#                             script_dir = script_dir
#                           ))
