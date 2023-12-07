'''
------------------------------------------------------------------------------------------------------------------------
This is gonna be my workflow for the estimation of speciation, extinction and migration of all vascular plant orders
species using the occurence records of Herbarium specimens combined with the phylogenetic tree of Smith & Brown 2018.

Lastly this script will fit the ESSE models from Tapestree to the data and estimate the parameters for each order.

The Idea is that if you download the repository then you can run the entire analysis using GWF on a cluster using slurm.
Ideally it should run the entire analysis as a pipeline where the output of the frist function is the input to the next
one and result in all the final data output.

Download the repository from Github into a folder for the project.
When you run the gwf run command, the pipeline should by itself download the necessary data from the internet, and perform
all the required steps in the necessary order. 

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
# https://www-personal.umich.edu/~eebsmith/big_seed_plant_datasets/trees/ # why didnt I just use this.....

def download_data(path_out,
                  smb_doi,
                  kew_doi,
                  gbif_doi,
                  output_smb,
                  output_kew,
                  output_gbif):
    """This function should download all the necessary files for the project"""
    inputs = []
    outputs = [path_out+output_kew,path_out+output_smb,path_out+output_gbif]
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

    # Checking if file has already been downloaded. If not download it at # https://api.gbif.org/v1/occurrence/download/request/0012129-230828120925497.zip
    # First im extracting the name of the downloaded file.
    filename_gbif=$(echo "https://api.gbif.org/v1/occurrence/download/request/0012129-230828120925497.zip" | awk -F "/" '{{print $NF}}')
    echo "Were trying to download $filename_gbif"

    #Then i am checking if file exists
    if [ -f $filename_gbif ]; then
        echo "$filename_gbif is already downloaded"
    else
        echo "starting download of GBIF data at: "
        date
        wget {gbif_doi}
        echo "Finished downloading GBIF data at: "
        date
    fi

    # Unzipping GBIF data
    if [ -f {path_out}{output_gbif} ]; then
        echo "Files from Gbif have already been unzipped"
    else
        echo "\n  starting to unzip at GBIF data at: "
        date
        unzip -o {path_out}$filename_gbif
        echo " Finished unzipping GBIF data at :"
        date
    fi


    '''.format(path_out = path_out, smb_doi = smb_doi, kew_doi = kew_doi, output_smb=output_smb, output_kew=output_kew, gbif_doi = gbif_doi, output_gbif=output_gbif)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



################################################################################################################################################################################
############################################--------------------- Working on occurrence data ---------------------##############################################################
################################################################################################################################################################################

##############################################################
############---- Selecting usefull columns ----###############
##############################################################
def Rm_cols(input_file, output_file, path_in,path_out, occurrence_dir):
    """Here I want to remove unimportant columns from the file in order to save some ram space and increase computing speed"""
    inputs = [path_in+input_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '50g',
        'account':"Trf_models",
        'walltime': "00:10:00"
    }

    spec = '''
    #Checking if output dir exists
    [ -d {occurrence_dir} ] && echo "{occurrence_dir} exist." || {{ echo "{occurrence_dir} does not exist."; mkdir {occurrence_dir}; }}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo Starting the column removal script
    
    date

    cd {path_in}

    # The numbers here are the columns which I want to keep in my data.
    # and they are 
    #"gbifID", 1
    #"establishmentMeans", 34
    #"habitat", 65
    #"eventRemarks" 71,
    #"decimalLatitude" 92
    #"decimalLongitude" 93
    #"acceptedNameUsageID" 137
    #"scientificName" 143
    #"family" 155
    #"genus" 157
    #"taxonRank" 164
    #"datasetKey" 171
    #"speciesKey 192
    #"species" 193


    cut -d "\t" -f 1,34,65,71,92,93,137,143,155,157,164,171,192,193 {input_file} > {output_file}

    echo ending the column removal script

    mv {output_file} {path_out}
    
    date

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, path_out = path_out, occurrence_dir = occurrence_dir)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
################---- Loading the dataset ----#################
##############################################################

def Load_data(input_file, output_file, path_in,path_out, script_dir):
    """Here I load the GBIF data and save it as an rds file for further analysis."""
    inputs = [path_in+input_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '50g',
        'account':"Trf_models",
        'walltime': "03:10:00"
    }

    spec = '''
     #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    cd {path_in}

    echo Starting the R script
    
    date

    Rscript --vanilla {script_dir}Downloading_data.r {input_file} {output_file}

    echo Ended the R script

    date
    
    mv {output_file} {path_out}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
################---- Taxon lookup on Gbif ----################
##############################################################
def taxon_look_up(input_file, output_file, path_in, script_dir, path_out):
    """Here the function looks up each unique species in the dataset and finds the taxonomic information for that species."""
    inputs = [path_in+input_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '20g',
        'account':"Trf_models",
        'walltime': "100:00:00"
    }

    spec = '''
    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    # Starting up conda env
    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Navigating to input data folder
    cd {path_in}

    echo Starting the R script
    
    date

    # Running the taxon look up script
    Rscript --vanilla {script_dir}Gbif_taxon_look_up.r {path_in}{input_file} {output_file}

    echo Ended the R script

    date

    mv {output_file} {path_out}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
##############---- Creating common format ----################
##############################################################
def create_common_format(input_file_occurrences,input_file_taxonomy, output_file, path_in, script_dir, path_out):
    """Here I want to create a common format between the file which needs names aligned to the WCVP and the WCVP."""
    inputs = [path_in+input_file_taxonomy, input_file_occurrences]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '100g',
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


    Rscript --vanilla {script_dir}Create_common_format.r {input_file_occurrences} {input_file_taxonomy} {output_file}

    echo Done with Rscript

    mv {output_file} {path_out}

    '''.format(input_file_taxonomy=input_file_taxonomy,input_file_occurrences=input_file_occurrences, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

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
##################---- Taxonomy Matcher ----##################
##############################################################
def taxon_match(input_file, output_file, path_in, script_dir, path_out, wcvp):
    """Here I match the taxa from the GBIF file to the WCVP."""
    inputs = [path_in+input_file, path_in+wcvp]
    outputs = [path_out+output_file]
    options = {
        'cores': 15,
        'memory': '75g',
        'account':"Trf_models",
        'walltime': "48:00:00"
    }

    spec = '''
    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    cd {path_in}

    echo Starting the Taxon_matcher script at: 
    date
    
    Rscript --vanilla {script_dir}Taxon_matcher.r {input_file} {wcvp} {output_file} 

    mv {output_file} {path_out}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, wcvp = wcvp)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



########################################################################################################################
####################################---- Taxon renamer ----#############################################################
########################################################################################################################

def Renamer(input_file, output_file, path_in, script_dir, path_out, wcvp, renaming_file):
    """This function renames all the species names in the GBIF datafile based on the taxon matcher.
    This should be done by looping through the GBIF data and looking up each species in the taxon matcher file.
    The name in the taxon matcher file would then be used to find the accepted_plant_name_id in the WCVP file.
    and get the correct name from that file. This name would then be used to replace the name in the GBIF data."""
    inputs = [input_file, wcvp, path_in+renaming_file]
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

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    cd {path_in}

    Rscript --vanilla {script_dir}renaming.r {input_file} {wcvp} {renaming_file} {output_file}


    mv {output_file} {path_out}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, wcvp = wcvp, renaming_file = renaming_file)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

################################################################################################################################################################################
############################################--------------------- Working on tree data ---------------------##############################################################
################################################################################################################################################################################

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
####---- Pruning the GBMB tree for tips not in WCVP ----######
##############################################################
def pruning_tree(wcp,tree, output_file, path_in, script_dir, path_out):
    """Here I am looping through the tips in the GBMB tree to see if they are in the WCVP.
    If there is no exact match for the tip, I will look for a fuzzy match allowing for 1 substitution, insertion or deletion.
    If the tip name is longer than 2 IE. a subsp or a variety I will then search for just the genus and species epithet.
    Lastly I check if this has introduced any duplicate species in the tree.
    If the duplicates are sister species I will then remove one of them at random and if they are not I will remove both of them"""
    inputs = [wcp, path_in+tree]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '20g',
        'account':"Trf_models",
        'walltime': "30:00:00" # should be atleast 24 hours if you dont load from RDS files
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


    Rscript --vanilla {script_dir}pruning.r {tree} {wcp} {output_file}

    echo Done with Rscript

    mv {output_file} {path_out}

    '''.format(wcp = wcp, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, tree = tree)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
#############---- Approach 1 finding the orders ----##############
##############################################################
def Slicing_trees(input_file, output_file, path_in,path_out, script_dir, wcvp_file, apg):
    """This function slices the GBMB tree into subtrees based on the order and families of the species in the tree. """
    inputs = [path_in+input_file, wcvp_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 5,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "03:00:00"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo Starting the Adding orders script
    date

    # Loading the input file which is the file containing the tips of the SmB tree
    Rscript --vanilla {script_dir}Slicing_tree.R {input_file} {output_file} {wcvp_file} {path_out} {apg}


    echo Ended the Adding orders script
    date


    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, input_file = input_file, output_file = output_file, wcvp_file = wcvp_file, apg = apg)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




##############################################################
###########---- Slicing the tree into orders ----############
##############################################################
def Forcing_orders(input_file_tree, output_file, path_in,path_out, script_dir, wcvp_file, apg):
    """This function searchers for the largest monophyletic clades in the orders which are not monophyletic in the GBMB tree."""
    inputs = [path_in+input_file_tree, wcvp_file]
    outputs = [path_out+output_file]
    options = {
        'cores': 10,
        'memory': '30g',
        'account':"Trf_models",
        'walltime': "12:00:00"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo Starting the Forcing monophyly script
    date

    # Running the R script
    Rscript --vanilla {script_dir}Finding_monophyletic_clades.R {input_file_tree} {wcvp_file} {output_file} {path_out} {apg}


    echo Ended the Forcing monophyly script
    date


    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, input_file_tree = input_file_tree, output_file = output_file, wcvp_file = wcvp_file, apg = apg)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
###########---- Downloading distribution data ----############
##############################################################
def Finding_areas_in_wcvp(input_file_tree, wcvp_file,path_out, output_file, path_in, order, script_dir):
    """This Function creates a states file for the tips in WCVP based on the climate column."""
    inputs = [path_in+input_file_tree]
    outputs = [path_out+output_file]
    options = {
        'cores': 2,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "00:10:00"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    echo Starting the script to find state data for the tips in the wcvp 
    date

    # Running the R script
    Rscript --vanilla {script_dir}Finding_monophyletic_clades.R {input_file_tree} {output_file} {wcvp_file} {path_out} {order}


    echo Ended the script to find state data for the tips in the wcvp
    date

    '''.format(path_out = path_out, output_file = output_file, wcvp_file = wcvp_file, order = order, input_file_tree = input_file_tree, path_in = path_in, script_dir = script_dir)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


########################################################################################################################
####################################---- Running pipeline ----##########################################################
########################################################################################################################

# Setting up some global variables with relative paths 
script_dir = os.path.dirname(getsourcefile(download_data)) + "/"
data_dir = os.path.normpath(os.path.join(script_dir, "../data/")) +"/"
workflow_dir = os.path.normpath(os.path.join(script_dir, "../workflow/")) +"/"



gwf.target_from_template (name = "Download_Data",
                          template=download_data(
                            path_out = data_dir,
                            smb_doi = "https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip",
                            kew_doi = "http://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip",
                            gbif_doi = "https://api.gbif.org/v1/occurrence/download/request/0012129-230828120925497.zip",
                            output_smb ="GBMB.tre",
                            output_kew = "wcvp_names.csv",
                            output_gbif ="occurrence.txt"
                          ))

#########################################################################################################################
####################################---- Working on the occurrence data ----#############################################
#########################################################################################################################

gwf.target_from_template(name = "Removing_cols",
                          template=Rm_cols(
                            input_file ="occurrence.txt",
                            output_file = "occurrence_select.txt",
                            path_in = data_dir,
                            occurrence_dir = workflow_dir+"01_distribution_data/",
                            path_out = workflow_dir+"01_distribution_data/01_rm_cols/",
                          ))

gwf.target_from_template(name = "Data_Parsing",
                             template = Load_data(
                                 input_file ="occurrence_select.txt",
                                 output_file = "gbif_parsed.rds",
                                 path_in = workflow_dir+"01_distribution_data/01_rm_cols/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/02_data_parsing/"
                             ))

gwf.target_from_template(name = "GBIF_lookup",
                             template = taxon_look_up(
                                 input_file ="gbif_parsed.rds",
                                 output_file = "gbif_parsed_taxon_data.rds",
                                 path_in = workflow_dir+"01_distribution_data/02_data_parsing/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/03_taxon_lookup/"
                             ))

gwf.target_from_template(name = "Creating_Common_Format",
                             template = create_common_format(
                                 input_file_taxonomy ="gbif_parsed_taxon_data.rds",
                                 input_file_occurrences =workflow_dir+"01_distribution_data/02_data_parsing/gbif_parsed.rds",
                                 output_file = "gbif_common_format.rds",
                                 path_in = workflow_dir+"01_distribution_data/03_taxon_lookup/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/04_common_format/"
                             ))

gwf.target_from_template(name = "APG_preb_occurrences",
                             template = apg_name_align(
                                 wcp = "wcvp_names.csv",
                                 apg ="apgweb_parsed.csv",
                                 output_file = "wcvp_names_apg_aligned.rds",
                                 path_in = data_dir,
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/04_common_format/"
                             ))

gwf.target_from_template(name = "Taxon_match",
                             template = taxon_match(
                                 input_file ="gbif_common_format.rds",
                                 wcvp = "wcvp_names_apg_aligned.rds",
                                 output_file = "gbif_taxon_matched.rds",
                                 path_in = workflow_dir+"01_distribution_data/04_common_format/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/05_Taxon_match/"
                             ))

gwf.target_from_template(name = "Renaming",
                                template = Renamer(
                                    input_file =workflow_dir+"01_distribution_data/04_common_format/gbif_common_format.rds",
                                    wcvp = workflow_dir+"01_distribution_data/04_common_format/wcvp_names_apg_aligned.rds",
                                    renaming_file = "gbif_taxon_matched.rds",
                                    output_file = "gbif_renamed.rds",
                                    path_in = workflow_dir+"01_distribution_data/05_Taxon_match/",
                                    script_dir = script_dir,
                                    path_out = workflow_dir+"01_distribution_data/06_Renamed"
                                ))

################################################################################################################################
############################------- Starting on the tree data -------###########################################################
################################################################################################################################

gwf.target_from_template(name = "Load_tree",
                          template=Load_tree(
                            input_file = "GBMB.tre",
                            output_file = "GBMB_tips.txt",
                            path_in = data_dir,
                            path_out = data_dir,
                            script_dir = script_dir
                          ))

gwf.target_from_template(name = "APG_preb_tree",
                             template = apg_name_align(
                                 wcp = "wcvp_names.csv",
                                 apg ="apgweb_parsed.csv",
                                 output_file = "wcvp_names_apg_aligned.rds",
                                 path_in = data_dir,
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"02_adding_orders/"
                             ))

gwf.target_from_template(name = "Pruning_tree",
                            template=pruning_tree(
                                wcp = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                tree = "GBMB.tre",
                                output_file = "GBMB_pruned.tre",
                                path_in = data_dir,
                                path_out = data_dir,
                                script_dir = script_dir
                            ))
    
gwf.target_from_template(name = "slicing_Trees_no_pruning",
                        template=Slicing_trees(
                            input_file = "GBMB.tre",
                            output_file = "GBMB_sp_per_orders_no_pruning.txt",
                            path_in = data_dir,
                            path_out = workflow_dir+"02_adding_orders/no_pruning/",
                            script_dir = script_dir,
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                            apg = script_dir+"apgweb_parsed.csv"
                            ))

gwf.target_from_template(name = "slicing_Trees_pruning",
                        template=Slicing_trees(
                            input_file = "GBMB_pruned.tre",
                            output_file = "GBMB_sp_per_orders_pruning.txt",
                            path_in = data_dir,
                            path_out = workflow_dir+"02_adding_orders/pruning/",
                            script_dir = script_dir,
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                            apg = script_dir+"apgweb_parsed.csv"
                            ))

gwf.target_from_template(name = "Finding_monophyletic_orders",
                        template=Forcing_orders(
                            input_file_tree = "GBMB_pruned.tre", 
                            output_file = "Orders_which_could_not_be_solved.txt", 
                            path_in = data_dir, 
                            path_out = workflow_dir+"02_adding_orders/pruning/", 
                            script_dir = script_dir, 
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds", 
                            apg = script_dir+"apgweb_parsed.csv" 
                            ))

order_trees = ["pruned_tree__order_Acorales_GBMB.txt", "pruned_tree__order_Crossosomatales_GBMB.txt", "pruned_tree__order_Gnetales_GBMB.txt", "pruned_tree__order_Piperales_GBMB.txt",
"pruned_tree__order_Alismatales_GBMB.txt", "pruned_tree__order_Cucurbitales_GBMB.txt", "pruned_tree__order_Gunnerales_GBMB.txt", "pruned_tree__order_Poales_GBMB.txt",
"pruned_tree__order_Amborellales_GBMB.txt", "pruned_tree__order_Cupressales_GBMB.txt", "pruned_tree__order_Huerteales_GBMB.txt",
"pruned_tree__order_Aquifoliales_GBMB.txt", "pruned_tree__order_Cycadales_GBMB.txt", "pruned_tree__order_Magnoliales_GBMB.txt", "pruned_tree__order_Ranunculales_GBMB.txt",
"pruned_tree__order_Arecales_GBMB.txt", "pruned_tree__order_Dilleniales_GBMB.txt", "pruned_tree__order_Malpighiales_GBMB.txt", "pruned_tree__order_Rosales_GBMB.txt",
"pruned_tree__order_Austrobaileyales_GBMB.txt", "pruned_tree__order_Dioscoreales_GBMB.txt", "pruned_tree__order_Malvales_GBMB.txt", "pruned_tree__order_Santalales_GBMB.txt",
"pruned_tree__order_Berberidopsidales_GBMB.txt", "pruned_tree__order_Dipsacales_GBMB.txt", "pruned_tree__order_Myrtales_GBMB.txt", "pruned_tree__order_Sapindales_GBMB.txt",
"pruned_tree__order_Bruniales_GBMB.txt", "pruned_tree__order_Equisetales_GBMB.txt", "pruned_tree__order_Nymphaeales_GBMB.txt", "pruned_tree__order_Solanales_GBMB.txt",
"pruned_tree__order_Buxales_GBMB.txt", "pruned_tree__order_Ericales_GBMB.txt", "pruned_tree__order_Osmundales_GBMB.txt", "pruned_tree__order_Trochodendrales_GBMB.txt",
"pruned_tree__order_Canellales_GBMB.txt", "pruned_tree__order_Escalloniales_GBMB.txt", "pruned_tree__order_Pandanales_GBMB.txt", "pruned_tree__order_Vahliales_GBMB.txt",
"pruned_tree__order_Celastrales_GBMB.txt", "pruned_tree__order_Fabales_GBMB.txt", "pruned_tree__order_Paracryphiales_GBMB.txt", "pruned_tree__order_Vitales_GBMB.txt",
"pruned_tree__order_Ceratophyllales_GBMB.txt", "pruned_tree__order_Garryales_GBMB.txt", "pruned_tree__order_Petrosaviales_GBMB.txt", "pruned_tree__order_Zingiberales_GBMB.txt",
"pruned_tree__order_Chloranthales_GBMB.txt", "pruned_tree__order_Ginkgoales_GBMB.txt", "pruned_tree__order_Picramniales_GBMB.txt", "pruned_tree__order_Zygophyllales_GBMB.txt",
"twice_pruned_tree_Liliales_GBMB.txt","twice_pruned_tree_Apiales_GBMB.txt","twice_pruned_tree_Metteniusales_GBMB.txt","twice_pruned_tree_Asterales_GBMB.txt",
"twice_pruned_tree_Boraginales_GBMB.txt", "twice_pruned_tree_Oxalidales_GBMB.txt","twice_pruned_tree_Brassicales_GBMB.txt","twice_pruned_tree_Pinales_GBMB.txt",
"twice_pruned_tree_Caryophyllales_GBMB.txt","twice_pruned_tree_Proteales_GBMB.txt","twice_pruned_tree_Cornales_GBMB.txt",
"twice_pruned_tree_Geraniales_GBMB.txt", "twice_pruned_tree_Saxifragales_GBMB.txt",
"twice_pruned_tree_Icacinales_GBMB.txt","twice_pruned_tree_Fagales_GBMB.txt", "twice_pruned_tree_Lamiales_GBMB.txt"
]

orders = ["Acorales", "Crossosomatales", "Gnetales", "Piperales","Alismatales", "Cucurbitales", "Gunnerales", "Poales","Amborellales", "Cupressales", "Huerteales",
"Aquifoliales", "Cycadales", "Magnoliales", "Ranunculales","Arecales", "Dilleniales", "Malpighiales", "Rosales","Austrobaileyales", "Dioscoreales", "Malvales", "Santalales","Berberidopsidales", "Dipsacales", "Myrtales", "Sapindales",
"Bruniales", "Equisetales", "Nymphaeales", "Solanales","Buxales", "Ericales", "Osmundales", "Trochodendrales","Canellales", "Escalloniales", "Pandanales", "Vahliales","Celastrales", "Fabales", "Paracryphiales", "Vitales",
"Ceratophyllales", "Garryales", "Petrosaviales", "Zingiberales","Chloranthales", "Ginkgoales", "Picramniales", "Zygophyllales","Liliales", "Apiales", "Metteniusales", "Asterales",
"Boraginales", "Oxalidales", "Brassicales", "Pinales","Caryophyllales", "Proteales", "Cornales",
 "Geraniales", "Saxifragales", "Icacinales","Fagales", "Lamiales"
]

# Problematic trees
# "Polypodiales-eupolypod_I", "pruned_tree__order_Polypodiales-eupolypod_I_GBMB.txt",
# "Alismatales","twice_pruned_tree_Alismatales_GBMB.txt"
# "Magnoliales","twice_pruned_tree_Magnoliales_GBMB.txt"
# "Aquifoliales" "twice_pruned_tree_Aquifoliales_GBMB.txt",
# "Crossosomatales", "twice_pruned_tree_Crossosomatales_GBMB.txt"
# "Cucurbitales", "twice_pruned_tree_Cucurbitales_GBMB.txt",
# "Gunnerales", "twice_pruned_tree_Gunnerales_GBMB.txt",
# "Cycadales", "twice_pruned_tree_Cycadales_GBMB.txt",
# "Ranunculales", "twice_pruned_tree_Ranunculales_GBMB.txt",
# "Arecales", "twice_pruned_tree_Arecales_GBMB.txt",
#  "Malpighiales", "twice_pruned_tree_Malpighiales_GBMB.txt",
# "Santalales", "twice_pruned_tree_Santalales_GBMB.txt",
# "Sapindales", "twice_pruned_tree_Sapindales_GBMB.txt",
#  "Nymphaeales",  "twice_pruned_tree_Nymphaeales_GBMB.txt",
# "Solanales",  "twice_pruned_tree_Solanales_GBMB.txt",
# "Buxales", "twice_pruned_tree_Buxales_GBMB.txt",
# "Canellales", "twice_pruned_tree_Canellales_GBMB.txt",
# "Pandanales", "twice_pruned_tree_Pandanales_GBMB.txt",
# "Celastrales", "twice_pruned_tree_Celastrales_GBMB.txt",
# "Zingiberales" "twice_pruned_tree_Zingiberales_GBMB.txt",
# "Chloranthales", "twice_pruned_tree_Chloranthales_GBMB.txt",
# "Picramniales", "twice_pruned_tree_Picramniales_GBMB.txt",
# "#",  "twice_pruned_tree_Zygophyllales_GBMB.txt"








for i in range(len(order_trees)):
    #### Running the script to find the environmental data for the tips in the trees
    gwf.target_from_template(name = orders[i]+"_distribution_data.",
                                                        template=Finding_areas_in_wcvp(
                                                        input_file_tree= order_trees[i],
                                                        path_in =  workflow_dir+"02_adding_orders/pruning/",
                                                        path_out = workflow_dir+"03_distribution_data",
                                                        output_file = orders[i]+"_distribution_data.txt",
                                                        wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                                        order = orders[i],
                                                        script_dir= script_dir
                                                        ))



# gwf.target_from_template(name = "Esse",
#                           template=Esse(
#                             path_in = data_dir,
#                             path_out = "/home/owrisberg/Trf_models/Esse_test",
#                             script_dir = script_dir
#                           ))