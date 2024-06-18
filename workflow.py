'''
------------------------------------------------------------------------------------------------------------------------
This is gonna be my workflow for the estimation of speciation, extinction and migration of all vascular plant orders
species using the occurence records of Herbarium specimens combined with the phylogenetic tree of Smith & Brown 2018.

The Idea is that if you download the repository then you can run the entire analysis using GWF on a cluster using slurm.
Ideally it should run the entire analysis as a pipeline where the output of the first function is the input to the next
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

###################################################################################################################################################################################
###################################################################################################################################################################################
#########################################################################---- Function for downloading data ----###################################################################
###################################################################################################################################################################################
###################################################################################################################################################################################

# Webpage for Smith and Brown combined trees
# https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip

# Webpage for trees for individual orders
# https://www-personal.umich.edu/~eebsmith/big_seed_plant_datasets/trees/ # why didnt I just use this.....

def download_data(path_out,
                  smb_doi,
                  kew_doi,
                  gbif_doi,
                  paleo_doi,
                  output_smb,
                  output_kew,
                  output_gbif,
                  output_paleo,
                  done_dir,
                  done):
    """This function should download all the necessary files for the project"""
    inputs = []
    outputs = [path_out+output_kew,path_out+output_smb,path_out+output_gbif, path_out+"paleo_clim/"+output_paleo, done_dir+done]
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

    # Webpage for paleoclimates
    # https://zenodo.org/records/6620748/files/All_CSV_files.zip?download=1

    # Webpage for GBIF data
    # https://api.gbif.org/v1/occurrence/download/request/0012129-230828120925497.zip

    # Webpage for Kew data
    # http://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip

    # Creating a done folder to keep track of which steps have been completed
    [ -d {done_dir} ] && echo "{done_dir} exist." || {{ echo "{done_dir} does not exist."; mkdir {done_dir}; }}

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

    #Checking if output dir exists
    [ -d {path_out}paleo_clim/ ] && echo "{path_out}paleo_clim/ exist." || {{ echo "{path_out}paleo_clim/ does not exist."; mkdir {path_out}paleo_clim/; }}

    # Going to subfolder
    cd {path_out}paleo_clim/
    
    # Extracting the name of the downloaded file without query parameters.
    filename_paleo=$(echo {paleo_doi} | awk -F "/" '{{split($NF, name, "?"); print name[1]}}')
    echo "We're trying to download $filename_paleo"


    #Then i am checking if file exists
    if [ -f $filename_paleo ]; then
        echo "$filename_paleo is already downloaded"
    else
        echo "Starting downlaod of paleoclimate data at: "
        date
        wget {paleo_doi}
        echo "Finished downloading paleoclimate data at: "
        date
    fi

    # If the file is downloaded but it has the wrong name, then rename it.
    if [ -f All_CSV_files.zip ]; then
        echo "All_CSV_files.zip is already renamed"
    else
    mv 'All_CSV_files.zip?download=1' All_CSV_files.zip
    fi

    # Unzipping paleo data
    if [ -f {path_out}paleo_clim/{output_paleo} ]; then
        echo "Files from 
        Paleo data has already been unzipped"
    else
        echo "\n  starting to unzip paleodata data at: "
        date
        unzip -o {path_out}paleo_clim/$filename_paleo
        echo " Finished unzipping paleodata data at :"
        date
        cd ${{filename_paleo%.zip}}
        mv * ../.
    fi


    touch {done_dir}{done}

    


    '''.format(path_out = path_out, smb_doi = smb_doi, kew_doi = kew_doi, output_smb=output_smb, output_kew=output_kew, gbif_doi = gbif_doi,
                output_gbif=output_gbif, paleo_doi = paleo_doi, output_paleo = output_paleo, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

###################################################################################################################################################################################
###################################################################################################################################################################################
############################################--------------------- Working on paleoclimatic data ---------------------##############################################################
###################################################################################################################################################################################
###################################################################################################################################################################################

def paleo_clim_area(output_file, data_dir, script_dir,done_dir, done):
    """Here I calculate the area of the Tropical rainforests through time."""
    inputs = [data_dir, done_dir+"Download_Data"]
    outputs = [data_dir+output_file,done_dir+done]
    options = {
        'cores': 5,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "00:10:00"
    }

    spec = '''

    package_name="libgdal-dev"

    if dpkg -s "$package_name" >/dev/null 2>&1; then
        echo "$package_name is installed."
        # Going to datadir
        cd "$data_dir"

        source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
        conda activate R_env

        echo Starting the R script
        date

        Rscript --vanilla {script_dir}calculating_paleoclim.r {data_dir} {output_file}

        echo Ended the R script

        touch {done_dir}{done}
    else
        echo "$package_name is not installed."
        echo output_file is : {output_file}
        echo script_dir is : {script_dir}
        if [ -e {script_dir}{output_file} ]; then
            echo "But {output_file} is in {script_dir}, so were copying it to {data_dir} instead, See ReadMe for more info"
            cp {script_dir}{output_file} {data_dir}
            touch {done_dir}{done}
        else
            echo "{output_file} is not in {script_dir}"
            echo "Please check the repository for solutions"
        fi
    fi

    '''.format(output_file=output_file, script_dir = script_dir, data_dir = data_dir, done=done, done_dir = done_dir)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

################################################################################################################################################################################
################################################################################################################################################################################
############################################--------------------- Working on occurrence data ---------------------##############################################################
################################################################################################################################################################################
################################################################################################################################################################################

##############################################################
############---- Selecting usefull columns ----###############
##############################################################
def Rm_cols(input_file, output_file, path_in,path_out, occurrence_dir, done, done_dir):
    """Here I want to remove unimportant columns from the file in order to save some ram space and increase computing speed"""
    inputs = [path_in+input_file, done_dir+"Download_Data"]
    outputs = [path_out+output_file,done_dir+done]
    options = {
        'cores': 5,
        'memory': '50g',
        'account':"Trf_models",
        'walltime': "00:20:00"
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

    touch {done_dir}{done}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, path_out = path_out, occurrence_dir = occurrence_dir, done = done, done_dir = done_dir)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
################---- Loading the dataset ----#################
##############################################################

def Load_data(input_file, output_file, path_in,path_out, script_dir, done_dir, done):
    """Here I load the GBIF data and save it as an rds file for further analysis."""
    inputs = [path_in+input_file, done_dir+"Removing_cols"]
    outputs = [path_out+output_file,done_dir+done]
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

    touch {done_dir}{done}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
################---- Taxon lookup on Gbif ----################
##############################################################
def taxon_look_up(input_file, output_file, path_in, script_dir, path_out, done_dir, done):
    """Here the function looks up each unique species in the dataset and finds the taxonomic information for that species."""
    inputs = [path_in+input_file, done_dir+"Data_Parsing"]
    outputs = [path_out+output_file,done_dir+done]
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

    touch {done_dir}{done}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
##############---- Creating common format ----################
##############################################################

def create_common_format(input_file_occurrences,input_file_taxonomy, output_file, path_in, script_dir, path_out, done_dir, done):
    """Here I want to create a common format between the file which needs names aligned to the WCVP and the WCVP."""
    inputs = [path_in+input_file_taxonomy, input_file_occurrences, done_dir+"GBIF_lookup"]
    outputs = [path_out+output_file,done_dir+done]
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

    touch {done_dir}{done}

    '''.format(input_file_taxonomy=input_file_taxonomy,input_file_occurrences=input_file_occurrences, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#################################################################
##############---- Preparing WCVP-names file ----################
#################################################################
def apg_name_align(apg,wcp, output_file, path_in, script_dir, path_out, done_dir, done):
    """Here I want to create a common format the wcvp and the previous file and update some family names to APGIV."""
    inputs = [script_dir+apg, path_in+wcp, done_dir+"Creating_Common_Format"]
    outputs = [path_out+output_file,done_dir+done]
    options = {
        'cores': 5,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    echo This is the script dir \n
    echo {script_dir}

    # Starting Conda environment
    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    #Navigating to folder
    cd {path_in}

    echo Starting the R script at:
    date
    
    Rscript --vanilla {script_dir}apg_name_aligner.r {script_dir}{apg} {wcp} {output_file}

    echo Done with Rscript at
    date

    mv {output_file} {path_out}

    touch {done_dir}{done}

    '''.format(apg=apg,wcp = wcp, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
##################---- Taxonomy Matcher ----##################
##############################################################
def taxon_match(input_file, output_file, path_in, script_dir, path_out, wcvp, done_dir, done):
    """Here I match the taxa from the GBIF file to the WCVP."""
    inputs = [path_in+input_file, path_in+wcvp, done_dir+"APG_preb_occurrences"]
    outputs = [path_out+output_file,done_dir+done]
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

    touch {done_dir}{done}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, wcvp = wcvp, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



###################################################################
######################---- Taxon renamer ----######################
###################################################################

def Renamer(input_file, output_file, path_in, script_dir, path_out, wcvp, renaming_file, done, done_dir, raw_occurrences):
    """This function renames all the species names in the GBIF datafile based on the taxon matcher.
    This should be done by looping through the GBIF data and looking up each species in the taxon matcher file.
    The name in the taxon matcher file would then be used to find the accepted_plant_name_id in the WCVP file.
    and get the correct name from that file. This name would then be used to replace the name in the GBIF data."""
    inputs = [input_file, wcvp, path_in+renaming_file, done_dir+"Taxon_match"]
    outputs = [path_out+output_file,done_dir+done]
    options = {
        'cores': 5,
        'memory': '100g',
        'account':"Trf_models",
        'walltime': "01:00:00"
    }
    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    cd {path_in}

    Rscript --vanilla {script_dir}renaming.r {input_file} {wcvp} {renaming_file} {output_file} {raw_occurrences}


    mv {output_file} {path_out}

    touch {done_dir}{done}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in, script_dir = script_dir,
                path_out = path_out, wcvp = wcvp, renaming_file = renaming_file, done = done,
                done_dir = done_dir, raw_occurrences = raw_occurrences)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

################################################################################################################################################################################
################################################################################################################################################################################
############################################--------------------- Working on tree data ---------------------####################################################################
################################################################################################################################################################################
################################################################################################################################################################################

##############################################################
#############---- Loading the SmB tree tips ----##############
##############################################################
def Load_tree(input_file, output_file, path_in,path_out, script_dir, done_dir, done):
    """Here I load the SmB tree and save it as an rds file for further analysis."""
    inputs = [path_in+input_file, done_dir+"Download_Data"]
    outputs = [path_out+output_file,done_dir+done]
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

    touch {done_dir}{done}

    '''.format(input_file=input_file, output_file=output_file, path_in = path_in,script_dir = script_dir, path_out = path_out, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
##############---- Preparing WCVP-names file ----################
##############################################################
def apg_name_align(apg,wcp, output_file, path_in, script_dir, path_out, done_dir, done):
    """Here I want to create a common format the wcvp and the previous file and update some family names to APGIV."""
    inputs = [script_dir+apg, path_in+wcp, done_dir+"Load_tree"]
    outputs = [path_out+output_file,done_dir+done]
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

    touch {done_dir}{done}

    '''.format(apg=apg,wcp = wcp, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



##############################################################
####---- Pruning the GBMB tree for tips not in WCVP ----######
##############################################################
def pruning_tree(wcp,tree, output_file, path_in, script_dir, path_out, done_dir, done):
    """Here I am looping through the tips in the GBMB tree to see if they are in the WCVP.
    If there is no exact match for the tip, I will look for a fuzzy match allowing for 1 substitution, insertion or deletion.
    If the tip name is longer than 2 IE. a subsp or a variety I will then search for just the genus and species epithet.
    Lastly I check if this has introduced any duplicate species in the tree.
    If the duplicates are sister species I will then remove one of them at random and if they are not I will remove both of them"""
    inputs = [wcp, path_in+tree, done_dir+"APG_preb_tree"]
    outputs = [path_out+output_file, done_dir+done]
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

    touch {done_dir}{done}

    '''.format(wcp = wcp, output_file=output_file, path_in = path_in, script_dir = script_dir, path_out = path_out, tree = tree, done_dir = done_dir, done = done)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
#############---- Approach 1 finding the orders ----##############
##############################################################
def Slicing_trees(input_file, output_file, path_in,path_out, script_dir, wcvp_file, apg, done_dir, done):
    """This function slices the GBMB tree into subtrees based on the order and families of the species in the tree. """
    inputs = [path_in+input_file, wcvp_file, done_dir+"Pruning_tree"]
    outputs = [path_out+output_file, done_dir+done]
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

    touch {done_dir}{done}

    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, input_file = input_file, output_file = output_file, wcvp_file = wcvp_file, apg = apg, done_dir = done_dir, done = done)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




##############################################################
###########---- Slicing the tree into orders ----############
##############################################################
def Forcing_orders(input_file_tree, output_file, path_in,path_out, script_dir, wcvp_file, apg, input_from_before, done_dir, done):
    """This function searchers for the largest monophyletic clades in the orders which are not monophyletic in the GBMB tree."""
    inputs = [path_in+input_file_tree, wcvp_file, input_from_before, done_dir+"slicing_Trees_pruning"]
    outputs = [path_out+output_file, done_dir+done]
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

    touch {done_dir}{done}
    
    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, input_file_tree = input_file_tree, output_file = output_file, wcvp_file = wcvp_file, apg = apg, done_dir = done_dir, done = done)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
##"##-- Slicing problematic orders into smaller clades --#####
##############################################################
def Creating_subclades(path_out, script_dir, done_dir, done, path_in, wcvp_file):
    """This Function cuts the problematic orders into smaller clades which ClaDs and Esse Should be able to handle"""
    inputs = [done_dir+"Finding_monophyletic_orders"]
    outputs = [done_dir+done]
    options = {
        'cores': 2,
        'memory': '5g',
        'account':"Trf_models",
        'walltime': "00:00:50"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    echo Starting the Creating subclades script
    date

    # Running the R script
    Rscript --vanilla {script_dir}Finding_monophylo_smaller_clades.r {path_in} {wcvp_file} {path_out}


    echo Ended the Creating subclades script
    date

    touch {done_dir}{done}
    
    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, wcvp_file = wcvp_file, done_dir = done_dir, done = done)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
####-- Slicing problematic families into smaller clades --####
##############################################################
def Creating_family_subclades(path_out, script_dir, done_dir, done, path_in):
    """This Function cuts the problematic orders into smaller clades which ClaDs and Esse Should be able to handle"""
    inputs = [done_dir+"Subdividing_problematic_orders"]
    outputs = [done_dir+done]
    options = {
        'cores': 2,
        'memory': '5g',
        'account':"Trf_models",
        'walltime': "00:00:50"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    echo Starting the Creating subclades script
    date

    # Running the R script
    Rscript --vanilla {script_dir}Creating_sub_phylogenies_for_big_families.r {path_in} {path_out}


    echo Ended the Creating subclades script
    date

    touch {done_dir}{done}
    
    '''.format(path_out = path_out, script_dir = script_dir, path_in = path_in, done_dir = done_dir, done = done)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

#########################################################################################################################################################################################
##########################################################################################################################################################################################
##################################################################---- Splitting the analysis into subtrees ----##########################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################

##############################################################
#####---- Finding the area of species using the wcvp  ----####
##############################################################

def Finding_areas_in_wcvp(input_file_tree, wcvp_file,path_out, output_file, path_in, order, script_dir, apg, done_dir, done, renamed_occurrences, koppen_biome):
    """This Function creates a states file for the tips in WCVP based on the climate column."""
    inputs = [done_dir+"Finding_monophyletic_orders", renamed_occurrences]
    outputs = [path_out+output_file]
    options = {
        'cores': 2,
        'memory': '100g',
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
    Rscript --vanilla {script_dir}wcvp_states.r {input_file_tree} {output_file} {wcvp_file} {path_out} {order} {apg} {renamed_occurrences} {koppen_biome}


    echo Ended the script to find state data for the tips in the wcvp
    date

    touch {done_dir}{done}

    '''.format(path_out = path_out, output_file = output_file, wcvp_file = wcvp_file, order = order,
     input_file_tree = input_file_tree, path_in = path_in, script_dir = script_dir, apg = apg, done_dir = done_dir, done = done,
     renamed_occurrences = renamed_occurrences, koppen_biome = koppen_biome)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
#########---- Finding the sampling frequency  ----############
##############################################################


def sampling_frequency(input_file_tree, wcvp_file,path_out, output_file, path_in, order, script_dir, apg, done_dir, done):
    """This function calculates the number of species sampled per genus in each subtree.
    This is then used by the ClaDs model to get a better result on speciation"""
    inputs = [done_dir+"Finding_monophyletic_orders"]
    outputs = [path_out+output_file, done_dir+done]
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

    echo Starting the script to find the sampling frequency for the tips in the wcvp 
    date

    # Running the R script
    Rscript --vanilla {script_dir}sampling_frequency.r {input_file_tree} {wcvp_file} {order} {apg} {path_out}


    echo Ended the script to find sampling frequency for the tips in the wcvp
    date

    touch {done_dir}{done}

    '''.format(path_out = path_out, output_file = output_file, wcvp_file = wcvp_file, order = order,
     input_file_tree = input_file_tree, path_in = path_in, script_dir = script_dir, apg = apg, done_dir = done_dir, done = done)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



##############################################################
#########---- Finding the sampling frequency per biome ----############
##############################################################


def sampling_frequency_per_biome(input_file_tree, wcvp_file,path_out,renamed_occurrences,koppen_biome, output_file, path_in, order, script_dir, apg, done_dir, done, percentage):
    """This function calculates the number of missing speciers missing from the tree per biome"""
    inputs = [done_dir+"Finding_monophyletic_orders", renamed_occurrences, koppen_biome]
    outputs = [path_out+output_file, done_dir+done]
    options = {
        'cores': 10,
        'memory': '100g',
        'account':"Trf_models",
        'walltime': "00:1:00"
    }

    spec = '''

    #Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    echo Starting the script to find the sampling frequency per biome 
    date

    # Running the R script
    Rscript --vanilla {script_dir}wcvp_states_esse.r {input_file_tree} {output_file} {wcvp_file} {path_out} {order} {apg} {renamed_occurrences} {koppen_biome} {percentage}


    echo Ended the script to find sampling frequency per biome
    date

    touch {done_dir}{done}

    '''.format(path_out = path_out, output_file = output_file, wcvp_file = wcvp_file, order = order, renamed_occurrences = renamed_occurrences, koppen_biome = koppen_biome,
     input_file_tree = input_file_tree, path_in = path_in, script_dir = script_dir, apg = apg, done_dir = done_dir, done = done, percentage = percentage)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


#########################################################################################
##########---- Removing tips from trees which have no distribution data  ----############
#########################################################################################


def rem_tips(input_file_tree, output_file, path_in, order, script_dir, done_dir, done, distribution_file):
    """This function finds tips in a tree which is not in the climate data and removes them from the tree."""
    inputs = [done_dir+order+"_Sampling_fraction"]
    outputs = [path_in+output_file, done_dir+done]
    options = {
        'cores': 1,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "00:10:00"
    }

    spec = '''

    # Loading Conda source folder
    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    echo Starting the script to remove tips from the tree which have no distribution data 
    date

    # Running the R script
    Rscript --vanilla {script_dir}tip_remover.r {input_file_tree} {distribution_file} {output_file}

    echo Ended the script
    date

    touch {done_dir}{done}

    '''.format(output_file = output_file, order = order, distribution_file = distribution_file,
     input_file_tree = input_file_tree, path_in = path_in, script_dir = script_dir, done_dir = done_dir, done = done)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

############################################################################
#########---- Finding the sampling frequency for Subclades  ----############
############################################################################

def sampling_frequency_subclades(input_file_tree, wcvp_file, path_out, output_file, path_in, script_dir, apg, done_dir, done, name):
    """This function calculates the number of species sampled per genus in each subtree.
    This is then used by the ClaDs model to get a better result on speciation"""
    inputs = [path_in + input_file_tree]
    outputs = [path_out + output_file, done_dir + done]
    options = {
        'cores': 2,
        'memory': '10g',
        'account': "Trf_models",
        'walltime': "00:10:00"
    }

    spec = '''

    # Checking if output dir exists
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    # Going to input folder
    cd {path_in}

    echo Starting the script to find the sampling frequency for the tips in the wcvp 
    date

    # Running the R script
    Rscript --vanilla {script_dir}sampling_frequency_subclades.r {input_file_tree} {wcvp_file} {apg} {path_out} {name}


    echo Ended the script to find sampling frequency for the tips in the wcvp
    date

    touch {done_dir}{done}

    '''.format(path_out=path_out, output_file=output_file, wcvp_file=wcvp_file,
               input_file_tree=input_file_tree, path_in=path_in, script_dir=script_dir, apg=apg,
               done_dir=done_dir, done=done, name = name)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
###########---- Runnning simple ClaDs models  ----############
##############################################################
def Clads(tree, done_file, path_in, output_file,wcvp_input, order, apg, script_dir, done_dir, sampling_frequency):
    """ """
    inputs = [path_in+tree,wcvp_input,apg,done_dir+"Finding_monophyletic_orders",done_dir+order+"_Sampling_fraction"]
    outputs = [done_dir+done_file, path_in+output_file]
    options = {
        'cores': 1,
        'memory': '1200g',
        'account':"Trf_models",
        'walltime': "168:00:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate Julia_env

    cd {path_in}

    echo Starting the Julia script at:
    date

    srun --unbuffered julia {script_dir}ClaDs.jl {path_in}{tree} {sampling_frequency} {output_file}

    echo Ended the Clads script at:
    date

    touch {done_dir}{done_file}

    '''.format(tree = tree, done_file = done_file, path_in = path_in, output_file = output_file, wcvp_input = wcvp_input, order = order,
                apg = apg, script_dir = script_dir, done_dir = done_dir, sampling_frequency = sampling_frequency)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
################---- Calculating Priors  ----#################
##############################################################
def Calculating_priors( done_file, path_in, output_file,script_dir, done_dir, input_folder):
    """ """
    inputs = [path_in + "Clads_output_Aquifoliales.R", path_in + "Clads_output_Canellales.R", path_in + "Clads_output_Cupressales.R",
                path_in + "Clads_output_Escalloniales.R", path_in + "Clads_output_Nymphaeales.R", path_in + "Clads_output_Zygophyllales.R",
                path_in + "Clads_output_Austrobaileyales.R", path_in + "Clads_output_Celastrales.R", path_in + "Clads_output_Cycadales.R",
                path_in + "Clads_output_Garryales.R", path_in + "Clads_output_Paracryphiales.R", path_in + "Clads_output_Berberidopsidales.R",
                path_in + "Clads_output_Chloranthales.R", path_in + "Clads_output_Dilleniales.R", path_in + "Clads_output_Gnetales.R", path_in + "Clads_output_Picramniales.R",
                path_in + "Clads_output_Bruniales.R", path_in + "Clads_output_Commelinales.R", path_in + "Clads_output_Dioscoreales.R", path_in + "Clads_output_Gunnerales.R",
                path_in + "Clads_output_Pinales.R", path_in + "Clads_output_Buxales.R", path_in + "Clads_output_Crossosomatales.R", path_in + "Clads_output_Dipsacales.R",
                path_in + "Clads_output_Huerteales.R", path_in + "Clads_output_Vitales.R"]
    outputs = [output_file, done_dir+done_file]
    options = {
        'cores': 1,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "00:30:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate R_env

    cd {path_in}

    echo Starting the R script at:
    date

    Rscript --vanilla {script_dir}prior_calculator.r {input_folder} {output_file}

    echo Ended the R script at:
    date

    touch {done_dir}{done_file}

    '''.format( done_file = done_file, path_in = path_in, output_file = output_file,script_dir = script_dir, done_dir = done_dir,input_folder = input_folder)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
#######---- Runnning ClaDs models with new priors  ----#######
##############################################################
def Clads_priors(tree, done_file, path_in, output_file,wcvp_input, order, apg, script_dir, done_dir, sampling_frequency, prior_file):
    """ """
    inputs = [path_in+tree,wcvp_input,apg,done_dir+"Finding_monophyletic_orders",done_dir+order+"_Sampling_fraction", prior_file]
    outputs = [done_dir+done_file, path_in+output_file]
    options = {
        'cores': 1,
        'memory': '1000g',
        'account':"Trf_models",
        'walltime': "168:00:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate Julia_env

    cd {path_in}

    echo Starting the Julia script at:
    date

    srun --unbuffered julia {script_dir}ClaDs_new_prior.jl {path_in}{tree} {sampling_frequency} {output_file} {prior_file}

    echo Ended the Clads script at:
    date

    touch {done_dir}{done_file}

    '''.format(tree = tree, done_file = done_file, path_in = path_in, output_file = output_file, wcvp_input = wcvp_input, order = order,
                apg = apg, script_dir = script_dir, done_dir = done_dir, sampling_frequency = sampling_frequency, prior_file = prior_file)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
############---- Runnning ClaDs on subclades  ----############
##############################################################
def Clads_subclades(tree, done_file, path_in, output_file,wcvp_input, order, script_dir, done_dir, sampling_frequency, prior_file):
    """ """
    inputs = [path_in+tree,done_dir+"Finding_monophyletic_orders",done_dir+order+"_Sampling_fraction", prior_file]
    outputs = [done_dir+done_file, path_in+"Clads_output_"+order+".Rdata"]
    options = {
        'cores': 10,
        'memory': '1000g',
        'account':"Trf_models",
        'walltime': "168:00:00"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate Julia_env

    cd {path_in}

    echo Starting the Julia script at:
    date

    srun --unbuffered julia {script_dir}ClaDs_new_prior_subclades.jl {path_in}{tree} {sampling_frequency} {output_file} {prior_file}

    echo Ended the Clads script at:
    date

    touch {done_dir}{done_file}

    '''.format(tree = tree, done_file = done_file, path_in = path_in, output_file = output_file, wcvp_input = wcvp_input, order = order,
                 script_dir = script_dir, done_dir = done_dir, sampling_frequency = sampling_frequency, prior_file = prior_file)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

##############################################################
######---- Calculating states files for Esse Model  ----######
##############################################################
def states_converter(path_in,tip_states_file, out_states_file, script_dir, done_dir, done, percentage_for_present):
    """ A small julia program which writes the necessary states file for a given proportion of occurrences present per biome."""
    inputs = [tip_states_file]
    outputs = [done_dir+done, out_states_file]
    options = {
        'cores': 1,
        'memory': '3g',
        'account':"Trf_models",
        'walltime': "00:00:600"
    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate Julia_env

    cd {path_in}

    echo tip_states_file {tip_states_file}
    echo out_states_file {out_states_file}
    echo percentage for present {percentage_for_present}

    echo Starting the Julia script at:
    date

    julia {script_dir}states_converter.jl {tip_states_file} {out_states_file} {percentage_for_present}

    echo Ended the Julia script at:
    date

    touch {done_dir}{done}

    '''.format(tip_states_file = tip_states_file, out_states_file = out_states_file, percentage_for_present = percentage_for_present, script_dir = script_dir,
                path_in = path_in, done_dir = done_dir, done = done)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


##############################################################
################---- Runnning ESSE model  ----################
##############################################################
def Esse(path_in, tree_file,tip_states_file,paleo_clim_file, out_states_file, out_file, hidden_states, script_dir, done_dir, done, save_file, output_folder, path_out):
    """ Function for running the ESSE model on the tree of each order. """
    inputs = [path_in+tree_file,tip_states_file,paleo_clim_file]
    outputs = [output_folder+out_file+".log", done_dir+done]
    options = {
        'cores': 10,
        'memory': '10g',
        'account':"Trf_models",
        'walltime': "168:00:00"

    }

    spec = '''

    source /home/owrisberg/miniconda3/etc/profile.d/conda.sh
    conda activate Julia_env

    # Cheking if path_out folder exists otherwise create it
    [ -d {path_out} ] && echo "{path_out} exist." || {{ echo "{path_out} does not exist."; mkdir {path_out}; }}
    
    # Cheking if output folder exists otherwise create it
    [ -d {output_folder} ] && echo "{output_folder} exist." || {{ echo "{output_folder} does not exist."; mkdir {output_folder}; }}

    cd {path_in}

    echo Starting the Julia script at:
    date
    echo using {processors} processors, {memory} gb-RAM and {hidden_states} hidden states.

    srun --unbuffered julia {script_dir}Esse.jl {processors} {tree_file} {tip_states_file} {paleo_clim_file} {out_states_file} {out_file} {save_file} {hidden_states} {output_folder}

    echo Ended the Julia script at:
    date

    touch {done_dir}{done}

    '''.format(processors = options['cores'], memory = options['memory'], tree_file = tree_file, tip_states_file = tip_states_file, paleo_clim_file = paleo_clim_file,
                out_states_file = out_states_file, out_file = out_file, hidden_states = hidden_states, script_dir = script_dir, path_in = path_in,
                done_dir = done_dir, done = done, save_file = save_file, output_folder = output_folder, path_out = path_out)


    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

################################################################################################################################################################################
################################################################################################################################################################################
######################################################################---- Running pipeline ----################################################################################
################################################################################################################################################################################
################################################################################################################################################################################

# Setting up some global variables with relative paths
script_dir = os.path.dirname(getsourcefile(download_data)) + "/"
data_dir = os.path.normpath(os.path.join(script_dir, "../data/")) +"/"
workflow_dir = os.path.normpath(os.path.join(script_dir, "../workflow/")) +"/"
done_dir = os.path.normpath(os.path.join(script_dir, "../done/")) +"/"
conda_executable_path = os.popen('which conda').read().strip()
conda_bin_dir = os.path.dirname(conda_executable_path)
conda_sh_path = os.path.normpath(os.path.join(conda_bin_dir, "../etc/profile.d/conda.sh"))

#########################################################################################################################
#############################################---- Downloading Data ----##################################################
#########################################################################################################################

gwf.target_from_template (name = "Download_Data",
                          template=download_data(
                            path_out = data_dir,
                            smb_doi = "https://github.com/FePhyFoFum/big_seed_plant_trees/releases/download/v0.1/v0.1.zip",
                            kew_doi = "http://sftp.kew.org/pub/data-repositories/WCVP/wcvp.zip",
                            gbif_doi = "https://api.gbif.org/v1/occurrence/download/request/0012129-230828120925497.zip",
                            paleo_doi = "https://zenodo.org/records/6620748/files/All_CSV_files.zip?download=1",
                            output_smb ="GBMB.tre",
                            output_kew = "wcvp_names.csv",
                            output_gbif ="occurrence.txt",
                            output_paleo="500Ma_Pohletal2022_DIB_PhaneroContinentalClimate.csv",
                            done = "Download_Data",
                            done_dir= done_dir
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
                            done = "Removing_cols",
                            done_dir = done_dir
                          ))

gwf.target_from_template(name = "Data_Parsing",
                             template = Load_data(
                                 input_file ="occurrence_select.txt",
                                 output_file = "gbif_parsed.rds",
                                 path_in = workflow_dir+"01_distribution_data/01_rm_cols/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/02_data_parsing/",
                                 done= "Data_Parsing",
                                 done_dir = done_dir
                             ))

gwf.target_from_template(name = "GBIF_lookup",
                             template = taxon_look_up(
                                 input_file ="gbif_parsed.rds",
                                 output_file = "gbif_parsed_taxon_data.rds",
                                 path_in = workflow_dir+"01_distribution_data/02_data_parsing/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/03_taxon_lookup/",
                                 done = "GBIF_lookup",
                                 done_dir = done_dir
                             ))

gwf.target_from_template(name = "Creating_Common_Format",
                             template = create_common_format(
                                 input_file_taxonomy ="gbif_parsed_taxon_data.rds",
                                 input_file_occurrences =workflow_dir+"01_distribution_data/02_data_parsing/gbif_parsed.rds",
                                 output_file = "gbif_common_format.rds",
                                 path_in = workflow_dir+"01_distribution_data/03_taxon_lookup/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/04_common_format/",
                                 done = "Creating_Common_Format",
                                 done_dir = done_dir
                             ))

gwf.target_from_template(name = "APG_preb_occurrences",
                             template = apg_name_align(
                                 wcp = "wcvp_names.csv",
                                 apg ="apgweb_parsed.csv",
                                 output_file = "wcvp_names_apg_aligned.rds",
                                 path_in = data_dir,
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/04_common_format/",
                                 done = "APG_preb_occurrences",
                                 done_dir = done_dir
                             ))

gwf.target_from_template(name = "Taxon_match",
                             template = taxon_match(
                                 input_file ="gbif_common_format.rds",
                                 wcvp = "wcvp_names_apg_aligned.rds",
                                 output_file = "gbif_taxon_matched.rds",
                                 path_in = workflow_dir+"01_distribution_data/04_common_format/",
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"01_distribution_data/05_Taxon_match/",
                                 done = "Taxon_match",
                                 done_dir = done_dir
                             ))

gwf.target_from_template(name = "Renaming",
                                template = Renamer(
                                    input_file =workflow_dir+"01_distribution_data/04_common_format/gbif_common_format.rds",
                                    wcvp = workflow_dir+"01_distribution_data/04_common_format/wcvp_names_apg_aligned.rds",
                                    renaming_file = "gbif_taxon_matched.rds",
                                    output_file = "gbif_renamed.rds",
                                    path_in = workflow_dir+"01_distribution_data/05_Taxon_match/",
                                    script_dir = script_dir,
                                    path_out = workflow_dir+"01_distribution_data/06_Renamed/",
                                    done = "Renaming",
                                    done_dir = done_dir,
                                    raw_occurrences=workflow_dir+"01_distribution_data/02_data_parsing/gbif_parsed.rds"
                                ))
################################################################################################################################
##########################------- Starting on the paleoclim data -------########################################################
################################################################################################################################

gwf.target_from_template(name = "Paleo_clim_area",
                                template = paleo_clim_area(
                                    output_file="paleoclim_area.txt",
                                    data_dir = data_dir,
                                    script_dir = script_dir,
                                    done_dir = done_dir,
                                    done = "Paleo_clim_area"
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
                            script_dir = script_dir,
                            done= "Load_tree",
                            done_dir= done_dir
                          ))

gwf.target_from_template(name = "APG_preb_tree",
                             template = apg_name_align(
                                 wcp = "wcvp_names.csv",
                                 apg ="apgweb_parsed.csv",
                                 output_file = "wcvp_names_apg_aligned.rds",
                                 path_in = data_dir,
                                 script_dir = script_dir,
                                 path_out = workflow_dir+"02_adding_orders/",
                                 done = "APG_preb_tree",
                                 done_dir= done_dir
                             ))

gwf.target_from_template(name = "Pruning_tree",
                            template=pruning_tree(
                                wcp = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                tree = "GBMB.tre",
                                output_file = "GBMB_pruned.tre",
                                path_in = data_dir,
                                path_out = data_dir,
                                script_dir = script_dir,
                                done= "Pruning_tree",
                                done_dir= done_dir
                            ))
    
gwf.target_from_template(name = "slicing_Trees_no_pruning",
                        template=Slicing_trees(
                            input_file = "GBMB.tre",
                            output_file = "GBMB_sp_per_orders_no_pruning.txt",
                            path_in = data_dir,
                            path_out = workflow_dir+"02_adding_orders/no_pruning/",
                            script_dir = script_dir,
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                            apg = script_dir+"apgweb_parsed.csv",
                            done = "slicing_Trees_no_pruning",
                            done_dir= done_dir
                            ))

gwf.target_from_template(name = "slicing_Trees_pruning",
                        template=Slicing_trees(
                            input_file = "GBMB_pruned.tre",
                            output_file = "GBMB_sp_per_orders_pruning.txt",
                            path_in = data_dir,
                            path_out = workflow_dir+"02_adding_orders/pruning/",
                            script_dir = script_dir,
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                            apg = script_dir+"apgweb_parsed.csv",
                            done = "slicing_Trees_pruning",
                            done_dir = done_dir
                            ))

gwf.target_from_template(name = "Finding_monophyletic_orders",
                        template=Forcing_orders(
                            input_file_tree = "GBMB_pruned.tre", 
                            output_file = "Orders_which_could_not_be_solved.txt", 
                            path_in = data_dir, 
                            path_out = workflow_dir+"02_adding_orders/pruning/", 
                            script_dir = script_dir, 
                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds", 
                            apg = script_dir+"apgweb_parsed.csv",
                            input_from_before=workflow_dir+"02_adding_orders/pruning/GBMB_sp_per_orders_pruning.txt",
                            done_dir= done_dir,
                            done = "Finding_monophyletic_orders"
                            ))

gwf.target_from_template(name = "Subdividing_problematic_orders",
                        template=Creating_subclades(
                            path_in = workflow_dir+"02_adding_orders/pruning/orders/", 
                            path_out = workflow_dir+"02_adding_orders/pruning/subset_of_orders/", 
                            script_dir = script_dir, 
                            wcvp_file = data_dir+"wcvp_names.csv",
                            done_dir= done_dir,
                            done = "Subdividing_problematic_orders"
                            ))

gwf.target_from_template(name = "Subdividing_problematic_families",
                        template=Creating_family_subclades(
                            path_in = workflow_dir+"02_adding_orders/pruning/subset_of_orders/", 
                            path_out = workflow_dir+"02_adding_orders/pruning/subset_of_orders/", 
                            script_dir = script_dir, 
                            done_dir= done_dir,
                            done = "Subdividing_problematic_families"
                            ))


# Remove impossible orders
impossible_orders = ["Acorales","Amborellales","Austrobaileyales","Cardiopteridales","Ceratophyllales","Cycadales","Desfontainiales","Dipsacales","Garryales","Ginkgoales","Oncothecales","Petrosaviales","Picramniales","Trochodendrales"]

# All Possible orders
orders = ["Alismatales", "Apiales", "Aquifoliales", "Arecales", "Asparagales", "Asterales","Berberidopsidales", "Boraginales", "Brassicales", "Bruniales", "Buxales", "Canellales", "Caryophyllales",
"Celastrales", "Chloranthales", "Commelinales", "Cornales", "Crossosomatales", "Cucurbitales", "Cupressales", "Dilleniales", "Dioscoreales", "Ericales", "Escalloniales", "Fabales", "Fagales",
"Gentianales", "Geraniales", "Gnetales", "Gunnerales", "Huerteales", "Icacinales", "Lamiales", "Laurales", "Liliales", "Magnoliales", "Malpighiales", "Malvales", "Metteniusales",
"Myrtales", "Nymphaeales", "Oxalidales", "Pandanales", "Paracryphiales", "Pinales", "Piperales", "Poales", "Proteales", "Ranunculales", "Rosales", "Santalales", "Sapindales",
"Saxifragales", "Solanales", "Vahliales", "Vitales", "Zingiberales", "Zygophyllales"
]

# Percentage of occurences required to be present in a biome.
percentages = ["0.33"]

#####################################################################################################################################################################
########################################################--- ClaDs on Orders with uniform prior  ---##################################################################
#####################################################################################################################################################################

# Orders that ran with the uniform prior 
orders_not_in_orders_new_prior = ["Aquifoliales", "Berberidopsidales", "Boraginales", "Bruniales", "Buxales",
                                "Canellales", "Celastrales", "Chloranthales", "Commelinales", "Cornales", "Crossosomatales", "Cucurbitales",
                                "Cupressales", "Dilleniales", "Dioscoreales", "Escalloniales", "Fagales", "Gunnerales", "Huerteales", "Icacinales",
                                "Liliales", "Magnoliales", "Metteniusales", "Nymphaeales", "Oxalidales", "Pandanales", "Paracryphiales",
                                "Pinales", "Piperales", "Proteales", "Santalales", "Vitales", "Zygophyllales"]         

for i in range(len(orders_not_in_orders_new_prior)):

    # Calculating the sampling fraction
    gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_samplingfraction",
                                      template = sampling_frequency(
                                          input_file_tree = orders_not_in_orders_new_prior[i],
                                          path_in = workflow_dir + "02_adding_orders/pruning/orders/",
                                          path_out = workflow_dir + "03_distribution_data/",
                                          output_file = orders_not_in_orders_new_prior[i] + "_sampling_fraction._ClaDs.txt",
                                          wcvp_file = workflow_dir + "02_adding_orders/wcvp_names_apg_aligned.rds",
                                          script_dir = script_dir,
                                          order = orders_not_in_orders_new_prior[i],
                                          apg = script_dir + "apgweb_parsed.csv",
                                          done_dir = done_dir,
                                          done = orders_not_in_orders_new_prior[i] + "_Sampling_fraction_ClaDs",
                                      ))

    # Running ClaDs model
    gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_ClaDs",
                                template= Clads(
                                tree= "pruned_tree_order_"+orders_not_in_orders_new_prior[i]+"_GBMB.tre",
                                wcvp_input = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                order = orders_not_in_orders_new_prior[i],
                                apg = script_dir+"apgweb_parsed.csv",
                                done_file = orders_not_in_orders_new_prior[i]+"_ClaDs",
                                path_in = workflow_dir+"02_adding_orders/pruning/orders/",
                                output_file = "Clads_output_"+orders_not_in_orders_new_prior[i]+".jld2",
                                script_dir=script_dir,
                                done_dir = done_dir,
                                sampling_frequency= workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_sampling_fraction.txt"
                             ))
    
    #Running the script to find the environmental data for the tips in the trees
    gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_distribution_data_ClaDs.",
                                template=Finding_areas_in_wcvp(
                                input_file_tree= "family_phylo_"+orders_not_in_orders_new_prior[i]+".tre", # 
                                path_in =  workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                path_out = workflow_dir+"03_distribution_data/",
                                output_file = orders_not_in_orders_new_prior[i]+"_distribution_data_ClaDs.txt",
                                wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                order = orders_not_in_orders_new_prior[i],
                                script_dir= script_dir,
                                apg = script_dir+"apgweb_parsed.csv",
                                done_dir= done_dir,
                                done= orders_not_in_orders_new_prior[i]+"_distribution_data",
                                renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                koppen_biome = script_dir+"koppen_geiger_0p01.tif"
                                ))
    
    # Running the script to find the biome of the tips in the tree
    for j in range(len(percentages)):
        gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_states_converter_ClaDs"+percentages[j],
                                template=states_converter(
                                path_in= workflow_dir+"03_distribution_data/",
                                tip_states_file= workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_distribution_data_ClaDs.txt",
                                out_states_file= workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_states_"+percentages[j]+"_ClaDs_.txt",
                                script_dir= script_dir,
                                done_dir= done_dir,
                                done= "States_converter_"+orders_not_in_orders_new_prior[i]+"_"+percentages[j]+"_ClaDs",
                                percentage_for_present= percentages[j]
                                ))
    
#####################################################################################################################################################################
########################################################--- ClaDs on Orders with calculated Prior ---#################################################################
#####################################################################################################################################################################

# Orders that need a to be run with a modified prior
orders_new_prior = ["Solanales"]
    
gwf.target_from_template(name = "Calculating_priors",
                            template=Calculating_priors(
                            done_file = "Calculating_priors",
                            path_in = workflow_dir + "02_adding_orders/pruning/orders/",
                            output_file = workflow_dir+"02_adding_orders/pruning/orders/priors.txt",
                            script_dir = script_dir,
                            done_dir = done_dir,
                            input_folder = workflow_dir + "02_adding_orders/pruning/orders"
                            ))

for i in range(len(orders_new_prior)):
    gwf.target_from_template(name = orders_new_prior[i]+"_ClaDs",
                                template= Clads_priors(
                                tree= "pruned_tree_order_"+orders_new_prior[i]+"_GBMB.tre",
                                wcvp_input = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                order = orders_new_prior[i],
                                apg = script_dir+"apgweb_parsed.csv",
                                done_file = orders_new_prior[i]+"_ClaDs",
                                path_in = workflow_dir+"02_adding_orders/pruning/orders/",
                                output_file = "Clads_output_"+orders_new_prior[i]+".jld2",
                                script_dir=script_dir,
                                done_dir = done_dir,
                                sampling_frequency= workflow_dir+"03_distribution_data/"+orders_new_prior[i]+"_sampling_fraction.txt",
                                prior_file = workflow_dir+"02_adding_orders/pruning/orders/priors.txt"
                             ))
    
#####################################################################################################################################################################
############################################################--- ClaDs on Order subclades ---#########################################################################
#####################################################################################################################################################################

# Here is the list of orders that I have had to split into smaller clades to be able to run the ClaDs on them.
# orders_split_for_clads = ["Zingiberales", "Laurales","Poales","Ranunculales","Rosales","Sapindales","Saxifragales","Myrtales","Malvales","Malpighiales","Lamiales","Gentianales","Fabales",
# "Ericales", "Apiales","Asterales","Asparagales","Caryophyllales","Arecales","Brassicales"]

# This is the list of Clades resulting from my personal splitting of the orders.
Clads_clades = ["Aizoaceae_Phytolaccaceae_Barbeuiaceae_Lophiocarpaceae_Gisekiaceae_Sarcobataceae", "Alzateaceae_Crypteroniaceae_Penaeaceae",
                "Araliaceae","Arecaceae","Asphodelaceae","Balsaminaceae_Marcgraviaceae_Tetrameristaceae","Berberidaceae",
                "Bignoniaceae","Cactaceae_Molluginaceae_Didiereaceae_Anacompserotaceae_Basellaceae_Montiaceae_Halophytaceae_Portulacaceae_Talinaceae", "Calyceraceae",
                "Campanulaceae_Rousseaceae","Cannabaceae","Capparaceae","Cercidiphyllaceae_Hamamelidaceae_Daphniphyllaceae_Altingiaceae_Paeoniaceae",
                "Chrysobalanaceae_Malpighiaceae_Caryocaraceae_Balanopaceae_Elatinaceae_Centroplacaceae_Dichapetalaceae_Putranjivaceae_Euphroniaceae_Lophopyxidaceae_Trigoniaceae","Cleomaceae",
                "Combretaceae","Costaceae","Crassulaceae_Aphanopetalaceae_Halograceae_Penthoraceae_Tetracarpaeaceae","Dipterocarpaceae_Bixaceae_Cistaceae_Sarcoleanaceae_Muntingiaceae_Sphaerosepalaceae",
                "Droseraceae_Ancistrocladaceae_Drosophyllaceae_Nepenthaceae_Dioncophyllaceae","Ebenaceae","Ericaceae_Clethraceae_Cyrillaceae","Gentianaceae",
                "Goodeniaceae","Heliconiaceae_Lowiaceae_Strelitziaceae","Juncaceae","Loganiaceae_Gelsemiaceae","Lythraceae_Onagraceae","Malvaceae","Marantaceae_Cannaceae",
                "Meliaceae","Menispermaceae","Menyanthaceae","Monimiaceae","Moraceae",
                "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae",
                "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae","Papaveraceae","Passifloraceae","Pentaphylacaceae_Sladeniaceae",
                "Pittosporaceae","Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae","Polemoniaceae_Lecythidaceae_Fouquieriaceae","Polygalaceae_Surianaceae",
                "Primulaceae","Resedaceae","Restionaceae","Rhamnaceae_Barbeyaceae_Dirachmaceae_Elaeagnaceae","Rutaceae","Salicaceae_Lacistemataceae","Sapindaceae",
                "Sapotaceae","Saxifragaceae_Iteaceae_Grossulariaceae","Scrophulariaceae","Simaroubaceae","Styracaceae_Diapensiaceae_Symplocaceae","Theaceae","Thymelaeaceae","Typhaceae","Ulmaceae","Urticaceae",
                "Verbenaceae_Schlegeliaceae_Lentibulariaceae_Thomandersiaceae","Violaceae_Goupiaceae","Xyridaceae_Eriocaulaceae","Zingiberaceae"]

for i in range(len(Clads_clades)):
    gwf.target_from_template(name = Clads_clades[i]+"_samplingfraction",
                                      template = sampling_frequency_subclades(
                                          input_file_tree = "family_phylo_" + Clads_clades[i] + ".tre",
                                          path_in = workflow_dir + "02_adding_orders/pruning/subset_of_orders/",
                                          path_out = workflow_dir + "03_distribution_data/",
                                          output_file = Clads_clades[i] + "_sampling_fraction.txt",
                                          wcvp_file = workflow_dir + "02_adding_orders/wcvp_names_apg_aligned.rds",
                                          script_dir = script_dir,
                                          apg = script_dir + "apgweb_parsed.csv",
                                          done_dir = done_dir,
                                          done = Clads_clades[i] + "_Sampling_fraction",
                                          name = Clads_clades[i]
                                      ))

    gwf.target_from_template(name = Clads_clades[i]+"_ClaDs",
                                template= Clads_subclades(
                                    tree= "family_phylo_"+Clads_clades[i]+".tre",
                                    wcvp_input = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                    order = Clads_clades[i],
                                    done_file = Clads_clades[i]+"_ClaDs",
                                    path_in = workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                    output_file = "Clads_output_"+Clads_clades[i]+".jld2",
                                    script_dir=script_dir,
                                    done_dir = done_dir,
                                    sampling_frequency= workflow_dir+"03_distribution_data/"+Clads_clades[i]+"_sampling_fraction.txt",
                                    prior_file = workflow_dir+"02_adding_orders/pruning/orders/priors.txt"
                             ))

    #Running the script to find the environmental data for the tips in the sub trees
    gwf.target_from_template(name = Clads_clades[i]+"_distribution_data.",
                                template=Finding_areas_in_wcvp(
                                input_file_tree= "family_phylo_"+Clads_clades[i]+".tre", # 
                                path_in =  workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                path_out = workflow_dir+"03_distribution_data/",
                                output_file = Clads_clades[i]+"_distribution_data.txt",
                                wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                order = Clads_clades[i],
                                script_dir= script_dir,
                                apg = script_dir+"apgweb_parsed.csv",
                                done_dir= done_dir,
                                done= Clads_clades[i]+"_distribution_data",
                                renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                koppen_biome = script_dir+"koppen_geiger_0p01.tif"
                                ))
    
    
    for j in range(len(percentages)):
        gwf.target_from_template(name = Clads_clades[i]+"_states_converter_"+percentages[j],
                                template=states_converter(
                                path_in= workflow_dir+"03_distribution_data/",
                                tip_states_file= workflow_dir+"03_distribution_data/"+Clads_clades[i]+"_distribution_data.txt",
                                out_states_file= workflow_dir+"03_distribution_data/"+Clads_clades[i]+"_states_"+percentages[j]+".txt",
                                script_dir= script_dir,
                                done_dir= done_dir,
                                done= "States_converter_"+Clads_clades[i]+"_"+percentages[j]+"",
                                percentage_for_present= percentages[j]
                                ))

#####################################################################################################################################################################
############################################################--- ClaDs on family  subclades ---#########################################################################
#####################################################################################################################################################################

# Families which have been divided further and therefore need to be removed from the Clads_clades list.
# Amaryllidaceae, Anacardiaceae_Burseraceae_Kirkiaceae, Apiaceae, Apocynaceae, Asparagaceae, Asteraceae,Brassicaceae, Euphorbiaceae, Fabaceae, Lamiaceae, Orchidaceae, Poaceae, Rubiaceae
# Rosaceae, Ranunculaceae, Plantaginaceae, Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae,Myrtaceae, Lauraceae, Iridaceae, Gesneriaceae_Calceolariaceae, Cyperaceae,
# Caryophyllaceae_Achatocarpaceae_Amaranthaceae, ,Bromeliaceae, Acanthaceae_Martyniaceae_Pedaliaceae, Melastomataceae,

sub_family_clades = ["sub_phylo_Acanthaceae_Martyniaceae_Pedaliaceae_1.tre","sub_phylo_Acanthaceae_Martyniaceae_Pedaliaceae_2.tre","sub_phylo_Acanthaceae_Martyniaceae_Pedaliaceae_3.tre",
                     "sub_phylo_Amaryllidaceae_1.tre","sub_phylo_Amaryllidaceae_2.tre","sub_phylo_Anacardiaceae_Burseraceae_Kirkiaceae_1.tre","sub_phylo_Anacardiaceae_Burseraceae_Kirkiaceae_2.tre",
                     "sub_phylo_Apiaceae_1.tre","sub_phylo_Apiaceae_2.tre","sub_phylo_Apiaceae_3.tre","sub_phylo_Apiaceae_4.tre","sub_phylo_Apiaceae_5.tre","sub_phylo_Apiaceae_6.tre",
                     "sub_phylo_Apocynaceae_1.tre","sub_phylo_Apocynaceae_2.tre","sub_phylo_Apocynaceae_3.tre","sub_phylo_Apocynaceae_4.tre","sub_phylo_Asparagaceae_1.tre","sub_phylo_Asparagaceae_2.tre",
                     "sub_phylo_Asparagaceae_3.tre","sub_phylo_Asteraceae_10.tre","sub_phylo_Asteraceae_11.tre","sub_phylo_Asteraceae_12.tre","sub_phylo_Asteraceae_13.tre","sub_phylo_Asteraceae_1.tre",
                     "sub_phylo_Asteraceae_2.tre","sub_phylo_Asteraceae_3.tre","sub_phylo_Asteraceae_4.tre","sub_phylo_Asteraceae_5.tre","sub_phylo_Asteraceae_6.tre","sub_phylo_Asteraceae_7.tre","sub_phylo_Asteraceae_8.tre",
                     "sub_phylo_Asteraceae_9.tre","sub_phylo_Brassicaceae_1.tre","sub_phylo_Brassicaceae_2.tre","sub_phylo_Brassicaceae_3.tre","sub_phylo_Brassicaceae_4.tre","sub_phylo_Bromeliaceae_1.tre","sub_phylo_Bromeliaceae_2.tre",
                     "sub_phylo_Caryophyllaceae_Achatocarpaceae_Amaranthaceae_1.tre","sub_phylo_Caryophyllaceae_Achatocarpaceae_Amaranthaceae_2.tre","sub_phylo_Cyperaceae_1.tre","sub_phylo_Cyperaceae_2.tre","sub_phylo_Cyperaceae_3.tre",
                     "sub_phylo_Euphorbiaceae_1.tre","sub_phylo_Euphorbiaceae_2.tre","sub_phylo_Euphorbiaceae_3.tre","sub_phylo_Fabaceae_10.tre","sub_phylo_Fabaceae_11.tre","sub_phylo_Fabaceae_12.tre","sub_phylo_Fabaceae_13.tre","sub_phylo_Fabaceae_1.tre",
                     "sub_phylo_Fabaceae_2.tre","sub_phylo_Fabaceae_3.tre","sub_phylo_Fabaceae_4.tre","sub_phylo_Fabaceae_5.tre","sub_phylo_Fabaceae_6.tre","sub_phylo_Fabaceae_7.tre","sub_phylo_Fabaceae_8.tre","sub_phylo_Fabaceae_9.tre",
                     "sub_phylo_Gesneriaceae_Calceolariaceae_1.tre","sub_phylo_Gesneriaceae_Calceolariaceae_2.tre","sub_phylo_Iridaceae_1.tre","sub_phylo_Iridaceae_2.tre","sub_phylo_Lamiaceae_1.tre","sub_phylo_Lamiaceae_2.tre","sub_phylo_Lamiaceae_3.tre",
                     "sub_phylo_Lamiaceae_4.tre","sub_phylo_Lamiaceae_5.tre","sub_phylo_Lamiaceae_6.tre","sub_phylo_Lamiaceae_7.tre","sub_phylo_Lauraceae_1.tre","sub_phylo_Lauraceae_2.tre","sub_phylo_Melastomataceae_1.tre",
                     "sub_phylo_Melastomataceae_2.tre","sub_phylo_Myrtaceae_1.tre","sub_phylo_Myrtaceae_2.tre","sub_phylo_Myrtaceae_3.tre","sub_phylo_Orchidaceae_10.tre","sub_phylo_Orchidaceae_11.tre","sub_phylo_Orchidaceae_12.tre",
                     "sub_phylo_Orchidaceae_13.tre","sub_phylo_Orchidaceae_14.tre","sub_phylo_Orchidaceae_15.tre","sub_phylo_Orchidaceae_16.tre","sub_phylo_Orchidaceae_1.tre","sub_phylo_Orchidaceae_2.tre","sub_phylo_Orchidaceae_3.tre",
                     "sub_phylo_Orchidaceae_4.tre","sub_phylo_Orchidaceae_5.tre","sub_phylo_Orchidaceae_6.tre","sub_phylo_Orchidaceae_7.tre","sub_phylo_Orchidaceae_8.tre","sub_phylo_Orchidaceae_9.tre",
                     "sub_phylo_Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_1.tre","sub_phylo_Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_2.tre","sub_phylo_Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_3.tre",
                     "sub_phylo_Plantaginaceae_1.tre","sub_phylo_Plantaginaceae_2.tre","sub_phylo_Plantaginaceae_3.tre","sub_phylo_Poaceae_1.tre","sub_phylo_Poaceae_2.tre","sub_phylo_Poaceae_3.tre","sub_phylo_Poaceae_4.tre",
                     "sub_phylo_Poaceae_5.tre","sub_phylo_Poaceae_6.tre","sub_phylo_Poaceae_7.tre","sub_phylo_Ranunculaceae_1.tre","sub_phylo_Ranunculaceae_2.tre","sub_phylo_Ranunculaceae_3.tre","sub_phylo_Rosaceae_1.tre","sub_phylo_Rosaceae_2.tre",
                     "sub_phylo_Rosaceae_3.tre","sub_phylo_Rosaceae_4.tre","sub_phylo_Rosaceae_5.tre","sub_phylo_Rubiaceae_1.tre","sub_phylo_Rubiaceae_2.tre","sub_phylo_Rubiaceae_3.tre",]


# Orders that are removed as they are too small: Alismatales, Amborellales, Petrosaviales, Trochodendrales, Vahliales

for k in range(len(sub_family_clades)):
    gwf.target_from_template(name = sub_family_clades[k]+"_samplingfraction",
                                      template = sampling_frequency_subclades(
                                          input_file_tree = sub_family_clades[k],
                                          path_in = workflow_dir + "02_adding_orders/pruning/subset_of_orders/",
                                          path_out = workflow_dir + "03_distribution_data/",
                                          output_file = sub_family_clades[k] + "_sampling_fraction.txt",
                                          wcvp_file = workflow_dir + "02_adding_orders/wcvp_names_apg_aligned.rds",
                                          script_dir = script_dir,
                                          apg = script_dir + "apgweb_parsed.csv",
                                          done_dir = done_dir,
                                          done = sub_family_clades[k] + "_Sampling_fraction",
                                          name = sub_family_clades[k]
                                      ))
    

    gwf.target_from_template(name = sub_family_clades[k]+"_ClaDs",
                                template= Clads_subclades(
                                    tree=sub_family_clades[k],
                                    wcvp_input = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                    order = sub_family_clades[k],
                                    done_file = sub_family_clades[k]+"_ClaDs",
                                    path_in = workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                    output_file = "Clads_output_"+sub_family_clades[k]+".jld2",
                                    script_dir=script_dir,
                                    done_dir = done_dir,
                                    sampling_frequency= workflow_dir+"03_distribution_data/"+sub_family_clades[k]+"_sampling_fraction.txt",
                                    prior_file = workflow_dir+"02_adding_orders/pruning/orders/priors.txt"
                             )) 

 #Running the script to find the environmental data for the tips in the sub trees
    gwf.target_from_template(name = sub_family_clades[k]+"_distribution_data.",
                                template=Finding_areas_in_wcvp(
                                input_file_tree= sub_family_clades[k], 
                                path_in =  workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                path_out = workflow_dir+"03_distribution_data/",
                                output_file = sub_family_clades[k]+"_distribution_data.txt",
                                wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                order = sub_family_clades[k],
                                script_dir= script_dir,
                                apg = script_dir+"apgweb_parsed.csv",
                                done_dir= done_dir,
                                done= sub_family_clades[k]+"_distribution_data",
                                renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                koppen_biome = script_dir+"koppen_geiger_0p01.tif"
                                ))
    
    
    for j in range(len(percentages)):
        gwf.target_from_template(name = sub_family_clades[k]+"_states_converter_"+percentages[j],
                                template=states_converter(
                                path_in= workflow_dir+"03_distribution_data/",
                                tip_states_file= workflow_dir+"03_distribution_data/"+sub_family_clades[k]+"_distribution_data.txt",
                                out_states_file= workflow_dir+"03_distribution_data/"+sub_family_clades[k]+"_states_"+percentages[j]+".txt",
                                script_dir= script_dir,
                                done_dir= done_dir,
                                done= "States_converter_"+sub_family_clades[k]+"_"+percentages[j]+"",
                                percentage_for_present= percentages[j]
                                ))

#####################################################################################################################################################################
############################################################--- ESSE on family  subclades ---########################################################################
#####################################################################################################################################################################
# All Possible orders
orders = ["Alismatales", "Apiales", "Aquifoliales", "Arecales", "Asparagales", "Asterales","Berberidopsidales", "Boraginales", "Brassicales", "Bruniales", "Buxales", "Canellales", "Caryophyllales",
"Celastrales", "Chloranthales", "Commelinales", "Cornales", "Crossosomatales", "Cucurbitales", "Cupressales", "Dilleniales", "Dioscoreales", "Ericales", "Escalloniales", "Fabales", "Fagales",
"Gentianales", "Geraniales", "Gnetales", "Gunnerales", "Huerteales", "Icacinales", "Lamiales", "Laurales", "Liliales", "Magnoliales", "Malpighiales", "Malvales", "Metteniusales",
"Myrtales", "Nymphaeales", "Oxalidales", "Pandanales", "Paracryphiales", "Pinales", "Piperales", "Poales", "Proteales", "Ranunculales", "Rosales", "Santalales", "Sapindales",
"Saxifragales", "Solanales", "Vahliales", "Vitales", "Zingiberales", "Zygophyllales"
]

###############################################################################################################################################################
#################################################################--- Esse  Runs ---############################################################################
###############################################################################################################################################################

###############################################################################################################################################################
#############################################---  Running ESSE on the orders which can finish in a week ---####################################################
###############################################################################################################################################################
# Orders that ran with the uniform prior 
orders_not_in_orders_new_prior = ["Arecales","Buxales","Canellales","Celastrales","Commelinales", "Cornales", "Crossosomatales",
                                "Cupressales", "Escalloniales", "Fagales","Geraniales","Gnetales","Huerteales","Icacinales","Malvales","Metteniusales", "Nymphaeales","Oxalidales", "Pandanales",
                                "Pinales", "Piperales", "Proteales", "Zygophyllales"]

for i in range(len(orders_not_in_orders_new_prior)):
            gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_distribution_data_Esse.",
                                template=Finding_areas_in_wcvp(
                                input_file_tree= "pruned_tree_order_"+orders_not_in_orders_new_prior[i]+"_GBMB.tre",
                                path_in =  workflow_dir+"02_adding_orders/pruning/orders/",
                                path_out = workflow_dir+"03_distribution_data/",
                                output_file = orders_not_in_orders_new_prior[i]+"_distribution_data_Esse.txt",
                                wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                order = orders_not_in_orders_new_prior[i],
                                script_dir= script_dir,
                                apg = script_dir+"apgweb_parsed.csv",
                                done_dir= done_dir,
                                done= orders_not_in_orders_new_prior[i]+"_distribution_data",
                                renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                koppen_biome = script_dir+"koppen_geiger_0p01.tif"
                                ))
            
            for j in range(len(percentages)):

                gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_states_converter_Esse_"+percentages[j],
                                        template=states_converter(
                                        path_in= workflow_dir+"03_distribution_data/",
                                        tip_states_file= workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_distribution_data_Esse.txt",
                                        out_states_file= workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_states_"+percentages[j]+"_Esse_.txt",
                                        script_dir= script_dir,
                                        done_dir= done_dir,
                                        done= "States_converter_"+orders_not_in_orders_new_prior[i]+"_"+percentages[j]+"_Esse",
                                        percentage_for_present= percentages[j]
                                        ))
        
                gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_samplingfraction_Esse",
                                     template = sampling_frequency(
                                            input_file_tree= "pruned_tree_order_"+orders_not_in_orders_new_prior[i]+"_GBMB.tre",
                                            path_in =  workflow_dir+"02_adding_orders/pruning/orders/",
                                            path_out = workflow_dir+"03_distribution_data/",
                                            output_file = orders_not_in_orders_new_prior[i]+"_sampling_fraction_Esse.txt",
                                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                            order = orders_not_in_orders_new_prior[i],
                                            script_dir= script_dir,
                                            apg = script_dir+"apgweb_parsed.csv",
                                            done_dir= done_dir,
                                            done= orders_not_in_orders_new_prior[i]+"_Sampling_fraction_Esse"
                                    ))
        
                gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_Tip_removal",
                                            template = rem_tips(
                                            input_file_tree = "pruned_tree_order_"+orders_not_in_orders_new_prior[i]+"_GBMB.tre",
                                            distribution_file= workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_states_"+percentages[j]+".txt",
                                            output_file = orders_not_in_orders_new_prior[i]+"_Esse_tree.tre",
                                            path_in = workflow_dir+"02_adding_orders/pruning/orders/",
                                            order = orders_not_in_orders_new_prior[i],
                                            script_dir= script_dir,
                                            done_dir= done_dir,
                                            done= "Tip_removal"+orders_not_in_orders_new_prior[i]
                                            ))
                
                gwf.target_from_template(name = orders[i]+"_Biome_sampling_fraction.",
                                            template=sampling_frequency_per_biome(
                                            input_file_tree= "pruned_tree_order_"+orders[i]+"_GBMB.tre",
                                            path_in =  workflow_dir+"02_adding_orders/pruning/orders/",
                                            path_out = workflow_dir+"03_distribution_data/",
                                            output_file = orders[i]+"_biome_sampling_fraction.txt",
                                            wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                            order = orders[i],
                                            script_dir= script_dir,
                                            apg = script_dir+"apgweb_parsed.csv",
                                            done_dir= done_dir,
                                            done= orders[i]+"_biome_sampling_fraction",
                                            renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                            koppen_biome = script_dir+"koppen_geiger_0p01.tif",
                                            percentage= percentages[j]
                                            ))

                gwf.target_from_template(name = orders_not_in_orders_new_prior[i]+"_Esse",
                                            template = Esse(
                                            tree_file = orders_not_in_orders_new_prior[i]+"_Esse_tree.tre", # Input tree
                                            tip_states_file = workflow_dir+"03_distribution_data/"+orders_not_in_orders_new_prior[i]+"_states_"+percentages[j]+".txt", 
                                            paleo_clim_file = data_dir+"paleoclim_area.txt", # File with paleoclimatic variables
                                            done = orders_not_in_orders_new_prior[i]+"_Esse",
                                            path_in = workflow_dir+"02_adding_orders/pruning/orders/",
                                            save_file = "Esse_output_"+orders_not_in_orders_new_prior[i]+"_"+percentages[j]+".jld2",
                                            script_dir=script_dir,
                                            done_dir = done_dir,
                                            out_states_file = "Esse_states_"+orders_not_in_orders_new_prior[i]+"_"+percentages[j], 
                                            hidden_states = 0,
                                            out_file = "Esse_output_"+orders_not_in_orders_new_prior[i]+"_hidden_states_"+percentages[j],
                                            output_folder = workflow_dir+"04_results/Esse_output/",
                                            path_out = workflow_dir+"04_results/"
                                         ))
        
###############################################################################################################################################################
#######################################################---  Running ESSE on the esse_clades ---################################################################
###############################################################################################################################################################

# This is the list of Clades resulting from my personal splitting of the orders. ( Here I need to remove)
esse_clades = ["Aizoaceae_Phytolaccaceae_Barbeuiaceae_Lophiocarpaceae_Gisekiaceae_Sarcobataceae", # Caryophyllales
                "Alzateaceae_Crypteroniaceae_Penaeaceae", # Myrtales
                "Araliaceae", # Apiales
                "Asphodelaceae", # Asparagales
                "Balsaminaceae_Marcgraviaceae_Tetrameristaceae", # Ericales
                "Berberidaceae", # Ranunculales
                "Bignoniaceae", # Lamiales
                "Cactaceae_Molluginaceae_Didiereaceae_Anacompserotaceae_Basellaceae_Montiaceae_Halophytaceae_Portulacaceae_Talinaceae", # Caryophyllales
                "Calyceraceae",  # Asterales
                "Campanulaceae_Rousseaceae", # Asterales
                "Cannabaceae", # Rosales
                "Capparaceae", # Brassicales
                "Cercidiphyllaceae_Hamamelidaceae_Daphniphyllaceae_Altingiaceae_Paeoniaceae", # Saxifragales
                "Chrysobalanaceae_Malpighiaceae_Caryocaraceae_Balanopaceae_Elatinaceae_Centroplacaceae_Dichapetalaceae_Putranjivaceae_Euphroniaceae_Lophopyxidaceae_Trigoniaceae", # Malpighiales
                "Cleomaceae", # Brassicales
                "Combretaceae", # Myrtales
                "Crassulaceae_Aphanopetalaceae_Halograceae_Penthoraceae_Tetracarpaeaceae", # Saxifragales
                "Dipterocarpaceae_Bixaceae_Cistaceae_Sarcoleanaceae_Muntingiaceae_Sphaerosepalaceae", # Malvales
                "Droseraceae_Ancistrocladaceae_Drosophyllaceae_Nepenthaceae_Dioncophyllaceae", # Caryophyllales
                "Ebenaceae", # Ericales
                "Ericaceae_Clethraceae_Cyrillaceae", # Ericales
                "Gentianaceae", # Gentianales
                "Goodeniaceae", # Asterales
                "Juncaceae", # Poales
                "Loganiaceae_Gelsemiaceae", # Gentianales
                "Lythraceae_Onagraceae", # Myrtales
                "Malvaceae", # Malvales
                "Meliaceae", # Sapindales
                "Menispermaceae", # Ranunculales
                "Menyanthaceae", # Asterales
                "Monimiaceae", # Laurales
                "Moraceae", # Rosales
                "Ochnaceae_Clusiaceae_Erythroxylaceae_Podostemaceae_Bonnetiaceae_Rhizophoraceae_Calophyllaceae_Hypericaceae_Ctenolophonaceae_Irvingiaceae_Pandaceae", # Malpighiales
                "Orobanchaceae_Phrymaceae_Mazaceae_Paulowniaceae", # Lamiales
                "Papaveraceae", # Ranunculales
                "Passifloraceae", # Malpighiales
                "Pentaphylacaceae_Sladeniaceae", # Ericales
                "Pittosporaceae", # Apiales
                "Plumbaginaceae_Polygonaceae_Frankeniaceae_Tamaricaceae", # Caryophyllales
                "Polemoniaceae_Lecythidaceae_Fouquieriaceae", # Ericales
                "Polygalaceae_Surianaceae", # Fabales
                "Primulaceae", # Ericales
                "Resedaceae", # Brassicales
                "Restionaceae", # Poales
                "Rhamnaceae_Barbeyaceae_Dirachmaceae_Elaeagnaceae", # Rosales
                "Rutaceae", # Sapindales
                "Salicaceae_Lacistemataceae", # Malpighiales
                "Sapindaceae", # Sapindales
                "Sapotaceae", # Ericales
                "Saxifragaceae_Iteaceae_Grossulariaceae", # Saxifragales
                "Scrophulariaceae", # Lamiales
                "Simaroubaceae", # Sapindales
                "Styracaceae_Diapensiaceae_Symplocaceae", # Ericales
                "Theaceae", # Ericales
                "Thymelaeaceae", # Malvales
                "Typhaceae", # Poales
                "Ulmaceae", # Rosales
                "Urticaceae", # Rosales
                "Verbenaceae_Schlegeliaceae_Lentibulariaceae_Thomandersiaceae", # Lamiales
                "Violaceae_Goupiaceae", # Malpighiales
                "Xyridaceae_Eriocaulaceae", # Poales
                ]

for i in range(len(esse_clades)):
                # gwf.target_from_template(name = esse_clades[i]+"_distribution_data.",
                #                             template=Finding_areas_in_wcvp(
                #                             input_file_tree= "pruned_tree_family_"+esse_clades[i]+"_GBMB.tre",
                #                             path_in =  workflow_dir+"02_adding_orders/pruning/families/",
                #                             path_out = workflow_dir+"03_distribution_data/",
                #                             output_file = esse_clades[i]+"_distribution_data.txt",
                #                             wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                #                             order = esse_clades[i],
                #                             script_dir= script_dir,
                #                             apg = script_dir+"apgweb_parsed.csv",
                #                             done_dir= done_dir,
                #                             done= esse_clades[i]+"_distribution_data",
                #                             renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                #                             koppen_biome = script_dir+"koppen_geiger_0p01.tif"
                #                             ))

                # for j in range(len(percentages)):
                #     gwf.target_from_template(name = esse_clades[i]+"_states_converter_"+percentages[j],
                #                             template=states_converter(
                #                             path_in= workflow_dir+"03_distribution_data/",
                #                             tip_states_file= workflow_dir+"03_distribution_data/"+esse_clades[i]+"_distribution_data.txt",
                #                             out_states_file= workflow_dir+"03_distribution_data/"+esse_clades[i]+"_states_"+percentages[j]+".txt",
                #                             script_dir= script_dir,
                #                             done_dir= done_dir,
                #                             done= "States_converter_"+esse_clades[i]+"_"+percentages[j]+"",
                #                             percentage_for_present= percentages[j]
                #                             ))

                #     gwf.target_from_template(name = esse_clades[i]+"_Sampling_fraction",
                #                          template = sampling_frequency(
                #                                 input_file_tree= "pruned_tree_family_"+esse_clades[i]+"_GBMB.tre",
                #                                 path_in =  workflow_dir+"02_adding_orders/pruning/families/",
                #                                 path_out = workflow_dir+"03_distribution_data/",
                #                                 output_file = esse_clades[i]+"_sampling_fraction.txt",
                #                                 wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                #                                 order = esse_clades[i],
                #                                 script_dir= script_dir,
                #                                 apg = script_dir+"apgweb_parsed.csv",
                #                                 done_dir= done_dir,
                #                                 done= esse_clades[i]+"_Sampling_fraction"
                #                          ))


                    gwf.target_from_template(name = esse_clades[i]+"_Tip_removal",
                                                template = rem_tips(
                                                input_file_tree = "pruned_tree_family_"+esse_clades[i]+"_GBMB.tre",
                                                distribution_file= workflow_dir+"03_distribution_data/"+esse_clades[i]+"_states_"+percentages[j]+".txt",
                                                output_file = esse_clades[i]+"_Esse_tree.tre",
                                                path_in = workflow_dir+"02_adding_orders/pruning/families/",
                                                order = esse_clades[i],
                                                script_dir= script_dir,
                                                done_dir= done_dir,
                                                done= "Tip_removal"+esse_clades[i]
                                                ))
                    
                    gwf.target_from_template(name = esse_clades[i]+"_Biome_sampling_fraction.",
                                                template=sampling_frequency_per_biome(
                                                input_file_tree= "pruned_tree_family_"+esse_clades[i]+"_GBMB.tre",
                                                path_in =  workflow_dir+"02_adding_orders/pruning/families/",
                                                path_out = workflow_dir+"03_distribution_data/",
                                                output_file = esse_clades[i]+"_biome_sampling_fraction.txt",
                                                wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                                order = esse_clades[i],
                                                script_dir= script_dir,
                                                apg = script_dir+"apgweb_parsed.csv",
                                                done_dir= done_dir,
                                                done= esse_clades[i]+"_biome_sampling_fraction",
                                                renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                                koppen_biome = script_dir+"koppen_geiger_0p01.tif",
                                                percentage= percentages[j]
                                                ))

                    gwf.target_from_template(name = esse_clades[i]+"_Esse",
                                                template = Esse(
                                                tree_file = esse_clades[i]+"_Esse_tree.tre", # Input tree
                                                tip_states_file = workflow_dir+"03_distribution_data/"+esse_clades[i]+"_states_"+percentages[j]+".txt", 
                                                paleo_clim_file = data_dir+"paleoclim_area.txt", # File with paleoclimatic variables
                                                done = esse_clades[i]+"_Esse",
                                                path_in = workflow_dir+"02_adding_orders/pruning/orders/",
                                                save_file = "Esse_output_"+esse_clades[i]+"_"+percentages[j]+".jld2",
                                                script_dir=script_dir,
                                                done_dir = done_dir,
                                                out_states_file = "Esse_states_"+esse_clades[i]+"_"+percentages[j], 
                                                hidden_states = 0,
                                                out_file = "Esse_output_"+esse_clades[i]+"_hidden_states_"+percentages[j],
                                                output_folder = workflow_dir+"04_results/Esse_output/",
                                                path_out = workflow_dir+"04_results/"
                                             ))


###############################################################################################################################################################
###################################################---  Running ESSE on the sub_family_clades ---##############################################################
###############################################################################################################################################################

# Families which have been divided further and therefore need to be removed from the Clads_clades list.
# Amaryllidaceae, Anacardiaceae_Burseraceae_Kirkiaceae, Apiaceae, Apocynaceae, Asparagaceae, Asteraceae,Brassicaceae, Euphorbiaceae, Fabaceae, Lamiaceae, Orchidaceae, Poaceae, Rubiaceae
# Rosaceae, Ranunculaceae, Plantaginaceae, Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae,Myrtaceae, Lauraceae, Iridaceae, Gesneriaceae_Calceolariaceae, Cyperaceae,
# Caryophyllaceae_Achatocarpaceae_Amaranthaceae, ,Bromeliaceae, Acanthaceae_Martyniaceae_Pedaliaceae, Melastomataceae,

sub_family_clades = ["Acanthaceae_Martyniaceae_Pedaliaceae_1",
                     "Acanthaceae_Martyniaceae_Pedaliaceae_2",
                     "Acanthaceae_Martyniaceae_Pedaliaceae_3",
                     "Amaryllidaceae_1",
                     "Amaryllidaceae_2",
                     "Anacardiaceae_Burseraceae_Kirkiaceae_1",
                     "Anacardiaceae_Burseraceae_Kirkiaceae_2",
                     "Apiaceae_1","Apiaceae_2",
                     "Apiaceae_3",
                     "Apiaceae_4",
                     "Apiaceae_5",
                     "Apiaceae_6",
                     "Apocynaceae_1",
                     "Apocynaceae_2",
                     "Apocynaceae_3",
                     "Apocynaceae_4",
                     "Asparagaceae_1",
                     "Asparagaceae_2",
                     "Asparagaceae_3",
                     "Asteraceae_10",
                     "Asteraceae_11",
                     "Asteraceae_12",
                     "Asteraceae_13",
                     "Asteraceae_1",
                     "Asteraceae_2",
                     "Asteraceae_3",
                     "Asteraceae_4",
                     "Asteraceae_5",
                     "Asteraceae_6",
                     "Asteraceae_7",
                     "Asteraceae_8",
                     "Asteraceae_9",
                     "Brassicaceae_1",
                     "Brassicaceae_2",
                     "Brassicaceae_3",
                     "Brassicaceae_4",
                     "Bromeliaceae_1",
                     "Bromeliaceae_2",
                     "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_1",
                     "Caryophyllaceae_Achatocarpaceae_Amaranthaceae_2",
                     "Cyperaceae_1",
                     "Cyperaceae_2",
                     "Cyperaceae_3",
                     "Euphorbiaceae_1",
                     "Euphorbiaceae_2",
                     "Euphorbiaceae_3",
                     "Fabaceae_10",
                     "Fabaceae_11",
                     "Fabaceae_12",
                     "Fabaceae_13",
                     "Fabaceae_1",
                     "Fabaceae_2",
                     "Fabaceae_3",
                     "Fabaceae_4",
                     "Fabaceae_5",
                     "Fabaceae_6",
                     "Fabaceae_7",
                     "Fabaceae_8",
                     "Fabaceae_9",
                     "Gesneriaceae_Calceolariaceae_1",
                     "Gesneriaceae_Calceolariaceae_2",
                     "Iridaceae_1",
                     "Iridaceae_2",
                     "Lamiaceae_1",
                     "Lamiaceae_2",
                     "Lamiaceae_3",
                     "Lamiaceae_4",
                     "Lamiaceae_5",
                     "Lamiaceae_6",
                     "Lamiaceae_7",
                     "Lauraceae_1",
                     "Lauraceae_2",
                     "Melastomataceae_1",
                     "Melastomataceae_2",
                     "Myrtaceae_1",
                     "Myrtaceae_2",
                     "Myrtaceae_3",
                     "Orchidaceae_10",
                     "Orchidaceae_11",
                     "Orchidaceae_12",
                     "Orchidaceae_13",
                     "Orchidaceae_14",
                     "Orchidaceae_15",
                     "Orchidaceae_16",
                     "Orchidaceae_1",
                     "Orchidaceae_2",
                     "Orchidaceae_3",
                     "Orchidaceae_4",
                     "Orchidaceae_5",
                     "Orchidaceae_6",
                     "Orchidaceae_7",
                     "Orchidaceae_8",
                     "Orchidaceae_9",
                     "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_1",
                     "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_2",
                     "Phyllanthaceae_Picodendraceae_Linaceae_Ixonanthaceae_3",
                     "Plantaginaceae_1",
                     "Plantaginaceae_2",
                     "Plantaginaceae_3",
                     "Poaceae_1",
                     "Poaceae_2",
                     "Poaceae_3",
                     "Poaceae_4",
                     "Poaceae_5",
                     "Poaceae_6",
                     "Poaceae_7",
                     "Ranunculaceae_1",
                     "Ranunculaceae_2",
                     "Ranunculaceae_3",
                     "Rosaceae_1",
                     "Rosaceae_2",
                     "Rosaceae_3",
                     "Rosaceae_4",
                     "Rosaceae_5",
                     "Rubiaceae_1",
                     "Rubiaceae_2",
                     "Rubiaceae_3",
                     ]

for i in range(len(sub_family_clades)):
            # gwf.target_from_template(name = sub_family_clades[i]+"_distribution_data.",
            #                     template=Finding_areas_in_wcvp(
            #                     input_file_tree= "sub_phylo_"+sub_family_clades[i]+".tre",
            #                     path_in =  workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
            #                     path_out = workflow_dir+"03_distribution_data/",
            #                     output_file = sub_family_clades[i]+"_distribution_data.txt",
            #                     wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
            #                     order = sub_family_clades[i],
            #                     script_dir= script_dir,
            #                     apg = script_dir+"apgweb_parsed.csv",
            #                     done_dir= done_dir,
            #                     done= sub_family_clades[i]+"_distribution_data",
            #                     renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
            #                     koppen_biome = script_dir+"koppen_geiger_0p01.tif"
            #                     ))
                        
            # for j in range(len(percentages)):
            #         gwf.target_from_template(name = sub_family_clades[i]+"_states_converter_"+percentages[j],
            #                                 template=states_converter(
            #                                 path_in= workflow_dir+"03_distribution_data/",
            #                                 tip_states_file= workflow_dir+"03_distribution_data/"+sub_family_clades[i]+"_distribution_data.txt",
            #                                 out_states_file= workflow_dir+"03_distribution_data/"+sub_family_clades[i]+"_states_"+percentages[j]+".txt",
            #                                 script_dir= script_dir,
            #                                 done_dir= done_dir,
            #                                 done= "States_converter_"+sub_family_clades[i]+"_"+percentages[j]+"",
            #                                 percentage_for_present= percentages[j]
            #                                 ))
                    
            #         gwf.target_from_template(name = sub_family_clades[i]+"_Sampling_fraction",
            #                              template = sampling_frequency(
            #                                     input_file_tree= "sub_phylo_"+sub_family_clades[i]+".tre",
            #                                     path_in =  workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
            #                                     path_out = workflow_dir+"03_distribution_data/",
            #                                     output_file = sub_family_clades[i]+"_sampling_fraction.txt",
            #                                     wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
            #                                     order = sub_family_clades[i],
            #                                     script_dir= script_dir,
            #                                     apg = script_dir+"apgweb_parsed.csv",
            #                                     done_dir= done_dir,
            #                                     done= sub_family_clades[i]+"_Sampling_fraction"
            #                              ))
                    
            
                    gwf.target_from_template(name = sub_family_clades[i]+"_Tip_removal",
                                                template = rem_tips(
                                                input_file_tree = "sub_phylo_"+sub_family_clades[i]+".tre",
                                                distribution_file= workflow_dir+"03_distribution_data/"+sub_family_clades[i]+"_states_"+percentages[j]+".txt",
                                                output_file = sub_family_clades[i]+"_Esse_tree.tre",
                                                path_in = workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                                order = sub_family_clades[i],
                                                script_dir= script_dir,
                                                done_dir= done_dir,
                                                done= "Tip_removal"+sub_family_clades[i]
                                                ))
                    
                    gwf.target_from_template(name = sub_family_clades[i]+"_Biome_sampling_fraction.",
                                                template=sampling_frequency_per_biome(
                                                input_file_tree= "sub_phylo_"+sub_family_clades[i]+".tre",
                                                path_in =  workflow_dir+"02_adding_orders/pruning/subset_of_orders/",
                                                path_out = workflow_dir+"03_distribution_data/",
                                                output_file = sub_family_clades[i]+"_biome_sampling_fraction.txt",
                                                wcvp_file = workflow_dir+"02_adding_orders/wcvp_names_apg_aligned.rds",
                                                order = sub_family_clades[i],
                                                script_dir= script_dir,
                                                apg = script_dir+"apgweb_parsed.csv",
                                                done_dir= done_dir,
                                                done= sub_family_clades[i]+"_biome_sampling_fraction",
                                                renamed_occurrences = workflow_dir+"01_distribution_data/06_Renamed/gbif_renamed.rds", 
                                                koppen_biome = script_dir+"koppen_geiger_0p01.tif",
                                                percentage= percentages[j]
                                                ))
                        
                    gwf.target_from_template(name = sub_family_clades[i]+"_Esse",
                                                template = Esse(
                                                tree_file = sub_family_clades[i]+"_Esse_tree.tre", # Input tree
                                                tip_states_file = workflow_dir+"03_distribution_data/"+sub_family_clades[i]+"_states_"+percentages[j]+".txt", 
                                                paleo_clim_file = data_dir+"paleoclim_area.txt", # File with paleoclimatic variables
                                                done = sub_family_clades[i]+"_Esse",
                                                path_in = workflow_dir+"02_adding_orders/pruning/orders/",
                                                save_file = "Esse_output_"+sub_family_clades[i]+"_"+percentages[j]+".jld2",
                                                script_dir=script_dir,
                                                done_dir = done_dir,
                                                out_states_file = "Esse_states_"+sub_family_clades[i]+"_"+percentages[j], 
                                                hidden_states = 0,
                                                out_file = "Esse_output_"+sub_family_clades[i]+"_hidden_states_"+percentages[j],
                                                output_folder = workflow_dir+"04_results/Esse_output/",
                                                path_out = workflow_dir+"04_results/"
                                             ))
            
            
#####################################################################################################################################################################
############################################ Random old parts of the code that I am saving for one reason or another. ################################################
#####################################################################################################################################################################


# If I need to find the list of trees again, or somehow change the scripts so they print all the good trees into a single folder where I can just run the script on all of them.
# order_trees=["pruned_tree_order_Alismatales_GBMB.tre", "pruned_tree_order_Crossosomatales_GBMB.tre", "pruned_tree_order_Gunnerales_GBMB.tre", "pruned_tree_order_Poales_GBMB.tre",
# "pruned_tree_order_Amborellales_GBMB.tre", "pruned_tree_order_Cucurbitales_GBMB.tre", "pruned_tree_order_Huerteales_GBMB.tre",
# "pruned_tree_order_Aquifoliales_GBMB.tre", "pruned_tree_order_Cupressales_GBMB.tre", "pruned_tree_order_Magnoliales_GBMB.tre", "pruned_tree_order_Ranunculales_GBMB.tre",
# "pruned_tree_order_Arecales_GBMB.tre", "pruned_tree_order_Cycadales_GBMB.tre", "pruned_tree_order_Malpighiales_GBMB.tre", "pruned_tree_order_Rosales_GBMB.tre",
# "pruned_tree_order_Austrobaileyales_GBMB.tre", "pruned_tree_order_Dilleniales_GBMB.tre", "pruned_tree_order_Malvales_GBMB.tre", "pruned_tree_order_Santalales_GBMB.tre",
# "pruned_tree_order_Berberidopsidales_GBMB.tre", "pruned_tree_order_Dioscoreales_GBMB.tre", "pruned_tree_order_Myrtales_GBMB.tre", "pruned_tree_order_Sapindales_GBMB.tre",
# "pruned_tree_order_Bruniales_GBMB.tre", "pruned_tree_order_Dipsacales_GBMB.tre", "pruned_tree_order_Nymphaeales_GBMB.tre", "pruned_tree_order_Solanales_GBMB.tre",
# "pruned_tree_order_Buxales_GBMB.tre", "pruned_tree_order_Ericales_GBMB.tre", "pruned_tree_order_Pandanales_GBMB.tre", "pruned_tree_order_Trochodendrales_GBMB.tre",
# "pruned_tree_order_Canellales_GBMB.tre", "pruned_tree_order_Escalloniales_GBMB.tre", "pruned_tree_order_Paracryphiales_GBMB.tre", "pruned_tree_order_Vahliales_GBMB.tre",
# "pruned_tree_order_Celastrales_GBMB.tre", "pruned_tree_order_Fabales_GBMB.tre", "pruned_tree_order_Petrosaviales_GBMB.tre", "pruned_tree_order_Vitales_GBMB.tre",
# "pruned_tree_order_Ceratophyllales_GBMB.tre", "pruned_tree_order_Garryales_GBMB.tre", "pruned_tree_order_Picramniales_GBMB.tre", "pruned_tree_order_Zingiberales_GBMB.tre",
# "pruned_tree_order_Chloranthales_GBMB.tre", "pruned_tree_order_Ginkgoales_GBMB.tre", "pruned_tree_order_Pinales_GBMB.tre", "pruned_tree_order_Zygophyllales_GBMB.tre",
# "pruned_tree_order_Commelinales_GBMB.tre", "pruned_tree_order_Gnetales_GBMB.tre", "pruned_tree_order_Piperales_GBMB.tre"
        


