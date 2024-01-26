# TroRaiMo ReadME

This page is currenrly under construction. Please check back soon for updates.


# TroRaiMo computational pipeline
This is a directory for the article "Speciation does not drive tropical rainforest biodiversity". 

The directory contains a pipeline which handles everything from the downloading of the raw data to the production of the figures in the article. 

This pipeline is created using the General workflow manager ([GWF](https://gwf.app/)) and [CONDA](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on a high performance computation cluster using [Slurm](https://slurm.schedmd.com/documentation.html). A basic understanding of these three software packages will help you understand how the pipeline works, but ideally it should not be necessary to run the pipeline. ***Be adviced, this workflow is meant for a high performance computational cluster and as such some of the jobs takes several days to run and requires several times the amount of computational power found on a standard computer. As such it is not feasible to run this workflow without the access to such a cluster***

The pipeline will on its own create all the necessary folders, download the data, clean the data, run the models and lastly produce all the figures.

## Set up
Your first step is to acess your high performance cluster and download the github repository.
```
git clone git@github.com:OscarWrisberg/TroRaiMo.git
```
Make sure conda is installed on your server.
This can by downloading it and running the installer as so:  
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh
chmod +x miniforge.sh
bash miniforge.sh -b
./miniforge3/bin/conda init bash
```

now navigate to the TroRaiMo folder and copy the conda environments saved in the repository, using the following commands:
```
conda env create -f gwf_env.yml
conda env create -f R_env.yml
conda env create -f julia_env.yml
```
You have now installed the required programming languages on the cluster for the pipeline to run. The pipeline will activate the required environments when it needs them. 

You can now activate the gwf_env using conda: 
 
```
conda activate gwf_env
```
now you need to set the backend for gwf which in this case is slurm:
```
gwf config set backend slurm
```

## Running the pipeline
you are now ready to run the pipeline. 
While in the gwf_env type in the terminal:
```
gwf run
```
All the jobs are now submitted to the queing system on the cluster. When all the dependencies of a job is met and slurm has the requested resources available your job will start running on the cluster.
if you write:
``` 
gwf status 
```
You will see all the status of all the jobs in the workflow.
in order to have a closer look at each individual job, you can write:
```
gwf logs "job_name"
```
Where you change "job_name" to the name of the job you want a closer look at.
If any of the jobs change status to "shouldrun" then something has gone wrong in the pipeline. You can check what has gone wrong with a job by looking at the error messages for a job:
```
gwf logs -e "job_name"
```
When a job reaches the status "completed" all the expected output files from the job have been produced, and the jobs which have these output files as their dependency will now start running automatically.

The pipeline might take several days to finish but you should keep track of it with ``` gwf status ``` at regular intervals to see if it progresses as expected. 

If you have any questions or problems regarding the pipeline dont hesitate to contact me at: oscar.wrisberg@bio.au.dk and I will do my best to help you. 


# Percularities
## Paleoenvironmental Data
Due too not being able to install libgdal-dev on the cluster used for this project, the paleoenvironmental data was downloaded and processed on a local machine. The processed data was then uploaded to the github repository. The pipeline will check if libgdal-dev is installed and if it is it will run the code producing the results. If libgdal-dev isent installed the code will fetch the results from the github repository. 

