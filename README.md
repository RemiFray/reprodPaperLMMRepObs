# Estimating population size simultaneous non-invasive sampling of the same individual liable to misidentification

This repository contains all the material needed to reproduce the results from the paper.
The instructions are given below.

## Opening the project

Once you have cloned or downloaded and extracted the folder, **open the project through the .Rproj file**. The working directory, will be automatically defined as the folder containing it. That way, code lines using file paths should work without modifying any code.


## Data

First generate the simulation data. In the folder "data" are 3 scripts to generate data according to the $M_{\lambda, \alpha}$ model for the different numbers of occasion.
Simply run the scripts in any order.

The otter data from the study *Non-Invasive Genetic Mark-Recapture as a Means to Study Population Sizes and Marking Behaviour of the Elusive Eurasian Otter (Lutra lutra)* ([paper DOI](https://doi.org/10.1371/journal.pone.0125684)), available in a pdf ([data DOI](https://doi.org/10.1371/journal.pone.0125684.s002)) and are also in the otter.csv in the data folder.


## Run the models

The scripts used to run the models are in the folder "**mainScripts**". There are two sub-folders for:

* the simulation study with simultaneous observations (simulations by the model $M_{\lambda, \alpha}$), 
  * the Mla_pop30 is the script for the study in the supplementary material 
  * the others are for the main study
* the otter study.

Most of the time, there is one script for one model and one number of occasion. For efficiency reasons, we ran each script divided in three (cheap parallelisation) but we regrouped them in here for clarity. The ones that were not regrouped are because the MCMC was not ran on the same number of iteration.  

The scripts will save the markov chains in the result folder.


## Extract the results

There are scripts to extract summaries of the MCMC available in the folder "plotResults". There is one script per model. At the end of one of these scripts per study, there is also a few lines of code to merge all summary tables from the different models of the study.

There is also one script per simulation study to generate the plots in the paper.


## Functions

This folder contains all the functions, models, distributions and samplers needed to run the models. Nothing need to be touched in there.


## Results

The summary tables of results from the simulation study and the otter study are directly available in the result folder.

