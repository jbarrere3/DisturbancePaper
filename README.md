# DisturbancePaper

R script to reproduce the analyses of the paper *Functional traits and climate drive interspecific differences in disturbance-induced tree mortality* by Julien Barrere, Björn Reineking, Thomas Cordonnier, Niko Kulha, Juha Honkaniemi, Kari T. Korhonen, Mikko Peltoniemi, Paloma Ruiz-Benito, Miguel A. Zavala and Georges Kunstler

The bayesian model presented in the paper can not be reproduced here as the National Forest Inventory data from Finland are not publicly accessible. This script reproduces the analyses testing the effects of functional traits and of climate on sensitivity. 

Before running the script, some data are required. The list of files needed and the structure that the "data" folder must have is visible in lines 42 to 49 of the ```_targets.R``` script. 
- Most files (data/sensitivity, data/climate, traits calculated from NFI data and shade tolerance) are stored in [Zenodo](https://doi.org/10.5281/zenodo.7603488)
- To get the file containing the traits extracted from TRY database, a requested must be sent to TRY [here](https://www.try-db.org/TryWeb/Prop0.php). The list of the traits and species to request is stored in the zenodo repository. 
- Data on wood density, P50, and root traits (to include in the data/traits directory) can respectively be downloaded in the following public databases: [Wood density database](https://datadryad.org/stash/dataset/doi:10.5061/dryad.234), [Sureau database](https://zenodo.org/record/854700) and [GROOT database](https://groot-database.github.io/GRooT/)

Package ```targets``` is needed to run the script. Once the package is installed, just run ```targets::tar_make()``` from R and the script should install automatically the packages missing, and run the analyses. The figures will be stored in the directory "output". 

For any question about the script or analyses, contact Julien BARRERE (julien.barrere@inrae.fr). 
