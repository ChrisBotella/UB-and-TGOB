# UB-and-TGOB
Script for the simulation experiment of Botella et al. (2019) Bias in presence-only niche models related to sampling effort and species niches: lessons for background point selection.

For reproducing the experiment
1) Open **make_species_table_ready.R** and change the variable **dir** to the desired local repository location.
2) Run **make_species_table_ready.R** with R.
3) Open **simu_fit_ready.R** and change the variable **dir** to the same location as in the previous script. 
4) Run **simu_fit_ready.R** with R. It will generate and save each fitted model in **dir** and finally saves PNG images files containing the results multiplots of all scenarios.  
5) For producing the three scenarios of "TG species density", set the variable **a** to successively to "flat", "thick" and "thin" and run the script.
