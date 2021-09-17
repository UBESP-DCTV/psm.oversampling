# Propensity score matching and Oversampling
This repository contains the R scripts to reproduce the results of the
Monte Carlo simulations of the study "__Oversampling and Replacement 
Strategies in Propensity Score Matching: A Critical Review Focused on 
Small Sample Size in Clinical Settings__".

Before running the scripts, please install the following R packages:

``` r
install.packages(
  c("tidyverse", "Matching", "cobalt", "glue", "furrr",
    "tictoc", "MASS", "doParallel", "see", "ggpubr", "knitr", "here",
    "assertive", "usethis"), 
    dependencies = TRUE
)
```

The project contains the following files:

- The file _psm-oversampling.Rproj_ is the file of the Rstudio project.

- The file _sim_functions.R_, that contains the functions used for the
  generation of the simulated datasets.
      
- The file _sim_setup.R_, that contains the code used to simulate
  the scenarios for the main Monte Carlo simulations.

- The file _analysis_functions.R_, that contains the functions 
  implemented for the analysis on the simulated datasets in main simulations.
      
- The file _analysis_mc_sim.R_, that contains the code used to run
  the analysis on the simulated datasets in main simulations.

- The file _sim_results.Rmd_, an Rmarkdown file that can be used to
  display the results of the main simulations.
  
- The file _sim_functions_linear.R_, that contains the functions used for the
  generation of the simulated datasets in secondary simulations.
      
- The file _sim_setup_linear.R_, that contains the code used to simulate
  the scenarios for the secondary Monte Carlo simulations.

- The file _analysis_functions_linear.R_, that contains the functions 
  implemented for the analysis on the simulated datasets in secondary simulations.
      
- The file _analysis_mc_sim_linear.R_, that contains the code used to run
  the analysis on the simulated datasets in secondary simulations.

- The file _sim_results_linear.Rmd_, an Rmarkdown file that can be used to
  display the results of the secondary simulations.


The results of the main simulations can be replicated as follows:

1. Run the script _sim_setup.R_, which saves the simulated scenarios
   in a file named __simulation_data_pop_10000000_n_1000.rda__.
   
2. Run the script _analysis_mc_sim.R_, which stores the results of
   the simulations in a file named __simulation_results_pop_1000000_n_1000.rda__.

3. Knitr the _sim_results.Rmd_ file to visualize the results of the
   simulations.
   
The results of the secondary simulations can be replicated as follows:

1. Run the script _sim_setup_linear.R_, which saves the simulated scenarios
   in a file named __simulation_data_pop_10000000_n_1000_linear.rda__.
   
2. Run the script _analysis_mc_sim_linear.R_, which stores the results of
   the simulations in a file named __simulation_results_pop_1000000_n_1000_linear.rda__.

3. Knitr the _sim_results_linear.Rmd_ file to visualize the results of the
   simulations.

