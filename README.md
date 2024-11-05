
# DES-HE-Schizophrenia

Discrete event simulation model for the economic evaluation in
schizophrenia

This folder contains the R codes for the discrete event simulation (DES)
model for the economic evaluation in schizophrenia and the necessary
data for the illustrative example of applying the model for such
evaluation in the UK. To run this model, you need to download all the
files (r script and csv files) in a specific folder, open the r script
"run_model.r" and specify the path of the folder where you store these
files at the beginning of the script, then follow the instruction in the
r script to run the script.

The file "run_model.r" includes the r codes of the functions to run the
DES model, and the r codes for running the DES model with an
illustrative example.

The file "dat_pats.csv" includes the data of 100 patients of the
synthetic cohort.

The file "dat_model.csv" includes the base-case input for the parameters
used in the model.

The file "dat_lifetable.csv" includes the input for the UK life table by
age and sex in 2020.

An R Shiny interface webpage has been developed for this model.
<https://livedataoxford.shinyapps.io/shiny_des_schizophrenia/>
