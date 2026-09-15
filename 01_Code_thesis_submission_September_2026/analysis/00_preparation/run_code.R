# Directory paths and loading of functions
# 
# In this file the project paths are defined and the scripts containing the main functions are seourced.
#
# no AI used

rm (list = ls())
path_wd = getwd() # or your path
path_preparation = "/04_code/00_preparation"
path_functions = "/04_code/functions"
path_stochastic_functions = "/04_code/functions/Simulations_stochastic"
path_deterministic_functions = "/04_code/functions/Simulations_deterministic"
path_data = "/04_code/data_preperation"

path_experiments = "/04_code/experiments"
path_start_intervention = "/Experiments_intervention_start_time"
path_influx = "/Experiments_influx"
path_age = "/Experiment_age"
path_results_plot = "/03_results/Experiment_plots"
path_results_data = "/03_results/Experiments_data"
path_finite_size = "/Experiment_finite_size"
path_all_exp = "/All_Experiments_ended_up_in_thesis"
# run basics
source(paste0(path_wd,path_preparation,"/read_packages.R"))

# run extra functions
source(paste0(path_wd, path_functions,"/compute_burden.R"))

# run deterministic functions
#source(paste0(path_wd,path_stochastic_functions,"/determine parameters.R"))
source(paste0(path_wd,path_deterministic_functions,"/deterministic_SIR.R"))
source(paste0(path_wd,path_deterministic_functions,"/det_SIR_with_influx.R"))
source(paste0(path_wd,path_deterministic_functions,"/deterministic_SIR_with_age.R"))

# run stochastic functions
source(paste0(path_wd,path_stochastic_functions,"/state_transition.R"))
source(paste0(path_wd,path_stochastic_functions,"/initiaization_initial_state.R"))
source(paste0(path_wd,path_stochastic_functions,"/intervention_functions.R"))
source(paste0(path_wd,path_stochastic_functions,"/rate_function_age.R"))
source(paste0(path_wd,path_stochastic_functions,"/Simulation_wrapper.R"))
source(paste0(path_wd,path_stochastic_functions,"/Multi_replicate_simulation.R"))

# Data 
source(paste0(path_wd,path_data,"/contact_matrix.R")) # das noch als funktion schreiben
#source(paste0(path_wd,path_experiments,"/setup_parameters_template.R")) # das noch als funktion schreiben

#Experiments
#source(paste0(path_wd,path_experiments,"/Experiments_get_to_know_the_Setting/Experiment6_no_intervention.R"))
#source(paste0(path_wd,path_experiments,"/Experiments_get_to_know_the_Setting/Experiment6_no_intervention_N100k.R"))
#source(paste0(path_wd,path_experiments,path_influx,"/Experiment4_influx_N100k.R"))


#to do 
