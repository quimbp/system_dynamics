# exp_enkf
# Run the lorenz model
# Get three sets of observations

ENSEMBLE=../bin/ensemble
ENKF=../bin/enkf

${ENSEMBLE} ensemble.namelist

${ENKF} enkf_02_010.namelist 
display_results.py lorenz_simulation.dat noise_02_sampl_010

