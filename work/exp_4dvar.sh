# exp_4dvar
# Run the lorenz model 4DVAR

FDVAR=../bin/fdvar

${FDVAR} fdvar_02_010.namelist 
display_results.py lorenz_simulation.dat noise_02_sampl_010

