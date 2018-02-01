# exp_interp
# From the observation sets, get the time interpolation

INTERPOLATION=../bin/interpolation

rm -f interpolation.log

${INTERPOLATION} ./obs_noise_02_sampl_010.dat \
                   interp_noise_02_sampl_010.dat >> interpolation.log

${INTERPOLATION} ./obs_noise_02_sampl_025.dat \
                   interp_noise_02_sampl_025.dat >> interpolation.log

${INTERPOLATION} ./obs_noise_02_sampl_050.dat \
                   interp_noise_02_sampl_050.dat >> interpolation.log

./display_results.py lorenz_simulation.dat noise_02_sampl_010
