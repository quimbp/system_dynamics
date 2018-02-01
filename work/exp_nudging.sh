# exp_nudging
# Run the lorenz model
# Perform Forward nudging

rm -f nudging.log

NUDGING=../bin/nudging

${NUDGING} nudging_02_010_100.namelist > nudging.log
${NUDGING} nudging_02_025_100.namelist >> nudging.log
${NUDGING} nudging_02_050_100.namelist >> nudging.log
${NUDGING} nudging_02_075_100.namelist >> nudging.log
display_results.py lorenz_simulation.dat noise_02_sampl_010

