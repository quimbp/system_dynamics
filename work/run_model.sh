# run_model.sh
# Run the lorenz model
# Gets four sets of observations

LORENZ_MODEL=../bin/lorenz
MAKEOBS=../bin/makeobs

# Erase log file
rm -f run_model.log 

${LORENZ_MODEL} -trajectory lorenz_simulation.dat \
                -end lorenz_end.dat   > run_model.log

python3 ../src/lorenz/lorenz.py lorenz_simulation.dat
mv lorenz.pdf plots/

${MAKEOBS} det_02_010.namelist   >> run_model.log
${MAKEOBS} det_02_025.namelist   >> run_model.log
${MAKEOBS} det_02_050.namelist   >> run_model.log
${MAKEOBS} det_02_075.namelist   >> run_model.log

python3 ../src/makeobs/display_obs.py ./lorenz_simulation.dat

echo 'run_model done'
