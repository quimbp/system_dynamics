# exp_4dvar
# Run the lorenz model 4DVAR

FDVAR=../bin/fdvar

rm -f fdvar.log

#${FDVAR} ff_00_10.namelist >> fdvar.log
${FDVAR} ff_10_20.namelist >> fdvar.log
#how lorenz_simulation.dat
#how ff_00_10.dat
#how ff_10_20.dat

