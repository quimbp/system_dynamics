# run_eof.sh
# Run the lorenz model
# Gets the covariance matrix

EOF=../bin/eof

# Erase log file
rm -f run_eof.log 

${EOF} 

echo 'run_eof done'
