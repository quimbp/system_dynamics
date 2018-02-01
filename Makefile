include make.inc

all: include lib bin work SYSDEM LORENZ EOF MAKEOBS INTERPOL NUDGING ENKF 4DVAR
	@echo "Done"

include:
	mkdir $@

lib:
	mkdir $@

bin:
	mkdir $@

work:
	mkdir $@
	cp python/*py work

SYSDEM:
	@echo 
	@echo "=============================================="
	@echo "Compiling the SYSDEM library"
	@echo "=============================================="
	@echo 
	(cd src/lib/; make all)

LORENZ:
	@echo 
	@echo "=============================================="
	@echo "Compiling LORENZ MODEL"
	@echo "=============================================="
	@echo 
	(cd src/lorenz/; make all)

EOF:
	@echo 
	@echo "=============================================="
	@echo "Compiling Empirical Orthogonal Functions"
	@echo "=============================================="
	@echo 
	(cd src/eof/; make all)

MAKEOBS:
	@echo 
	@echo "=============================================="
	@echo "Compiling MAKEOBS UTILITY"
	@echo "=============================================="
	@echo 
	(cd src/makeobs/; make all)

INTERPOL:
	@echo 
	@echo "=============================================="
	@echo "Compiling INTERPOLATION UTILITY"
	@echo "=============================================="
	@echo 
	(cd src/interpolation/; make all)

NUDGING:
	@echo 
	@echo "=============================================="
	@echo "Compiling NUDGING METHOD"
	@echo "=============================================="
	@echo 
	(cd src/nudging/; make all)

ENKF:
	@echo 
	@echo "=============================================="
	@echo "Compiling Ensemble KF"
	@echo "=============================================="
	@echo 
	(cd src/enkf/; make all)

4DVAR:
	@echo 
	@echo "=============================================="
	@echo "Compiling 4DVAR METHOD"
	@echo "=============================================="
	@echo 
	(cd src/4dvar/; make all)

clean:
	(cd src/lib; make clean)
	(cd src/lorenz; make clean)
	(cd src/eof; make clean)
	(cd src/makeobs; make clean)
	(cd src/interpolation; make clean)
	(cd src/nudging; make clean)
	(cd src/enkf; make clean)
	(cd src/4dvar; make clean)

