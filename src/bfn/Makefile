include ../../make.inc

TARGET = bfn
OBJ = bfn.o

FFLAGS=-fbounds-check

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -I../../include

all: $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ) -L../../lib -lsysdyn
	$(INSTALL) $(TARGET) ../../bin


clean:
	rm -f *.o *.exe *.mod $(TARGET) ../../bin/$(TARGET)

