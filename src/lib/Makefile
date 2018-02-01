include ../../make.inc

TARGET = $(LIBDIR)/libsysdyn.a
INC = $(CDFINC) -I$(INCDIR)

all: sysdyn
OBJ = types.o constants.o utils.o lineargs.o help.o dates.o \
      runge_kutta.o math.o lorenz.o obs.o minimization.o    \
      sysdyn.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(INC)

clean:
	rm -f *.mod *.o $(TARGET)
	rm -f $(INCDIR)/*
	rm -f $(LIBDIR)/*

sysdyn: $(OBJ)
	/bin/rm -f $(TARGET)
	$(FC) $(FFLAGS) -c sysdyn.f90 $(INC)
	$(INSTALL) *.mod $(INCDIR)
	$(AR) $(TARGET) *.o

