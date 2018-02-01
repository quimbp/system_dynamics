include ../../make.inc

TARGET = enkf ensemble
OBJ = etools.o etkf.o

%.o: %.f90 $(LIBDIR)/libsysdyn.a
	$(FC) $(FFLAGS) -c $< -I$(INCDIR)

all: $(TARGET) 

enkf: enkf.f90 $(OBJ)
	$(FC) $(FFLAGS) -o $@ $@.f90 $(OBJ) -I$(INCDIR) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $@ $(BINDIR)

ensemble: ensemble.f90 etools.o
	$(FC) $(FFLAGS) -o $@ $@.f90 etools.o -I$(INCDIR) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $@ $(BINDIR)

clean:
	rm -f $(OBJ) *.exe *.mod $(TARGET)
	cd $(BINDIR); rm -f $(TARGET)
