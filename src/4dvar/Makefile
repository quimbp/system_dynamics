include ../../make.inc

TARGET = test_ltm test_adj test_cost fdvar
OBJ = adjoint.o cost.o

#FFLAGS = -fbounds-check -Wall

%.o: %.f90 $(LIBDIR)/libsysdyn.a
	$(FC) $(FFLAGS) -c $< -I$(INCDIR)

all: $(TARGET) 

fdvar: fdvar.f90 $(OBJ)
	$(FC) $(FFLAGS) -o $@ $@.f90 $(OBJ) -I$(INCDIR) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $@ $(BINDIR)

test_ltm: test_ltm.f90 
	$(FC) $(FFLAGS) -o $@ $@.f90 -I$(INCDIR) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $@ $(BINDIR)

test_adj: test_adj.f90 
	$(FC) $(FFLAGS) -o $@ $@.f90 -I$(INCDIR) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $@ $(BINDIR)

test_cost: test_cost.f90 $(OBJ)
	$(FC) $(FFLAGS) -o $@ $@.f90 $(OBJ) -I$(INCDIR) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $@ $(BINDIR)

clean:
	rm -f $(OBJ) *.exe *.mod $(TARGET)
	cd $(BINDIR); rm -f $(TARGET)
