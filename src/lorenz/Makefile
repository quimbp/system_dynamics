include ../../make.inc

TARGET = lorenz
OBJ = main.o

%.o: %.f90 $(LIBDIR)/libsysdyn.a
	$(FC) $(FFLAGS) -c $< -I$(INCDIR) $(CDFINC)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ) -L$(LIBDIR) -lsysdyn $(CDFLIB)
	$(INSTALL) $(TARGET) $(BINDIR)

clean:
	rm -f $(OBJ) *.mod $(TARGET) $(BINDIR)/$(TARGET)
