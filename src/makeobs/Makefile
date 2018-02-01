include ../../make.inc

TARGET = makeobs
OBJ = makeobs.o

%.o: %.f90 $(LIBDIR)/libsysdyn.a
	$(FC) $(FFLAGS) -c $< -I$(INCDIR)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ) -L$(LIBDIR) -lsysdyn
	$(INSTALL) $(TARGET) $(BINDIR)

clean:
	rm -f *.o *.exe *.mod $(TARGET) $(BINDIR)/$(TARGET)
