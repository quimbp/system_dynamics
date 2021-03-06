# SYSTEM DYNAMICS PROJECT
# Quim Ballabrera, November 2017
#
# Compiler options and libraries
# The code has been created and tested using gfortran.
# The NF90_LIB and NF90_INC variables had been fed from the
# bash environment.

# Modify according operative system and library paths:
#
# SYSDYN = the project root folder
# FC = Fortran95 compiler
# FFLAGS = Compiler options
# INSTALL = Program used to copy files setting theis mode bits, etc.
# AR = Library archive manager 
# NF90_INC = Path to the  NetCDF Fortran 90 modules
# NF90_LIB = Path to the  NetCDF Fortran 90 libraries
# 

SYSDYN = $(HOME)/system_dynamics

FC = gfortran
FFLAGS = -O4 -fno-range-check -finit-real=nan -finit-integer=0 \
         -finit-character=0 -finit-logical=false 
INSTALL = install
AR = ar cq

#NF90_INC = /usr/lib64/gfortran/modules/
#NF90_LIB = /usr/lib64/

# Folders where files are going to be installed
#
BINDIR = $(SYSDYN)/bin
LIBDIR = $(SYSDYN)/lib
INCDIR = $(SYSDYN)/include

#CDFLIB = -L$(NF90_LIB) -lnetcdff
#CDFINC = -I$(NF90_INC)
CDFLIB = 
CDFINC = 
