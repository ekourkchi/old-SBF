FC= gfortran
CC= gcc
CFLAGS= -O  
# for ultrix
# FFLAGS= -O -static
# For SUNOS
# FFLAGS= -O -Bstatic
# For Solaris
FFLAGS= -O -fno-automatic -Wno-implicit-function-declaration
PLIBS= -lmongo
XLIBS= -lX11

match: match.o
	$(FC) $(FFLAGS) match.o sort.o lnb.o \
	$(FLIBS) -o match

linear: linear.o invert.o
	$(FC) $(FFLAGS) linear.o invert.o $(FLIBS) -o linear
