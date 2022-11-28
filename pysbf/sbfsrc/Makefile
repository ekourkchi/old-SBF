# for Mac OSX Snow Leopard and Lion (25 Mar 2013)
FC= gfortran
CC= gcc
CFLAGS= -O  -fPIC -Wno-implicit-function-declaration  -c 
FFLAGS= -O -fno-automatic  -fPIC -c
FLIBS=
#XLIBS= -L /usr/X11R6/lib -lX11
ULOC= ../util


all: likenew6.so

likenew6.so: likenew6_ $(ULOC)/mini.o $(ULOC)/invert.o \
	$(ULOC)/rwfits.o $(ULOC)/lnb.o
	$(FC) -shared likenew6_.o $(ULOC)/mini.o $(ULOC)/invert.o \
	$(ULOC)/rwfits.o $(ULOC)/lnb.o -o likenew6.so

likenew6_: likenew6_.o $(ULOC)/mini.o $(ULOC)/invert.o \
	$(ULOC)/rwfits.o $(ULOC)/lnb.o 
	$(FC) $(FFLAGS) likenew6_.o $(ULOC)/mini.o $(ULOC)/invert.o \
	$(ULOC)/rwfits.o $(ULOC)/lnb.o $(FLIBS) -o likenew6_

clean:
	rm -rf *o likenew60
	rm -rf $(ULOC)/*o 
