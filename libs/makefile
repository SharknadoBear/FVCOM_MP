# $Id$
# $Name$
# $Revision$
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh

MYINSTALLDIR = $(INSTALLDIR)

#PACKAGES =	proj	fproj	julian	metis	netcdf	
PACKAGES =	fproj	julian	metis


all:
	for item in $(PACKAGES); do (./untar.sh $$item ) || exit 1; done
#	cd proj && ./configure CC=$(CC) CFLAGS=-O3 CXX=$(CC) CXXFLAGS=-O3 F77=$(FC) FFLAGS=-O3 --prefix=$(MYINSTALLDIR)
#	cd proj && make install
	cd fproj && ./configure CPPFLAGS='$(COMPILER)' CC=$(CC) CFLAGS=-O3 CXX=$(CXX) CXXFLAGS="-O3 -xhost" FC=$(FC) FFLAGS="-O3 -xhost" --with-proj4=/share/apps/proj4/4.9.3/ --prefix=$(MYINSTALLDIR)
	cd fproj && make install
#	cd netcdf && ./configure CC=$(CC) CFLAGS=-O3 CXX=$(CC) CXXFLAGS=-O3 F77=$(FC) F90=$(FC) FFLAGS=-O3 --prefix=$(MYINSTALLDIR) --build=$(MACHTYPE)
#	cd netcdf && make install
	cd metis && make install
	cd julian && make install



clean:
	for item in $(PACKAGES); do (cd $$item && make clean) || exit 1; done

clobber:
	for item in $(PACKAGES); do (rm -rf $$item) || exit 1; done
	rm -rf $(MYINSTALLDIR)/*

includes::
include ../make.inc
