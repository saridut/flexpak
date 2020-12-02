#Author: Sarit Dutta
#Date: 07/2017

#See below for header file dependencies
#http://scottmcpeak.com/autodepend/autodepend.html
#http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/

# This Makefile resides in 'build' directory. The source files reside in
# 'src', object files in 'build', and executables in 'exec'.

SHELL = /bin/bash
SRC_DIR = . # ../src
BUILD_DIR = .
EXEC_DIR = . # ../exec

vpath %.c   $(SRC_DIR)
vpath %.h   $(SRC_DIR)
vpath $(EXECUTABLE) $(EXEC_DIR)

#C compiler
CC = $(mpgcc)

CFLAGS = -c -O0 -Wall -std=c11 -m64\
		 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans

#Define name of the executable
EXECUTABLE = els

HDF5_ROOT = /Users/saritdutta/soft/hdf5/hdf5-1.10.0
HDF5_INCLUDE = $(HDF5_ROOT)/include
HDF5_LIB = $(HDF5_ROOT)/lib
H5CC = $(HDF5_ROOT)/bin/h5cc

GSL_ROOT = /Users/saritdutta/soft/gsl/gsl-2.4
GSL_INCLUDE = $(GSL_ROOT)/include
GSL_LIB = $(GSL_ROOT)/lib

FFTW_ROOT = /Users/saritdutta/soft/fftw/fftw3-3.3.6-pl1/fftw
FFTW_INCLUDE = $(FFTW_ROOT)/include
FFTW_LIB = $(FFTW_ROOT)/lib

INCS = \
	-I$(HDF5_INCLUDE) 	\
	-I$(GSL_INCLUDE) 	\
	-I$(FFTW_INCLUDE)	\
	-I$(MKLROOT)/include

LIBS = \
	-L$(HDF5_LIB) $(HDF5_LIB)/libhdf5_hl.a $(HDF5_LIB)/libhdf5.a -lz -ldl -lm

#		-L$(GSL_LIB) -lgsl \
#		-L$(FFTW_LIB) -lfftw3 \
#		-L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib \
#		-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm # -ldl

OBJECTS = main.o \
		  mesh.o \
 		  io.o
#
#move : $(EXECUTABLE)
#	if [ -e "$(EXECUTABLE)" ]; then mv $(EXECUTABLE) $(EXECDIR);fi 

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(INCS) $(LIBS) $(OBJECTS) -o $(EXECUTABLE) 

main.o: main.c
	$(CC) $(CFLAGS) $(INCS) $<

mesh.o: mesh.c 
	$(CC) $(CFLAGS) $(INCS) $<

io.o: io.c 
	$(CC) $(CFLAGS) $(INCS) $<


.PHONY : clean clobber
clean:
	-rm -r *.o
clobber:
	-rm -r $(EXECDIR)/$(EXECUTABLE) *.o



# DO NOT DELETE
