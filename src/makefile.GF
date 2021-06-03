#-----------------------------------------------------------------
#
# AML
#
#-----------------------------------------------------------------
#
# $Header: /home/cvs/phsht/AML/src/makefile.IC,v 1.2 2007/09/25 10:30:33 phsht Exp $
#
# $Log: makefile.IC,v $
# Revision 1.2  2007/09/25 10:30:33  phsht
# rewrote to make "look" nicer and removed an ERROR for NEvals and Energy; now makefile has proper flags to check for these things automatically
#
# Revision 1.1  2007/09/20 16:53:39  phrfar
# previous files from a project of Rudo's. To be used as templates.
#
#-----------------------------------------------------------------#

F77 = gfortran -ffixed-form -std=legacy -fimplicit-none
F90 = gfortran
#F90 = gfortran

F77FLAGS = -fbounds-check -Warray-bounds  # For debugging
#F77FLAGS = -fno-range-check -Warray-bounds  # For optimizing
#F77FLAGS = # For profiling

F90FLAGS = -g -std=legacy -fbounds-check -fno-range-check -fimplicit-none -Warray-bounds #-Wall # For debugging
#F90FLAGS = -ffree-line-length-none -fbounds-check -mcmodel=large -fdefault-integer-8 -fno-range-check -fimplicit-none -fcray-pointer -fexternal-blas # For optimizing
#F90FLAGS = -fno-range-check # for profling

CC = gcc
CCFLAGS = -O2

LIBDIR = $(HOME)/f77/lib
LIBFLAGS  = -L../lib/INT64NGNU -ljadamilu -lblas -llapack #-lg2c
#LIBFLAGS  = -llapack -lblas 

PRINT1 = a2ps -g -o-
PRINT2 = > AMLdiag.ps # | fleet -odouble

CommonModules.o:      CommonModules.f90
	$(F90) $(F90FLAGS) -c CommonModules.f90

AMLModules.o:     AMLModules.f90
	$(F90) $(F90FLAGS) -c AMLModules.f90

random.o:     random.f90
	$(F90) $(F90FLAGS) -c random.f90

main.o:     main.f90
	$(F90) $(F90FLAGS) -c main.f90

util.o:     util.f90
	$(F90) $(F90FLAGS) -c util.f90

inout.o:     inout.f90
	$(F90) $(F90FLAGS) -c inout.f90

etime.o:     etime.f90
	$(F90) $(F90FLAGS) -c etime.f90

AMLdiag.GF:    CommonModules.o AMLModules.o random.o main.o util.o inout.o etime.o makefile.GF
	$(F90) $(F90FLAGS) -o AMLdiag.GF CommonModules.o AMLModules.o\
                util.o main.o inout.o random.o etime.o $(LIBFLAGS)

tar:	AMLdiag.tar	
AMLdiag.tar:   makefile.* main.f90 util.f90 inout.f90 \
                random.f90 AMLdiag.inp
	tar -cvf AMLdiag.tar makefile.* main.f90 util.f90 \
                inout.f90 random.f90 AMLdiag.inp

print:  
	$(PRINT1) AMLdiag.inp main.f90 util.f90 inout.f90 random.f90 $(PRINT2)

clean:
	rm -f core *.mod *.o
