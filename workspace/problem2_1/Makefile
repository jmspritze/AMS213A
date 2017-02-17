# ============================================================================
# Name        : Makefile
# Author      : John Spritzer
# Version     :
# Copyright   : no copywrite
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran

all:
	$(FORTRAN_COMPILER) LinAl.f90 driver1.f90 -o Chol_LS 	
	
    
clean:
	rm -f  *.mod *.o
