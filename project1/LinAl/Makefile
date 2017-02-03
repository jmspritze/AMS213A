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
	$(FORTRAN_COMPILER) LinAl.f90 writeToScreen.f90 GaussElim.f90 LU.f90 -o LinAl 	
	
    
clean:
	rm -f  *.mod *.o
