#Name - Sheel Nidhan
#Date - 8 January 2020
 
 
programs = main.f90 read_restart.f90 readgrid.f90 centering.f90 truncation.f90 gaussian_filter_parameters.f90 filtering.f90
files = main.o read_restart.o readgrid.o centering.o truncation.o gaussian_filter_parameters.o filtering.o
 
OPT = -O3 
#DEBUG = -g -traceback
EXTRA = -heap-arrays 100000
#INCL= -I/opt/fftw-3.3.8/api
#LIB1 = -L/usr/local/lib -lfftw3 -lm 
#LIB2 = -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0 -lgfortran 
#LIB1 = -L/usr/local/lib -lfftw3 -lm 
#LIB2 =  -L/home/sheel/packages/lapack-3.8.0 -llapack -L/home/sheel/packages/BLAS-3.8.0 -lblas
 
 
fcomps= ifort $(OPT) $(EXTRA) $(DEBUG)
 
 # Making the Makefile
 
output: $(files) 
	$(fcomps) $(INCL) -o output $(files)
 
read_restart.o: read_restart.f90 
	$(fcomps) -c read_restart.f90

readgrid.o: readgrid.f90 
	$(fcomps) -c readgrid.f90

centering.o: centering.f90 
	$(fcomps) -c centering.f90

truncation.o: truncation.f90 
	$(fcomps) -c truncation.f90

gaussian_filter_parameters.o: gaussian_filter_parameters.f90
	$(fcomps) -c gaussian_filter_parameters.f90

filtering.o: filtering.f90 
	$(fcomps) -c filtering.f90

main.o: main.f90 
	$(fcomps) -c main.f90
 
 
 ## cleaning
clean: 
	rm -f *.o *.mod *.bin output                                     
