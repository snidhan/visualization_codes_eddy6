#!/bin/bash
#ulimit -s unlimited
ulimit -s 999999999
#rm a.out *.vtk *.dat
#gfortran -O3 -mcmodel=medium -o a.out 3d_vtk.f90

# problem with -r8 -132 (can't see in paraview)
#ifort -O3 -r8 -132 -mcmodel=large -o a.out 3d_vtk.f90

#ifort -O3 -mcmodel=large -o a.out 3d_vtk.f90 -L/usr/local/LAPACK_KARU  -llapack -lblas
ifort -O3 -mcmodel=medium -heap-arrays 100000 -o a.out 3d_vtk_v2_1e4.f90 -L/opt/lapack-3.5.0 -llapack -L/opt/BLAS-3.6.0 -lblas /home/sheel/Work/system_utilities/tecplot10/lib/tecio64.a -lstdc++
#ifort -O3 -mcmodel=large  -o a.out 3d_vtk_v2_1e4.f90 -L/opt/lapack-3.6.1  -llapack -lblas -lm 
./a.out
