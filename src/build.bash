#

rm *.mod
rm *.o
rm *.vtk
gfortran -O3 -c  write_vtk2d.f90
gfortran -O3 -c  basic.f90
gfortran -O3 -c  main.f90
gfortran -O3 -o  lbm main.o basic.o write_vtk2d.o
./lbm
bash plot.bash
