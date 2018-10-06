SET PATH=C:\Users\hk\Documents\gfortran\bin;C:\Users\hk\Documents\gfortran\lib\gcc\i586-pc-mingw32\4.4.0;%PATH%
SET LIBRARY_PATH=C:\Users\hk\Documents\gfortran\lib;C:\Users\hk\Documents\gfortran\lib\gcc\i586-pc-mingw32\4.4.0

del *.o *.mod

gfortran -O2 -c modules.f90
gfortran -O2 -c watches.f90
gfortran -O2 -c artificial.f90
gfortran -O2 -c collectcells.f90
gfortran -O2 -c checkout.f90
gfortran -O2 -c diffschemes.f90
gfortran -O2 -c gradients.f90
gfortran -O2 -c initialfield.f90
gfortran -O2 -c lapacks.f
gfortran -O2 -c particles.f90
gfortran -O2 -c opendx.f90
gfortran -O2 -c patches.f90
gfortran -O2 -c patchscalars.f90
gfortran -O0 -c readcontrolfile1.f90
gfortran -O0 -c readcontrolfile2.f90
gfortran -O2 -c saverestart.f90
gfortran -O2 -c solverinterface.f90
gfortran -O2 -c solver_sparsekit2.f
gfortran -O2 -c tecplt.f90
gfortran -O2 -c gmv.f90
gfortran -O2 -c gmsh.f90
gfortran -O2 -c tools.f90
gfortran -O2 -c user.f90
gfortran -O2 -c vtk.f90

gfortran -O2 -o dolfyn.exe dolfyn.f90 *.o
gfortran -O2 -o doldge.exe preprocessor.f90
