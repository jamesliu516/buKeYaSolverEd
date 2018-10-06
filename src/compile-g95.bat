SET PATH=C:\Users\hk\Documents\g95\bin;C:\Users\hk\Documents\g95\lib\gcc-lib\i686-pc-mingw32\4.0.3;%PATH%
SET LIBRARY_PATH=C:\Users\hk\Documents\g95\lib;C:\Users\hk\Documents\g95\lib\gcc-lib\i686-pc-mingw32\4.0.3

del *.o *.mod

g95 -O2 -c modules.f90
g95 -O2 -c watches.f90
g95 -O2 -c artificial.f90
g95 -O2 -c collectcells.f90
g95 -O2 -c checkout.f90
g95 -O2 -c diffschemes.f90
g95 -O2 -c gradients.f90
g95 -O2 -c initialfield.f90
g95 -O2 -c lapacks.f
g95 -O2 -c particles.f90
g95 -O2 -c opendx.f90
g95 -O2 -c patches.f90
g95 -O2 -c patchscalars.f90
g95 -O0 -c readcontrolfile1.f90
g95 -O0 -c readcontrolfile2.f90
g95 -O2 -c saverestart.f90
g95 -O2 -c solverinterface.f90
g95 -O2 -c solver_sparsekit2.f
g95 -O2 -c tecplt.f90
g95 -O2 -c gmv.f90
g95 -O2 -c gmsh.f90
g95 -O2 -c tools.f90
g95 -O2 -c user.f90
g95 -O2 -c vtk.f90

g95 -O2 -o dolfyn.exe dolfyn.f90 *.o
g95 -O2 -o doldge.exe preprocessor.f90
