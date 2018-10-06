SHELL = /bin/sh

#F90 = ifc
#F90 = g95    
#F90 = ifort

F90 = gfortran
F90FLAGS = -O2   # -ffast-math -funroll-loops -O3

F77FILES = \
  solver_sparsekit2.f\
  lapacks.f  
# solver_bcg2.f
# solver_ma27.f

F90FILES=\
 modules.f90\
 watches.f90\
 patches.f90\
 initialfield.f90\
 saverestart.f90\
 artificial.f90\
 solverinterface.f90\
 readcontrolfile1.f90\
 readcontrolfile2.f90\
 patchscalars.f90\
 collectcells.f90\
 particles.f90\
 user.f90\
 tools.f90\
 opendx.f90\
 vtk.f90\
 tecplt.f90\
 gmv.f90\
 gmsh.f90\
 checkout.f90\
 diffschemes.f90\
 gradients.f90\
 dolfyn.f90

F77OBJS = ${F77FILES:.f=.o}
F90OBJS = ${F90FILES:.f90=.o}

%.o : %.mod 

##################################################################
# Targets
##################################################################

all:	dolgeo gmsh2dolfyn dolfyn 

dolgeo: preprocessor.f90
	@echo  "Compiling preprocessor" $@ "... "
	${F90} -o ../bin/dolgeo ${F90FLAGS} preprocessor.f90
	@rm *.mod

gmsh2dolfyn: gmsh2dolfyn.f90
	@echo  "Compiling gmsh2dolfyn" $@ "... "
	${F90} -o ../bin/gmsh2dolfyn ${F90FLAGS} gmsh2dolfyn.f90

dolfyn:	${F90OBJS} ${F77OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} ${F77OBJS} ${LFLAGS} ${INCS} -o ../bin/dolfyn 

clean:
	@rm *.o *.mod

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 .f 

.f90.o:
#	${F90} -O2 -c -Wuninitialized  ${INCS} ${@:.o=.f90}
#	${F90} -O2  -c -check bounds  -warn usage     ${INCS} ${@:.o=.f90}
#	${F90} -O0 -debug all  -check bounds -CB -traceback -c ${INCS} ${@:.o=.f90}
#	${F90} -O0 -g -c ${INCS} ${@:.o=.f90}
	${F90} -c ${F90FLAGS} ${INCS}  ${@:.o=.f90}
.f.o:
	${F90} -O2 -c ${INCS} ${@:.o=.f}

