!
! Copyright 2003-2008 Henk Krus, Cyclone Fluid Dynamics BV
! All Rights Reserved.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
! http://www.dolfyn.net/license.html
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an
! "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied. See the License for the specific
! language governing permissions and limitations under the License.
!
! With contributions by B. Tuinstra
! (www.home.zonnet.nl/bouke_1/dolfyn)
!
subroutine CheckOutTests(Divergence)

   use constants
   use geometry
   use variables

   logical, intent(IN) :: Divergence
   logical             :: Failed, Passed, WallsPresent, QTransferPresent

   real, dimension(3)  :: ShearForce 

   if( Divergence )then
     write(IOdef,*) '*** Creating *.DIV file'
     call openfile(IOchk,casename,'.DIV','FORMATTED','APPEND','UNKNOWN',debug)
     write(IOchk,'(A)')'DIV'
     close(IOchk)
     return
   endif

   Passed = .true.   

   !
   ! now only volume averaged tests, future todo mass averaged
   !
   volume = sum(cell(1:Ncel)%vol)  
   
   Failed = .false.
   if( CheckU )then
     umin = minval( U(1:Ncel) )
     umax = maxval( U(1:Ncel) )
     if( umin < CheckLowU ) Failed = .true. 
     if( CheckUpU < umax  ) Failed = .true. 
     if( CheckReprtU )then
       write(IOdef,*) '*** Test ',umin,' <= U <= ',umax
       write(IOdbg,*) '*** Test ',umin,' <= U <= ',umax
       write(IOrun,*) '*** Test ',umin,' <= U <= ',umax
     endif

     if( CheckAverageBandU > -0.5 )then
       usum = sum( U(1:Ncel)*Cell(1:Ncel)%vol )
       uavr = usum/volume
       if( uavr < CheckAverageU-CheckAverageBandU .or. &
           uavr > CheckAverageU+CheckAverageBandU ) Failed = .true.
       if( CheckReprtU )then
         write(IOdef,*) '*** Test Average U ',uavr
         write(IOdbg,*) '*** Test Average U ',uavr
         write(IOrun,*) '*** Test Average U ',uavr
       endif
     endif
     
     if( Failed )then
       write(IOdef,*)'*** Test Variable U FAILED *** '
       write(IOdbg,*)'*** Test Variable U FAILED *** '
       write(IOrun,*)'*** Test Variable U FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.
   
   Failed = .false.
   if( CheckV )then
     vmin = minval( V(1:Ncel) )
     vmax = maxval( V(1:Ncel) )
     if( vmin < CheckLowV ) Failed = .true. 
     if( CheckUpV < vmax  ) Failed = .true. 
     if( CheckReprtV )then
       write(IOdef,*) '*** Test ',vmin,' <= V <= ',vmax
       write(IOdbg,*) '*** Test ',vmin,' <= V <= ',vmax
       write(IOrun,*) '*** Test ',vmin,' <= V <= ',vmax
     endif

     if( CheckAverageBandV > -0.5 )then
       vsum = sum( V(1:Ncel)*Cell(1:Ncel)%vol )
       vavr = vsum/volume
       if( vavr < CheckAverageV-CheckAverageBandV .or. &
           vavr > CheckAverageV+CheckAverageBandV ) Failed = .true.
       if( CheckReprtV )then
         write(IOdef,*) '*** Test Average V ',vavr
         write(IOdbg,*) '*** Test Average V ',vavr
         write(IOrun,*) '*** Test Average V ',vavr
       endif
     endif

     if( Failed )then
       write(IOdef,*)'*** Test Variable V FAILED *** '
       write(IOdbg,*)'*** Test Variable V FAILED *** '
       write(IOrun,*)'*** Test Variable V FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.
   
   Failed = .false.
   if( CheckW )then
     wmin = minval( W(1:Ncel) )
     wmax = maxval( W(1:Ncel) )
     if( wmin < CheckLowW ) Failed = .true. 
     if( CheckUpW < wmax  ) Failed = .true. 
     if( CheckReprtW )then
       write(IOdef,*) '*** Test ',wmin,' <= W <= ',wmax
       write(IOdbg,*) '*** Test ',wmin,' <= W <= ',wmax
       write(IOrun,*) '*** Test ',wmin,' <= W <= ',wmax
     endif

     if( CheckAverageBandW > -0.5 )then
       wsum = sum( W(1:Ncel)*Cell(1:Ncel)%vol )
       wavr = wsum/volume
       if( wavr < CheckAverageW-CheckAverageBandW .or. &
           wavr > CheckAverageW+CheckAverageBandW ) Failed = .true.
       if( CheckReprtW )then
         write(IOdef,*) '*** Test Average W ',wavr
         write(IOdbg,*) '*** Test Average W ',wavr
         write(IOrun,*) '*** Test Average W ',wavr
       endif
     endif

     if( Failed )then
       write(IOdef,*)'*** Test Variable W FAILED *** '
       write(IOdbg,*)'*** Test Variable W FAILED *** '
       write(IOrun,*)'*** Test Variable W FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.
   
   Failed = .false.
   if( CheckP )then
     pmin = minval( P(1:Ncel) )
     pmax = maxval( P(1:Ncel) )
     if( pmin < CheckLowP ) Failed = .true. 
     if( CheckUpP < pmax  ) Failed = .true. 
     if( CheckReprtP )then
       write(IOdef,*) '*** Test ',pmin,' <= P <= ',pmax
       write(IOdbg,*) '*** Test ',pmin,' <= P <= ',pmax
       write(IOrun,*) '*** Test ',pmin,' <= P <= ',pmax
     endif

     if( CheckAverageBandP > -0.5 )then
       psum = sum( P(1:Ncel)*Cell(1:Ncel)%vol )
       pavr = psum/volume
       if( pavr < CheckAverageP-CheckAverageBandP .or. &
           pavr > CheckAverageP+CheckAverageBandP ) Failed = .true.
       if( CheckReprtP )then
         write(IOdef,*) '*** Test Average P ',pavr
         write(IOdbg,*) '*** Test Average P ',pavr
         write(IOrun,*) '*** Test Average P ',pavr
       endif
     endif

     if( Failed )then
       write(IOdef,*)'*** Test Variable P FAILED *** '
       write(IOdbg,*)'*** Test Variable P FAILED *** '
       write(IOrun,*)'*** Test Variable P FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.
   
   Failed = .false.
   if( CheckTE )then
     temin = minval( TE(1:Ncel) )
     temax = maxval( TE(1:Ncel) )
     if( temin < CheckLowTE ) Failed = .true. 
     if( CheckUpTE < temax  ) Failed = .true. 
     if( CheckReprtTE )then
       write(IOdef,*) '*** Test ',temin,' <= TE <= ',temax
       write(IOdbg,*) '*** Test ',temin,' <= TE <= ',temax
       write(IOrun,*) '*** Test ',temin,' <= TE <= ',temax
     endif

     if( CheckAverageBandTE > -0.5 )then
       tesum = sum( TE(1:Ncel)*Cell(1:Ncel)%vol )
       teavr = tesum/volume
       if( teavr < CheckAverageTE-CheckAverageBandTE .or. &
           teavr > CheckAverageTE+CheckAverageBandTE ) Failed = .true.
       if( CheckReprtTE )then
         write(IOdef,*) '*** Test Average TE ',teavr
         write(IOdbg,*) '*** Test Average TE ',teavr
         write(IOrun,*) '*** Test Average TE ',teavr
       endif
     endif

     if( Failed )then
       write(IOdef,*)'*** Test Variable TE FAILED *** '
       write(IOdbg,*)'*** Test Variable TE FAILED *** '
       write(IOrun,*)'*** Test Variable TE FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.

   Failed = .false.
   if( CheckED )then
     edmin = minval( ED(1:Ncel) )
     edmax = maxval( ED(1:Ncel) )
     if( edmin < CheckLowED ) Failed = .true. 
     if( CheckUpED < edmax  ) Failed = .true. 
     if( CheckReprtED )then
       write(IOdef,*) '*** Test ',edmin,' <= ED <= ',edmax
       write(IOdbg,*) '*** Test ',edmin,' <= ED <= ',edmax
       write(IOrun,*) '*** Test ',edmin,' <= ED <= ',edmax
     endif

     if( CheckAverageBandED > -0.5 )then
       edsum = sum( ED(1:Ncel)*Cell(1:Ncel)%vol )
       edavr = edsum/volume
       if( edavr < CheckAverageED-CheckAverageBandED .or. &
           edavr > CheckAverageED+CheckAverageBandED ) Failed = .true.
       if( CheckReprtED )then
         write(IOdef,*) '*** Test Average ED ',edavr
         write(IOdbg,*) '*** Test Average ED ',edavr
         write(IOrun,*) '*** Test Average ED ',edavr
       endif
     endif

     if( Failed )then
       write(IOdef,*)'*** Test Variable ED FAILED *** '
       write(IOdbg,*)'*** Test Variable ED FAILED *** '
       write(IOrun,*)'*** Test Variable ED FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.

   Failed = .false.
   if( CheckT )then
     tmin = minval( T(1:Ncel) )
     tmax = maxval( T(1:Ncel) )
     if( tmin < CheckLowT ) Failed = .true. 
     if( CheckUpT < tmax  ) Failed = .true. 
     if( CheckReprtT )then
       write(IOdef,*) '*** Test ',tmin,' <= T <= ',tmax
       write(IOdbg,*) '*** Test ',tmin,' <= T <= ',tmax
       write(IOrun,*) '*** Test ',tmin,' <= T <= ',tmax
     endif

     if( CheckAverageBandT > -0.5 )then
       tsum = sum( T(1:Ncel)*Cell(1:Ncel)%vol )
       tavr = tsum/volume
       if( tavr < CheckAverageT-CheckAverageBandT .or. &
           tavr > CheckAverageT+CheckAverageBandT ) Failed = .true.
       if( CheckReprtT )then
         write(IOdef,*) '*** Test Average T ',tavr
         write(IOdbg,*) '*** Test Average T ',tavr
         write(IOrun,*) '*** Test Average T ',tavr
       endif
     endif

     if( Failed )then
       write(IOdef,*)'*** Test Variable T FAILED *** '
       write(IOdbg,*)'*** Test Variable T FAILED *** '
       write(IOrun,*)'*** Test Variable T FAILED *** '
     endif
   endif
   if( Failed ) Passed = .false.

   !
   ! loop over wall boundaries and report overall results
   !
   QtransferIn  = 0.0
   QtransferOut = 0.0
   ShearForce   = 0.0
   WallsPresent     = .false.
   QTransferPresent = .false.
   
   do ib=1,Nbnd
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RWall )then
       WallsPresent = .true.
       ShearForce = ShearForce + Bnd(ib)%shear 

       if( SolveEnthalpy )then
         QTransferPresent = .true.
         i = Bnd(ib)%face
         if( Bnd(ib)%q > 0.0 )then
           QtransferIn  = QtransferIn  + Bnd(ib)%q * Face(i)%area
         else
           QtransferOut = QtransferOut + Bnd(ib)%q * Face(i)%area
         endif
       endif
     endif
   end do

   if( WallsPresent )then
     write(IOdef,1) ' *** Test Shear Components  ',ShearForce
     write(IOdbg,1) ' *** Test Shear Components  ',ShearForce
     write(IOrun,1) ' *** Test Shear Components  ',ShearForce
     if( QTransferPresent )then
       write(IOdef,1) ' *** Test Heat Transfer IN  ',QtransferIn 
       write(IOdef,1) ' *** Test Heat Transfer OUT ',QtransferOUT
       write(IOdbg,1) ' *** Test Heat Transfer IN  ',QtransferIn 
       write(IOdbg,1) ' *** Test Heat Transfer OUT ',QtransferOUT
       write(IOrun,1) ' *** Test Heat Transfer IN  ',QtransferIn 
       write(IOrun,1) ' *** Test Heat Transfer OUT ',QtransferOUT
     endif
   endif
 1 format(A,3(1x,f12.6))
 
   if( Passed )then
     call openfile(IOchk,casename,'.OK','FORMATTED','APPEND','UNKNOWN',debug)
     write(IOchk,'(A)')'PASSED'
     close(IOchk)
   else
     call openfile(IOchk,casename,'.FAIL','FORMATTED','APPEND','UNKNOWN',debug)
     write(IOchk,'(A)')'FAILED'
     close(IOchk)     
   endif

end subroutine CheckOutTests
