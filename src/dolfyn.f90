!
! Copyright 2003-2010 Henk Krus, Cyclone Fluid Dynamics BV
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
program dolfyn

   use constants
   use geometry
   use variables
   use scalars
   use watches

   real    :: TimeElapsed, TimeNew
   logical :: Ready, Exists, Divergence

   Exists      = .false.
   Divergence  = .false.

   call watch_setup(40)
   call cpu_time( TimeElapsed )
   TimeNew = -1.0

   call ReadCase
   !
   ! first read geometry (needed to set certain constants)
   !
   call ReadGeometry

   call ReadControlFile

   if( UsePatches ) call PatchesSetUp(1)

   !
   ! initialisation and calculate geometry data
   !
   call InitialiseVariables

   !call InitializeProperties

   if( UsePatches   ) call PatchesSetUp(2)  ! stage 2, variables allocated
   if( UseParticles ) call ParticleSetUp
   if( UseSensors   ) call SensorsSetUp

   call SetBoundaryConditions

   call PrettyPrint(IOdef)
   call PrettyPrint(IOrun)
   call PrettyPrint(IOdbg)

   if( count( Solver(:) == SparseKit2 ) > 0 )then
     write(IOdef,*)'Allocating SparsKit2 Work array'
     allocate( Work(NNZ,8),stat=istat)
     call TrackMemory(istat,(NNZ)*8*2,'Work array allocated')
   endif

   IterStart = 0
   Time      = 0.0
   dppmax    = 0.0
   
 !  if ( BoolShearRortex>0 ) then
 !      write (*,*) CoefShearRortex(1), CoefShearRortex(2), CoefShearRortex(2)
 !  endif

   if( Restart > 0 )then

     write(IOdef,*)'Reading restart data, ',Restart
     write(IOdbg,*)'Reading restart data'
     write(IOrun,*)'Reading restart data'
     call ReadRestartField(IterStart)

     if( Restart == 2 .or. Restart == 3 )then
       write(IOdef,*) 'Overriding boundary conditions'
       write(IOdbg,*) 'Overriding boundary conditions'
       write(IOrun,*) 'Overriding boundary conditions'

       IterStart = 0
       Time      = 0.0
       call SetBoundaryConditions
       call Set_Normalisation_Factors
       if( Debug > 0 )write(IOdef,*) 'Opening residuals file (1)'
       call openfile(IOres,Casename,'.res','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
       call Set_Normalisation_Factors

       call InitialField(dppmax)

       if( Debug > 0 )write(IOdef,*) 'Opening monitor file (1)'
       call openfile(IOmon,Casename,'.mon','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)

     else if( Restart == 4 .or.  Restart == 5 )then
       write(IOdef,*) 'Ignoring stored boundary conditions'
       write(IOdbg,*) 'Ignoring stored boundary conditions'
       write(IOrun,*) 'Ignoring stored boundary conditions'

       write(IOdef,*) 'Restart ',Restart
       write(IOdbg,*) 'Restart ',Restart
       write(IOrun,*) 'Restart ',Restart

       IterStart = 0
       Time      = 0.0
       call SetBoundaryConditions
       call Set_Normalisation_Factors

       call openfile(IOres,Casename,'.res','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
       call openfile(IOmon,Casename,'.mon','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)

     else
       if( Debug > 0 )write(IOdef,*) 'Opening residuals file (2)'
       call openfile(IOres,Casename,'.res','FORMATTED','APPEND','UNKNOWN',Debug)
       call Set_Normalisation_Factors

       if( Debug > 0 )write(IOdef,*) 'Opening monitor file (2)'
       call openfile(IOmon,Casename,'.mon','FORMATTED','APPEND','UNKNOWN',Debug)
     endif

   else
     if( Debug > 0 )write(IOdef,*) 'Opening residuals file (3) ',Restart
     call openfile(IOres,Casename,'.res','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
     call Set_Normalisation_Factors

     call InitialField(dppmax)

     if( Debug > 0 )write(IOdef,*) 'Opening monitor file (3)'
     call openfile(IOmon,Casename,'.mon','FORMATTED','SEQUENTIAL','UNKNOWN',Debug)
   endif

   !
   ! the big time loop
   !
   if( Debug > -1 )then
     write(IOdef,'(//,'' Go! '',/)')
     write(IOrun,'(//,'' Go! '',/)')
   endif

   TimeSteps: do iter=IterStart+1,IterStart+Niter

     !
     ! reset the array of konsole output flags
     !
     Flags(1:NFlags) = ' '

     Time = Time + dt

     if( Debug > 1 )then
       write(IOdef,'(1x,A)')     '****************************'
       if( Transient )then
         write(IOdef,'(1x,A,i5,1x,f12.5)')  '*** Time step:  ',iter,time
       else
         write(IOdef,'(1x,A,i5)')  '*** Iteration:  ',iter
       endif
       write(IOdef,'(1x,A)')     '****************************'
     endif
     write(IOdbg,'(1x,A)')   '**********************************************'
     write(IOdbg,'(1x,A,i5,1x,f12.5)')  '*** Iteration/Time step:  ',iter,time
     write(IOdbg,'(1x,A)')   '**********************************************'

     !
     ! store old solutions etc
     !
     if( Transient )then
       if( Debug > 2 )write(IOdef,*) 'Storing old values'
       if( QuadTime )then
         Uold2 = Uold
         Vold2 = Vold
         Wold2 = Wold
       endif
       Uold = U
       Vold = V
       Wold = W
       if( SolveEnthalpy .and. QuadTime ) Told2 = Told
       if( SolveEnthalpy ) Told = T
       if( SolveTurb )then
         if( QuadTime )then
           TEold2 = TEold
           EDold2 = EDold
         endif
         TEold = TE
         EDold = ED
       endif
       if( SolveScalars )then
         if( QuadTime ) SCold2 = SCold
         SCold = SC
       endif
       if( MovingGrid )then
         write(IOdef,*)'Not implemented yet'
       endif
       call Set_Normalisation_Factors  ! inlet might be changed
     endif

     iouter   =  0
     ready    = .false.
     OuterSteps: do while( .not. ready )

       iouter = iouter + 1
       !
       ! velocity gradients
       !
       call GradientPhi(VarU,U,dUdX)      !
       call GradientPhi(VarV,V,dVdX)      ! supporting stuff
       call GradientPhi(VarW,W,dWdX)      ! needed at various places
       call GradientPhi(VarP,P,dPdX)      !

       if( SolveUVW) call CalculateUVW

       if( SolveP )  call CalculatePressure(dppmax)

       call UpdateBC

       if( SolveTurbEnergy ) call CalculateScalar(VarTE,TE,dpdx,TEold,TEold2)
       if( SolveTurbDiss )   call CalculateScalar(VarED,ED,dpdx,EDold,EDold2)
       if( SolveVisc )       call CalculateViscosity

       if( SolveEnthalpy )   call CalculateScalar(VarT,T,dpdx,Told,Told2)

       if( SolveScalars )then
         do is=1,Nscal
           if( ScProp(is)%solve ) &
           call CalculateScalar(VarS(is),SC(:,is),dpdx,SCold(:,is),SCold2(:,is))
         end do
       endif

       !
       ! fluid properties (den,cp,vis) will be set here
       !

       !if( UseParticles )    call ParticleInterface ! only when coupled
       !if( UseSensors   )    call EvaluateSensors   ! only when switched on

       call PrintSummary(IOdbg,iouter)

       !
       ! stopping criteria
       !
       if( iouter >= MaxOUTER .or. (.not. Transient) ) ready = .true.

       tmpres = Residual(6)   ! do not use the epsilon-residual
       Residual(6) = Small    ! because its nature is arbitrary

       if( maxval(Residual(1:8)) < ResMax )then
         if( Debug > 1 .or. .not. Transient )then
           write(IOdef,'(/)')
           write(IOdef,*)'*** CONVERGENCE ***',maxval(Residual(1:8))
           write(IOrun,'(/)')
           write(IOrun,*)'*** CONVERGENCE ***',maxval(Residual(1:8))
         endif
         Residual(6) = tmpres
         if( Transient )then
           exit OuterSteps
         else
           exit TimeSteps
         endif
       endif

       if( maxval(Residual(1:8)) > 1.e+18 )then
         write(IOdef,'(/)')
         write(IOdef,*)'*** DIVERGENCE ***',maxval(Residual(1:8))
         write(IOrun,'(/)')
         write(IOdbg,*)'*** DIVERGENCE ***',maxval(Residual(1:8))
         write(IOdbg,'(/)')
         write(IOrun,*)'*** DIVERGENCE ***',maxval(Residual(1:8))
         Divergence = .true.
         exit TimeSteps
       endif

       Residual(6) = tmpres

     end do OuterSteps

     call PrintSummary(IOdef,iouter)
     call PrintSummary(IOrun,iouter)

     !write(IOext,*) TIME, Qtransfer

     if( Debug > 2 ) write(IOdef,*)'Maximum change in PP:', &
                                   dppmax,maxval(Residual(1:8)),ResMax

     !if( dppmax < Small .and. iter > 20 .and. .not. Transient) EXIT timesteps

     !
     ! intermediate opendx debug dump
     !
     !if( mod(Iter,DXdump) == 0 .and. Iter /= IterStart+Niter )then
     !
     !  call dolfyn2opendx( Divergence )
     !
     !endif
     !
     ! intermediate output data options
     !
     if( mod(Iter,NOutput) == 0 .and. Iter /= IterStart+Niter )then
       write(IOdef,*)'Writing postprocessing file(s)'
       if( UseOpenDX  ) call dolfyn2opendx( Divergence )
       if( UseVTK     ) call dolfyn2vtk
       if( UseTECPLOT ) call dolfyn2tecplt(0)
       if( UseGMV     ) call dolfyn2gmv
       if( UseGMSH    ) call dolfyn2gmsh
     endif

     if( Iter == IOutput)then
       write(IOdef,*)'Writing postprocessing file(s) at iteration ',Iter
       if( UseOpenDX  ) call dolfyn2opendx( Divergence )
       if( UseVTK     ) call dolfyn2vtk
       if( UseTECPLOT ) call dolfyn2tecplt(0)
       if( UseGMV     ) call dolfyn2gmv
       if( UseGMSH    ) call dolfyn2gmsh
     endif

     if( TOutput /= -1. .and. Time >= TOutput .and. .not. TOutputDone )then
       write(IOdef,*)'Writing postprocessing file(s) at time ',Time
       if( UseOpenDX  ) call dolfyn2opendx( Divergence )
       if( UseVTK     ) call dolfyn2vtk
       if( UseTECPLOT ) call dolfyn2tecplt(0)
       if( UseGMV     ) call dolfyn2gmv
       if( UseGMSH    ) call dolfyn2gmsh
       TOutputDone = .true.
     endif
     !
     ! saving restart data options
     !
     if( mod(Iter,NSave) == 0 .and. Iter /= IterStart+Niter )then

       call WriteRestartField

     endif

     if( Iter == ISave )then
       write(IOdef,*)'Saving restart data at iteration ',Iter
       call WriteRestartField
     endif

     if( TSave > 0.0 .and. Time >= Tsave .and. .not. TSaveDone )then
       write(IOdef,*)'Saving restart data at time ',Time
       call WriteRestartField
       TSaveDone = .true.
     endif

     if( TCPU > 0.0 )then
       if( TimeNew < 0.0 ) TimeNew = TCPU              ! if not set, do it
       call cpu_time( TimeElapsed )                    ! get current time

       if( TimeElapsed >= TimeNew )then
         if( CPUstop )then
           write(IOdef,*)'*** CPU time exceeded ***'
           write(IOdbg,*)'*** CPU time exceeded ***'
           write(IOrun,*)'*** CPU time exceeded ***'
           exit timesteps
         else
           TimeNew = TimeNew + TCPU                      ! set new save time
           write(IOdef,*)'Saving restart data CPU time ',TimeElapsed
           call WriteRestartField
         endif
       endif
     endif
     !
     ! trick to stop a running job (all cases in directory)
     !
     inquire(file='STOP',exist=Exists)
     if( Exists )then
       write(IOdef,*)'*** OK stop ***'
       write(IOdbg,*)'*** OK stop ***'
       write(IOrun,*)'*** OK stop ***'
       exit timesteps
     endif

     inquire(file=casename(1:lens(casename))//'.STOP',exist=Exists)
     if( Exists )then
       write(IOdef,*)'*** OK stop ***'
       write(IOdbg,*)'*** OK stop ***'
       write(IOrun,*)'*** OK stop ***'
       exit timesteps
     endif
     !
     ! option to fix values for ABL
     !
     if( UseFixABL )then
      !write(IOdef,*) 'Fixing values in ABL citd ',IdFixABL
       write(IOdbg,*) 'Fixing values in ABL citd ',IdFixABL
       write(IOdbg,*) '=>',UFixABL,VFixABL,WFixABL,TeFixABL,EdFixABL
       do i=1,Ncel
         if( cell(i)%ctid == IdFixABL )then
           U(i)  = UFixABL
           V(i)  = VFixABL
           W(i)  = WFixABL
           TE(i) = TeFixABL
           ED(i) = EdFixABL
         endif
       end do
     endif
     !
     ! warning: 'flush' is not a standard call
     !
     call flush(IOdbg)
     call flush(IOrun)

   end do timesteps

   if( Iter == IterStart+Niter+1 )then
     iter = iter - 1
     write(IOdef,'(/)')
     write(IOdef,*)'Number of requested steps done. ',Niter
     write(IOdef,'(/)')
     write(IOdbg,'(/)')
     write(IOdbg,*)'Number of requested steps done. ',Niter
     write(IOdbg,'(/)')
     write(IOrun,'(/)')
     write(IOrun,*)'Number of requested steps done. ',Niter
     write(IOrun,'(/)')
   endif

   if( .not. Divergence ) call WriteRestartField

   if( PrintCellVar .or. PrintWallVar ) call ShowCells

   !
   ! clean up a bit before the external interfaces are called
   !
   if( allocated( Work ) )then
     deallocate( Work )
     call TrackMemory(istat,-NNZ*8*2,'SparsKit2 work array deallocated')
   endif

   if( UseParticles ) call ParticleInterface      ! only passive particles!
   if( UseSensors   ) call EvaluateSensors


   if( UseOpenDX  ) call dolfyn2opendx( Divergence )
   if( UseVTK     ) call dolfyn2vtk
   if( UseTECPLOT ) call dolfyn2tecplt(1)
   if( UseGMV     ) call dolfyn2gmv
   if( UseGMSH    ) call dolfyn2gmsh

   if( Debug > -1 ) write(IOdef,*)'Maximum change in PP:', &
                                  dppmax,maxval(Residual(1:8)),ResMax

   !
   ! checkout tests
   !
   if( CheckOut ) call CheckOutTests( Divergence )

   call CleanUp

   call watch_leave
   call watch_print(IOdef)
   call watch_print(IOdbg)
   call watch_print(IOrun)
   call watch_destroy

   write(IOdef,*)'Done ',casename(1:lens(casename))
   write(IOdbg,*)'Done ',casename(1:lens(casename))
   write(IOrun,*)'Done ',casename(1:lens(casename))

end program Dolfyn
subroutine ReadGeometry

   use constants
   use geometry
   use watches

   real, dimension(3)     :: Xp, Xn, X
   real, dimension(3)     :: X1, X2, X3, Xcg
   real, dimension(3)     :: Xpn, Xnp, Xnn, Xpac, Xnac, Xpac2, Xnac2
   real, dimension(3)     :: Xnorm, Xface

   integer, allocatable, dimension(:) :: TmpCellNodes

   integer, parameter     :: Ncg = 30                 ! Waar was dit voor???
   real, dimension(Ncg)   :: TetVol
   real, dimension(Ncg,3) :: TetSurf, TetNorm, TetCG

   real :: pi = 3.1415927

   character(len=32) string
   integer      geoversion
   logical      bin

   logical      warning

   logical      found 
   
   call watch_enter('ReadGeometry')
   !
   ! extra information will be dumped in this file:
   !
   call openfile(IOdbg,casename,'.dbg','FORMATTED','APPEND','UNKNOWN',debug)

   call Disclaimer

   !
   ! inquire for the type of dolfyn geometry file
   !
   call feelfile(IOgeo,IFORM,casename,Ext(IOgeo),'FORMATTED','SEQUENTIAL','OLD',debug)

   if( IFORM == 0 )then
     bin = .true.
   else
     bin = .false.
   endif
   !
   ! read dolfyn files
   !
   write(IOdef,*) 'Opening geometry file'
   write(IOdbg,*) 'Opening geometry file'
   write(IOrun,*) 'Opening geometry file'

   if( bin )then
     call openfile(IOgeo,casename,Ext(IOgeo),'UNFORMATTED','SEQUENTIAL','OLD',debug)
   else
     call openfile(IOgeo,casename,Ext(IOgeo),'FORMATTED','SEQUENTIAL','OLD',debug)
   endif

   if( bin ) then
     !
     ! binary
     !
     read(IOgeo) string(1:12)
     if( string(1:12) /= 'dolfyn bin g' )then
       write(IOdef,*)'*** Error: Not a dolfyn geometry file'
       stop
     endif
     read(IOgeo) geoversion
     if( geoversion /= version )then
       write(IOdef,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOdbg,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOrun,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
     endif
     read(IOgeo) ScaleFactor
     write(IOdef,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOdbg,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOrun,*) 'Using ScaleFactor = ',ScaleFactor

     read(IOgeo) string(1:4)
     if( string(1:3) /= 'ce:' )then
       write(IOdef,*)'*** Error: Cell definitions expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     read(IOgeo) icel
     write(IOdef,*)'Reading cells ',icel,' BINARY'
     write(IOdbg,*)'Reading cells ',icel,' BINARY'
     write(IOrun,*)'Reading cells ',icel,' BINARY'

     Ncel = icel
     allocate(Cell(Ncel),stat=istat)
     call TrackMemory(istat,Ncel,'Cell array allocated')

     do i=1,Ncel
       read(IOgeo) j,cell(j)%ctid,cell(j)%x,cell(j)%vol
       if( j /= i ) write(IOdef,*)'*** Warning: cell array corrupted'
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a12)') string(1:12)
     if( string(1:12) /= 'dolfyn asc g' )then
       write(IOdef,*)'*** Error: Not a dolfyn geometry file'
       stop
     endif
     read(IOgeo,'(8x,i6)') geoversion
     if( geoversion /= version )then
       write(IOdef,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOdbg,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOrun,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
     endif
     read(IOgeo,'(6x,e12.5)') ScaleFactor
     write(IOdef,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOdbg,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOrun,*) 'Using ScaleFactor = ',ScaleFactor

     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'ce:' )then
       write(IOdef,*)'*** Error: Cell definitions expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     read(IOgeo,*) icel
     write(IOdef,*)'Reading cells ',icel
     write(IOdbg,*)'Reading cells ',icel
     write(IOrun,*)'Reading cells ',icel

     Ncel = icel
     allocate(Cell(Ncel),stat=istat)
     call TrackMemory(istat,Ncel,'Cell array allocated')

     do i=1,Ncel
       read(IOgeo,'(2(i8,1x),4(1pe16.9,1x))') &
                                    j,cell(j)%ctid,cell(j)%x,cell(j)%vol
       if( j /= i ) write(IOdef,*)'*** Warning: cell array corrupted'
     end do
   endif

   !
   ! apply scalefactor
   !
   Cell(1:Ncel)%x(1) = ScaleFactor * Cell(1:Ncel)%x(1)
   Cell(1:Ncel)%x(2) = ScaleFactor * Cell(1:Ncel)%x(2)
   Cell(1:Ncel)%x(3) = ScaleFactor * Cell(1:Ncel)%x(3)

   ScaleFactor3 = ScaleFactor**3

   Cell(1:Ncel)%vol = ScaleFactor3 * Cell(1:Ncel)%vol

   !
   ! some output
   !
   write(Iodef,'(1x,A,i2,A,i2)') 'Cell types from ', &
                       minval(cell(:)%ctid),' to ',maxval(cell(:)%ctid)
   write(IOdbg,'(1x,A,i2,A,i2)') 'Cell types from ', &
                       minval(cell(:)%ctid),' to ',maxval(cell(:)%ctid)
   write(IOrun,'(1x,A,i2,A,i2)') 'Cell types from ', &
                       minval(cell(:)%ctid),' to ',maxval(cell(:)%ctid)

   if( debug > 0 )then
     write(IOdef,'(1x,A)') 'Cell centers are:'
     write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
       'Dimension ',minval(cell(:)%x(1)),' <  x  < ',maxval(cell(:)%x(1))
     write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
       'Dimension ',minval(cell(:)%x(2)),' <  y  < ',maxval(cell(:)%x(2))
     write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
       'Dimension ',minval(cell(:)%x(3)),' <  z  < ',maxval(cell(:)%x(3))
     write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
       'Volumes   ',minval(cell(:)%vol), ' < Vol < ',maxval(cell(:)%vol)
     write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
       'Total volume: ',sum(cell(:)%vol)
   endif
   write(IOdbg,'(1x,A)') 'Cell centers are:'
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(cell(:)%x(1)),' <  x  < ',maxval(cell(:)%x(1))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(cell(:)%x(2)),' <  y  < ',maxval(cell(:)%x(2))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(cell(:)%x(3)),' <  z  < ',maxval(cell(:)%x(3))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Volumes   ',minval(cell(:)%vol), ' < Vol < ',maxval(cell(:)%vol)
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Total volume: ',sum(cell(:)%vol)

   !
   ! list of cell faces
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'cf: ' )then
       write(IOdef,*)'*** Error: List of cell faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading list of cell faces'

     read(IOgeo) maxfaces
     if( maxfaces > 6 )then
       write(IOdef,*)'*** Warning: Unsupported feature'
     endif
     allocate(NFaces(Ncel),stat=istat)
     call TrackMemory(istat,Ncel,'NFaces array allocated')

     allocate(CFace(Ncel,12),stat=istat)
     call TrackMemory(istat,Ncel*12,'CFace array allocated')

     do i=1,Ncel
       read(IOgeo) j,NFaces(j)
       do k=1,Nfaces(j)
         read(IOgeo) CFace(j,k)
       end do
       if( j /= i ) write(IOdef,*)'*** Warning: cell face array corrupted'
     end do

   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'cf:' )then
       write(IOdef,*)'*** Error: List of cell faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading list of cell faces'

     read(IOgeo,*) maxfaces
     if( maxfaces > 6 )then
       write(IOdef,*)'*** Warning: Unsupported feature'
     endif
     allocate(NFaces(Ncel),stat=istat)
     call TrackMemory(istat,Ncel,'NFaces array allocated')

     allocate(CFace(Ncel,12),stat=istat)
     call TrackMemory(istat,Ncel*12,'CFace array allocated')

     do i=1,Ncel
       read(IOgeo,'(12(i8,1x))') j,NFaces(j),(CFace(j,k),k=1,NFaces(j))
       if( j /= i ) write(IOdef,*)'*** Warning: cell face array corrupted'
     end do
   endif

   !
   ! the faces
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:3) /= 'fc: ' )then
       write(IOdef,*)'*** Error: Cell faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif

     read(IOgeo) Nfac
     allocate(Face(Nfac),stat=istat)
     call TrackMemory(istat,Nfac*16,'Face array allocated')

     allocate(RFace(Nfac,2),stat=istat)
     call TrackMemory(istat,Nfac*2,'RFace array allocated')

     write(IOdef,*)'Reading cell faces ',Nfac
     write(IOdbg,*)'Reading cell faces ',Nfac
     write(IOrun,*)'Reading cell faces ',Nfac

     !
     ! read explicitely because in dolfyn the face type is
     ! now augmented by extra elements, Face%area has been
     ! taken out. later the vertices will have to go into
     ! a compressed array of their own
     !
     do i=1,Nfac
       read(IOgeo) j, Face(j)%bnd,     &
                      Face(j)%cell1,   &
                      Face(j)%cell2,   &
                      Face(j)%vertices,&
                      Face(j)%n,       &
                      Face(j)%x,       &
                      Face(j)%lambda
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'fc:' )then
       write(IOdef,*)'*** Error: Cell faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif

     read(IOgeo,*) Nfac
     allocate(Face(Nfac),stat=istat)
     call TrackMemory(istat,Nfac*16,'Face array allocated')

     allocate(RFace(Nfac,2),stat=istat)
     call TrackMemory(istat,Nfac*2,'RFace array allocated')

     write(IOdef,*)'Reading cell faces ',Nfac
     write(IOdbg,*)'Reading cell faces ',Nfac
     write(IOrun,*)'Reading cell faces ',Nfac

     !
     ! read explicitely because in dolfyn the face type is
     ! now augmented by extra elements, Face%area has been
     ! taken out. later the verices will have to go into
     ! a compressed array of their own
     !
     do i=1,Nfac
       read(IOgeo,'(8(i8,1x),8(1pe16.9,1x))') &
                   j, Face(j)%bnd,     &
                      Face(j)%cell1,   &
                      Face(j)%cell2,   &
                      Face(j)%vertices,&
                      Face(j)%n,       &
                      Face(j)%x,       &
                      Face(j)%lambda
     end do
   endif

   write(IOdef,*)'Read cell faces ',Nfac
   write(IOdbg,*)'Read cell faces ',Nfac
   write(IOrun,*)'Read cell faces ',Nfac
   !
   ! rapidly calculate the face area
   !
   do i=1,Nfac
     Face(i)%area = sqrt( dot_product(Face(i)%n,Face(i)%n) )
   end do

   write(IOdef,*)'Cell face area',Nfac,ScaleFactor
   write(IOdbg,*)'Cell face area',Nfac,ScaleFactor
   write(IOrun,*)'Cell face area',Nfac,ScaleFactor
   !
   ! apply scalefactor
   !
   ScaleFactor2 = ScaleFactor**2

write(*,*) '..1'
   Face(1:Nfac)%area = ScaleFactor2 * Face(1:Nfac)%area

write(*,*) '..2'
   Face(1:Nfac)%n(1) = ScaleFactor2 * Face(1:Nfac)%n(1)
   Face(1:Nfac)%n(2) = ScaleFactor2 * Face(1:Nfac)%n(2)
   Face(1:Nfac)%n(3) = ScaleFactor2 * Face(1:Nfac)%n(3)

write(*,*) '..3'
   Face(1:Nfac)%x(1) = ScaleFactor * Face(1:Nfac)%x(1)
   Face(1:Nfac)%x(2) = ScaleFactor * Face(1:Nfac)%x(2)
   Face(1:Nfac)%x(3) = ScaleFactor * Face(1:Nfac)%x(3)

   write(IOdef,*)'Scaled ',ScaleFactor,ScaleFactor2
   write(IOdbg,*)'Scaled ',ScaleFactor,ScaleFactor2
   write(IOrun,*)'Scaled ',ScaleFactor,ScaleFactor2
   !
   ! boundaries
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'bn: ' )then
       write(IOdef,*)'*** Error: Boundary faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif

     read(IOgeo) Nbnd                           
     allocate(Bnd(Nbnd),stat=istat)
     call TrackMemory(istat,Nbnd,'Boundary array allocated')

     write(IOdef,*)'Reading boundaries ',Nbnd
     write(IOdbg,*)'Reading boundaries ',Nbnd
     write(IOrun,*)'Reading boundaries ',Nbnd

     do i=1,Nbnd
       read(IOgeo) j,bnd(j)%face, &
                     bnd(j)%vertices, bnd(j)%rid, bnd(j)%distance
       bnd(j)%yplus =  0.0
       bnd(j)%uplus =  0.0
       bnd(j)%shear =  0.0
       bnd(j)%h     =  0.0
       bnd(j)%q     =  0.0
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'bn:' )then
       write(IOdef,*)'*** Error: Boundary faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif

     read(IOgeo,*) Nbnd                          
     allocate(Bnd(Nbnd),stat=istat)
     call TrackMemory(istat,Nbnd,'Boundary array allocated')

     write(IOdef,*)'Reading boundaries ',Nbnd
     write(IOdbg,*)'Reading boundaries ',Nbnd
     write(IOrun,*)'Reading boundaries ',Nbnd

     do i=1,Nbnd
       read(IOgeo,'(7(i8,1x),1(e12.5,1x))') j,bnd(j)%face, &
                     bnd(j)%vertices, bnd(j)%rid, bnd(j)%distance
       bnd(j)%yplus =  0.0
       bnd(j)%uplus =  0.0
       bnd(j)%shear =  0.0
       bnd(j)%h     =  0.0
       bnd(j)%q     =  0.0
     end do
   endif
   !
   ! apply scalefactor
   !
   Bnd(1:Nbnd)%distance = ScaleFactor * Bnd(1:Nbnd)%distance

   Nreg = maxval( Bnd(:)%rid )
   write(IOdef,*) 'Maximum region number found: ',Nreg
   write(IOdbg,*) 'Maximum region number found: ',Nreg
   write(IOrun,*) 'Maximum region number found: ',Nreg
   allocate(Reg(0:Nreg),stat=istat)                        ! note: 0 <= Nreg !!!
   call TrackMemory(istat,(1+Nreg)*17, &
     'Boundary region definitions array allocated')

   !
   ! finally the vertices
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'vr: ' )then
       write(IOdef,*)'*** Error: Vertices expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif

     read(IOgeo) Nvrt
     allocate(Vert(Nvrt,3),stat=istat)
     call TrackMemory(istat,Nvrt,'Vertices array allocated')

     write(IOdef,*)'Reading Vertices ',Nvrt
     write(IOdbg,*)'Reading Vertices ',Nvrt
     write(IOrun,*)'Reading Vertices ',Nvrt

     do i=1,Nvrt
       read(IOgeo) j,(vert(j,k),k=1,3)
     end do
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'vr:' )then
       write(IOdef,*)'*** Error: Vertices expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif

     read(IOgeo,*) Nvrt
     allocate(Vert(Nvrt,3),stat=istat)
     call TrackMemory(istat,Nvrt,'Vertices array allocated')

     write(IOdef,*)'Reading Vertices ',Nvrt
     write(IOdbg,*)'Reading Vertices ',Nvrt
     write(IOrun,*)'Reading Vertices ',Nvrt

     do i=1,Nvrt
       read(IOgeo,'(1(i8,1x),3(1pe16.9,1x))') j,(vert(j,k),k=1,3)
     end do
   endif
   !
   ! apply scalefactor
   !
   Vert = ScaleFactor * Vert

   !
   ! stored region names
   !
   if( bin )then
     read(IOgeo) string(1:4)
     if( string(1:4) /= 'rd: ' )then
       write(IOdef,*)'*** Error: Region Data expected'
       write(IOdef,*)'    Found: ',string(1:3)
     else
     write(IOdbg,*)'Reading Region Data'

     read(IOgeo) k
     do i=0,k
       read(IOgeo) j,Reg(j)%name
     end do
     endif
   else
     !
     ! ascii
     !
     read(IOgeo,'(a3)') string
     if( string(1:3) /= 'rd:' )then
       write(IOdef,*)'*** Error: Region Data expected'
       write(IOdef,*)'    Found: ',string(1:3)
     else
     write(IOdbg,*)'Reading Region Data'
     read(IOgeo,*) k
     do i=0,k
       read(IOgeo,'(i8,1x,A32)') j,Reg(j)%name
     end do
     endif
   endif

   !
   ! *** consistency checks ***
   !
   write(IOdef,*) 'Geometry checks'
   write(IOdbg,*) 'Geometry checks'
   write(IOrun,*) 'Geometry checks'
   !
   ! check volume using Gauss (see eq. 8.42)
   !
   vol1 = 0.0
   vol2 = 0.0
   vol3 = 0.0
   do i=1,Ncel

     Xp = Cell(i)%x
     v1 = 0.0
     v2 = 0.0
     v3 = 0.0

     do j=1,NFaces(i)
       k  = CFace(i,j)
       ip = Face(k)%cell1
       in = Face(k)%cell2

       if( ip == i )then
         v1 = v1 + Face(k)%x(1) * Face(k)%n(1)
         v2 = v2 + Face(k)%x(2) * Face(k)%n(2)
         v3 = v3 + Face(k)%x(3) * Face(k)%n(3)
       else if( in == i )then
         v1 = v1 - Face(k)%x(1) * Face(k)%n(1)
         v2 = v2 - Face(k)%x(2) * Face(k)%n(2)
         v3 = v3 - Face(k)%x(3) * Face(k)%n(3)
       else
         write(IOdef,*)'Error in check volume consistency'
       endif
     end do
     vol1 = vol1 + v1
     vol2 = vol2 + v2
     vol3 = vol3 + v3
   end do
   if( abs(vol1-vol2) > 1.e-3 * vol1 .or. &
       abs(vol1-vol3) > 1.e-3 * vol1 .or. &
       abs(vol2-vol3) > 1.e-3 * vol1 )then
     write(IOdef,*) 'WARNING: Check volumes'
     write(IOdbg,*) 'WARNING: Check volumes'
     write(IOrun,*) 'WARNING: Check volumes'
   endif

   write(IOdbg,*)'Volume based on xi:',vol1
   write(IOdbg,*)'Volume based on yj:',vol2
   write(IOdbg,*)'Volume based on zk:',vol3
   !
   ! check non-orthogonality
   !
   anglemx  =  0.0
   anglemn  = 90.0
   anglelm  = 70.0
   iangle   =   0

   do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )then
       ! internal face
       Xp = Cell(ip)%x   ! cell centre face owner
       Xn = Cell(in)%x   ! cell centre neigbor

       Xp = Xn - Xp      ! vector from P to N
       Xn = Face(i)%n    ! surface vector

       call normalise(Xp)
       call normalise(Xn)

       dotp  = dot_product( Xp , Xn )
       ! BT, Visual Fortran does not like acos(1.0)
       if(dotp > 0.99999 ) dotp = 0.99999
       angle = acos( dotp )*180./pi
     else
       Xp = Cell(ip)%x   ! cell centre face owner
       Xn = Face(i)%x    ! face centre

       Xp = Xn - Xp      ! vector from P to N
       Xn = Face(i)%n    ! surface vector

       call normalise(Xp)
       call normalise(Xn)

       dotp  = dot_product( Xp , Xn )
       ! BT, Visual Fortran does not like acos(1.0)
       if(dotp > 0.99999 ) dotp = 0.99999
       angle = acos( dotp )*180./pi
     endif
     anglemx =  max(anglemx,angle)
     anglemn =  min(anglemn,angle)
     if( angle > anglelm ) iangle = iangle + 1
   end do

   write(IOdef,*) 'Angles:',anglemn,anglemx,iangle
   write(IOdbg,*) 'Angles:',anglemn,anglemx,iangle
   write(IOrun,*) 'Angles:',anglemn,anglemx,iangle

   write(IOdef,'(1x,A)') 'Vertices are:'
   write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,1)),' <  x  < ',maxval(vert(:,1))
   write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,2)),' <  y  < ',maxval(vert(:,2))
   write(IOdef,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,3)),' <  z  < ',maxval(vert(:,3))

   write(IOdbg,'(1x,A)') 'Vertices are:'
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,1)),' <  x  < ',maxval(vert(:,1))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,2)),' <  y  < ',maxval(vert(:,2))
   write(IOdbg,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,3)),' <  z  < ',maxval(vert(:,3))

   write(IOrun,'(1x,A)') 'Vertices are:'
   write(IOrun,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,1)),' <  x  < ',maxval(vert(:,1))
   write(IOrun,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,2)),' <  y  < ',maxval(vert(:,2))
   write(IOrun,'(1x,A,1pe10.3,A,e10.3)') &
     'Dimension ',minval(vert(:,3)),' <  z  < ',maxval(vert(:,3))
   !
   ! final stuff
   !
   close(IOgeo)

   write(Iodef,'(1x,A)') 'End reading geometry'

   write(IOdbg,*)'End reading geometry data'
   write(IOdbg,*)'Current allocated memory (words): ',memory
   write(IOdbg,*)'Status:'
   write(IOdbg,'(1x,A,i8)') 'Number of cells      =',Ncel
   write(IOdbg,'(1x,A,i8)') '            Specials =',count( Nfaces > 6 )
   write(IOdbg,'(1x,A,i8)') '               Hexas =',count( Nfaces == 6 )
   write(IOdbg,'(1x,A,i8)') '    Wedges, Pyramids =',count( Nfaces == 5 )
   write(IOdbg,'(1x,A,i8)') '              Tetras =',count( Nfaces == 4 )
   write(IOdbg,'(1x,A,i8)') 'Number of vertices   =',Nvrt
   write(IOdbg,'(1x,A,i8)') 'Number of faces      =',Nfac
   Nint = count( Face(:)%cell2 /= 0)
   write(IOdbg,'(1x,A,i8)') '      Internal faces =',Nint
   if( (Nfac-Nint) /= Nbnd ) write(IOdbg,*)'Error: inconsistent geometry data'
   write(IOdbg,'(1x,A,i8)') 'Number of boundaries =',Nbnd
   do i=0,Nreg
     write(IOdbg,'(1x,A,i2,A,i8,A)') &
               '           Region ',i,' =',count( Bnd(:)%rid == i)
   end do

   write(IOrun,*)'End reading geometry data'
   write(IOrun,*)'Status:'
   write(IOrun,'(1x,A,i8)') 'Number of cells      =',Ncel
   write(IOrun,'(1x,A,i8)') '            Specials =',count( Nfaces > 6 )
   write(IOrun,'(1x,A,i8)') '               Hexas =',count( Nfaces == 6 )
   write(IOrun,'(1x,A,i8)') '    Wedges, Pyramids =',count( Nfaces == 5 )
   write(IOrun,'(1x,A,i8)') '              Tetras =',count( Nfaces == 4 )
   write(IOrun,'(1x,A,i8)') 'Number of vertices   =',Nvrt
   write(IOrun,'(1x,A,i8)') 'Number of faces      =',Nfac
   Nint = count( Face(:)%cell2 /= 0)
   write(IOrun,'(1x,A,i8)') '      Internal faces =',Nint
   if( (Nfac-Nint) /= Nbnd ) write(IOdbg,*)'Error: inconsistent geometry data'
   write(IOrun,'(1x,A,i8)') 'Number of boundaries =',Nbnd
   do i=0,Nreg
     Reg(i)%n = count( Bnd(:)%rid == i)
     write(IOrun,'(1x,A,i2,A,i8,A)') '           Region ',i,' =',Reg(i)%n
   end do

   !
   ! now everything is read in calculate extra geometric data
   ! (regions are known, faces etc as well so fill in some arrays)
   !
   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     Reg(ir)%area = Reg(ir)%area + Face(i)%area
   end do

   !
   ! these calculations appear relatively often and originaly
   ! were in the main loops of FluxMass and FluxMass2
   !
   do i=1,Nfac
     ip = face(i)%cell1
     in = face(i)%cell2
     if( in > 0 )then
       !
       ! internal cell face which points to two cells
       !
       facn = Face(i)%lambda
       facp = 1.0 - facn

       Xpn   = Cell(in)%x - Cell(ip)%x
        
       Xnorm = Face(i)%n
       call normalise(Xnorm)
       
       Xface = Face(i)%x
       !
       ! find auxiliary nodes  (eq. 8.53)
       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  P    e       e   P
       !
       Xpac = Xface - dot_product(Xface-Cell(ip)%x,Xnorm)*Xnorm

       Xnp  = Xpac - Xface

       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  E    e       e   E
       !
       Xnac = Xface - dot_product(Xface-Cell(in)%x,Xnorm)*Xnorm

       Xnn  = Xnac - Xface

       !
       ! use the shortest
       !
       dp   = vector_length(Xnp)
       dn   = vector_length(Xnn)

       if( dn >= dp )then
         Xpac2 = Xface + Xnp
         Xnac2 = Xface - Xnp
       else
         Xpac2 = Xface - Xnn
         Xnac2 = Xface + Xnn
       endif

     else
       !
       ! boundary face
       !   
       
     endif
     !
     ! store the auxiliary nodes in the face
     !
     Face(i)%Xpac = Xpac2
     Face(i)%Xnac = Xnac2
     
   end do

   !
   ! the diffusion term uses the over relaxed approach reeds, see e.g.
   ! www.cfd-online.com/Wiki/Discretization_of_the_diffusion_term
   !
   do i=1,Nfac
     ip = face(i)%cell1
     in = face(i)%cell2
     if( in > 0 )then
       Xpn = Cell(in)%x - Cell(ip)%x
     else
       Xpn = Face(i)%x - Cell(ip)%x     
     endif
     
     Face(i)%rlencos = Face(i)%area / vector_length(Xpn) / &
                                      vector_cosangle(Face(i)%n,Xpn)
   end do
      

   !
   ! list arrays of cell faces and nodes in very compact form. 
   ! example cfaces(i,j)-array becomes iface(i)-array. suppose cell 1
   ! has faces 1,2,3 and 4, and cell 2 has faces 3,4,5,6,7 and 8 and 
   ! finally cell 3 has 6,7,8 and 9. thus a compact list is of the 
   ! form: 1,2,3,4, 3,4,5,6,7,8,  6,7,8,9 +
   ! idx:  1        5            11       15
   ! nr:   5-1=4    11-5=6       15-11=4 
   !
   ! so the NFaces-array is actually superfluous and can be replaced by 
   ! subtracting index(ip+1)-index(ip). 
   !
   ! The Compressed Row Storage (CRS) type of format is for sparse matrices
   ! but there is no 'column' here. So one array can be saved. 
   ! Another point is that there are no real values involved only
   ! integers.
   ! 
   icnt = sum(NFaces(1:Ncel))
   
   allocate(iFaces%i(1:Ncel+1),stat=istat)
   call TrackMemory(istat,Ncel+1,'Cell face index array allocated')

   allocate(iFaces%list(1:icnt),stat=istat)
   call TrackMemory(istat,icnt,'Cell face list array allocated')

   iFaces%i    = -1
   iFaces%list = -1
   
   idx = 1
   do i=1,Ncel
     iFaces%i(i) = idx
     do j=1,NFaces(i)
       iFaces%list(idx) = CFace(i,j)
       idx              = idx + 1
     end do
   end do

   if( idx == icnt + 1 )then
     iFaces%i(Ncel+1) = idx
   else
     write(IOdef,*)'+++ internal error: Cell face array corrupted'
     stop
   endif 
 
   write(IOdbg,'(/,''    Cell    :    Faces ...'')')
   do i=1,min(10,Ncel)
     write(IOdbg,'(i8,1x,i2,'' : '',$)') i,NItemsFromList(iFaces,i)
     do j=1,NItemsFromList(iFaces,i)
       write(IOdbg,'(i8,1x,$)') ItemFromList(iFaces,i,j)
     end do
     write(IOdbg,'('' '')')
   end do
   if( Ncel > 10 )then
     write(IOdbg,'('' ...'')')
     do i=Ncel-10,Ncel
       write(IOdbg,'(i8,1x,i2,'' : '',$)') i,NItemsFromList(iFaces,i)
       do j=1,NItemsFromList(iFaces,i)
         write(IOdbg,'(i8,1x,$)') ItemFromList(iFaces,i,j)
       end do
       write(IOdbg,'('' '')')
     end do
   endif
   write(IOdbg,'('' '')')
   
   !
   ! in principle the CFace array can now be deleted 
   !
   !deallocate(CFace)
   !call TrackMemory(istat,icnt,...
   !
   ! same idea for the nodes of a cell 
   ! use a temporary array...
   !
   allocate(NNodes(Ncel),stat=istat)
   call TrackMemory(istat,Ncel,'NNodes array allocated')
  
   NNodes(:) = 0
  
   maxfac = maxval(NFaces(1:Ncel))
   maxnod = maxfac * MaxFaceNodes
      
   allocate(TmpCellNodes(maxnod),stat=istat)
   call TrackMemory(istat,maxnod,'Temporary TmpCellNodes array allocated')

   write(IOdbg,*)'Maximum faces of a cell: ',maxfac
   !
   ! see also subroutine CollectCellNodes (which can be deleted...)
   !
   do i=1,Ncel
     TmpCellNodes = 0
     Nnod = 0 
     do j=1,NItemsFromList(iFaces,i)
       k = ItemFromList(iFaces,i,j)
       do iv=1,MaxFaceNodes
         if( Face(k)%vertices(iv) < 0 ) cycle
         Found = .false.
         do m=1,Nnod
           if( TmpCellNodes(m) == Face(k)%vertices(iv) )then
             Found = .true.
             exit
           endif
         end do
         if( .not. Found )then
           Nnod = Nnod + 1
           TmpCellNodes(Nnod) = Face(k)%vertices(iv)
         endif
       end do
     end do
     NNodes(i) = Nnod
   end do

   write(IOdbg,*)'Maximum nodes of a cell: ',maxval(NNodes(1:Ncel))
   !
   ! now the compressed sizes are known
   !    
   icnt = sum(NNodes(1:Ncel))
   
   allocate(iNodes%i(1:Ncel+1),stat=istat)
   call TrackMemory(istat,Ncel+1,'Cell node index array allocated')

   allocate(iNodes%list(1:icnt),stat=istat)
   call TrackMemory(istat,icnt,'Cell node list array allocated')

   iNodes%i    = -1
   iNodes%list = -1

   idx = 1
   do i=1,Ncel
     iNodes%i(i)  = idx
     TmpCellNodes = 0

     Nnod = 0 
     do j=1,NItemsFromList(iFaces,i)
       k = ItemFromList(iFaces,i,j)
       do iv=1,MaxFaceNodes
         if( Face(k)%vertices(iv) < 0 ) cycle
         
         Found = .false.
         do m=1,Nnod
           if( TmpCellNodes(m) == Face(k)%vertices(iv) )then
             Found = .true.
             exit
           endif
         end do
         
         if( .not. Found )then
           Nnod = Nnod + 1
           TmpCellNodes(Nnod) = Face(k)%vertices(iv)
           iNodes%list(idx)   = Face(k)%vertices(iv)
           idx = idx + 1
         endif
       end do
     end do
   end do
   
   if( idx == icnt + 1 )then
     iNodes%i(Ncel+1) = idx
   else
     write(IOdef,*)'+++ internal error: Cell node array corrupted'
     stop
   endif 

   !
   ! todo: reshuffle the iNodes%list? faster access?
   !

   write(IOdbg,'(/,''    Cell    :    Nodes ...'')')
   do i=1,min(10,Ncel)
     write(IOdbg,'(i8,1x,i2,'' : '',$)') i,NItemsFromList(iNodes,i)
     do j=1,NItemsFromList(iNodes,i)
       write(IOdbg,'(i8,1x,$)') ItemFromList(iNodes,i,j)
     end do
     write(IOdbg,'('' '')')
   end do
   if( Ncel > 10 )then
     write(IOdbg,'('' ...'')')
     do i=Ncel-10,Ncel
       write(IOdbg,'(i8,1x,i2,'' : '',$)') i,NItemsFromList(iNodes,i)
       do j=1,NItemsFromList(iNodes,i)
         write(IOdbg,'(i8,1x,$)') ItemFromList(iNodes,i,j)
       end do
       write(IOdbg,'('' '')')
     end do
   endif
   write(IOdbg,'('' '')')
   
   deallocate(TmpCellNodes)
   call TrackMemory(0,-maxnod,'Temporary TmpCellNodes array deallocated')
     
   call watch_leave('ReadGeometry')
   
end subroutine ReadGeometry
integer function ItemFromList(list,i,j)

   use geometry
   
   type(StoredList), intent(IN) :: list
   integer, intent(IN)          :: i, j
   
   idx = list%i(i)
   ItemFromList = list%list(idx+j-1)

end function ItemFromList
integer function NItemsFromList(list,i)

   use geometry
   
   type(StoredList), intent(IN) :: list
   integer, intent(IN)          :: i
   
   idx1 = list%i(i)
   idx2 = list%i(i+1)
   NItemsFromList = idx2 - idx1

end function NItemsFromList
subroutine InitialiseVariables
!========================================================================

   use constants
   use geometry
   use variables
   use scalars
   use watches
   
  !call watch_enter('InitialiseVariables')
   
   Small = epsilon(U(1)) * 0.1
   Large = huge(large) * 1.e-3

   SMALL = 1.e-12

   if( Debug > -1 )then
     write(IOdbg,*)'Constants small/large:',Small,Large
     write(IOdbg,*)'Tiny     :',tiny(small)
     write(IOdbg,*)'Epsilon  :',epsilon(U(1))
     write(IOdbg,*)'Precision:',precision(small),digits(small), &
                            minexponent(small),maxexponent(small)
   endif

   i = 0
   allocate( U(Ncel+Nbnd),stat=istat)
   i = i + istat
   allocate( V(Ncel+Nbnd),stat=istat )
   i = i + istat
   allocate( W(Ncel+Nbnd),stat=istat )
   istat = i + istat
   call TrackMemory(istat,(Ncel+Nbnd)*3,&
                   'Velocity components array allocated')

   U  = Guess(VarU)
   V  = Guess(VarV)
   W  = Guess(VarW)

   i = 0
   allocate( P(Ncel+Nbnd),stat=istat)
   i = i + istat
   allocate( PP(Ncel+Nbnd),stat=istat)
   istat = i + istat
   call TrackMemory(istat,(Ncel+Nbnd)*2,'Pressure array allocated')

   P  = Guess(VarP)
   PP = 0.0

   i = 0
   allocate( dudx(Ncel+Nbnd,3),stat=istat )
   i = i + istat
   allocate( dvdx(Ncel+Nbnd,3),stat=istat )
   i = i + istat
   allocate( dwdx(Ncel+Nbnd,3),stat=istat )
   i = i + istat
   allocate( dpdx(Ncel+Nbnd,3),stat=istat )
   istat = i + istat
   call TrackMemory(istat,(Ncel+Nbnd)*12,'Gradient arrays allocated')

   dudx = 0.0
   dvdx = 0.0
   dwdx = 0.0

   dpdx = 0.0

   allocate( Den(Ncel+Nbnd),stat=istat)
   call TrackMemory(istat,Ncel+Nbnd,'Density array allocated')

   Den = DensRef                                                  !<=== let op material 1

   allocate( Cp(Ncel+Nbnd),stat=istat)
   call TrackMemory(istat,Ncel+Nbnd,'Cp array allocated')

   Cp  = CpStd                                                    !<=== let op material 1

   allocate( VisEff(Ncel+Nbnd),stat=istat)
   call TrackMemory(istat,Ncel+Nbnd,'Effective viscosity array allocated')

   VisEff  = VisLam                                               !<=== let op material 1

   if( SolveEnthalpy )then
     allocate( T(Ncel+Nbnd),stat=istat)
     call TrackMemory(istat,Ncel+Nbnd,'Temperature array allocated')

     T = Tref          !Tref wordt 0! we werken straks relatief en per materiaal!

   endif

   if( SolveTurb )then
     i = 0
     allocate( TE(Ncel+Nbnd),stat=istat )
     i = i + istat
     allocate( TurbP(Ncel+Nbnd),stat=istat )
     i = i + istat
     allocate( ED(Ncel+Nbnd),stat=istat )
     i = i + istat
     call TrackMemory(istat,(Ncel+Nbnd)*3,'K-eps model arrays allocated')

     TE    = Guess(VarTE)
     TurbP = 0.0

     ED    = Guess(VarED)

   endif

   if( SolveScalars )then

     write(IOdef,*)'Allocating ',Nscal,' scalar arrays'
     allocate( SC(Ncel+Nbnd,Nscal),stat=istat)
     call TrackMemory(istat,(Ncel+Nbnd)*Nscal,&
                     'Scalar array(s) allocated')

     SC = 0.0

   endif

   if( DXdebugdata )then
     allocate( DXdebug(Ncel),stat=istat)
     dxdebug = 0.0

     allocate( DXgrad(Ncel,3),stat=istat)
     call TrackMemory(istat,Ncel,'DX debug arrays allocated')
     dxgrad = 0.0
   endif

   allocate( dPhidXo(Ncel+Nbnd,3),stat=istat)  ! vraag? moet Nbnd???
   call TrackMemory(istat,Ncel*6,'Old gradient dPhi/dX array allocated')
   dPhidXo = 0.0

   allocate( MassFlux(Nfac),stat=istat)
   call TrackMemory(istat,Nfac,'Mass fluxes array allocated')
   MassFlux = 0.0

   !
   ! linear solver section
   !
   i = 0
   allocate( Ar(Ncel),stat=istat)
   i = i + istat
   allocate( Au(Ncel),stat=istat)
   i = i + istat
   allocate( Av(Ncel),stat=istat)
   i = i + istat
   allocate( Aw(Ncel),stat=istat)
   i = i + istat
   allocate( Su(Ncel),stat=istat)
   i = i + istat
   allocate( Sv(Ncel),stat=istat)
   i = i + istat
   allocate( Sw(Ncel),stat=istat)
   i = i + istat
   allocate( Res(Ncel),stat=istat)
   i = i + istat
   call TrackMemory(istat,Ncel*8,'Cell matrices allocated')

   Ar  = 0.0
   Au  = 0.0
   Su  = 0.0
   Av  = 0.0
   Sv  = 0.0
   Aw  = 0.0
   Sw  = 0.0
   Res = 0.0

   NNZ = Ncel + 2*Nint
   i = 0
   allocate( Acoo(NNZ),stat=istat)     ! double prec.
   i = i + istat
   allocate( Arow(NNZ),stat=istat)
   i = i + istat
   allocate( Acol(NNZ),stat=istat)
   i = i + istat
   allocate( Acsr(NNZ),stat=istat)     ! double prec.
   i = i + istat
   allocate( Arwc(Ncel+1),stat=istat)
   i = i + istat
   allocate( Aclc(NNZ),stat=istat)
   istat = i + istat

   allocate( RHS(Ncel),stat=istat)     ! righthand side in dble. prec.
   istat = i + istat
   allocate( SOL(Ncel),stat=istat)     ! solution in dble. prec.
   istat = i + istat

   call TrackMemory(istat,(NNZ)*10,'Matrix A in COO/CSR format allocated')
   if( debug > 2 ) write(IOdef,*)'NNZ: ',NNZ

   !
   ! Preconditioner arrays
   !
   ! allocate( ALU(NNZ),stat=istat)     ! double prec.
   ! istat = i + istat
   ! allocate( JLU(NNZ),stat=istat)
   ! istat = i + istat
   ! allocate( JU(Ncel),stat=istat)
   ! istat = i + istat
   ! allocate( JW(Ncel*2),stat=istat)
   ! istat = i + istat

   !
   ! store arrays for transient simulations
   !
   if( Transient )then
     i = 0
     allocate( Uold(Ncel+Nbnd),stat=istat)
     i = i + istat
     allocate( Vold(Ncel+Nbnd),stat=istat )
     i = i + istat
     allocate( Wold(Ncel+Nbnd),stat=istat )
     istat = i + istat
     call TrackMemory(istat,(Ncel+Nbnd)*3,&
                     'Old velocity components arrays allocated')

     Uold = Guess(VarU)
     Vold = Guess(VarV)
     Wold = Guess(VarW)

     if( QuadTime )then
       i = 0
       allocate( Uold2(Ncel+Nbnd),stat=istat)
       i = i + istat
       allocate( Vold2(Ncel+Nbnd),stat=istat )
       i = i + istat
       allocate( Wold2(Ncel+Nbnd),stat=istat )
       istat = i + istat
       call TrackMemory(istat,(Ncel+Nbnd)*3,&
                       'Old2 velocity components arrays allocated')

       Uold2 = Guess(VarU)
       Vold2 = Guess(VarV)
       Wold2 = Guess(VarW)
     endif

     if( SolveEnthalpy )then
       allocate( Told(Ncel+Nbnd),stat=istat)
       call TrackMemory(istat,(Ncel+Nbnd)*3,&
                       'Old enthalpy array allocated')
       Told = Guess(VarT)

       if( QuadTime )then
         allocate( Told2(Ncel+Nbnd),stat=istat)
         call TrackMemory(istat,(Ncel+Nbnd)*3,&
                       'Old2 enthalpy array allocated')
         Told2 = Guess(VarT)
       endif
     endif

     if( SolveTurb )then
       i = 0
       allocate( TEold(Ncel+Nbnd),stat=istat)
       i = i + istat
       allocate( EDold(Ncel+Nbnd),stat=istat )
       istat = i + istat
       call TrackMemory(istat,(Ncel+Nbnd)*3,&
                     'Old turbulence model arrays allocated')
       TEold = Guess(VarTE)
       EDold = Guess(VarED)

       if( QuadTime )then
         i = 0
         allocate( TEold2(Ncel+Nbnd),stat=istat)
         i = i + istat
         allocate( EDold2(Ncel+Nbnd),stat=istat )
         istat = i + istat
         call TrackMemory(istat,(Ncel+Nbnd)*3,&
                     'Old2 turbulence model arrays allocated')
         TEold2 = Guess(VarTE)
         EDold2 = Guess(VarED)
       endif

     endif

     if( SolveScalars )then

       allocate( SCold(Ncel+Nbnd,Nscal),stat=istat)
       call TrackMemory(istat,(Ncel+Nbnd)*Nscal,&
                       'Old scalar array(s) allocated')
       SCold = 0.0

       if( QuadTime )then
         allocate( SCold2(Ncel+Nbnd,Nscal),stat=istat)
         call TrackMemory(istat,(Ncel+Nbnd)*Nscal,&
                       'Old2 scalar array(s) allocated')
         SCold2 = 0.0
       endif
     endif

   endif

  !call watch_leave('InitialiseVariables')
   
end subroutine InitialiseVariables
subroutine SetBoundaryConditions
!========================================================================
!
!   this routine calculates the convective and diffusive fluxes
!   through the cell faces.
!

   use constants
   use geometry
   use variables
   use scalars

   real Xn(3)

   if( Debug > 3 ) write(IOdef,*)'*** SetBoundaryConditions'
   !
   ! first collect areas and inlet fluxes
   !
   areain  = 0.0
   fluxin  = 0.0
   areaout = 0.0
   icntin  = 0
   icntout = 0
   do ib=1,Nbnd

     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RInlet )then
       icntin = icntin + 1
       areain = areain + Face(i)%area
       if( Reg(ir)%User )then

         ip = Face(i)%cell1

         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den

         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area

         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)

         fluxin = fluxin +   Dtmp * ( Utmp*Face(i)%n(1)   &
                + Vtmp*Face(i)%n(2) + Wtmp*Face(i)%n(3) )

       else
         fluxin = fluxin + Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )
       endif
     else if( it == ROutlet )then
       icntout = icntout + 1
       areaout = areaout + Face(i)%area
     endif
   end do
   !
   ! check if valid
   !
   if( icntin > 0 .and. icntout == 0 )then
     write(IOdef,*) 'Invalid case: inflow without outflow. Stop'
     stop
   endif
   if( icntin == 0 .and. icntout > 0 )then
     write(IOdef,*) 'Invalid case: outflow without inflow. Stop'
     stop
   endif
   if( icntin == 0 .and. icntout == 0 )then
     write(IOdef,*) 'Case without in and outflow.'
     write(IOdbg,*) 'Case without in and outflow.'
     write(IOrun,*) 'Case without in and outflow.'
   endif

   if( icntin /= 0 )then
     rate = -1.0 * fluxin/(areain+Small)        ! mf < 0 bij inlet
     if( Debug > 1 )then
       write(IOdef,*)'Mass flow in:',fluxin,fluxin/areain
       write(IOdef,*)'     Area in:',areain
       write(IOdef,*)'    Area out:',areaout
     endif
     write(IOdbg,*)'Mass flow in:',fluxin,fluxin/areain
     write(IOdbg,*)'     Area in:',areain
     write(IOdbg,*)'    Area out:',areaout

     rate = areain / (areaout+Small) * rate     ! corr. voor opp. verschil
     if( Debug > 1 ) write(IOdef,*)' Rate in/out:',rate
     write(IOdbg,*)' Rate in/out:',rate
     if( rate <= 0.0 )then
       write(IOdef,*)''
       write(IOdef,*)'+++ Warning: Check inlet conditions'
       write(IOdef,*)''
       rate = Small
     endif
   else
     write(IOdef,*)'Mass flow in: absent'
     write(IOdbg,*)'Mass flow in: absent'
     write(IOrun,*)'Mass flow in: absent'
     rate = 0.0
   endif

   do ib=1,Nbnd

     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in /= 0 )write(IOdef,*)'+ internal error: corrupt bc ',ib,i,ip,in

     if( it == RInlet )then
       U(Ncel+ib)       = Reg(ir)%uvw(1)
       V(Ncel+ib)       = Reg(ir)%uvw(2)
       W(Ncel+ib)       = Reg(ir)%uvw(3)
       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Reg(ir)%den
       if( SolveEnthalpy) T(Ncel+ib) = Reg(ir)%T

       if( SolveTurbEnergy ) TE(Ncel+ib) = Reg(ir)%k
       if( SolveTurbDiss )   ED(Ncel+ib) = Reg(ir)%e
       if( SolveTurb )then
         VisEff(Ncel+ib) = VisLam + &
           Den(Ncel+ib)*TMCmu*TE(Ncel+ib)/(ED(Ncel+ib)+Small)
       else
         VisEff(Ncel+ib) = VisLam
       endif

       if( Restart == 0 )then
         U(ip)  = U(Ncel+ib)
         V(ip)  = V(Ncel+ib)                     ! initial guess
         W(ip)  = W(Ncel+ib)
       endif
       !
       ! extra
       !
       MassFlux(i) = Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )
       !
       if( Reg(ir)%user )then

         ip = Face(i)%cell1

         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den

         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area

         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)

         MassFlux(i) = Dtmp * ( Utmp*Face(i)%n(1)   &
                              + Vtmp*Face(i)%n(2)   &
                              + Wtmp*Face(i)%n(3) )

         U(Ncel+ib)       = Utmp
         V(Ncel+ib)       = Vtmp
         W(Ncel+ib)       = Wtmp
         P(Ncel+ib)       = P(ip)
         if( SolveTurbEnergy ) TE(Ncel+ib) = TEtmp
         if( SolveTurbDiss )   ED(Ncel+ib) = EDtmp
         if( SolveEnthalpy)    T(Ncel+ib) = Ttmp
         Den(Ncel+ib)     = Dtmp

         if( SolveTurb )then
           VisEff(Ncel+ib) = VisLam + Dtmp*TMCmu*TEtmp/(EDtmp+Small)
         else
           VisEff(Ncel+ib) = VisLam
         endif

         if( Restart == 0 )then
           U(ip)  = U(Ncel+ib)
           V(ip)  = V(Ncel+ib)                  ! initial guess
           W(ip)  = W(Ncel+ib)
         endif

       endif

     else if( it == ROutlet .and. (Restart == 0) )then

       Xn     = Face(i)%n                      ! de surface vector
       call normalise(Xn)                      ! normaal vector met lengte 1

       U(Ncel+ib) = Xn(1) * rate / DensRef !
       V(Ncel+ib) = Xn(2) * rate / DensRef !
       W(Ncel+ib) = Xn(3) * rate / DensRef !

       U(ip)  = U(Ncel+ib)
       V(ip)  = V(Ncel+ib)                     ! initial guess
       W(ip)  = W(Ncel+ib)

       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Den(ip)
       if( SolveEnthalpy) T(Ncel+ib) = T(ip)
       VisEff(Ncel+ib)  = VisEff(ip)

     else if( it == RSymp .and. (Restart == 0 ) )then
       U(Ncel+ib)       = U(ip)
       V(Ncel+ib)       = V(ip)
       W(Ncel+ib)       = W(ip)
       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Den(ip)
       if( SolveEnthalpy) T(Ncel+ib) = T(ip)
       VisEff(Ncel+ib)  = VisEff(ip)
     else if( it == RWall .and. (Restart == 0 ) )then
       U(Ncel+ib)       = Reg(ir)%uvw(1)
       V(Ncel+ib)       = Reg(ir)%uvw(2)
       W(Ncel+ib)       = Reg(ir)%uvw(3)
       P(Ncel+ib)       = P(ip)
       Den(Ncel+ib)     = Den(ip)
       if( SolveEnthalpy) T(Ncel+ib) = Reg(ir)%T
       VisEff(Ncel+ib)  = VisEff(ip)
     else if( Restart == 0 )then
       write(IOdef,*)'SetBoundaryConditions, unknown Boundary Condition'
       write(IOdef,*)'i,ir,it: ',i,ir,it
     endif

   end do

   if( Debug > 3 ) write(IOdef,*)'=== SetBoundaryConditions'

end subroutine SetBoundaryConditions
subroutine CalculateScalar(ivar,Phi,dPhidX,PhiOld,PhiOld2)
!========================================================================

   use constants
   use geometry
   use variables
   use scalars
   use watches
   
   !use hypre_mod, only: solve_type                        ! === Hypre ===

   integer, intent(in)          :: iVar      ! variable/scalar id

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd)   :: PhiOld
   real, dimension(Ncel+Nbnd)   :: PhiOld2
   real, dimension(Ncel+Nbnd,3) :: dPhidX

   logical :: Warning = .false.

   real, dimension(3)           :: ds

   call watch_enter('CalculateScalar')
   
itmpdebug = Debug
!Debug = 0

   if( Debug > 3 ) write(IOdef,*)'*** CalculateScalar: ',Variable(ivar)
   !
   ! calculate gradients of phi
   !
!debug=8
   call GradientPhi(iVar,Phi,dPhidX)
!debug=0

   Au  = 0.0  !
   Su  = 0.0  ! use the arrays for the u-velocity component
   Acoo= 0.0  !

   !
   ! calculate fluxes through all inner faces
   !
   call FluxScalar(ivar,Phi,dPhidX)

   !
   ! unsteady term
   !
   if( Transient )then
     rdt = 1.0 / dt
     if( Debug > 2 ) write(IOdef,*)'Instat. Scalar term, rdt =',rdt
     if( Euler )then
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt
         Su(i) = Su(i) + f * PhiOld(i)
         Au(i) = Au(i) + f
       end do
     else
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt
         Su(i) = Su(i) + f *( (1.0 + GammaTime )*PhiOld(i) - 0.5*GammaTime*PhiOld2(i) )
         Au(i) = Au(i) + f *( 1.0 + 0.5*GammaTime )
       end do
     endif
   endif

   !
   ! de term P*div(Vol) ?
   !

   !
   ! moleculaire dissipatie
   !
   if( SolveEnthalpy .and. iVar == VarT )then
     !write(IOdef,*)'Molecular dissipation term'
     !do i=1,Ncel
     !  Dissip = GenP * Vism * vol
     !  Su(i) = Su(i) + Dissip
     !end do
     !QMdis(mat) = QMdis(mat) + Dissip(1:Ncel)
   endif

   !
   ! source and wall terms for turbulence models
   !
   if( iVar == VarTE .or. ivar == VarED )then

     call TurbulenceModels(iVar,Phi,dPhidX)

   endif

   if( UsePatches ) call PatchesSource(iVar, Phi, Au, Su)    !<====!!!!!!

   !
   ! pick underrelaxation factor
   !
   if( iVar > NVar )then
     URFScalar  = URF(VarSC)
   else
     URFScalar  = URF(iVar)
   endif
   
   call SetUpMatrixA( iVar,URFScalar, Au,Phi,Su  )

   !solve_type = BiCGstab2                                 ! === Hypre ===

   call SolveMatrixA( iVar, ResS0,ResS1, Au,Phi,Su )

   Residual(iVar)    = ResS1

   write(IOdbg,*)'S:',minval(phi(:)),'< ',Variable(ivar),' <',maxval(phi(:))
   write(IOdbg,*)'Residu:',ResS0,' -> ',ResS1
   !
   ! final feasibility check!
   !
   if( iVar == VarT .and. LimitLowEnthalpy )then
     ic = 0
     write(IOdbg,*)'Laagste waarde T:',minval(phi(1:Ncel))
     do ip=1,Ncel
       if( phi(ip) < LowerLimitT )then
         ic = ic + 1
         phi(ip) = LowerLimitT
       endif
     end do
     write(IOdbg,*)'Ingrepen low in T:',ic,LowerLimitT
   endif
   if( iVar == VarT .and. LimitUpEnthalpy )then
     ic = 0
     write(IOdbg,*)'hoogste waarde T:',maxval(phi(1:Ncel))
     do ip=1,Ncel
       if( phi(ip) > UpperLimitT )then
         ic = ic + 1
         phi(ip) = UpperLimitT
       endif
     end do
     write(IOdbg,*)'Ingrepen upper in T:',ic,UpperLimitT
   endif
   if( iVar == NVar+1 )then
     if( LimitLowScalar(iVar-Nvar) )then
       ic = 0
       write(IOdbg,*)'Laagste waarde Sc:',iVar-Nvar,minval(phi(1:Ncel))
       do ip=1,Ncel
         if( phi(ip) < LowerLimitSC(iVar-Nvar) )then
           ic = ic + 1
           phi(ip) = LowerLimitSC(iVar-Nvar)
         endif
       end do
       write(IOdbg,*)'Ingrepen low in Scalar :',iVar-Nvar,':',ic, &
                      LowerLimitSC(iVar-Nvar)
     endif
     if( LimitUpScalar(iVar-Nvar) )then
       ic = 0
       write(IOdbg,*)'Hoogste waarde Sc:',iVar-Nvar,maxval(phi(1:Ncel))
       do ip=1,Ncel
         if( phi(ip) > UpperLimitSC(iVar-Nvar) )then
           ic = ic + 1
           phi(ip) = UpperLimitSC(iVar-Nvar)
         endif
       end do
       write(IOdbg,*)'Ingrepen upper in Scalar :',iVar-Nvar,':',ic, &
                     UpperLimitSC(iVar-Nvar)
     endif
   endif



   if( iVar == VarTE .or. ivar == VarED )then

     Warning = .false.
     icnt    = 0

     do ip=1,Ncel
       if( Phi(ip) <= 0.0 )then
         icnt = icnt + 1
         if( icnt <= 10 ) write(IOdbg,*)'+ neg. value ',Phi(ip), &
                                       ' for cell ',ip,Variable(iVar)
         Warning = .true.

         PhiAv = 0.0
         icells = 0
         do j=1,NFaces(ip)
           k   = CFace(ip,j)
           ipp = Face(k)%cell1
           inn = Face(k)%cell2
           if( inn > 0 )then
             ! internal
             icells = icells + 1
             if( ipp == ip )then
               PhiAv = PhiAv + Phi(inn)
             elseif( inn == ip )then
               PhiAv = PhiAv + Phi(ipp)
             endif
           endif
         end do
         if( icells > 0 )then
           PhiAv = PhiAv / float(icells )
           !Phi(ip) = max(PhiAv,Small)            ! must be positive
           if( iVar == VarTE )then
             Phi(ip) = max(PhiAv,1.e-12)           ! Small is too small!
           elseif( iVar == VarED )then
             Phi(ip) = max(PhiAv,1.e-12)
           else
             Phi(ip) = max(PhiAv,Small)
           endif
           if( icnt <= 10 ) write(IOdbg,*) '+ replaced by:',Phi(ip)
         else
           Phi(ip) = Small
         endif
       endif
     end do

     if( Warning )then
       write(IOdbg,*) '+ neg. values found for ',icnt,' cells'
       Flags(IFlagTE) = 't'
     endif

     if( iVar == VarTE .and. LimitLowTurbEnergy )then
       ic = 0
       write(IOdbg,*)'laagste waarde TE:',minval(phi(1:Ncel))
       do ip=1,Ncel
         if( phi(ip) < LowerLimitTE )then
           ic = ic + 1
           phi(ip) = LowerLimitTE
         endif
       end do
       write(IOdbg,*)'Ingrepen low in TE:',ic,LowerLimitTE
     endif
     if( iVar == VarTE .and. LimitUpTurbEnergy )then
       ic = 0
       write(IOdbg,*)'hoogste waarde TE:',maxval(phi(1:Ncel))
       do ip=1,Ncel
         if( phi(ip) > UpperLimitTE )then
           ic = ic + 1
           phi(ip) = UpperLimitTE
         endif
       end do
       write(IOdbg,*)'Ingrepen upper in TE:',ic,UpperLimitTE
     endif
     if( iVar == VarED .and. LimitLowTurbDiss )then
       ic = 0
       write(IOdbg,*)'laagste waarde ED:',minval(phi(1:Ncel))
       do ip=1,Ncel
         if( phi(ip) < LowerLimitED )then
           ic = ic + 1
           phi(ip) = LowerLimitED
         endif
       end do
       write(IOdbg,*)'Ingrepen low in ED:',ic,LowerLimitED
     endif
     if( iVar == VarED .and. LimitUpTurbDiss )then
       ic = 0
       write(IOdbg,*)'hoogste waarde TE:',maxval(phi(1:Ncel))
       do ip=1,Ncel
         if( phi(ip) > UpperLimitED )then
           ic = ic + 1
           phi(ip) = UpperLimitED
         endif
       end do
       write(IOdbg,*)'Ingrepen upper in ED:',ic,UpperLimitED
     endif

   endif

   !
   ! update sym. boundaries
   !
   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RSymp .or. it == ROutlet)then

       ds = Face(i)%x - Cell(ip)%x

       Phi(Ncel+ib) = Phi(ip) + dot_product( dPhidX(ip,:) , ds )

     endif
   end do
   !call Update_Scalars_at_boundaries2(ivar,Phi)

   if( Debug > 3 ) write(IOdef,*)'=== CalculateScalar'

 Debug = itmpdebug
 
   call watch_leave('CalculateScalar')
   
end subroutine CalculateScalar
subroutine FluxScalar(iVar,Phi,dPhidX)
!========================================================================
!
!   this routine calculates the convective and diffusive scalar fluxes
!   through the cell faces.
!
   use constants
   use geometry
   use variables
   use scalars
   use watches
   
   integer, intent(in)          :: iVar      ! variable/scalar id

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX

   real Xac(3), Xpn(3), dPhidXac(3)
   real d(3), d1(3), d2(3), d3(3), ds(3)
   real Xf(3), Xp(3), Xn(3), Xnorm(3), Xpac(3), Xnac(3), delp(3), deln(3)
 
   if( Debug > 3 ) write(IOdef,*)'*** FluxScalar: ',Variable(ivar),iVar

   call watch_enter('FluxScalar')
   
   pe0 =  Large
   pe1 = -Large

   GammaBlend = Gamma(iVar)
   IScheme    = Scheme(ivar)

   QtransferIn  =   0.0
   QtransferOut =   0.0
   Atot         =   0.0
   qmin         =  Large
   qmax         = -Large
   qtot         =   0.0
   hmin         =  Large
   hmax         = -Large
   htot         =   0.0

   ! reset counters
   call DiffSchemeCounter(0,0,0,iVar)
   
icntface=0

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2

     if( in > 0 )then
       !
       ! internal cell face which points to two cells
       ! see eq. 4.13 and 4.14, note lambda is from Xf to Xp
       !
       facn = Face(i)%lambda
       facp = 1.0 - facn
       !
       ! Viscosity Uncorrected! to be done?
       !
       Xac      =   Cell(in)%x * facn +   Cell(ip)%x * facp
       Phiac    =      Phi(in) * facn +      Phi(ip) * facp != central diff. scheme!
       Visac    =   VisEff(in) * facn +   VisEff(ip) * facp !<= moet alg. diff. coef worden

       !
       ! dif. coef. is dependent on the kind of simulation
       ! (turb. on or off) and the variable
       !
       ! recall: VisEff = VisLam + VisTurb
       !
       if( SolveTurb )then
         !
         ! turbulent diffusion
         !
         Visac = Visac - VisLam

         if( ivar == VarT )then
           Visac = ( VisLam + Visac / Sigma_T )/Prandtl
         else if( ivar == VarTE )then
           Visac = VisLam + Visac / Sigma_k
         else if( ivar == VarED )then
           !
           ! of is het alleen: Visac / Sigma_e ????
           ! nee want anders kloppen de extreme diss. sim. niet
           !
           Visac = VisLam + Visac / Sigma_e
         else
           Visac = ( VisLam + Visac / Sigma_s )/Schmidt
         endif
       else
         !
         ! only laminar scalars (incl. T)
         !
         ! note: thermal dif. is defined by: lambda/Cp
         !       using the definition of the Prandtl-number
         !       Pr = (vislam Cp)/lambda one can use
         !       vislam/Pr instead.
         !
         if( ivar == VarT )then
           Visac  = Visac / Prandtl
         else
           Visac  = Visac / Schmidt
         endif
       endif

       dPhidXac = dPhidX(in,:) * facn + dPhidX(ip,:) * facp
       delta    = dot_product( dPhidXac , Face(i)%x - Xac )  ! correction

       Xpn      = Cell(in)%x - Cell(ip)%x
      !VisFace  = Visac * Face(i)%area/vector_length(Xpn)
       VisFace  = Visac * Face(i)%rlencos

       call SelectDiffSchemeScalar(i,iScheme,ip,in, &
                                   Phi,dPhidX,PhiFace)

       !
       ! explicit higher order fluxes
       !
       fce = MassFlux(i) * PhiFace

       fde1= Visac * dot_product( dPhidXac , Face(i)%n )    ! eq. 8.17 

       !
       ! Jasak's over-relaxed alternative
       !
       d1  = dot_product( Xpn , Face(i)%n )
       s2  = Face(i)%area * Face(i)%area 
       
       d2  = Xpn * s2/d1
       d3  = Face(i)%n - d2
        
       fde2= Visac * dot_product( d3 , dPhidXac )

       !
       ! implicit lower order (simple upwind)
       ! convective and diffusive fluxes
       !
       fci = min( MassFlux(i) , 0.0 ) * Phi(in) &
           + max( MassFlux(i) , 0.0 ) * Phi(ip)

       fdi = VisFace * dot_product( dPhidXac , Xpn )

       !
       ! convective coefficients with deferred correction with
       ! gamma as the blending factor (0.0 <= gamma <= 1.0)
       !
       !      low            high    low  OLD
       ! F = F    + gamma ( F     - F    )
       !     ----   -------------------------
       !      |                  |
       !  implicit           explicit (dump into source term)
       !
       Rface(i,1) = -VisFace - max( MassFlux(i) , 0.0 )  ! Ap, Ae
       Rface(i,2) = -VisFace + min( MassFlux(i) , 0.0 )  ! An, Aw

       blend = GammaBlend*( fce - fci )

       !
       ! assemble the two source terms
       !
       Su(ip) = Su(ip) - blend + fde1 - fdi    
       Su(in) = Su(in) + blend - fde1 + fdi

      !Su(ip) = Su(ip) - blend + fde2  
      !Su(in) = Su(in) + blend - fde2  

       !
       ! Pe_f = (rho U)_f dXpn / Vis_f = (A rho U)_f dXpn /(A_f Vis_f)
       !      = MassFlux dXpn /(A_f Vis_f)
       !
       peclet = MassFlux(i)/Face(i)%area*vector_length(Xpn)/(Visac+Small)

       pe0 = min( pe0 , peclet )
       pe1 = max( pe1 , peclet )

     else
       !
       ! boundary face
       !
       ib = Face(i)%bnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       ip = Face(i)%cell1

       if( it == RInlet )then

         dPhidXac = dPhidX(ip,:)
         Xac      = Face(i)%x

         if( Reg(ir)%user )then
           Utmp  = Reg(ir)%uvw(1)
           Vtmp  = Reg(ir)%uvw(2)
           Wtmp  = Reg(ir)%uvw(3)
           TEtmp = Reg(ir)%k
           EDtmp = Reg(ir)%e
           Ttmp  = Reg(ir)%T
           Dtmp  = Reg(ir)%den

           Xtmp  = Face(i)%x(1)
           Ytmp  = Face(i)%x(2)
           Ztmp  = Face(i)%x(3)
           Atmp  = Face(i)%area

           call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                          Pref,Tref,Iter,Time,                  &
                          Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                          Reg(ir)%name1,Reg(ir)%name2)

           if( ivar == VarT )then
             PhiFace  = Ttmp
           else if( ivar == VarTE )then
             PhiFace  = TEtmp
           else if( ivar == VarED )then
             PhiFace  = EDtmp
           !else
           !  PhiFace  = ScReg(ir,(iVar-Nvar))%value
           endif
         else
           if( ivar == VarT )then
             PhiFace  = Reg(ir)%T
           else if( ivar == VarTE )then
             PhiFace  = Reg(ir)%k
           else if( ivar == VarED )then
             PhiFace  = Reg(ir)%e
           else if( ivar > NVar )then
             PhiFace  = ScReg(ir,(iVar-Nvar))%value
           else
             write(IOdef,*)'+++ internal error FluxScalar: inlet for ',&
                       Variable(ivar),ivar
           endif
         endif

         Visac    = VisEff(Ncel+ib)

         if( SolveTurb )then
           !
           ! turbulent diffusion
           !
           Visac = Visac - VisLam

           if( ivar == VarT )then
             Visac = ( VisLam + Visac / Sigma_T )/Prandtl
           else if( ivar == VarTE )then
             Visac = VisLam + Visac / Sigma_k
           else if( ivar == VarED )then
             Visac = VisLam + Visac / Sigma_e
           else
             Visac = ( VisLam + Visac / Sigma_s )/Schmidt
           endif
         else
           !
           ! only laminar scalars (incl. T)
           !
           if( ivar == VarT )then
             Visac  = Visac / Prandtl
           else
             Visac  = Visac / Schmidt
           endif
         endif

         Xpn      = Xac - Cell(ip)%x
        !VisFace  = Visac * Face(i)%area/vector_length(Xpn)
         VisFace  = Visac * Face(i)%rlencos
         !
	 ! explicit part
	 !
	 fde = Visac * dot_product( dPhidXac , Face(i)%n )

         fce = MassFlux(i) * PhiFace
         fde = fde

         !
	 ! implicit part
	 !
         fci = min( MassFlux(i) , 0.0 ) * PhiFace &
             + max( MassFlux(i) , 0.0 ) * Phi(ip)   ! voor het geval van neg. inlet

         fdi = VisFace * dot_product( dPhidXac , Xpn )

         f   = -VisFace + min( MassFlux(i) , 0.0 )

         Au(ip) = Au(ip) - f
         Su(ip) = Su(ip) + fde - fdi - f*PhiFace
         Phi(Ncel+ib) = PhiFace

       else if( it == ROutlet )then

         dPhidXac = dPhidX(ip,:)

         Xac      = Face(i)%x
         Visac    = VisEff(ip)
         if( SolveTurb )then
           !
           ! turbulent diffusion
           !
           Visac = Visac - VisLam

           if( ivar == VarT )then
             Visac = ( VisLam + Visac / Sigma_T )/Prandtl
           else if( ivar == VarTE )then
             Visac = VisLam + Visac / Sigma_k
           else if( ivar == VarED )then
             Visac = VisLam + Visac / Sigma_e
           else
             Visac = ( VisLam + Visac / Sigma_s )/Schmidt
           endif
         else
           !
           ! only laminar scalars (incl. T)
           !
           if( ivar == VarT )then
             Visac  = Visac / Prandtl
           else
             Visac  = Visac / Schmidt
           endif
         endif

         Xpn      = Xac - Cell(ip)%x
         PhiFace  = Phi(ip) + dot_product( dPhidX(ip,:) , Xpn )
        !VisFace  = Visac * Face(i)%area/vector_length(Xpn)
         VisFace  = Visac * Face(i)%rlencos

         fce = MassFlux(i) * PhiFace
         fde = Visac * dot_product( dPhidXac , Face(i)%n )

         fci = MassFlux(i) * Phi(ip)                    ! per. def naar buiten
         fdi = VisFace * dot_product( dPhidXac , Xpn )

         Su(ip) = Su(ip)  + fde - fdi
         Phi(Ncel+ib) = PhiFace

       else if( it == RSymp )then
         !
         ! do nothing just set PhiFace equal to Phi(ip)
         ! (te simpel maar ok nu)
         !
         ds = Face(i)%x - Cell(ip)%x
         Phi(Ncel+ib) = Phi(ip) + dot_product( dPhidX(ip,:) , ds )

       else if( it == RWall )then
         !write(IOdef,*)'Scalar wand: ',SolveEnthalpy,Reg(ir)%adiab,Reg(ir)%flux
         if( SolveEnthalpy .and. ivar == VarT )then
           if( Reg(ir)%adiab )then
             !
             ! adiabatic wall, kan beter nu even simpel
             !
             Phi(Ncel+ib) = Phi(ip)
             Bnd(ib)%q    = 0.0                     ! safety first
           else
             !
             ! wall temp or heat flux
             !
             if( Reg(ir)%flux )then
               PhiFlux      = Reg(ir)%T
             else
               PhiFace      = Reg(ir)%T
               Phi(Ncel+ib) = PhiFace
             endif

             Visac  = VisEff(Ncel+ib)               ! = VisLam + VisTurb
             dn     = Bnd(ib)%distance              ! dist. wall center to P
             Resist = Reg(ir)%R                     ! extra thermal resistance

             if( .not. SolveTurb )then
               !
               ! laminar case; note lambda/Cp = VisLam/Pr
               !
               VisFace = VisLam / Prandtl / dn

               Hcoef   = 1.0/( 1.0/VisFace + Resist*Cp(Ncel+ib) ) * Face(i)%area
	     else
               !
               ! turbulent case
               !
               ! the standard Jayatilleke extra sublayer resistance factor
               !
               SLres = 9.24*( (Prandtl/Sigma_T)**(0.75) - 1.0 ) *  &
                      (1.0 + 0.28*exp(-0.007*Prandtl/Sigma_T) )

               Cmu25 = TMCmu**0.25

               Tplus = Sigma_T * ( Bnd(ib)%uplus + SLres )
               utau  = Cmu25 * sqrt( TE(ip) )

               if( Bnd(ib)%yplus < Reg(ir)%ylog  )then
                 !
                 ! in the viscous sublayer
                 !
                 VisFace = VisLam / Prandtl / dn

                 Hcoef   = 1.0/( 1.0/VisFace + Resist*Cp(Ncel+ib) ) * Face(i)%area
               else
                 VisFace = Den(ip)*Utau/(Tplus + Small)

                 Hcoef   = 1.0/( 1.0/VisFace + Resist*Cp(Ncel+ib) ) * Face(i)%area
               endif
             endif

             if( Reg(ir)%flux )then
               PhiFace = Phi(ip) + PhiFlux / (Hcoef*Cp(Ncel+ib)/Face(i)%area)
               Phi(Ncel+ib) = PhiFace
             endif

!    write(*,*)'>>',ib,':',visface,hcoef,phiface
    
             Au(ip) = Au(ip) + Hcoef
             Su(ip) = Su(ip) + Hcoef * PhiFace

             Tdif      = ( PhiFace - Phi(ip) )     ! temperature difference
             Bnd(ib)%h = Hcoef * Cp(Ncel+ib) / &   ! take cp out of Hcoef
                                 Face(i)%area   ! DENK AAN rlencos !!!!
             Bnd(ib)%q = bnd(ib)%h  * Tdif         ! flux q in W/m2 (incl. R!)
             Bnd(ib)%h = VisFace * Cp(Ncel+ib)     ! take out cp and resistance
             Bnd(ib)%T = PhiFace                   ! just record the wall temp.

             if( Bnd(ib)%q > 0.0 )then
               QtransferIn  = QtransferIn  + Bnd(ib)%q * Face(i)%area
             else
               QtransferOut = QtransferOut + Bnd(ib)%q * Face(i)%area
             endif

             qmin = min(qmin,Bnd(ib)%q)
             qmax = max(qmax,Bnd(ib)%q)
             hmin = min(hmin,Bnd(ib)%h)
             hmax = max(hmax,Bnd(ib)%h)

             Atot = Atot + Face(i)%area
             qtot = qtot + Bnd(ib)%q * Face(i)%area
             htot = htot + Bnd(ib)%h * Face(i)%area

           endif

         elseif( SolveTurb )then
           !
           ! doe even niets
           !
         elseif( SolveScalars .and. iVar > Nvar )then
           !
           ! 'adiabatic' scalar only now
           !
           Phi(Ncel+ib) = Phi(ip)

         else

           write(IOdef,*)'+ internal error FluxScalar (KAN NOG NIET)'
           write(IOdef,*)'  Variable:',ivar,Variable(ivar)

         endif
       endif

     endif

   end do faceloop

   Qtransfer = QtransferIn + QtransferOut

   ! counters to IOdbg
   call DiffSchemeCounter(2,0,0,iVar)

   if( Debug > 0 .and. ivar == VarT ) &
     write(IOdef,*)'Total heat transferred:',Qtransfer,' W'

   if( ivar == VarT )then
     write(IOdbg,*)'Heat transferred IN    :',QtransferIN,' W'
     write(IOdbg,*)'Heat transferred OUT   :',QtransferOUT,' W'
     write(IOdbg,*)'Total heat transferred :',Qtransfer,' W'
     if( atot > 0.0 )then
       write(IOdbg,*) hmin,'<  h  <',hmax,' '
       write(IOdbg,*) qmin,'< q_w <',qmax,' W/m2'
       write(IOdbg,*)'Total heat area        :',atot,' m2'
       write(IOdbg,*)'Average heat flux      :',qtot/atot,' W/m2'
       write(IOdbg,*)'Average heat tran.coef.:',htot/atot,' W/m2K'
     endif
   endif

   !if( Debug > 0 )then
   !    write(IOdef,*) 'Scl: ',pe0,' < Peclet < ',pe1,'(',Variable(iVAR),')'
       write(IOdbg,*) 'Scl: ',pe0,' < Peclet < ',pe1,'(',Variable(iVAR),')'
   !endif

   call watch_leave('FluxScalar')
   if( Debug > 3 ) write(IOdef,*)'=== FluxScalar'

end subroutine FluxScalar
subroutine CalculateUVW
!========================================================================

   use constants
   use geometry
   use variables
   use watches

   !use hypre_mod, only: solve_type                        ! === Hypre ===

   logical :: copied = .false.   ! is the Ar array copied?

   if( Debug > 3 ) write(IOdef,*)'*** CalculateUVW'

   call watch_enter('CalculateUVW')

   copied = .false.

   Ar  = 0.0      ! reciprocal of one of the A-matrices

   Au  = 0.0      ! A-matrix for U-component
   Su  = 0.0      ! source term

   Av  = 0.0      ! V-component
   Sv  = 0.0

   Aw  = 0.0      ! W-component
   Sw  = 0.0

   ResU0 = 0.0
   ResU1 = 0.0
   ResV0 = 0.0
   ResV1 = 0.0
   ResW0 = 0.0
   ResW1 = 0.0

   !
   ! calculate fluxes through all inner faces
   ! (including processing of bc's)
   !
   call FluxUVW

   !
   ! bouyancy
   ! (later ook = (Den(i)-DenRef(imat))*Cell(i)%vol)
   ! algemeen: bodyforce(1-3,i)*Cell(i)%vol
   !
   if( SolveEnthalpy )then
     if( Debug > 2 )write(IOdef,*)'BodyForce',Gravity
     do i=1,Ncel
       BodyForce  = -Beta*Den(i)*Cell(i)%vol*( T(i) - Tref )

       Su(i) = Su(i) + Gravity(1) * BodyForce
       Sv(i) = Sv(i) + Gravity(2) * BodyForce
       Sw(i) = Sw(i) + Gravity(3) * BodyForce
     end do
   endif

   if( UsePatches )then
    call PatchesSource(VarU, U, Au, Su)
    call PatchesSource(VarV, V, Av, Sv)                      !<====!!!!!!
    call PatchesSource(VarW, W, Aw, Sw)
   endif
   !
   ! pressure force
   !
   if( Debug > 2 )write(IOdef,*)'PressureForce'
   do i=1,Ncel
     Su(i) = Su(i) - dPdX(i,1)*Cell(i)%vol
     Sv(i) = Sv(i) - dPdX(i,2)*Cell(i)%vol
     Sw(i) = Sw(i) - dPdX(i,3)*Cell(i)%vol
   end do
   !
   ! unsteady term
   !
   if( Transient )then
     rdt = 1.0 / dt
     if( Debug > 2 ) write(IOdef,*)'CalculateUVW: Instat. term, rdt =',rdt
     if( Euler )then
       !
       ! see section 6.3.2 Implicit Euler Method
       !
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt

         Su(i) = Su(i) + f * Uold(i)
         Sv(i) = Sv(i) + f * Vold(i)   ! op een later tijdstip
         Sw(i) = Sw(i) + f * Wold(i)   ! ook 2 oude tijden

         Au(i) = Au(i) + f
         Av(i) = Av(i) + f
         Aw(i) = Aw(i) + f
       end do
     else
       !
       ! see chapter 6.3.2 the Three Time Level Method
       !
       do i=1,Ncel
         f = Den(i) * Cell(i)%vol * rdt

         Su(i) = Su(i) + f *( (1.0 + GammaTime )*Uold(i) - 0.5*GammaTime*Uold2(i) )
         Sv(i) = Sv(i) + f *( (1.0 + GammaTime )*Vold(i) - 0.5*GammaTime*Vold2(i) )
         Sw(i) = Sw(i) + f *( (1.0 + GammaTime )*Wold(i) - 0.5*GammaTime*Wold2(i) )

         Au(i) = Au(i) + f *( 1.0 + 0.5*GammaTime )
         Av(i) = Av(i) + f *( 1.0 + 0.5*GammaTime )
         Aw(i) = Aw(i) + f *( 1.0 + 0.5*GammaTime )
       end do
     endif
   endif

   ! ====================
   ! Velocity component U
   ! ====================
   if( SolveU )then

     !
     ! set under relaxation factor
     !
     URFactor = URF(VarU)

     !
     ! According to Ferziger/Peric chap. 8.8
     ! "since A_P^{u_i} is the same for all vel. comp. in a
     ! given CV (exept near some boundaries), one can replace
     ! A_P^{v_n} by this quantity."
     !
     ! see also 8.10.3 (wall /symm. treatment)

     call SetUpMatrixA( VarU,URFactor, Au,U,Su  )

     !
     ! later anders ivm solid/fluid en poreuze media  ?!?!
     !
     if( .not. copied )then
       where( Au /= 0.0 )
         Ar = 1.0 / Au
       elsewhere
         Ar = 0.0
       end where
       copied = .true.
       if( Debug > 3 ) write(IOdef,*)'Calculate UVW (U): Ar based on Au'
     endif

     !solve_type = GMRES                                   ! === Hypre ===

     call SolveMatrixA( VarU, ResU0,ResU1, Au,U,Su )

   endif
   ! ====================
   ! Velocity component V
   ! ====================
   if( SolveV )then

     URFactor = URF(VarV)

     call SetUpMatrixA( VarV,URFactor,  Av,V,Sv  )

     if( .not. copied )then
       where( Av /= 0.0 )
         Ar = 1.0 / Av
       elsewhere
         Ar = 0.0
       end where
       copied = .true.
       if( Debug > 3 ) write(IOdef,*)'Calculate UVW (V): Ar based on Av'
     endif

     !solve_type = GMRES                                   ! === Hypre ===

     call SolveMatrixA( VarV, ResV0,ResV1, Av,V,Sv )

   endif
   ! ====================
   ! Velocity component W
   ! ====================
   if( SolveW )then

     URFactor = URF(VarW)

     call SetUpMatrixA( VarW,URFactor,  Aw,W,Sw  )

     if( .not. copied )then
       where( Aw /= 0.0 )
         Ar = 1.0 / Aw
       elsewhere
         Ar = 0.0
       end where
       copied = .true.
       if( Debug > 3 ) write(IOdef,*)'Calculate UVW (W): Ar based on Aw'
     endif

     !solve_type = GMRES                                   ! === Hypre ===

     call SolveMatrixA( VarW, ResW0,ResW1, Aw,W,Sw )

  endif

   if( Debug > 2 )then
     write(IOdef,*)'Overzicht residuen:'
     write(IOdef,*)'U:',ResU0,ResU1,ResU1/(ResU0+Small)
     write(IOdef,*)'V:',ResV0,ResV1,ResV1/(ResV0+Small)
     write(IOdef,*)'W:',ResW0,ResW1,ResW1/(ResW0+Small)
   endif

   if( SolveU )then
     if( LimitLowU )then
       ic = 0
       write(IOdbg,*)'Lowest value U:',minval(U(1:Ncel))
       do ip=1,Ncel
	 if( U(ip) < LowerLimitU )then
           ic = ic + 1
           U(ip) = LowerLimitU
	 endif
       end do
       write(IOdbg,*)'Number of cells U changed:',ic,LowerLimitU
     endif
     if( LimitUpU )then
       ic = 0
       write(IOdbg,*)'Highest value U:',maxval(U(1:Ncel))
       do ip=1,Ncel
	 if( U(ip) > UpperLimitU )then
           ic = ic + 1
           U(ip) = UpperLimitU
	 endif
       end do
       write(IOdbg,*)'Number of cells U changed:',ic,UpperLimitU
     endif
   endif

   if( SolveV )then
     if( LimitLowV )then
       ic = 0
       write(IOdbg,*)'Lowest value V:',minval(V(1:Ncel))
       do ip=1,Ncel
	 if( V(ip) < LowerLimitV )then
           ic = ic + 1
           V(ip) = LowerLimitV
	 endif
       end do
       write(IOdbg,*)'Number of cells V changed:',ic,LowerLimitV
     endif
     if( LimitUpV )then
       ic = 0
       write(IOdbg,*)'Highest value V:',maxval(V(1:Ncel))
       do ip=1,Ncel
	 if( V(ip) > UpperLimitV )then
           ic = ic + 1
           V(ip) = UpperLimitV
	 endif
       end do
       write(IOdbg,*)'Number of cells V changed:',ic,UpperLimitV
     endif
   endif

   if( SolveW )then
     if( LimitLowW )then
       ic = 0
       write(IOdbg,*)'Lowest value W:',minval(W(1:Ncel))
       do ip=1,Ncel
	 if( W(ip) < LowerLimitW )then
           ic = ic + 1
           W(ip) = LowerLimitW
	 endif
       end do
       write(IOdbg,*)'Number of cells W changed:',ic,LowerLimitW
     endif
     if( LimitUpW )then
       ic = 0
       write(IOdbg,*)'Highest value W:',maxval(W(1:Ncel))
       do ip=1,Ncel
	 if( W(ip) > UpperLimitW )then
           ic = ic + 1
           W(ip) = UpperLimitW
	 endif
       end do
       write(IOdbg,*)'Number of cells W changed:',ic,UpperLimitW
     endif
   endif

   if( SolveU .or. SolveV .or. SolveW )then
     if( LimitLowVMag )then
       ic  =  0
       dum = 1.e9
       do ip=1,Ncel
	 velmag = sqrt( U(ip)*U(ip) + V(ip)*V(ip) + W(ip)*W(ip) )
	 dum = min(dum,velmag)
	 if( velmag < LowerLimitVMag .and. &
	     velmag > 0.0                  )then
           ic = ic + 1
           fact  = LowerLimitVMag / velmag
	   U(ip) = fact * U(ip)
	   V(ip) = fact * V(ip)
	   W(ip) = fact * W(ip)
	 endif
       end do
       write(IOdbg,*)'Lowest value Velocity Magnitude:',dum
       write(IOdbg,*)'Number of cells magnitude changed:',ic,LowerLimitVMag
     endif
     if( LimitUpVMag )then
       ic  =   0
       dum = -1.e9
       do ip=1,Ncel
	 velmag = sqrt( U(ip)*U(ip) + V(ip)*V(ip) + W(ip)*W(ip) )
	 dum = max(dum,velmag)
	 if( velmag > UpperLimitVMag .and. &
	     velmag > 0.0                  )then
           ic = ic + 1
           fact  = UpperLimitVMag / velmag
	   U(ip) = fact * U(ip)
	   V(ip) = fact * V(ip)
	   W(ip) = fact * W(ip)
	 endif
       end do
       write(IOdbg,*)'Highest value Velocity Magnitude:',dum
       write(IOdbg,*)'Number of cells magnitude changed:',ic,UpperLimitVMag
     endif
   endif

  
   call PrintSummary(IOdbg,iter)

   call watch_leave('CalculateUVW')
  
   if( Debug > 3 ) write(IOdef,*)'=== CalculateUVW'

end subroutine CalculateUVW
subroutine FluxUVW
!========================================================================
!
!   this routine calculates the convective and diffusive fluxes
!   through the cell faces.
!

   use constants
   use geometry
   use variables
   use watches

   real Xac(3), Xpn(3), ds(3),   Up(3), Ub(3), Un(3), Ut(3), &
        Force(3), TotalForce(3), Uw(3), Us(3), Xn(3), Xp(3), Xt(3)

   real d(3), d1(3), d2(3)
   real Xf(3), Xpa(3), Xna(3), dsp(3), dsn(3), dX(3)
   real Tau_nn(3), Tau_nt(3)
   real dUdXac(3), dVdXac(3), dWdXac(3)

   real dUdXP(3),  dVdXP(3),  dWdXP(3), &
        dUdXN(3),  dVdXN(3),  dWdXN(3), &
        dUdXC(3),  dVdXC(3),  dWdXC(3), &
        dUdXD(3),  dVdXD(3),  dWdXD(3), &
        dUdXe(3),  dVdXe(3),  dWdXe(3), &
		vorticityVec(3), shearVec(3), rortexVec(3) , &
		shearHalf(3)
  real  shearM, rortexM, ushror,vshror,wshror
  
  	  real :: Amat(3,3)
	  real :: Bmat(3,3)
	  real :: dltUmat(3,3)
	  real :: omgMethod, retmp1, retmp2, ushtmp,vshtmp, omgBar
	  integer :: boolRort
	  
		
!vorticityVec, shearVec, rortexVec, shearM, rortexM 		

   if( Debug > 3 ) write(IOdef,*)'*** FluxUVW'

   call watch_enter('FluxUVW')

   GammaBlend = Gamma(VarU)
   IScheme    = Scheme(VarU)

   if( Debug > 3 ) write(IOdef,*)'Variable ',VarU,GammaBlend,IScheme

   ! reset counters
   call DiffSchemeCounter(0,0,0,VarU)

   pe0 =  9999.
   pe1 = -9999.
   TotalForce = 0.0

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     it = -1

	   vorticityVec=0.0
	   shearVec=0.0
	   rortexVec=0.0
	   shearM=0.0
	   rortexM=0.0
	   ushror=0.0
	   vshror=0.0
	   wshror=0.0	
	   ushtmp=0.0
	   vshtmp=0.0
	 
     if( in > 0 )then
       !
       ! internal cell face which points to two cells
       !
       facn = Face(i)%lambda
       facp = 1.0 - facn
       !
       ! Viscosity Uncorrected! to be done?
       !
       Xac      =   Cell(in)%x * facn +   Cell(ip)%x * facp

       Uac      =      U(in)   * facn +      U(ip)   * facp
       Vac      =      V(in)   * facn +      V(ip)   * facp
       Wac      =      W(in)   * facn +      W(ip)   * facp

       dUdXac   =   dUdX(in,:) * facn +   dUdX(ip,:) * facp
       dVdXac   =   dVdX(in,:) * facn +   dVdX(ip,:) * facp
       dWdXac   =   dWdX(in,:) * facn +   dWdX(ip,:) * facp

       Visac    =   VisEff(in) * facn +   VisEff(ip) * facp

       Xpn      = Cell(in)%x - Cell(ip)%x

      !VisFace  = Visac * Face(i)%area/vector_length(Xpn)
       VisFace  = Visac * Face(i)%rlencos

       call SelectDiffSchemeVector(i,iScheme,iP,iN,     &
                                   U,V,W,               &
				   dUdX,dVdX,dWdX,      &
				   UFace,VFace,Wface)
       !
       ! explicit higher order convective flux (see eg. eq. 8.16)
       !
       fuce = MassFlux(i) * UFace
       fvce = MassFlux(i) * VFace
       fwce = MassFlux(i) * WFace

       sx = Face(i)%n(1)
       sy = Face(i)%n(2)
       sz = Face(i)%n(3)

       !
       ! explicit higher order diffusive flux based on simple uncorrected
       ! interpolated cell centred gradients(see eg. eq. 8.19)
       !
       fude1 = (dUdXac(1)+dUdXac(1))*sx + (dUdXac(2)+dVdXac(1))*sy + (dUdXac(3)+dWdXac(1))*sz
       fvde1 = (dUdXac(2)+dVdXac(1))*sx + (dVdXac(2)+dVdXac(2))*sy + (dVdXac(3)+dWdXac(2))*sz
       fwde1 = (dUdXac(3)+dWdXac(1))*sx + (dWdXac(2)+dVdXac(3))*sy + (dWdXac(3)+dWdXac(3))*sz

       !
       ! expliciete diffusieve flux/gecorrigeerde diffusie
       ! Jasaks over-relaxed approach (thesis).
       !
       !d     = Cell(in)%x - Cell(ip)%x
       !dXpn  = vector_length( d )
       !S2    = Face(i)%area**2

       !d1    = d *S2/dot_product( d , Face(i)%n )  ! orthogonaal deel
       !d2    = Face(i)%n - d1                      ! niet-orthogonaal deel

       !dUdXe = vector_length(d1)*( UNs - UPs )/dXpn + &
       !        dot_product( d2 , dUdXac )

       !dVdXe = vector_length(d1)*( VNs - VPs )/dXpn + &
       !        dot_product( d2 , dVdXac )

       !dWdXe = vector_length(d1)*( WNs - WPs )/dXpn + &
       !        dot_product( d2 , dWdXac )

       !fude2 = (dUdXe(1)+dUdXe(1)) + (dUdXe(2)+dVdXe(1)) + (dUdXe(3)+dWdXe(1))
       !fvde2 = (dUdXe(2)+dVdXe(1)) + (dVdXe(2)+dVdXe(2)) + (dVdXe(3)+dWdXe(2))
       !fwde2 = (dUdXe(3)+dWdXe(1)) + (dWdXe(2)+dVdXe(3)) + (dWdXe(3)+dWdXe(3))
	   
!subroutine calcShRortex(dudxyz,dvdxyz,dwdxyz, vorticityVec, shearVec, rortexVec, shearM, rortexM )	 
	   

	   if (BoolShearRortex > 0) then
	    call calcShRortex(dUdXac,dVdXac,dWdXac, vorticityVec, shearVec, rortexVec, shearM, rortexM, boolRort )	
		   shearHalf=0.5*shearVec
		   ushror= coefShearRortex(1) * (shearHalf(2)*sz- shearHalf(3)*sy) 
		   vshror= coefShearRortex(2) * (shearHalf(3)*sx- shearHalf(1)*sz) 
		   wshror= coefShearRortex(3) * (shearHalf(1)*sy- shearHalf(2)*sx) 	   
		   
		   dltUmat(1,1:3)= dUdXac
		   dltUmat(2,1:3)= dVdXac
		   dltUmat(3,1:3)= dWdXac
		   retmp1=0.0
		   retmp2=0.0
		   
		   Amat=0.5*(dltUmat + TRANSPOSE(dltUmat) );
		   Bmat=0.5*(dltUmat - TRANSPOSE(dltUmat) );
		   
		   do iii=1,3
		   do jjj=1,3
		      retmp1=retmp1+ Amat(iii,jjj) ** 2
			  retmp2=retmp2+ Bmat(iii,jjj) ** 2
		   enddo
		   enddo
		   
		   omgMethod=retmp2/(retmp1+retmp2+1.0e-6)	
		   omgBar=(abs(omgMethod-0.5)+ (omgMethod-0.5))/(abs(2* omgMethod-1.0) + 1.0e-6)*omgMethod
		   
		   
		   ushror = omgBar*ushror
		   vshror = omgBar*vshror
		   wshror = omgBar*wshror
		   
		   
		!   write(*,*)  ushror, vshror, wshror
		   
		   
	     !ushtmp=0.5*sqrt((dVdXac(2)-dUdXac(1))**2 + (dVdXac(1)+dUdXac(2))**2) 
	     !vshtmp=ushtmp
	     !ushror=omgBar*coefShearRortex(1)*ushtmp*(sx+sy+sz)
		 !vshror=omgBar*coefShearRortex(1)*ushtmp*(sx+sy+sz) 
		 !wshror=0.0		   
	   
		   
	   
		!  if (omgMethod < 0.52) then
		!!   if (BoolRort==0) then
		!     ushror=0.0 
		!     vshror=0.0 
		!     wshror=0.0			 
		!	 !ushtmp=0.0
		!	 !vshtmp=0.0
		!   endif
		   
		   !ushror=coefShearRortex(1)*shearVec(1)*sx + coefShearRortex(1)*shearVec(2)*sy +  &
		   !      coefShearRortex(1)*shearVec(3)*sz
		   !vshror=coefShearRortex(2)*shearVec(1)*sx + coefShearRortex(2)*shearVec(2)*sy +  &
		   !      coefShearRortex(2)*shearVec(3)*sz	
		   !wshror=coefShearRortex(3)*shearVec(1)*sx + coefShearRortex(3)*shearVec(2)*sy +  &
		   !      coefShearRortex(3)*shearVec(3)*sz	
			!	 write(*,*) BoolShearRortex, coefShearRortex(1), shearVec(1), shearVec(2), shearVec(3), sz
		   
	   endif
	   
       fude = Visac * fude1      + ushror
       fvde = Visac * fvde1      + vshror
       fwde = Visac * fwde1      + wshror
	   

	   
       !fude = Visac * fude1      + coefShearRortex(1) *ushtmp *(sx+sy+sz)
       !fvde = Visac * fvde1      + coefShearRortex(1) * vshtmp *(sx+sy+sz)
       !fwde = Visac * fwde1     	   

       !
       ! implicit lower order (simple upwind)
       ! convective and diffusive fluxes
       !
       fmin = min(MassFlux(i),0.0)
       fmax = max(MassFlux(i),0.0)

       fuci = fmin * U(in) + fmax * U(ip) 
       fvci = fmin * V(in) + fmax * V(ip)
       fwci = fmin * W(in) + fmax * W(ip)

       fudi = VisFace * dot_product( dUdXac , Xpn )
       fvdi = VisFace * dot_product( dVdXac , Xpn )
       fwdi = VisFace * dot_product( dWdXac , Xpn )
       !
       ! convective coefficients with deferred correction with
       ! gamma as the blending factor (0.0 <= gamma <= 1.0)
       !
       !      low            high    low  OLD
       ! F = F    + gamma ( F     - F    )
       !     ----   -------------------------
       !      |                  |
       !  implicit           explicit (dump into source term)
       !
       !            diffusion       convection
       !                v               v
       Rface(i,1) = -VisFace - max( MassFlux(i) , 0.0 )  ! P (e)
       Rface(i,2) = -VisFace + min( MassFlux(i) , 0.0 )  ! N (w)

       blend_u = GammaBlend*( fuce - fuci )
       blend_v = GammaBlend*( fvce - fvci )
       blend_w = GammaBlend*( fwce - fwci )
       !
       ! assemble the two source terms
       !
       Su(ip) = Su(ip) - blend_u + fude - fudi 
       Su(in) = Su(in) + blend_u - fude + fudi

       Sv(ip) = Sv(ip) - blend_v + fvde - fvdi
       Sv(in) = Sv(in) + blend_v - fvde + fvdi

       Sw(ip) = Sw(ip) - blend_w + fwde - fwdi
       Sw(in) = Sw(in) + blend_w - fwde + fwdi


       peclet = MassFlux(i)/Face(i)%area*vector_length(Xpn)/(Visac+Small)
       pe0 = min( pe0 , peclet )
       pe1 = max( pe1 , peclet )

     else
       !
       ! boundary face               Ohne Gewaehr!  kritisch blijven !!!!
       !                             ====================================
       ib = Face(i)%bnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       ip = Face(i)%cell1

       if( it == RInlet )then
         dUdXac   = dUdX(ip,:)
         dVdXac   = dVdX(ip,:)
         dWdXac   = dWdX(ip,:)

         Xac      = Face(i)%x
         if( Reg(ir)%user )then
           Utmp  = Reg(ir)%uvw(1)
           Vtmp  = Reg(ir)%uvw(2)
           Wtmp  = Reg(ir)%uvw(3)
           TEtmp = Reg(ir)%k
           EDtmp = Reg(ir)%e
           Ttmp  = Reg(ir)%T
           Dtmp  = Reg(ir)%den

           Xtmp  = Face(i)%x(1)
           Ytmp  = Face(i)%x(2)
           Ztmp  = Face(i)%x(3)
           Atmp  = Face(i)%area

           call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                          Pref,Tref,Iter,Time,                  &
                          Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                          Reg(ir)%name1,Reg(ir)%name2)

           UFace    = Utmp
           VFace    = Vtmp
           WFace    = Wtmp

         else
           UFace    = Reg(ir)%uvw(1)
           VFace    = Reg(ir)%uvw(2)
           WFace    = Reg(ir)%uvw(3)
         endif

         Visac    = VisEff(Ncel+ib)

         Xpn      = Xac - Cell(ip)%x
        !VisFace  = Visac * Face(i)%area/vector_length(Xpn)
         VisFace  = Visac * Face(i)%rlencos

         fuce = MassFlux(i) * UFace
         fvce = MassFlux(i) * VFace
         fwce = MassFlux(i) * WFace

         sx = Face(i)%n(1)
         sy = Face(i)%n(2)
         sz = Face(i)%n(3)

         fude = (dUdXac(1)+dUdXac(1))*sx + (dUdXac(2)+dVdXac(1))*sy + (dUdXac(3)+dWdXac(1))*sz
         fvde = (dUdXac(2)+dVdXac(1))*sx + (dVdXac(2)+dVdXac(2))*sy + (dVdXac(3)+dWdXac(2))*sz
         fwde = (dUdXac(3)+dWdXac(1))*sx + (dWdXac(2)+dVdXac(3))*sy + (dWdXac(3)+dWdXac(3))*sz

         fude = Visac * fude
         fvde = Visac * fvde
         fwde = Visac * fwde

         fmin = min(MassFlux(i),0.0)   ! in inlaat stroomt het naar binnen mf < 0
         fmax = max(MassFlux(i),0.0)   ! fmin /= 0.0 en fmax = 0.0 !

         fuci = fmin * UFace + fmax * U(ip)     !
         fvci = fmin * VFace + fmax * V(ip)     ! voor het geval van neg. inlet ?
         fwci = fmin * WFace + fmax * W(ip)     ! nog sjekken of ok!

         fudi = VisFace * dot_product( dUdXac , Xpn )
         fvdi = VisFace * dot_product( dVdXac , Xpn )
         fwdi = VisFace * dot_product( dWdXac , Xpn )
         !
         ! by definition points a boundary normal outwards
         ! therefore an inlet results in a mass flux < 0.0
         !
         f   = -VisFace + min( MassFlux(i) , 0.0 )

         Au(ip) = Au(ip) - f
         Su(ip) = Su(ip) - f * UFace + fude - fudi
         U(Ncel+ib) = UFace

         Av(ip) = Av(ip) - f
         Sv(ip) = Sv(ip) - f * VFace + fvde - fvdi
         V(Ncel+ib) = VFace

         Aw(ip) = Aw(ip) - f
         Sw(ip) = Sw(ip) - f * WFace + fwde - fwdi
         W(Ncel+ib) = WFace

       else if( it == ROutlet )then

         dUdXac   = dUdX(ip,:)
         dVdXac   = dVdX(ip,:)
         dWdXac   = dWdX(ip,:)

         Xac      = Face(i)%x
         Visac    = VisEff(ip)

         Xpn      = Xac - Cell(ip)%x

        !UFace    = U(ip) + dot_product( dUdX(ip,:) , Xpn )
        !VFace    = V(ip) + dot_product( dVdX(ip,:) , Xpn )
        !WFace    = W(ip) + dot_product( dWdX(ip,:) , Xpn )
        !
        ! modification by Harry, discard the gradient
        !
         UFace    = U(ip)
         VFace    = V(ip)
         WFace    = W(ip)

        !VisFace  = Visac * Face(i)%area/vector_length(Xpn)
         VisFace  = Visac * Face(i)%rlencos

         fuce = MassFlux(i) * UFace
         fvce = MassFlux(i) * VFace
         fwce = MassFlux(i) * WFace

         sx = Face(i)%n(1)
         sy = Face(i)%n(2)
         sz = Face(i)%n(3)

         fude = (dUdXac(1)+dUdXac(1))*sx + (dUdXac(2)+dVdXac(1))*sy + (dUdXac(3)+dWdXac(1))*sz
         fvde = (dVdXac(1)+dUdXac(2))*sx + (dVdXac(2)+dVdXac(2))*sy + (dVdXac(3)+dWdXac(2))*sz
         fwde = (dWdXac(1)+dUdXac(3))*sx + (dWdXac(2)+dVdXac(3))*sy + (dWdXac(3)+dWdXac(3))*sz

         fude = Visac * fude
         fvde = Visac * fvde
         fwde = Visac * fwde

         fmin = min(MassFlux(i),0.0)   ! in uitlaat stroomt het naar buiten mf >=0
         fmax = max(MassFlux(i),0.0)   ! fmin = 0.0 en fmax > 0 !

         fuci = fmin * UFace + fmax * U(ip)
         fvci = fmin * VFace + fmax * V(ip)
         fwci = fmin * WFace + fmax * W(ip)

         fudi = VisFace * dot_product( dUdXac , Xpn )
         fvdi = VisFace * dot_product( dVdXac , Xpn )
         fwdi = VisFace * dot_product( dWdXac , Xpn )
         !
         ! by definition points a boundary normal outwards
         ! therefore an outlet results in a mass flux >= 0.0
         !
         if( MassFlux(i) < 0.0 )then
           write(IOdef,*)'+ internal error: neg. massflux in outlet',massflux(i)
           MassFlux(i) = Small
         endif
         f    = -VisFace + min( MassFlux(i) , 0.0 )

         Au(ip) = Au(ip) - f
         Su(ip) = Su(ip) - f * UFace + fude - fudi
         U(Ncel+ib) = UFace

         Av(ip) = Av(ip) - f
         Sv(ip) = Sv(ip) - f * VFace + fvde - fvdi
         V(Ncel+ib) = VFace

         Aw(ip) = Aw(ip) - f
         Sw(ip) = Sw(ip) - f * WFace + fwde - fwdi
         W(Ncel+ib) = WFace

        ! if( massflux(i) == 0.0 )then
        !   write(IOdef,*)'fluxuvw:',i,MassFlux(i),f
        !   write(IOdef,*)'    uvw:',Uface,Vface,Wface
        !   write(IOdef,*)'  uvw_P:',U(ip),V(ip),W(ip)
        !   write(IOdef,*)'     sx:',sx,sy,sz
        ! endif

       else if( it == RSymp )then
         !
         ! er is alleen een schuifspanning _|_ op het vlak
         ! zie hoofdstuk 8.10.4
         !
         Xac     = Face(i)%x                 ! face centre
         Xpn     = Xac - Cell(ip)%x          ! dist. Face to P

         dUdXac  = dUdX(ip,:)                !
         dVdXac  = dVdX(ip,:)                ! use values in symm.
         dWdXac  = dWdX(ip,:)                ! cell center

         Visac   = VisEff(ip)                !

        !VisFace = Visac * Face(i)%area/vector_length(Xpn)
         VisFace = Visac * Face(i)%rlencos

        !
        ! modification by Harry, discard the gradient
        !
         dU      = dot_product( dUdXac , Xpn )
         dV      = dot_product( dVdXac , Xpn )
         dW      = dot_product( dWdXac , Xpn )

         Us(1)   = U(ip) + dU
         Us(2)   = V(ip) + dV
         Us(3)   = W(ip) + dW

        !Us(1)   = U(ip)
        !Us(2)   = V(ip)
        !Us(3)   = W(ip)

         Xn = -Face(i)%n                     ! flip the normal inwards
         call normalise(Xn)                  ! normaal vector met lengte 1

         dp = dot_product( Us , Xn )
         Un = dp * Xn                        ! normale snelh. component

         dn     = Bnd(ib)%distance

         Tau_nn = 2.0 * Visac * Un / dn      ! eq. (8.75)
         Force  = Tau_nn * Face(i)%area !&   ! tau_nn zorgt voor de richting
                                        !* abs(Xn) ! dus |Xn| ontbindt juist

         !if( ip == 1 )then
         !  write(IOdef,*)' '
         !  write(IOdef,*)'S ip xac:',ip,Xac
         !  write(IOdef,*)'S ir xn :',ir,Xn
         !  write(IOdef,*)'S ib un :',ib,un
         !  write(IOdef,*)'S i  tau:',i,tau_nn
         !  write(IOdef,'(A,i3,5(1x1pe10.3))')' S force:',ip,Force(1:3),Bnd(ib)%yplus,dn
         !endif
         !
         !
         !
         if( .not. Initialisation )then
           Au(ip) = Au(ip) + VisFace
           Av(ip) = Av(ip) + VisFace
           Aw(ip) = Aw(ip) + VisFace

           Su(ip) = Su(ip) + VisFace*U(ip) - Force(1)
           Sv(ip) = Sv(ip) + VisFace*V(ip) - Force(2)
           Sw(ip) = Sw(ip) + VisFace*W(ip) - Force(3)
         endif

         U(Ncel+ib) = Us(1)
         V(Ncel+ib) = Us(2)
         W(Ncel+ib) = Us(3)

       else if( it == RWall )then
         !
         ! schuifspanning in het vlak; hoofdstuk 8.10.3
         !
         Xac      = Face(i)%x                 ! face centre
         Uw(1)    = Reg(ir)%uvw(1)            !
         Uw(2)    = Reg(ir)%uvw(2)            ! wall bc's
         Uw(3)    = Reg(ir)%uvw(3)            !

         dUdXac   = dUdX(ip,:)                !
         dVdXac   = dVdX(ip,:)                ! use values in wall
         dWdXac   = dWdX(ip,:)                ! cell center

         Visac    = VisEff(Ncel+ib)           !

         Xpn      = Face(i)%x - Cell(ip)%x    ! dist. wall center to P

         coef = Visac * Face(i)%area / vector_length(Xpn)

         Xn = Face(i)%n                       ! de surface vector
         call normalise(Xn)                   ! normaal vector met lengte 1

         Up(1) = U(ip)
         Up(2) = V(ip)
         Up(3) = W(ip)

         Up = Up - Uw                         ! snelh. verschil

         dp = dot_product( Up , Xn )
         Un = dp * Xn                         ! normale snelh. component
         Ut = Up - Un                         ! tang. snelh. component

         Uvel = abs(Ut(1)) + abs(Ut(2)) + abs(Ut(3))
         if( Uvel > Small )then               ! simpele test if |Ut| > 0

           dn     = Bnd(ib)%distance
           Tau_nt = Visac * Ut / dn
           Force  = Tau_nt * Face(i)%area     ! tau_nt zorgt voor de richting

           Bnd(ib)%shear = Force              ! absolute force (in N)

         else

           Force = 0.0                        ! voor het geval dat...
           Bnd(ib)%shear = Force

         endif

         if( .not. Initialisation )then

           TotalForce = TotalForce + Force

           !               standard
           !               implicit
           !                  V
           Au(ip) = Au(ip) + coef       !
           Av(ip) = Av(ip) + coef       ! deze 3 MOETEN gelijk
           Aw(ip) = Aw(ip) + coef       ! zijn ivm druk iteratie
           !
           !                    corr.     expliciet
           !                  impliciet
           !                     V           V
           Su(ip) = Su(ip) + coef*U(ip) - Force(1)
           Sv(ip) = Sv(ip) + coef*V(ip) - Force(2)
           Sw(ip) = Sw(ip) + coef*W(ip) - Force(3)

         endif

         U(Ncel+ib) = Uw(1)
         V(Ncel+ib) = Uw(2)
         W(Ncel+ib) = Uw(3)


       endif

     endif

   end do faceloop

   ! counters to IOdbg
   call DiffSchemeCounter(2,0,0,VarU)

   !if( Debug > 0 )then
   !  write(IOdef,*) 'UVW: ',pe0,' < Peclet < ',pe1
   !  write(IOdef,*) 'Total shear force:',TotalForce
   !endif
   write(IOdbg,*) 'UVW: ',pe0,' < Peclet < ',pe1

   do ir=0,Nreg
    !Reg(ir)%ForceNormal(:) = 0.0
     Reg(ir)%ForceTangnt(:) = 0.0
   end do

   do ib=1,Nbnd
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
     if( it == Rwall) &
       Reg(ir)%ForceTangnt = Reg(ir)%ForceTangnt + Bnd(ib)%shear
   end do

   write(IOdbg,*) 'Total Shear force:',TotalForce
   do ir=0,Nreg
     if( Reg(ir)%typ == Rwall ) &
       write(IOdbg,*) 'Shear Force, reg ',ir,'=',Reg(ir)%ForceTangnt(:)
   end do

   call watch_leave('FluxUVW')

   if( Debug > 3 ) write(IOdef,*)'=== FluxUVW'

end subroutine FluxUVW
subroutine CalculatePressure(dppmax)
!========================================================================
!     This routine assembles and solves the pressure-correction
!     equation using colocated grid. SIMPLE algorithm with one
!     or more corrector steps (non-orthogonality effects taken
!     into account as described in Sect. 8.8.
!
   use constants
   use geometry
   use variables
   use watches

   !use hypre_mod, only: new_matrix, solve_type           ! === Hypre ===
   use artificial_compressibility                         ! === AC ===

   real, dimension(3) :: PressureForce, ds, Xn

   logical, save :: ready  = .false., break = .false.
   integer, save :: idumpA = 0
   integer, save :: iprint = 0

   integer ipar(16)
   real*4  fpar(16) ! hint by Shibo

   double precision, dimension(:), allocatable :: Alu
   integer, dimension(:), allocatable :: Jlu, Ju, Jw

   if( Debug > 3 ) write(IOdef,*)'*** CalculatePressure'

   call watch_enter('CalculatePressure')

   if( UseArtificialComp )then                            ! === AC ===
     if( Debug > 3 ) write(IOdef,*)'Allocating art. compr. arrays'
     allocate( art_snd_speed(Ncel) )
     allocate( delta_t(Ncel) )
     allocate( density_change(Ncel) )
   endif

   dppmax = Large

   iloop = 0
 1 continue
   iloop = iloop + 1

   Au  = 0.0         ! re-use of the U-velocity component arrays
   Su  = 0.0         !
   MassFlux = 0.0
   PressureForce = 0.0

   Rface = 0.0       ! safety first; op nul zetten kan straks weg!
   !
   ! calculate fluxes through all inner faces
   ! in this routine the coef for the p'-equation are set
   ! Au and Su are set here
   !
   call FluxMass

   !
   ! assemble p' equations
   !
   ready = .false.
   break = .false.
   PP(:) = 0.0
   ia    = 0                     ! index for A-matrix

   reference = 0.0
   if( UseArtificialComp )then                            ! === AC ===

      call artificial_comp()

    endif

   do i=1,Ncel

     app = 0.0                   ! <= *** LET OP *** AANGEPAST

     do j=1,NFaces(i)
       k  = CFace(i,j)
       ip = Face(k)%cell1
       in = Face(k)%cell2
       if( in > 0 )then

         ! internal faces only
         if( ip == i )then
           aface   = RFace(k,2)
           ia = ia + 1
           Acoo(ia) = DBLE( aface )
           Arow(ia) = i
           Acol(ia) = in

           Su(i) = Su(i) - MassFlux(k)

         else if( in == i )then
           aface = RFace(k,1)
           ia = ia + 1
           Acoo(ia) = DBLE( aface )
           Arow(ia) = i
           Acol(ia) = ip

           Su(i) = Su(i) + MassFlux(k)

         else
           write(IOdef,*)'+ internal error: assembly A-matrix. in=',in
         endif

         app = app - aface

       endif
     end do

     !
     ! andere solver nog voor echt perfect symm. A-matrix?
     ! lijkt niet per se nodig; scheelt zeker in cpu-tijd
     !
     if( UseArtificialComp )then                          ! === AC ===
       Au(i)  = app + (Cell(i)%vol)/(art_snd_speed(i)**2*delta_t(i))
     else
       Au(i)  = app
     endif

     ia = ia + 1
     Acoo(ia) = DBLE( Au(i) )
     Arow(ia) = i
     Acol(ia) = i

     !write(31,'(i6,1x,i6,1x,e12.5)') i,i,app

   end do

   if( Debug > 2 )then
     sumsu = sum( Su(:) )
     write(IOdef,*)'Starting residual P:',sumsu
   endif

   if( Transient .and. SolveDensity )then
     ! iets met de dichtheid
   endif

   ! opendx dump
   !close(31)

   if( ia /= NNZ ) write(IOdef,*)'+ error: Calculate P: NNZ =',ia,' =/=',NNZ
   !
   ! solving equation system Ax=B with x=Phi
   ! even makkelijk omgooien later direct in CSR-formaat
   !
   call coocsr(Ncel,NNZ,Acoo,Arow,Acol,Acsr,Aclc,Arwc)

   !new_matrix = .true.                                   ! === Hypre ===

   ipcor = 0

   do while( .not. ready )

     ipcor = ipcor + 1

     if( Debug > 3 )then
       write(IOdef,*) 'Ervoor:'
       write(IOdef,*) minval(PP(1:Ncel)),'< PP <',maxval(PP(1:Ncel))
     endif

     RHS(1:Ncel) = DBLE( Su(1:Ncel) )
     SOL(1:Ncel) = DBLE( PP(1:Ncel) )

     if( Solver(VarP) == SparseKit2 )then

       ipar  =  0
       fpar  = 0.0

       ipar(2) = 0           ! 0 == no preconditioning
       ipar(3) = 0
       ipar(4) = NNZ*8       ! workspace for BGStab
       ipar(6) = 800         ! maximum number of matrix-vector multiplies

       fpar(1) = RTOL(VarP)  ! relative tolerance, must be between (0, 1)
       fpar(2) = ATOL(VarP)  ! absolute tolerance, must be positive

       write(IOdbg,*)'SolveMatrixA: SparsKit2 pressure'

  allocate( Alu(1) ) 
  allocate( Jlu(1) )
  allocate( Ju(1)  )
  allocate( Jw(1)  )

      !call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
      !                       WORK,Acsr,Aclc,Arwc)
       call solve_spars(debug,IOdbg,Ncel,RHS,SOL,ipar,fpar,&
                              WORK, Acsr,Aclc,Arwc,        &
	                      NNZ,  Alu,Jlu,Ju,Jw          )
 
  deallocate( Alu  )
  deallocate( Jlu  )
  deallocate(  Ju  )
  deallocate(  Jw  )

       if( ipar(1) ==  -3 ) break = .true.

       if( ipar(1) ==  -1 ) Flags(IFlagP) = 'p'

     else if( Solver(VarP) == SolverHypre )then           ! === Hypre ===

      !if( ipcor == 1 )then
      !  call solve_hypre(VarP,RHS,SOL,0,SolverResidual)
      !else
      !  call solve_hypre(VarP,RHS,SOL,1,SolverResidual)
      !endif

     else if( Solver(VarP) == SolverUser )then

       write(IOdef,*)'+ User solver for pressure'
       !call solve_user(debug,IOdbg,Ncel,NNZ,Acoo,Arow,Acol,S,Phi)

     else

       write(IOdef,*)'+ internal error: solver not set'

     endif

     PP(1:Ncel) = REAL( SOL(1:Ncel) )

     if( Debug > 2 )then
       write(IOdef,*) minval(PP(1:Ncel)),'< PP <',maxval(PP(1:Ncel))
       write(IOdef,*) 'Using tol.:',RTOL(VarP),ATOL(VarP)
     endif
     dppmax = maxval( abs(pp) )
     if( ipcor == 1 ) dppref = dppmax
     iprint = iprint + 1

     !
     ! check if residual get's bigger
     !
     if( .not. break )then
       if( ipcor == 2 )then
         reference = fpar(6)
         if( Debug > 2 ) write(IOdef,*)'Using Residual:',reference
       else if( ipcor > 2 )then
         if( fpar(6) > reference )then
           break = .true.
         else
           reference = fpar(6)
         if( Debug > 2 ) write(IOdef,*)'Using new Residual:',reference
         endif
       endif
     endif
     !
     ! update P' (PP) along boundaries, needed for gradient
     !
     call Update_P_at_boundaries(VarPP,PP)
     call GradientPhi(VarPP,PP,dPdX)

     do i=1,Nfac
       ip = Face(i)%cell1
       in = Face(i)%cell2
       if( in > 0 )then

         MassFlux(i) = MassFlux(i) + RFace(i,1)*( PP(in) - PP(ip) )

       endif
     end do

     PPref = pp(IPref)             ! <= IPref voor material 1
     URFp  = URF(VarP)

     if( UseArtificialComp )then                          ! === AC ===
       do i=1,NCel
         P(i) = P(i) + URFp*( PP(i) - PPref )
         density_change(i) = PP(i)/art_snd_speed(i)

         fact = Cell(i)%vol * Ar(i)
         if( SolveU ) U(i) = U(i) - dPdX(i,1) * fact + &
                            (density_change(i) * U(i))/Densref
         if( SolveV ) V(i) = V(i) - dPdX(i,2) * fact + &
                            (density_change(i) * V(i))/Densref
         if( SolveW ) W(i) = W(i) - dPdX(i,3) * fact + &
                            (density_change(i) * W(i))/Densref
       end do
     else
       do i=1,NCel
         P(i) = P(i) + URFp*( PP(i) - PPref )

         fact = Cell(i)%vol * Ar(i)
         if( SolveU ) U(i) = U(i) - dPdX(i,1) * fact
         if( SolveV ) V(i) = V(i) - dPdX(i,2) * fact
         if( SolveW ) W(i) = W(i) - dPdX(i,3) * fact

       end do
     endif
     !
     ! een residu bouwen
     !
     Res(:) = 0.0
     do i=1,Ncel
       do j=1,NFaces(i)
         k  = CFace(i,j)
         ip = Face(k)%cell1
         in = Face(k)%cell2
         if( in > 0 )then
           ! internal faces
           if( ip == i )then
             Res(i) = Res(i) - MassFlux(k)
           else if( in == i )then
             Res(i) = Res(i) + MassFlux(k)
           endif
         else
           ! boundary face
           Res(i) = Res(i) - MassFlux(k)
         endif
       end do
     end do
     Residual(VarP) = sum(abs(Res)) * ResiNorm(VarP)

     Su = 0.0
     P(IPref) = 0.0               ! <= let op mat(1) en in input deck!

     !
     ! update the source term Su for the second step
     !
     call FluxMass2

     PP = 0.0                     ! reset p' (=p''), is this really needed?

     sumsu = sum( Su(:) )

     if( debug > 2 )write(IOdef,*)'P: ',ipcor,reference,Residual(VarP),sumsu,dppmax

     !
     ! stopping criteria
     !
     if( ipcor >= MaxPCOR .or. dppmax <= FactDPP*dppref .or. break ) ready = .true.

   end do

   write(IOdbg,*)'Number of pressure corrections:',ipcor,dppmax,' of',maxpcor,FactDPP 
   if( Debug > 1 )then
     write(IOdef,*)'Number of pressure corrections:',ipcor,dppmax
   endif

   if( dppmax > Small )then
     dppmax = dppref
   else
     dppmax = 0.0
   endif
   !
   ! update pressure P along boundaries
   !
   call Update_P_at_boundaries(VarP,P)
   call GradientPhi(VarP,P,dPdX)

   call UpdateP(P,dPdX)


   if( UseArtificialComp )then                            ! === AC ===

     deallocate( art_snd_speed )
     deallocate( delta_t )
     deallocate( density_change )

   endif
!  if( SolverHypre )then                                 ! === Hypre ===
!
!    call hypre_destroy()
!
!  endif
   !
   ! extract pressure force on solid walls
   !
   ! real    :: ForceNormal(3) = (/0.0,0.0,0.0/) 
   ! real    :: ForceTangnt(3) = (/0.0,0.0,0.0/) 
   
   PressureForce = 0.0

   do ir=0,Nreg
     Reg(ir)%ForceNormal(:) = 0.0
    !Reg(ir)%ForceTangnt(:) = 0.0
   end do
   
   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == Rwall)then

       ds    = Face(i)%x - Cell(ip)%x                     ! distance
       Pwall = P(ip) + dot_product( dPdX(ip,:) , ds )     ! extrapolation
       Xn    = Face(i)%n                                  ! surface vector

      !Bnd(ib)%normal = Pwall * Xn
       Reg(ir)%ForceNormal = Reg(ir)%ForceNormal + Pwall * Xn
       
       PressureForce = PressureForce + Pwall * Xn

     endif
   end do

   if( (PressureForce(1)+PressureForce(2)+PressureForce(3)) > Small )then

     write(IOdbg,*) 'Total pressure force:',PressureForce
     do ir=0,Nreg
       if( Reg(ir)%typ == Rwall ) &
         write(IOdbg,*) 'Normal Pressure Force, reg ',ir,'=',Reg(ir)%ForceNormal(:)
     end do

   endif

   if( SolveP )then
     if( LimitLowP )then
       ic = 0
       write(IOdbg,*)'Lowest value P:',minval(P(1:Ncel))
       do ip=1,Ncel
	 if( P(ip) < LowerLimitP )then
           ic = ic + 1
           P(ip) = LowerLimitP
	 endif
       end do
       write(IOdbg,*)'Number of cells P changed:',ic,LowerLimitP
     endif
     if( LimitUpP )then
       ic = 0
       write(IOdbg,*)'Highest value P:',maxval(P(1:Ncel))
       do ip=1,Ncel
	 if( P(ip) > UpperLimitP )then
           ic = ic + 1
           P(ip) = UpperLimitP
	 endif
       end do
       write(IOdbg,*)'Number of cells P changed:',ic,UpperLimitP
     endif
   endif


   call watch_leave('CalculatePressure')

   if( Debug > 3 ) write(IOdef,*)'=== CalculatePressure'

end subroutine CalculatePressure
subroutine FluxMass
!========================================================================

   use constants
   use geometry
   use variables
   use scalars
   use watches

   real :: facn, facp
   real :: Xac(3), Xpn(3), Xn(3), Xface(3), Xnorm(3)
   real :: Xpac(3), Xnac(3), Xnp(3), Xnn(3)
   real :: deln(3), delp(3), Xpn2(3)
   real :: dUdXac(3), dVdXac(3), dWdXac(3) ,delta(3)
   real :: Uin(3), Uout(3), Uf(3)

  !real :: FlowRegion(30), FlowFact(30)   ! MOET BETER
   real, allocatable, dimension(:) ::  FlowRegion, FlowFact

   logical :: WarningOutlet

   if( Debug > 3 )write(IOdef,*)'*** FluxMass'

   call watch_enter('FluxMass')

   allocate(FlowRegion(0:Nreg),stat=istat)
   call TrackMemory(istat,Nreg+1,'FlowRegion array allocated')
   allocate(FlowFact(0:Nreg),stat=istat)
   call TrackMemory(istat,Nreg+1,'FlowFact array allocated')

   WarningOutlet = .false.

   icinl = 0
   icout = 0
   icsym = 0
   icwal = 0

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2

     if( in > 0 )then

       facn = Face(i)%lambda
       facp = 1.0 - facn

       dUdXac = dUdX(in,:) * facn + dUdX(ip,:) * facp
       dVdXac = dVdX(in,:) * facn + dVdX(ip,:) * facp
       dWdXac = dWdX(in,:) * facn + dWdX(ip,:) * facp

       !
       ! density here with a CD-like scheme;
       ! to be done introduce a blend with UD (small amount, 5-10%)
       !
       Denf   = Den(in)* facn + Den(ip)* facp

       Xac    = Cell(in)%x * facn + Cell(ip)%x * facp

       Xface  = Face(i)%x
       delta  = Xface - Xac

       Uface  = U(in)*facn + U(ip)*facp + dot_product( dUdXac , delta )
       Vface  = V(in)*facn + V(ip)*facp + dot_product( dVdXac , delta )
       Wface  = W(in)*facn + W(ip)*facp + dot_product( dWdXac , delta )

       MassFlux(i) = Denf * ( Uface*Face(i)%n(1) + &   !
                              Vface*Face(i)%n(2) + &   ! = rho.v_n
                              Wface*Face(i)%n(3) )     !

       !
       ! this was the easy part
       !

       Xnorm = Face(i)%n
       call normalise(Xnorm)

       !
       ! find auxiliary nodes  (eq. 8.53)
       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  P    e       e   P
       !
      !Xpac = Xface - dot_product(Xface-Cell(ip)%x,Xnorm)*Xnorm
      !
      !Xnp  = Xpac - Xface
      !
       ! _    _       _   _   _  _
       ! r ,= r  - [( r - r ).n] n
       !  E    e       e   E
       !
      !Xnac = Xface - dot_product(Xface-Cell(in)%x,Xnorm)*Xnorm
      !
      !Xnn  = Xnac - Xface
      !
       !
       ! use the shortest
       !
      !dp   = vector_length(Xnp)
      !dn   = vector_length(Xnn)
      !
      !if( dn >= dp )then
      !  Xpac2 = Xface + Xnp
      !  Xnac2 = Xface - Xnp
      !  dl = 2 * dp
      !else
      !  Xpac2 = Xface - Xnn
      !  Xnac2 = Xface + Xnn
      !  dl = 2 * dn
      !endif

       !
       ! preset in subroutine ReadGeometry
       !
       Xpac = Face(i)%Xpac     
       Xnac = Face(i)%Xnac     
       
       !
       ! set p ,and p , (see eqn. (8.54))
       !      P      E
       !
       delp = Xpac - Cell(ip)%x
       Pip  = P(ip) + dot_product( dPdX(ip,:) , delp )

       deln = Xnac - Cell(in)%x
       Pin  = P(in) + dot_product( dPdX(in,:) , deln )

       Xpn  = Xnac - Xpac
       Xpn2 = Cell(in)%x - Cell(ip)%x

       apv1 = Den(ip)*Ar(ip)         
       apv2 = Den(in)*Ar(in)          
      !ApV  = 0.5 *( apv1 + apv2 )                               ! simple averaging 1/ApV
       ApV  = apv2 * facn + apv1 * facp                          ! average 1/ApV in 8.56

       FactV= Cell(in)%vol * facn + Cell(ip)%vol* facp           ! volume Omega_e in 8.56

       !
       ! eq. 8.59 starts with (rho dOmega S)_e * (1/Apv)_e
       ! the term (rho 1/Apv)_e is ApV and weighted
       ! the volume term is weighted in FactV
       ! the area S divided by (rE-rP).n is calculated here
       !
       ApV  = ApV*Face(i)%area * FactV/dot_product(Xpn2,Xnorm)

       dpx  = ( dPdX(in,1) * facn + dPdX(ip,1) * facp ) * Xpn(1) !
       dpy  = ( dPdX(in,2) * facn + dPdX(ip,2) * facp ) * Xpn(2) ! (gradP).Xpn term
       dpz  = ( dPdX(in,3) * facn + dPdX(ip,3) * facp ) * Xpn(3) !

       fact = ApV

       Rface(i,1) = -fact                              ! the coefs for the
       Rface(i,2) = -fact                              ! P'-equation see 8.56


       !
       ! See eq. 8.56 FP 3rd edition
       !
       !            first part           the term in square brackets
       !
       !               m*           /delta Omega     1
       !              v          - --------------- (--)[ .. ]
       !               n           ( r   - r  ).n   Ap
       !                              E     P
       !
       MassFlux(i) = MassFlux(i) - fact*( (Pin-Pip) - dpx - dpy - dpz )

     else
       !
       ! boundary face
       !
       ib = Face(i)%bnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       ip = Face(i)%cell1

       if( it == RInlet )then
         icinl = icinl + 1
         !
         ! inlet simply has prescribed velocities
         ! get them from the input (or use user data)
         !

         Xac    = Face(i)%x

         if( Reg(ir)%user )then
           Utmp  = Reg(ir)%uvw(1)
           Vtmp  = Reg(ir)%uvw(2)
           Wtmp  = Reg(ir)%uvw(3)
           TEtmp = Reg(ir)%k
           EDtmp = Reg(ir)%e
           Ttmp  = Reg(ir)%T
           Dtmp  = Reg(ir)%den

           Xtmp  = Face(i)%x(1)
           Ytmp  = Face(i)%x(2)
           Ztmp  = Face(i)%x(3)
           Atmp  = Face(i)%area

           call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                          Pref,Tref,Iter,Time,                  &
                          Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                          Reg(ir)%name1,Reg(ir)%name2)

           Uin(1) = Utmp
           Uin(2) = Vtmp
           Uin(3) = Wtmp

           dens   = Dtmp
         else
           Uin(1) = Reg(ir)%uvw(1)
           Uin(2) = Reg(ir)%uvw(2)
           Uin(3) = Reg(ir)%uvw(3)
           dens   = Reg(ir)%den
         endif

         ! piet's stukje
         !Uin(1) = Face(i)%x(2)*(1.0-Face(i)%x(2))
         !
         ! stuwpunt
         !Uin(1)    =  Face(i)%x(1)
         !Uin(2)    = -Face(i)%x(2)
         !
         ! standaard test
         !Xac(3) = 0.0
         !Umag   = vector_length( Xac )
         !Uin(1) = Umag * 0.866025
         !Uin(2) = Umag * 0.5
         !write(IOdef,*)'Uin:',umag,Uin

         MassFlux(i) = dens * dot_product( Uin , Face(i)%n )

         su(ip) = su(ip) - MassFlux(i)    ! mf < 0 => source

         !write(IOdef,*)'inlet:',ip,Uin,MassFlux(i)

       else if( it == ROutlet )then
         icout = icout + 1
         !
         ! outlet simply looks at node P, correct using gradient
         !
         delta  = Face(i)%x - Cell(ip)%x

        !
        ! modification by Harry, discard the gradient
        !
        !Uface  = U(ip) + dot_product( dUdX(ip,:) , delta )
        !Vface  = V(ip) + dot_product( dVdX(ip,:) , delta )
        !Wface  = W(ip) + dot_product( dWdX(ip,:) , delta )
        !
         Uface  = U(ip)
         Vface  = V(ip)
         Wface  = W(ip)

         Denf   = Den(ip)                            ! NIET GECORR.
         Den(Ncel+ib) = Denf

         MassFlux(i) = Denf * ( Uface*Face(i)%n(1) + &   !
                                Vface*Face(i)%n(2) + &   ! = rho.v_n
                                Wface*Face(i)%n(3) )     !
         !
         ! For an outlet MassFlux must be 0.0 or positive
         !
         if( MassFlux(i) < 0.0 )then
           if( .not. WarningOutlet )then
             write(IOdbg,*)'+ Warning: Inflow detected on outflow boundary'
             Flags(IFlagInflow) = 'i'
             WarningOutlet      = .true.
           endif
           MassFlux(i) = Small

           !
           ! to be sure reset add. variables too
           !
           if( SolveTurbEnergy ) TE(Ncel+ib) = TE(ip)
           if( SolveTurbDiss   ) ED(Ncel+ib) = ED(ip)
           if( SolveVisc       ) VisEff(Ncel+ib) = VisEff(ip)
           if( SolveEnthalpy   ) T(Ncel+ib) = T(ip)
           if( SolveScalars    ) SC(Ncel+ib,1:Nscal) = SC(ip,1:Nscal)
         endif
         !
         ! set Su below after the mass flux into the domain is known
         !
       else if( it == RSymp )then
         icsym = icsym + 1
         MassFlux(i) = 0.0
       else if( it == RWall )then
         icwal = icwal + 1
         MassFlux(i) = 0.0
       else
         write(IOdef,*)'MassFlux: unknown bc, type ',it
       endif
     endif
   end do faceloop
   if( Debug > 2 )then
     write(IOdef,*)'MassFlux: bc inlet :',icinl
     write(IOdef,*)'          bc outlet:',icout
     write(IOdef,*)'          bc symp. :',icsym
     write(IOdef,*)'          bc wall  :',icwal
   endif

   if( Debug > 2 )write(IOdef,*)'Checking global mass flux'
   i1 = count( Reg%typ == RInlet  )
   i2 = count( Reg%typ == ROutlet )

   if( i1 > 0 .and. i2 == 0 )write(IOdef,*)'FOUT IN BC'

   if( i1 > 0 )then                        ! flow with in/out
     FlowRegion = 0.0                      ! reset to zero

     flowin = 0.0
     do ib=1,Nbnd
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == RInlet )then
         i  = Bnd(ib)%face
         flowin = flowin + MassFlux(i)
       endif
     end do
     write(IOdbg,*)'Mass flow rate in: ',FlowIn

     flowout = 0.0
     do ib=1,Nbnd                         ! flowsplit moet nog!
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == ROutlet )then
         i  = Bnd(ib)%face
         FlowOut = FlowOut + MassFlux(i)
         FlowRegion(ir) = FlowRegion(ir) + MassFlux(i)
       endif
     end do
     write(IOdbg,*)'Mass flow rate out:',FlowOut
     do ir=1,Nreg
       it = Reg(ir)%typ
       if( it == ROutlet ) &
         write(IOdbg,*)'Region ',ir,FlowRegion(ir),FlowRegion(ir)/FlowOut
     end do
     do ir=1,Nreg
       it = Reg(ir)%typ
       if( it == ROutlet )then
         split = Reg(ir)%splvl
         FlowFact(ir) = -(split * FlowIn) / FlowRegion(ir)
         write(IOdbg,*)'Fact Region ',ir,split,FlowFact(ir),split * FlowIn
       else
         FlowFact(ir) = 0.0
       endif
     end do

     if( flowout < Small )then
       ! nothing set yet, thus an initial guess,
       ! start with calculating area
       areaout = 0.0
       do ib=1,Nbnd
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         if( it == ROutlet )then
           i  = Bnd(ib)%face
           areaout = areaout + Face(i)%area
         endif
       end do
       ratearea = -flowin/areaout
       write(IOdbg,*)'Flowrate per area:',areaout,' m2',ratearea,' kg/s/m2'

       flowout = 0.0
       do ib=1,Nbnd
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         if( it == ROutlet )then
           i     = Bnd(ib)%face
           split = Reg(ir)%splvl
           MassFlux(i) = ratearea*Face(i)%area * split   ! flowsplit
           FaceFlux    = MassFlux(i)/Den(Ncel+ib)/Face(i)%area

           Xnorm = Face(i)%n
           call normalise(Xnorm)
           
           U(Ncel+ib) = FaceFlux*Xnorm(1)
           V(Ncel+ib) = FaceFlux*Xnorm(2)
           W(Ncel+ib) = FaceFlux*Xnorm(3)

           flowout = flowout + MassFlux(i)
         endif
       end do
       write(IOdbg,*)'Initial flow rate out: ',flowout
     endif

     fact = -flowin /( flowout + Small )         ! amount of unbalance
     write(IOdbg,*)'Correct flow rate:',fact

     if( abs(1.0-fact) > 0.001 ) Flags(IFlagMass) = 'm'

     flowout2 = 0.0
     do ib=1,Nbnd                                 
       ir = Bnd(ib)%rid
       it = Reg(ir)%typ
       if( it == ROutlet )then                   
         i     = Bnd(ib)%face

         MassFlux(i) = MassFlux(i) * FlowFact(ir)
         flowout2 = flowout2 + MassFlux(i)

         if( SolveU ) U(Ncel+ib) = U(Ncel+ib) * fact
         if( SolveV ) V(Ncel+ib) = V(Ncel+ib) * fact
         if( SolveW ) W(Ncel+ib) = W(Ncel+ib) * fact

         ip = Face(i)%cell1
         su(ip) = su(ip) - MassFlux(i)

         !write(IOdef,*)'outlet:',ip,U(Ncel+ib),V(Ncel+ib),W(Ncel+ib),MassFlux(i)

       endif
     end do
     if( Debug > 2 )write(IOdef,*)'Adapted mass flow rate out:',flowout2

   endif

   deallocate(FlowRegion,stat=istat)
   call TrackMemory(istat,-Nreg-1,'FlowRegion array deallocated')
   deallocate(FlowFact,stat=istat)
   call TrackMemory(istat,-Nreg-1,'FlowFact array deallocated')

   call watch_leave('FluxMass')

   if( Debug > 3 )write(IOdef,*)'=== FluxMass'

end subroutine FluxMass
subroutine FluxMass2
!========================================================================
!
! this routine should be consistent with FluxMass
! internal faces only
!
   use constants
   use geometry
   use variables
   use watches

   real :: facn, facp
   real :: Xface(3), Xnorm(3), Xp(3), Xn(3), Xpn(3) ,delta(3), Xpa(3), Xna(3)
   real :: Xpac(3), Xnac(3), Xnp(3), Xnn(3)

   if( Debug > 3 )write(IOdef,*)'*** FluxMass2'

   call watch_enter('FluxMass2')

   faceloop: do i=1,Nfac
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )then

       facn = Face(i)%lambda
       facp = 1.0 - facn

       Xface = Face(i)%x

       !
       ! get auxiliary nodes etc  (eq. 8.53)
       ! preset in subroutine ReadGeometry
       !
       Xpac = Face(i)%Xpac     
       Xnac = Face(i)%Xnac     

       Xn = Xnac - Cell(in)%x
       Xp = Xpac - Cell(ip)%x

       !
       ! the coeff. is the same as in FluxMass
       !
       fact = Rface(i,1)

       !          ,
       !   (grad p ).(r ,-r )
       !          E    E   E
       !
       dpx = dPdX(in,1)*Xn(1) - dpdx(ip,1)*Xp(1)
       dpy = dPdX(in,2)*Xn(2) - dpdx(ip,2)*Xp(2)
       dpz = dPdX(in,3)*Xn(3) - dpdx(ip,3)*Xp(3)

       ! hint from Max... under relax PP
       !
       fc  = fact*( dpx + dpy + dpz )*URF(VarPP)

       MassFlux(i) = MassFlux(i) + fc

       Su(ip) = Su(ip) - fc
       Su(in) = Su(in) + fc

     endif
   end do faceloop

   call watch_leave('FluxMass2')

   if( Debug > 3 )write(IOdef,*)'=== FluxMass2'

end subroutine FluxMass2
subroutine UpdateP(Phi,dPhiDx)
!========================================================================

   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX

   real :: ds(3)

   if( Debug > 3 )write(IOdef,*)'*** UpdateP'

  !call watch_enter('UpdateP')

   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == ROutlet )then
       !
       ! the pressure gradient in the outlet should be 0
       !
       Phi(Ncel+ib) = Phi(ip)

     else
       !
       ! modification by Harry, discard the gradient
       !
       ds = Face(i)%x - Cell(ip)%x
       Phi(Ncel+ib) = Phi(ip) + dot_product( dPhiDx(ip,:) , ds )

       !Phi(Ncel+ib) = Phi(ip)

     endif

   end do

  !call watch_leave('UpdateP')

   if( Debug > 3 )write(IOdef,*)'=== UpdateP'

end subroutine UpdateP
subroutine Update_P_at_boundaries(ivar,Phi)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(3) :: GradPhi, Xp, dX, ds

   if( Debug > 3 ) write(IOdef,*)'*** Update_P_at_boundaries',variable(ivar)

   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     Xp = Cell(ip)%x
     GradPhi = 0.0
     icnt  = 0
     icnt1 = 0
     icnt2 = 0
     icnt3 = 0

     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2

       if( ipn > 0 )then
         ! internal face
         if( ipp == ip )then
           dPhi = Phi(ipn) - Phi(ipp)
           dX   = Cell(ipn)%x - Cell(ipp)%x
         else
           dPhi = Phi(ipp) - Phi(ipn)
           dX   = Cell(ipp)%x - Cell(ipn)%x
         endif

         icnt = icnt + 1
         if( dX(1) /= 0.0 )then
           icnt1      = icnt1 + 1
           GradPhi(1) = GradPhi(1) + dPhi/dX(1)
         endif
         if( dX(2) /= 0.0 )then
           icnt2      = icnt2 + 1
           GradPhi(2) = GradPhi(2) + dPhi/dX(2)
         endif
         if( dX(3) /= 0.0 )then
           icnt3      = icnt3 + 1
           GradPhi(3) = GradPhi(3) + dPhi/dX(3)
         endif

       endif
     end do

     if( icnt > 0 ) GradPhi = GradPhi / float(icnt)

     ds = Face(i)%x - Xp

     !Phi(Ncel+ib) = Phi(ip) + dot_product( ds , GradPhi )

     Phi(Ncel+ib) = Phi(ip)

   end do

   if( Debug > 3 )write(IOdef,*)'=== Update_P_at_boundaries'
end subroutine Update_P_at_boundaries
subroutine UpdateBC
!========================================================================

   use constants
   use geometry
   use variables
   use watches

   real :: ds(3)

   if( Debug > 3 )write(IOdef,*)'*** UpdateBC (update symm. bc''s)'

  !call watch_enter('UpdateBC')

   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     if( it == RSymp )then

       ds = Face(i)%x - Cell(ip)%x

       if( SolveU ) U(Ncel+ib) = U(ip) + dot_product( dUdX(ip,:) , ds )
       if( SolveV ) V(Ncel+ib) = V(ip) + dot_product( dVdX(ip,:) , ds )
       if( SolveW ) W(Ncel+ib) = W(ip) + dot_product( dWdX(ip,:) , ds )
       ! (P already done)

     endif
   end do

  !call watch_leave('UpdateBC')

   if( Debug > 3 ) write(IOdef,*)'=== UpdateBC'
end subroutine UpdateBC
subroutine ShowMass(btype)
!========================================================================

   use constants
   use geometry
   use variables

   integer btype

   sumtot = 0.0
   do ib=1,Nbnd
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ
     if( it == btype )then
       i   = Bnd(ib)%face
       ic  = Face(i)%cell1
       sum = 0.0
       do j=1,NFaces(ic)
         k = CFace(ic,j)
         ip = Face(k)%cell1
         in = Face(k)%cell2
         if( ip == ic )then
           sum = sum - MassFlux(k)
         else
           sum = sum + MassFlux(k)
         endif
       end do
       sumtot = sumtot + sum
       !write(IOdef,*)'Cell:',ip,sum,sumtot
     endif
   end do

end subroutine ShowMass
subroutine ShowCells
!========================================================================

   use constants
   use geometry
   use variables

   IOloc = IOdef
   if( PrintCellVar .or. PrintWallVar )then
     if( PrintFile(1:1) /= '-' )then
       call openfile(IOprt,PrintFile,'.dat','FORMATTED','SEQUENTIAL', &
                     'UNKNOWN',Debug)
       IOloc = IOprt
       write(IOdef,*) 'Printing to ',PrintFile(1:lens(PrintFile)),'.dat'
     endif
   endif

   if( PrintCellVar )then
     if( PrintCellVarUser )then

       call UserPrintCell(IOloc)

     else
       write(IOloc,'(/,'' Cell data'')')
       write(IOloc,'(A,A)')'    Cell      U          V          W          P    ',&
                       '      k       epsilon    T_rel      Vis_t'

       do i=IPrintCellVarStart,min(IPrintCellVarEnd,NCel),IPrintCellVarInc

         TEtmp = 0.0
         EDtmp = 0.0
         Ttmp  = Tref

         if( SolveTurbEnergy )TEtmp = TE(i)
         if( SolveTurbDiss )  EDtmp = ED(i)
         if( SolveEnthalpy )  Ttmp  = T(i)

         write(IOloc,'(i8,8(1x,1pe10.3))') i,u(i),v(i),w(i),p(i), &
                                     TEtmp,EDtmp,Ttmp-Tref,VisEff(i)-Vislam

       end do
     endif
   endif

   if( PrintWallVar )then
     if( PrintWallVarUser )then

       call UserPrintWall(IOloc)

     else
       write(IOloc,'(/,'' Wall data'')')
       write(IOloc,'(A,A)')'    Cell   Reg  X-shear    Y-shear    Z-shear',&
        '      Y+         U+         dn        T_w        h_w        q_w'

       do ib=IPrintWallVarStart,min(IPrintWallVarEnd,NBnd),IPrintWallVarInc

         i  = Bnd(ib)%face
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ

         ip = Face(i)%cell1
         in = Face(i)%cell2
         if( in /= 0 )write(IOdef,*)'+ internal error: corrupt bc ',ib,i,ip,in

         if( it == RWall )then

           Ttmp = Tref
           if( SolveEnthalpy )  Ttmp  = T(Ncel+ib)

           write(IOloc,'(i8,1x,i4,9(1x,1pe10.3))') ip,ir,bnd(ib)%shear, bnd(ib)%yplus,&
                      bnd(ib)%uplus,bnd(ib)%distance,Ttmp-Tref,bnd(ib)%h,bnd(ib)%q

         endif
       end do
     endif
   endif
   !write(IOdef,'(/,'' Face data'')')
   !write(IOdef,'(A,A)')'    Face    ip       in   massflux'

   !do i=1,Nfac
   !  write(IOdef,'(i8,2(1x,i4),9(1x,1pe10.3))') i, &
   !           Face(i)%cell1,Face(i)%cell2, MassFlux(i)
   !end do

   if( IOprt /= IOdef ) close(IOprt)

end subroutine ShowCells
subroutine Update_P_at_boundaries2(ivar,Phi)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(3)   :: GradPhi, Xp, dX, ds

   real, dimension(3,3) :: A
   real, dimension(3)   :: RHS_A

   real, dimension(2,2) :: B
   real, dimension(2)   :: RHS_B

   real :: p1, p2, p3, p4, p5, p6

   integer, dimension(3):: IPIV  ! LAPACK

   logical :: UseX, UseY, UseZ

   if( Debug > -1 )write(IOdef,*)'*** Update_P_at_boundaries2: ',variable(ivar)

   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     Xp = Cell(ip)%x
     GradPhi = 0.0
     A       = 0.0       ! reset A
     B       = 0.0       ! reset B
     RHS_A   = 0.0
     RHS_B   = 0.0

     icnt = 0
     UseX = .false.
     UseY = .false.
     UseZ = .false.

     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2

       if( ipn > 0 )then
         ! internal faces only
         icnt = icnt + 1

         if( ipp == ip )then
           dPhi = Phi(ipn) - Phi(ipp)
           dX   = Cell(ipn)%x - Cell(ipp)%x
         else
           dPhi = Phi(ipp) - Phi(ipn)
           dX   = Cell(ipp)%x - Cell(ipn)%x
         endif

         A(1,1) = A(1,1) + dX(1)*dX(1)
         A(2,1) = A(2,1) + dX(1)*dX(2)
         A(3,1) = A(3,1) + dX(1)*dX(3)

         A(1,2) = A(2,1)
         A(2,2) = A(2,2) + dX(2)*dX(2)
         A(3,2) = A(3,2) + dX(2)*dX(3)

         A(1,3) = A(3,1)
         A(2,3) = A(3,2)
         A(3,3) = A(3,3) + dX(3)*dX(3)

         RHS_A(1) = RHS_A(1) + dX(1)*dPhi
         RHS_A(2) = RHS_A(2) + dX(2)*dPhi
         RHS_A(3) = RHS_A(3) + dX(3)*dPhi

         !if( ib == 16 )then
         !  write(IOdef,*)'dphi,dx:',dphi,dx
         !  write(IOdef,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
         !  write(IOdef,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
         !  write(IOdef,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
         !endif

       endif
     end do

     if( icnt <= 0 )then
       write(IOdef,*)'Foutje!'
     !else
     !  write(IOdef,*)'Aantal relaties:',ip,j,icnt
     endif

     if( UseLapack )then

       call SGESV( 3, 1, A, 3, IPIV, RHS_A, 3, INFO )

       if( info /= 0 ) write(IOdef,*)'Lapack (sgesv) info(1):',INFO

       !if( info /= 0 .and. ib ==16)then
       !  write(IOdef,*)'Lapack (sgesv) info(2):',INFO
       !  write(IOdef,*)'ib,ip,ir,it:',ib,ip,ir,it
       !  write(IOdef,*)'Xp:',Xp
       !  write(IOdef,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
       !  write(IOdef,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
       !  write(IOdef,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
       !endif

       GradPhi(1) = RHS_A(1)
       GradPhi(2) = RHS_A(2)
       GradPhi(3) = RHS_A(3)

     else

       irang = 0
       if( A(1,1) > Small )then  ! altijd >= 0
         UseX = .true.
         irang = irang + 1
       else
         GradPhi(1) = 0.0
       endif
       if( A(2,2) > Small )then  ! altijd >= 0
         UseY = .true.
         irang = irang + 1
       else
         GradPhi(2) = 0.0
       endif
       if( A(3,3) > Small )then  ! altijd >= 0
         UseZ = .true.
         irang = irang + 1
       else
         GradPhi(3) = 0.0
       endif

       if( irang == 1 )then
         !
         ! too simple, do nothing
         !

         GradPhi = 0.0

       else if( irang == 2 )then
         !
         ! only two components can be calculated
         ! assume the third to be zero
         !
         if( UseX .and. UseY )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,2)
           B(2,1) = A(2,1)
           B(2,2) = A(2,2)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(2)
         else if( UseX .and. UseZ )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,3)
           B(2,1) = A(3,1)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(3)
         else if( UseY .and. UseZ )then
           B(1,1) = A(2,2)
           B(1,2) = A(2,3)
           B(2,1) = A(3,2)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(2)
           RHS_B(2) = RHS_A(3)
         endif
         !
         ! B inverteren
         !
         if( B(1,1) /= 0.0 )then
           p1 = - B(2,1)/B(1,1)
           B(2,1) = B(2,1) + p1 * B(1,1)
           B(2,2) = B(2,2) + p1 * B(1,2)

           RHS_B(2) = RHS_B(2) + p1 * RHS_B(1)

         else
           write(IOdef,1) 1
         endif
         if( B(2,2) /= 0.0 )then
           p2 = - B(1,2)/B(2,2)
           B(1,1) = B(1,1) + p2 * B(2,1)
           B(1,2) = B(1,2) + p2 * B(2,2)

           RHS_B(1) = RHS_B(1) + p2 * RHS_B(2)

         else
           write(IOdef,1) 2
         endif
         !
         ! normeren
         !
         if( B(1,1) /= 0.0 )then
           p3 = 1.0/B(1,1)
           B(1,1) = p3 * B(1,1)
           B(1,2) = p3 * B(1,2)

           RHS_B(1) = p3 * RHS_B(1)

         else
           write(IOdef,1) 3,B(1,1)
         endif
         if( B(2,2) /= 0.0 )then
           p4 = 1.0/B(2,2)
           B(2,1) = p4 * B(2,1)
           B(2,2) = p4 * B(2,2)

           RHS_B(2) = p4 * RHS_B(2)

         else
           write(IOdef,1) 4,B(2,2)
         endif
         !
         ! de gradient
         !
         if( UseX .and. UseY )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = RHS_B(2)
           GradPhi(3) = 0.0
         else if( UseX .and. UseZ )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = 0.0
           GradPhi(3) = RHS_B(2)
         else if( UseY .and. UseZ )then
           GradPhi(1) = 0.0
           GradPhi(2) = RHS_B(1)
           GradPhi(3) = RHS_B(2)
         endif

       else if( irang == 3 )then
         !
         ! the case were all three components
         ! can be calculated
         !
         if( A(1,1) /= 0.0 )then
           p1 = - A(2,1)/A(1,1)
           A(2,1)   = A(2,1) + p1 * A(1,1)
           A(2,2)   = A(2,2) + p1 * A(1,2)
           A(2,3)   = A(2,3) + p1 * A(1,3)
           RHS_A(2) = RHS_A(2) + p1 * RHS_A(1)

           p2 = - A(3,1)/A(1,1)
           A(3,1)   = A(3,1) + p2 * A(1,1)
           A(3,2)   = A(3,2) + p2 * A(1,2)
           A(3,3)   = A(3,3) + p2 * A(1,3)
           RHS_A(3) = RHS_A(3) + p2 * RHS_A(1)
         else
           write(IOdef,1) 5
         endif
         if( A(2,2) /= 0.0 )then
           p3 = - A(1,2)/A(2,2)
           A(1,1)   = A(1,1) + p3 * A(2,1)
           A(1,2)   = A(1,2) + p3 * A(2,2)
           A(1,3)   = A(1,3) + p3 * A(2,3)
           RHS_A(1) = RHS_A(1) + p3 * RHS_A(2)

           p4 = - A(3,2)/A(2,2)
           A(3,1)   = A(3,1) + p4 * A(2,1)
           A(3,2)   = A(3,2) + p4 * A(2,2)
           A(3,3)   = A(3,3) + p4 * A(2,3)
           RHS_A(3) = RHS_A(3) + p4 * RHS_A(2)
         else
           write(IOdef,1) 6
         endif
         if( A(3,3) /= 0.0 )then
           p5 = - A(1,3)/A(3,3)
           A(1,1)   = A(1,1) + p5 * A(3,1)
           A(1,2)   = A(1,2) + p5 * A(3,2)
           A(1,3)   = A(1,3) + p5 * A(3,3)
           RHS_A(1) = RHS_A(1) + p5 * RHS_A(3)

           p6 = - A(2,3)/A(3,3)
           A(2,1)   = A(2,1) + p6 * A(3,1)
           A(2,2)   = A(2,2) + p6 * A(3,2)
           A(2,3)   = A(2,3) + p6 * A(3,3)
           RHS_A(2) = RHS_A(2) + p6 * RHS_A(3)
         else
           write(IOdef,1) 7
           !write(IOdef,*)'c : ',ip,Nfaces(ip)
           !write(IOdef,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
           !write(IOdef,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
           !write(IOdef,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
           !write(IOdef,*)'Use:',UseX,UseY,UseZ,irang
         endif
         !
         ! normeren
         !
         if( abs(A(1,1)) > Small )then
           p7 = 1.0/A(1,1)
           A(1,1)   = p7 * A(1,1)
           A(1,2)   = p7 * A(1,2)
           A(1,3)   = p7 * A(1,3)
           RHS_A(1) = p7 * RHS_A(1)
         else
           write(IOdef,1) 8
         endif
         if( abs(A(2,2)) > Small )then
           p8 = 1.0/A(2,2)
           A(2,1)   = p8 * A(2,1)
           A(2,2)   = p8 * A(2,2)
           A(2,3)   = p8 * A(2,3)
           RHS_A(2) = p8 * RHS_A(2)
         else
           write(IOdef,1) 9
         endif
         if( abs(A(3,3)) > Small )then
           p9 = 1.0/A(3,3)
           A(3,1)   = p9 * A(3,1)
           A(3,2)   = p9 * A(3,2)
           A(3,3)   = p9 * A(3,3)
           RHS_A(3) = p9 * RHS_A(3)
         else
           write(IOdef,1) 10
           write(IOdef,*)'c : ',ip,Nfaces(ip)
           write(IOdef,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_a(1)
           write(IOdef,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_a(2)
           write(IOdef,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_a(3)
           write(IOdef,*)'p13:',p1,p2,p3
           write(IOdef,*)'p46:',p4,p5,p6
           write(IOdef,*)'p79:',p7,p8,p9
           write(IOdef,*)'Use:',UseX,UseY,UseZ,irang

           RHS_A(3) = 0.0
         endif
         !
         ! de gradient
         !

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

       else
         write(IOdef,*)'+ internal error catch 22'
       endif

     2 continue

       !write(IOdef,*)'achter de rug'
       !write(IOdef,*)'B1:',B(1,1),B(1,2),'|',RHS_B(1)
       !write(IOdef,*)'B2:',B(2,1),B(2,2),'|',RHS_B(2)
       !write(IOdef,*)'---'
       !write(IOdef,*)'Cel:',ip,GradPhi

     endif

     ds = Face(i)%x - Xp

     Phi(Ncel+ib) = Phi(ip) + dot_product( ds , GradPhi )

   end do

 1 format('+ internal error: Update_P_at_boundaries2 ',i2,e14.4)
   if( Debug > 3 )write(IOdef,*)'=== Update_P_at_boundaries2'

end subroutine Update_P_at_boundaries2
subroutine Update_Scalars_at_boundaries2(ivar,Phi)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(3)   :: GradPhi, Xp, dX, ds
   real, dimension(3,3) :: A
   real, dimension(2,2) :: B
   real, dimension(3)   :: RHS_A
   real, dimension(2)   :: RHS_B

   logical :: UseX, UseY, UseZ

   if( Debug > 3 )write(IOdef,*)'*** Update_Scalars_at_boundaries2: ',variable(ivar)

   do ib=1,Nbnd
     i  = Bnd(ib)%face
     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     !
     ! LET OP: alleen voor k en epsilon ook de wanden
     !         bij T wordt de wandtemp. elders gezet
     !         en mag niet overschreven worden.
     !
     if( it == RSymp .or. it == ROutlet )then

       Xp = Cell(ip)%x
       GradPhi = 0.0
       A     = 0.0           ! reset A
       B     = 0.0           ! reset B
       RHS_A = 0.0
       RHS_B = 0.0

       icnt = 0
       UseX = .false.
       UseY = .false.
       UseZ = .false.

       do j=1,Nfaces(ip)
         k  = CFace(ip,j)
         ipp = Face(k)%cell1
         ipn = Face(k)%cell2

         if( ipn > 0 )then
           ! internal faces only
           icnt = icnt + 1

           if( ipp == ip )then
             dPhi = Phi(ipn) - Phi(ipp)
             dX   = Cell(ipn)%x - Cell(ipp)%x
           else
             dPhi = Phi(ipp) - Phi(ipn)
             dX   = Cell(ipp)%x - Cell(ipn)%x
           endif

           A(1,1) = A(1,1) + dX(1)*dX(1)
           A(2,1) = A(2,1) + dX(1)*dX(2)
           A(3,1) = A(3,1) + dX(1)*dX(3)

           A(1,2) = A(2,1)
           A(2,2) = A(2,2) + dX(2)*dX(2)
           A(3,2) = A(3,2) + dX(2)*dX(3)

           A(1,3) = A(3,1)
           A(2,3) = A(3,2)
           A(3,3) = A(3,3) + dX(3)*dX(3)

           RHS_A(1) = RHS_A(1) + dX(1)*dPhi
           RHS_A(2) = RHS_A(2) + dX(2)*dPhi
           RHS_A(3) = RHS_A(3) + dX(3)*dPhi

         endif
       end do

       if( icnt <= 0 )then
         write(IOdef,*)'Foutje!'
       !else
       !  write(IOdef,*)'Aantal relaties:',ip,j,icnt
       endif


       irang = 0
       if( A(1,1) > Small )then  ! altijd >= 0
         UseX = .true.
         irang = irang + 1
       else
         A(1,1) = 0.0
         GradPhi(1) = 0.0
       endif
       if( A(2,2) > Small )then  ! altijd >= 0
         UseY = .true.
         irang = irang + 1
       else
         A(2,2) = 0.0
         GradPhi(2) = 0.0
       endif
       if( A(3,3) > Small )then  ! altijd >= 0
         UseZ = .true.
         irang = irang + 1
       else
         A(3,3) = 0.0
         GradPhi(3) = 0.0
       endif

       if( irang == 1 )then
         !
         ! too simple, do nothing
         !

         GradPhi = 0.0

       else if( irang == 2 )then
         !
         ! only two components can be calculated
         ! assume the third to be zero
         !
         if( UseX .and. UseY )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,2)
           B(2,1) = A(2,1)
           B(2,2) = A(2,2)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(2)
         else if( UseX .and. UseZ )then
           B(1,1) = A(1,1)
           B(1,2) = A(1,3)
           B(2,1) = A(3,1)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(1)
           RHS_B(2) = RHS_A(3)
         else if( UseY .and. UseZ )then
           B(1,1) = A(2,2)
           B(1,2) = A(2,3)
           B(2,1) = A(3,2)
           B(2,2) = A(3,3)

           RHS_B(1) = RHS_A(2)
           RHS_B(2) = RHS_A(3)
         endif
         !
         ! B inverteren
         !
         if( B(1,1) /= 0.0 )then
           p1 = - B(2,1)/B(1,1)
           B(2,1) = B(2,1) + p1 * B(1,1)
           B(2,2) = B(2,2) + p1 * B(1,2)

           RHS_B(2) = RHS_B(2) + p1 * RHS_B(1)

         else
           write(IOdef,1) 1
         endif
         if( B(2,2) /= 0.0 )then
           p2 = - B(1,2)/B(2,2)
           B(1,1) = B(1,1) + p2 * B(1,2)
           B(1,2) = B(1,2) + p2 * B(2,2)

           RHS_B(1) = RHS_B(1) + p2 * RHS_B(2)

         else
           write(IOdef,1) 2
         endif
         !
         ! normeren
         !
         if( abs(B(1,1)) > Small )then
           p3 = 1.0/B(1,1)
           B(1,1) = p3 * B(1,1)
           B(1,2) = p3 * B(1,2)

           RHS_B(1) = p3 * RHS_B(1)

         else
           write(IOdef,1) 3
         endif
         if( abs(B(2,2)) > Small )then
           p4 = 1.0/B(2,2)
           B(2,1) = p4 * B(2,1)
           B(2,2) = p4 * B(2,2)

           RHS_B(2) = p4 * RHS_B(2)

         else
           write(IOdef,1) 4
         endif
         !
         ! de gradient
         !
         if( UseX .and. UseY )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = RHS_B(2)
           GradPhi(3) = 0.0
         else if( UseX .and. UseZ )then
           GradPhi(1) = RHS_B(1)
           GradPhi(2) = 0.0
           GradPhi(3) = RHS_B(2)
         else if( UseY .and. UseZ )then
           GradPhi(1) = 0.0
           GradPhi(2) = RHS_B(1)
           GradPhi(3) = RHS_B(2)
         endif

       else if( irang == 3 )then
         !
         ! the case were all three components
         ! can be calculated
         !
         if( A(1,1) /= 0.0 )then
           p1 = - A(2,1)/A(1,1)
           A(2,1) = A(2,1) + p1 * A(1,1)
           A(2,2) = A(2,2) + p1 * A(1,2)
           A(2,3) = A(2,3) + p1 * A(1,3)
           RHS_A(2) = RHS_A(2) + p1 * RHS_A(1)

           p2 = - A(3,1)/A(1,1)
           A(3,1) = A(3,1) + p2 * A(1,1)
           A(3,2) = A(3,2) + p2 * A(1,2)
           A(3,3) = A(3,3) + p2 * A(1,3)
           RHS_A(3) = RHS_A(3) + p2 * RHS_A(1)
         else
           write(IOdef,1) 5
         endif
         if( A(2,2) /= 0.0 )then
           p3 = - A(1,2)/A(2,2)
           A(1,1) = A(1,1) + p3 * A(2,1)
           A(1,2) = A(1,2) + p3 * A(2,2)
           A(1,3) = A(1,3) + p3 * A(2,3)
           RHS_A(1) = RHS_A(1) + p3 * RHS_A(2)

           p4 = - A(3,2)/A(2,2)
           A(3,1) = A(3,1) + p4 * A(2,1)
           A(3,2) = A(3,2) + p4 * A(2,2)
           A(3,3) = A(3,3) + p4 * A(2,3)
           RHS_A(3) = RHS_A(3) + p4 * RHS_A(2)
         else
           write(IOdef,1) 6
         endif
         if( A(3,3) /= 0.0 )then
           p5 = - A(1,3)/A(3,3)
           A(1,1) = A(1,1) + p5 * A(3,1)
           A(1,2) = A(1,2) + p5 * A(3,2)
           A(1,3) = A(1,3) + p5 * A(3,3)
           RHS_A(1) = RHS_A(1) + p5 * RHS_A(3)

           p6 = - A(2,3)/A(3,3)
           A(2,1) = A(2,1) + p6 * A(3,1)
           A(2,2) = A(2,2) + p6 * A(3,2)
           A(2,3) = A(2,3) + p6 * A(3,3)
           RHS_A(2) = RHS_A(2) + p6 * RHS_A(3)
         else
           write(IOdef,1) 7
           !write(IOdef,*)'c : ',ip,Nfaces(ip)
           !write(IOdef,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
           !write(IOdef,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
           !write(IOdef,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
           !write(IOdef,*)'Use:',UseX,UseY,UseZ,irang
         endif
         !
         ! normeren
         !
         if( abs(A(1,1)) > Small )then
           p7 = 1.0/A(1,1)
           A(1,1) = p7 * A(1,1)
           A(1,2) = p7 * A(1,2)
           A(1,3) = p7 * A(1,3)
           RHS_A(1) = p7 * RHS_A(1)
         else
           write(IOdef,1) 8
         endif
         if( abs(A(2,2)) > Small )then
           p8 = 1.0/A(2,2)
           A(2,1) = p8 * A(2,1)
           A(2,2) = p8 * A(2,2)
           A(2,3) = p8 * A(2,3)
           RHS_A(2) = p8 * RHS_A(2)
         else
           write(IOdef,1) 9
         endif
         if( abs(A(3,3)) > Small )then
           p9 = 1.0/A(3,3)
           A(3,1) = p9 * A(3,1)
           A(3,2) = p9 * A(3,2)
           A(3,3) = p9 * A(3,3)
           RHS_A(3) = p9 * RHS_A(3)
         else
           write(IOdef,1) 10
           write(IOdef,*)'c : ',ip,Nfaces(ip)
           write(IOdef,*)'a1: ',a(1,1),a(1,2),a(1,3),'|',rhs_A(1)
           write(IOdef,*)'a2: ',a(2,1),a(2,2),a(2,3),'|',rhs_A(2)
           write(IOdef,*)'a3: ',a(3,1),a(3,2),a(3,3),'|',rhs_A(3)
           write(IOdef,*)'Use:',UseX,UseY,UseZ,irang

           rhs_A(3) = 0.0
         endif
         !
         ! de gradient
         !

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

       else
         write(IOdef,*)'+ internal error catch 22'
       endif

       !write(IOdef,*)'achter de rug'
       !write(IOdef,*)'B1:',B(1,1),B(1,2),'|',RHS_B(1)
       !write(IOdef,*)'B2:',B(2,1),B(2,2),'|',RHS_B(2)
       !write(IOdef,*)'---'
       !write(IOdef,*)'Cel:',ip,GradPhi

       ds = Face(i)%x - Xp

       Phi(Ncel+ib) = Phi(ip) + dot_product( ds , GradPhi )

     endif
   end do

 1 format('+ internal error: Update_Scalars_at_boundaries2 ',i2,e14.4)
   if( debug > 3 )write(IOdef,*)'=== Update_Scalars_at_boundaries2'

end subroutine Update_Scalars_at_boundaries2
subroutine Set_Normalisation_Factors
!========================================================================

   use constants
   use geometry
   use variables
   use scalars

   if( Debug > 3 ) write(IOdef,*)'*** Set_Normalisation_Factors'

   areain  = 0.0
   fluxin  = 0.0

   sum0 = 0.0
   sum1 = 0.0
   sum2 = 0.0
   sum3 = 0.0

   areaout = 0.0

   icntin  = 0
   icntout = 0

   do ib=1,Nbnd

     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     if( it == RInlet )then
       icntin = icntin + 1
       areain = areain + Face(i)%area

       if( Reg(ir)%user )then

         ip = Face(i)%cell1

         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den

         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area

         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)

         fluxin =  Dtmp * ( Utmp*Face(i)%n(1)   &
                          + Vtmp*Face(i)%n(2)   &
                          + Wtmp*Face(i)%n(3) )
         sum0   = sum0 + fluxin

         sum1   = sum1 + &
                  fluxin *sqrt( Utmp*Utmp + Vtmp*Vtmp + Wtmp*Wtmp )

         sum2   = sum2 + &
                  fluxin * ( Utmp*Utmp + Vtmp*Vtmp + Wtmp*Wtmp )

         if( SolveEnthalpy )then

           sum3 = sum3 + &
                     Cp(Ncel+ib) * Ttmp * Dtmp * &
                       ( Utmp*Face(i)%n(1)       &
                       + Vtmp*Face(i)%n(2)       &
                       + Wtmp*Face(i)%n(3) )

         endif
       else

         fluxin = Reg(ir)%den * dot_product( Reg(ir)%uvw , Face(i)%n )

         sum0   = sum0 + fluxin

         sum1   = sum1 + &
                  fluxin *sqrt( dot_product( Reg(ir)%uvw , Reg(ir)%uvw ) )

         sum2   = sum2 + &
                  fluxin * dot_product( Reg(ir)%uvw , Reg(ir)%uvw )

         if( SolveEnthalpy )then

           sum3 = sum3 + &
                     Cp(Ncel+ib) * Reg(ir)%T * Reg(ir)%den * &
                     dot_product( Reg(ir)%uvw , Face(i)%n )

         endif
       endif

     else if( it == ROutlet )then
       icntout = icntout + 1
       areaout = areaout + Face(i)%area
     endif
   end do

   write(IOdbg,*)'Residuals with:'
   write(IOdbg,*)'   fluxin: ',fluxin
   write(IOdbg,*)'   areain: ',areain
   write(IOdbg,*)'  areaout: ',areaout
   write(IOdbg,*)'     sum0: ',sum0
   write(IOdbg,*)'     sum1: ',sum1
   write(IOdbg,*)'     sum2: ',sum2
   if( SolveEnthalpy ) write(IOdbg,*)'  sum3: ',sum3

   ResiNorm(:) = 0.0

   if( icntin /= 0 )then

     if( SolveU ) ResiNorm(VarU) = 1.0 / (sum1 + Small)
     if( SolveV ) ResiNorm(VarV) = 1.0 / (sum1 + Small)
     if( SolveW ) ResiNorm(VarW) = 1.0 / (sum1 + Small)

     if( SolveP ) ResiNorm(VarP) = 1.0 / (sum0 + Small)

     if( SolveTurbEnergy ) ResiNorm(VarTE) = 1.0 / (sum2 + Small)

     if( SolveTurbDiss )then
        Uin = fluxin / areain              ! mean velocity
        dL  = 1.0                          ! some length scale
        ResiNorm(VarED) = 1.0 / (sum2 + Small) * dL / Uin
     endif

     if( SolveEnthalpy ) ResiNorm(VarT) = 1.0 / (sum3 + Small)

     if( SolveScalars ) ResiNorm(VarSC) = 1.0   ! not used!              !<= MOET BETER!

     ResiNorm(:) = abs( ResiNorm )

   else

     ResiNorm(:) = 1.0

     if( SolveEnthalpy .and. sum3 > Small )then
       ResiNorm(VarT) = 1.0/sum3
     endif

   endif

   write(IOdbg,*)'Using ResiNorm 1-3:',ResiNorm(1:3)
   write(IOdbg,*)'      ResiNorm 4-6:',ResiNorm(4:6)
   write(IOdbg,*)'      ResiNorm 7-8:',ResiNorm(7:8)

   if( SolveScalars ) call Set_Scalar_Norm()


   if( Debug > 3 ) write(IOdef,*)'=== Set_Normalisation_Factors'

end subroutine Set_Normalisation_Factors
subroutine CalculateViscosity
!========================================================================
!                                     2
!  Eddy viscosity: VisTurb = rho Cmu k / epsilon (eq. 9.43)
!
!  Effective viscosity: VisEff = VisLam + VisTurb
!
!========================================================================

   use constants
   use geometry
   use variables
   use watches

   logical warning

   if( Debug > 3 ) write(IOdef,*)'*** CalculateViscosity'

   call watch_enter('CalculateViscosity')

   warning = .false.
   Cmu    = TMCmu
   Cmu25  = Cmu**(0.25)
   VisURF = 1.0   ! 0.99

   do ip=1,Ncel
     VisOld = VisEff(ip)
     if( ED(ip) > Small )then
       VisNew = VisLam + Cmu*Den(ip)*TE(ip)**2/ED(ip)
     else
       VisNew = VisLam
     endif
     VisEff(ip) = VisOld + VisURF*( VisNew - VisOld )
   end do

   !
   ! that was the easy part... now modify the viscosities
   ! at the boundaries
   !

   do ir=1,Nreg

     if( Reg(ir)%typ == RWall .and. Reg(ir)%std )then
       !
       ! standard smooth wall: u+ = 1/k ln(Ey+)
       !
       yplus = 11.
       elog  = Reg(ir)%elog

       i = 0
       do while( abs( yplus - yplustmp ) > 0.001 .and. i < 50 )
         yplustmp = yplus
         yplus    = log( elog * yplustmp )/kappa
         i = i + 1
       end do

       Reg(ir)%ylog = yplus

       write(IOdbg,'(1x,''Region'',i2,'' y+m ='',f5.1,i3)') ir,yplus,i

     else if( Reg(ir)%typ == RWall .and. .not. Reg(ir)%std )then
       !
       ! rough wall:  u+ = B_rough + 1/k ln( (y-do)/ks )
       ! y+m = 0.0 by definition
       !
       Reg(ir)%ylog = 0.0              

       write(IOdbg,'(1x,''Region'',i2,'' y+m ='',f5.1,'' rough'')') ir,yplus

     endif
   end do

   YplusMin = Large
   YplusMax = Small
   SplusMin = Large
   SplusMax = Small
   SdistMin = Large
   SdistMax = Small
   
   ii1 = 0
   ii2 = 0
   ii3 = 0
   
   do ib=1,Nbnd

     i  = Bnd(ib)%face
     ir = Bnd(ib)%rid
     it = Reg(ir)%typ

     ip = Face(i)%cell1
     in = Face(i)%cell2
     if( in /= 0 )write(IOdef,*)'+ internal error: corrupt bc ',ib,i,ip,in

     if( it == RInlet )then

       if( Reg(ir)%user )then

         ip = Face(i)%cell1

         Utmp  = Reg(ir)%uvw(1)
         Vtmp  = Reg(ir)%uvw(2)
         Wtmp  = Reg(ir)%uvw(3)
         TEtmp = Reg(ir)%k
         EDtmp = Reg(ir)%e
         Ttmp  = Reg(ir)%T
         Dtmp  = Reg(ir)%den

         Xtmp  = Face(i)%x(1)
         Ytmp  = Face(i)%x(2)
         Ztmp  = Face(i)%x(3)
         Atmp  = Face(i)%area

         call UserInlet(ir,ib,ip,Xtmp,Ytmp,Ztmp,Atmp,         &
                        Pref,Tref,Iter,Time,                  &
                        Utmp,Vtmp,Wtmp,TEtmp,EDtmp,Ttmp,Dtmp, &
                        Reg(ir)%name1,Reg(ir)%name2)

         VisEff(Ncel+ib)  = VisLam + &
           Dtmp * Cmu * TEtmp**2 / (EDtmp + Small)

       else

         VisEff(Ncel+ib)  = VisLam + &
           Reg(ir)%den * Cmu * Reg(ir)%k**2 / (Reg(ir)%e + Small)

       endif
     else if( it == ROutlet )then
       !
       ! simpel
       !
       VisEff(Ncel+ib)  = VisEff(ip)

     else if( it == RSymp )then

       VisEff(Ncel+ib)  = VisEff(ip)                ! uitvoeriger zoals bij Rwall?

     else if( it == RWall )then

       dist = Bnd(ib)%distance                      ! wall distance
       turb = TE(ip)                                ! k of wall cel

       if( turb < Zero )then
         turb = Zero                                ! don't let k become negative
         write(IOdef,*)'+ CalculateViscosity: Found negative k, set to zero.'
       endif

       Utau  = Cmu25*sqrt( turb )                   ! schuifspanningssnelheid (9.48)

       Yplus = dist * Utau * Den(ip) /        &     ! y+ (9.47)
                               (VisLam + Small)
       Bnd(ib)%yplus = Yplus

       if( Reg(ir)%std  )then      !<<===== LET OPPPPP!!!!
         !
         ! standard wall function based on:
         !
         !   u+ = 1/kappa ln( E y+ )
         !
         ! test if y+ < Y+m; then in viscous layer => u+ = y+
         !
         ! the product of E y+ may not become less or equal than 1
         ! and the natural logarithm becomes negative
         ! (which means that the central node submerges into the
         !  sandrougness). Use a slightly higher value of say 1.1
         !
         if( Yplus < Reg(ir)%ylog  )then

           Uplus = Yplus

           if( Yplus < Small )then
             Yplus = 0.0
             Uplus = Small
           endif
         else
           tmp = Reg(ir)%elog * Bnd(ib)%yplus
           if( tmp < 1.1 )then
             write(IOdef,*)'+ CalculateViscosity: Warning: E y+ < 1.1 ',tmp
             tmp = 1.1
           endif
           !
           !  u+ = 1/kappa ln( E y+ )
           !

           Uplus = log( tmp ) / kappa

         endif

         Bnd(ib)%uplus = Uplus

         YplusMin = min( YplusMin, Yplus )
         YplusMax = max( YplusMax, Yplus )
         !
         ! adhere to the workflow of sub. FluxUVW in which
         !
         !       Tau_(wall) = VisEff * U_(tan) / n
         !
         ! now Tau_(wall) is known and therefore:
         !
         !     VisEff = Tau_(wall) / ( U_(tan) / n )
         !
         !                        2
         !            = rho * Utau / ( U_(tan) / n )
         !
         !            = ( rho * Utau * n / lamvisc )*
         !                               ( lamvisc / ( U_(tan) / Utau )
         !
         !            = Y+ * lamvisc / U+
         !
         !            = Y+/U+ * lamvisc     (ie a factor for the viscosity)
         !
         VisEff(Ncel+ib) = &                       ! a viscosity lower than
               max( 1.0 , Yplus/Uplus ) * VisLam   ! VisLam is quite unlikely

       else
         !
         ! rough wall: based on u+ = 1/kappa ln(y/s) + B_rough + Delta B
         !
         ! for a very rough wall Brough = 8.5 (see Schlichting) and a
         ! sandroughness s. by definition is y+m = 0 (no viscous sublayer)
         ! 
         ! with z0 = s/exp(kappa Brough) = s/exp(0.4 8.5) = s/30.
         !
         s0 = 30.0 * Reg(ir)%z0
         Cs = Reg(ir)%cs
         ! 
         !                             y+ - d+
         ! more general:  u+ = 1/k ln( ------- ) + Br
         !                               r+

         d0 = 0.0                                  ! displacement thickness

         Dplus  = d0 / dist * Yplus
         Splus  = s0 / dist * Yplus
                  
         if( Splus < 2.25 )then
         
           ii1 = ii1 + 1
           DeltaB = 1.0                            ! smooth
         
         elseif( Splus <= 90 )then
         
           ii2 = ii2 + 1
           tmp1   = (Splus-2.25)/87.75 + Cs*Splus
           tmp2   = sin(0.4258*(log(Splus)-0.811)) ! transition
           DeltaB = tmp1**tmp2                        
         
         else
         
           ii3 = ii3 + 1
           DeltaB = 1.0 + Cs*Splus                 ! fully rough
         
         endif
         
        !Uplus = log(Reg(ir)%elog*(Bnd(ib)%yplus-Dplus)/DeltaB)/kappa
         Uplus = log(Reg(ir)%elog*Bnd(ib)%yplus/DeltaB)/kappa
         
         if( Uplus < 1.0 )then
           write(IOdef,*)'+ CalculateViscosity: Warning: u+ < 1.0 ',Uplus
           Uplus = 1.0
         endif

         u11 = log(9.0*11.3/DeltaB)/kappa
         if( Yplus <= 11.3 )then
           Uplus = u11/11.3*Yplus
         endif

         Bnd(ib)%uplus = Uplus

         if( Face(i)%x(1) > 1190. ) write(*,'(A,8(1x,f12.4),A,4(1x,f7.2))') &
                 'y+,u+,k+:',Yplus,Uplus,Splus,CS

         VisEff(Ncel+ib) = &
               max( 1.0 , Yplus/Uplus ) * VisLam

         SplusMin = min( SplusMin, Splus  )
         SplusMax = max( SplusMax, Splus  )
         SdistMin = min( SdistMin, S/dist )
         SdistMax = max( SdistMax, S/dist )
       endif
     else
       write(IOdef,*)'+ CalculateViscosity: Unknown boundary condition'
       write(IOdef,*)'i,ir,it: ',i,ir,it
     endif

   end do
   if( YplusMax > 0.001 ) write(IOdbg,*) YplusMin,'< y+  <',YplusMax
   if( SplusMax > 0.001 ) write(IOdbg,*) SplusMin,'< k+  <',SplusMax
   if( SdistMax > 0.001 ) write(IOdbg,*) SdistMin,'< s/d <',SdistMax
   if( SplusMax > 0.001 ) write(IOdbg,*) 'Smooth, trans, rough: ',ii1,ii2,ii3

   if( SdistMin < 1000. ) &
     write(IOdbg,*)'Warning: Minimum rougness s/d > 1.0 everywhere'
   if( SdistMax > 1.0 ) &
     write(IOdbg,*)'Warning: Maximum rougness s/d > 1.0 '

   !
   ! at the start of a turb. simulation k and epsilon may produce
   ! curious results and as a result viseff may increase to unrealistically
   ! high values. this is a kind of limiter (6 orders is sufficient?)
   !
   ! was       1000000
   !
   FactVisEff =  10000000.0
   VisEffMax = maxval( VisEff )
   if( VisEffMax > FactVisEff * VisLam ) Flags(IFlagVis) = 'v'

  !do ip=1,Ncel 
   do ip=1,Ncel+Nbnd
     VisEff(ip) = min( VisEff(ip) , FactVisEff * VisLam )
   end do

   call watch_leave('CalculateViscosity')

   if( Debug > 3 ) write(IOdef,*)'=== CalculateViscosity'

end subroutine CalculateViscosity
subroutine TurbulenceModels(iVar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX
   real, dimension(3)           :: dUdXp, dVdXp, dWdXp
   real, dimension(3)           :: Uw, Up, Xpn, Xn, Un, Ut

   if( Debug > 3 ) write(IOdef,*)'*** TurbulenceModels',TurbModel,Variable(iVar)

   call watch_enter('TurbulenceModels')


   if( TurbModel == TMkeps .or. TurbModel == TMrng )then
     !
     ! standard k-epsilon model
     !
     if( Debug > 3 ) write(IOdef,*)'k-epsilon model'
     if( iVar == VarTE )then
       if( Debug > 3 ) write(IOdef,*)'Turb. kinetic energy'
       !
       ! LHS = VisT * Production - rho * epsilon
       !     - compr. effects + non lin. term
       !
       ! the diffusion term is absorbed in the transp. eq.
       !
       ! source term
       !
       do ip=1,Ncel

         dUdXp = dUdX(ip,:)
         dVdXp = dVdX(ip,:)
         dWdXp = dWdX(ip,:)
         !
         ! rate of production of turbulent energy (eq. 9.40)
         !
         s1   = (dUdXp(1)+dUdXp(1))*dUdXp(1) + (dUdXp(2)+dVdXp(1))*dUdXp(2) + (dUdXp(3)+dWdXp(1))*dUdXp(3)
         s2   = (dVdXp(1)+dUdXp(2))*dVdXp(1) + (dVdXp(2)+dVdXp(2))*dVdXp(2) + (dVdXp(3)+dWdXp(2))*dVdXp(3)
         s3   = (dWdXp(1)+dUdXp(3))*dWdXp(1) + (dWdXp(2)+dVdXp(3))*dWdXp(2) + (dWdXp(3)+dWdXp(3))*dWdXp(3)

         VisT = VisEff(ip) - VisLam

         Pk   = VisT * ( s1 + s2 + s3 )
         !
         ! bouyancy production term: - Gi/(sigma_h,t rho) drho/dx
         !
         !Pbouy =

         TurbP(ip) = Pk

         !
         ! compres. amplification term
         !
         !Comp =

         !
         ! dissipation
         !

         Dis = Den(ip) * ED(ip)

         !
         ! finally the coefficients
         !
         ! pull -rho.eps on the RHS to the LHS by
         ! multiplying with k/k (=1.0) and using
         ! rho.eps/k as the coefficient.
         !
         ! if one uses LowerLimitTE (which is >= Small)
         ! Small can be removed here.

         Su(ip) = Su(ip) + TurbP(ip) * Cell(ip)%vol
         Au(ip) = Au(ip) + Dis /(TE(ip) + Small) * Cell(ip)%vol

       end do
       !
       ! solid walls
       !
       Cmu   = TMCmu
       Cmu25 = Cmu**0.25

       do ib=1,Nbnd
         i  = Bnd(ib)%face
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         ip = Face(i)%cell1
         in = Face(i)%cell2
         if( in /= 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

         if( it == RWall .and. .not. Initialisation )then

           !
           ! in the preceding loop the source term for ip has been
           ! calculated; it needs to be replaced next to a solid
           ! wall; thus first subtract the 'standard' production term
           !
           Su(ip) = Su(ip) - TurbP(ip) * Cell(ip)%vol

           !pk_oud = TurbP(ip)

           Uw(1)    = Reg(ir)%uvw(1)            !
           Uw(2)    = Reg(ir)%uvw(2)            ! wall bc's
           Uw(3)    = Reg(ir)%uvw(3)            !

           Up(1) = U(ip)
           Up(2) = V(ip)
           Up(3) = W(ip)

           Visac    =  VisEff(Ncel+ib)

           Xpn      = Face(i)%x - Cell(ip)%x    ! dist. wall center to P

           !coef = Visac * Face(i)%area / vector_length(Xpn)

           Xn = Face(i)%n                       ! de surface vector
           call normalise(Xn)                   ! normaal vector met lengte 1

           Up = Up - Uw                         ! snelh. verschil

           dp = dot_product( Up , Xn )
           Un = dp * Xn                         ! normale snelh. component
           Ut = Up - Un                         ! tang. snelh. component

          !Uvel  = vector_length( Ut )          !
           Uvel  = sqrt(dot_product(Ut,Ut))     ! inlined
           dn    = Bnd(ib)%distance             ! alleen de grootte nodig
           Tau_w = Visac * Uvel / dn            ! niet de richting

           if( dn <= Zero )then
             write(IOdef,*)'*** Error: Unexpected wall distance ',dn
             dn  = 1.e-6
           endif
           if( Tau_w < Zero )then
             write(IOdef,*)'*** Error: Wall shear stress < zero ',Tau_w
           endif
           !
           ! the production in the wall region is (eq. 9.50)
           !
           !    Pk = Tau_w d( Ut )/dn
           !
           ! for a node in the law of the wall region is (eq. 9.51)
           !
           !   d( Ut )/dn = u_tau/( kappa dn )
           !
           ! note that Uvel is already known based on the the law-of-the-wall
           ! or any other law so there is no need to use eq. 9.49.
           !
           ! 1/(kappa*dn) is used several times so set it in advance
           !
           rkapdn = 1.0/( kappa * dn )

           if( Bnd(ib)%yplus > Reg(ir)%ylog  )then

             turb = TE(ip)
             if( turb < Zero )then              ! reformulate using LowerLimitTE???
                turb = 0.001                    ! don't let k become negative
                write(IOdef,*)'+ TurbulenceModels: Found negative k, set to zero.'
             endif

             Utau      = Cmu25*sqrt(turb)                 ! (eq. 9.48)
             Tau_w     = Visac * Uvel / dn

             TurbP(iP) = Tau_w * Utau * rkapdn            ! (eq. 9.50/9.51)

           else

             Tau_w     = Visac * Uvel / dn
             Utau      = sqrt( Tau_w / Den(ip) )

             TurbP(iP) = Tau_w * Utau * rkapdn

           endif

           TE(Ncel+ib) = TE(ip)
           
           Su(ip) = Su(ip) + TurbP(ip) * Cell(ip)%vol

           !
           ! the dissipation term based on eq. 9.52:
           !
           !   dis_P = Cmu75 * k**(3/2) /( kappa * dn )
           !
           ! as above we need: -rho*dis*k (times volume) so
           ! use k**(3/2) = k * sqrt(k) and pull it to the RHS.
           !
           DisP   = Cmu75*sqrt(TE(ip))*rkapdn
           Au(ip) = Au(ip) + Den(ip) * DisP  * Cell(ip)%vol

         else if( it == RWall .and. Initialisation )then

           TE(Ncel+ib) = TE(ip)

         endif
       end do

     elseif( iVar == VarED )then
       if( Debug > 3 ) write(IOdef,*)'Turb. dissipation'

       !
       ! see eq. 9.42 for:
       !
       ! RHS = Ce1 Prod eps/k - rho Ce2 eps^2/k
       !
       !     = Ce1 Prod eps/k - rho Ce2 eps/k * eps
       !                        -------------
       !                              v
       !                   in the diagonal A-term * eps
       ! source term
       !
       do ip=1,Ncel

         fact = ED(ip)/(TE(ip)+Small) * Cell(ip)%vol

         Su(ip) = Su(ip) + TMCeps1 * fact * TurbP(ip)
         Au(ip) = Au(ip) + TMCeps2 * fact * Den(ip)

       end do

       !
       ! solid walls
       !
       Cmu   = TMCmu
       Cmu75 = Cmu**0.75

       do ib=1,Nbnd
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ
         if( it == RWall .and. .not. Initialisation )then
           i  = Bnd(ib)%face
           ip = Face(i)%cell1
           in = Face(i)%cell2
           if( in > 0 )write(IOdef,*)'*** Error: Bnd array corrupt'

           !
           ! due to the requirement Prod = Diss just
           ! force/fix the value in the wall cell (eq. 9.52)
           !

           turb = TE(ip)
           if( turb < Zero )then                ! reformulate using LowerLimitTE???
             turb = 0.001                       ! don't let k become negative
             write(IOdef,*)'+ TurbulenceModels: Found negative k, set to zero (2).'
           endif

           dn  = Bnd(ib)%distance

           Dis = Cmu75 * turb**1.5 /( dn * kappa )

           !
           ! to brutally fix a value enforce Ap * phi = Su
           ! (more efficient than bumping Su and Au?)
           !
           do j=1,NFaces(ip)
             k   = CFace(ip,j)
             ipf = Face(k)%cell1
             inf = Face(k)%cell2
             if( inf > 0 )then
               ! internal
               if( ipf == ip )then
                 RFace(k,2) = 0.0
               elseif( inf == ip )then
                 RFace(k,1) = 0.0
               else
                 write(IOdef,*)'+ internal error: No fix'
               endif
             endif
           end do

           if( .not. Initialisation )then
             ED(ip) = Dis
             Su(ip) = Dis
             Au(ip) = 1.0
           endif

           ED(Ncel+ib) = ED(ip)

         else if( it == RWall .and. Initialisation )then

           ED(Ncel+ib) = ED(ip)

         endif
       end do
     else
       write(IOdef,*)'+ internal error: invalid turb. model scalar:',Variable(iVar)
     endif
   !elseif( TurbModel == TMrng )then
   !  !
   !  ! RNG k-epsilon model
   !  !
   !  if( debug > 0 ) write(IOdef,*)'RNG k-epsilon model'

   else
     !
     ! unknown model
     !
     write(IOdef,*) 'Unknown turb. model / not implemented'
   endif

   call watch_leave('TurbulenceModels')

   if( Debug > 3 ) write(IOdef,*)'=== TurbulenceModels ',Variable(iVar)

end subroutine TurbulenceModels

