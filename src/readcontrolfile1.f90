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
subroutine ScanControlFile

   use constants
   use geometry
   use variables
   use scalars
   use particles
   
   character(len=120) string

   character(len=32)  :: key, key2, key3

   integer, parameter :: MaxKeys = 20
   character(len=32)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer, save      :: MathMode = 0

   write(IOdbg,*)'Start scan of input deck'

   call openfile(IOinp,casename,'.din','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)

   MaxSC  = 0
   MaxMat = 1
   
   do 
     read(IOinp,'(A)',end=20) string   ! read string

     if( string(1:1) == '#' ) cycle    ! skip if comment line
     if( string(1:4) == 'titl' ) cycle ! skip if comment line

     indx = index(string,'#')
     if( indx >= 2 )then               ! strip trailing comments
       do i=indx,len(string)
         string(i:i) = ' '
       end do
     endif

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     if( Nkeys == 0 ) cycle            ! skip blank line

     key = keys(1)
     call lowercase(key)

     select case (key) 

       case ('dimension')
         key = keys(2)
         call lowercase(key)
         select case (key(1:3))
           case ('scalars')
             idum = keyi(3)
             if( idum < 1 )then
               write(*,*)'+ error number of scalars should be > 0 :',idum
             endif

             MaxSC = max(MaxSC,abs(idum))
             write(IOdbg,'(1x,A,i3)')'+ scalars: MaxSC=',MaxSC 
           
           case ('materials')
             idum = keyi(3)
             if( idum < 1 )then
               write(*,*)'+ error number of materials should be > 0 :',idum
             endif

             MaxMat = max(MaxMat,abs(idum))
             write(IOdbg,'(1x,A,i3)')'+ materials: MaxMat=',MaxMat
           
           case default
             write(*,*)'+ error unknown dimension command: ',string
           
         end select
         
       case ('use')
         key = keys(2)
         call lowercase(key)
         select case (key(1:7))
           case ('patches')
             idum = keyi(3)
             if( idum < 0 )then
               write(*,*)'+ error number of scalars should be >= 0 :',idum
             endif 

             MaxSC = max(MaxSC,idum)
             write(IOdbg,'(1x,A,i3)')'+ scalars: MaxSC=',MaxSC 
           case ('particl')
             idum = keyi(3)
             if( idum < 0 )then
               write(*,*)'+ error number of particles should be >= 0 :',idum
             endif 

             Npart = idum
             write(IOdbg,'(1x,A,i3)')'+ particles: Npart=',Npart 
           case ('sensors')
             idum = keyi(3)
             if( idum < 0 )then
               write(*,*)'+ error number of sensors should be >= 0 :',idum
             endif 

             Nsens = idum
             write(IOdbg,'(1x,A,i3)')'+ sensors: Nsens=',Nsens 
         end select

       case ('set')
         key = keys(2)
         call lowercase(key)
         !write(*,*)'set>',key
         r1 = keyr(3)
         r2 = keyr(4)

         call setvariable(key(1:8),r1,r2)

       case ('math')
         key2 = keys(2)
         call lowercase(key2)

         if( key2(1:3) == 'deg' )then
           write(IOdbg,*)'Mathmode set to degrees'
           MathMode = 0
         else if( key2(1:3) == 'rad' )then
           write(IOdbg,*)'Mathmode set to radians'
           MathMode = 1         
         else
         
           key3 = keys(3)
           call lowercase(key3)

           r1  = keyr(4)
           r2  = keyr(5)
           r3  = keyr(6)
           r4  = keyr(7)

           call mathvariable(key3,MathMode,r1,r2,r3,r4)

           r2  = 0.0        
           call setvariable(key2(1:8),r1,r2)

         endif

       case default
         !write(*,*) '+ skipped keyword: ',key(1:lens(key))

     end select

   end do
20 continue
   write(IOdef,*)'End scan of input deck'
   write(IOdbg,*)'End scan of input deck'

   close(IOinp)
   rewind(IOinp)

   !
   ! dimension array's !!!
   !
   ! scalars
   !
   MaxVar = NVar + MaxSC

   write(IOdbg,*)'MaxSc  = ',MaxSC
   write(IOdbg,*)'MaxVar = ',MaxVar
   
   allocate(Variable(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Variable array allocated')

   allocate(Residual(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Residual array allocated')

   allocate(ResiNorm(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'ResiNorm array allocated')

   Residual = 0.0
   ResiNorm = 1.0

   allocate(Solver(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Solver array allocated')

   Solver   = SparseKit2

   allocate(Gamma(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Gamma array allocated')

   Gamma    = 0.0

   allocate(Scheme(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Scheme array allocated')

   Scheme   = DSud

   allocate(Solve(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Solve array allocated')

   Solve(:)        = .false.
   Solve(1:4)      = .true. 

   allocate(Store(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'Store array allocated')

   Store(:)        = .false.
   Store(1:4)      = .true. 

   allocate(PostC(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'PostC array allocated')

   PostC(:)        = .false. 
   PostC(1:4)      = .true. 

   allocate(PostV(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar,'PostV array allocated')

   PostV(:)        = .false.

   allocate(UseSlopeLimiter(MaxVar),stat=istat)
   allocate(URFSlopeLimiter(MaxVar),stat=istat)
   allocate(LimiterSlope(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar*3,'Slope limiter arrays allocated')

  !UseSlopeLimiter(:) = .false.
  !URFSlopeLimiter(:) = 1.0
  !LimiterSlope(:)    = SlopeLimiterOFF

   UseSlopeLimiter(:) = .true.
   URFSlopeLimiter(:) = 1.0
   LimiterSlope(:)    = SlopeLimiterVNf
  !LimiterSlope(VarPP)= SlopeLimiterOff

   allocate(GradVar(MaxVar),stat=istat)
   allocate(nGradVar(MaxVar),stat=istat)
   call TrackMemory(istat,MaxVar*2,'Variable Gradient arrays allocated')

   GradVar(:)         = GradGauss 
   nGradVar(:)        = Ngradient

   if( MaxSC > 0 )then
     allocate(Scalar(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'Scalar array allocated')

     allocate(SolveScalar(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'SolveScalar array allocated')

     allocate(StoreScalar(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'StoreScalar array allocated')

     SolveScalar = .false.
     StoreScalar = .false.

     allocate(VarS(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'VarS array allocated')

     allocate(LimitLowScalar(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'LimitLowScalar array allocated')
     allocate(LimitUpScalar(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'LimitUpScalar array allocated')
     
     LimitLowScalar = .false.
     LimitUpScalar  = .false.
     
     allocate(LowerLimitSC(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'LowerLimitSC array allocated')
     allocate(UpperLimitSC(MaxSC),stat=istat)
     call TrackMemory(istat,MaxSC,'UpperLimitSC array allocated')

     LowerLimitSC =  0.0
     UpperLimitSC =  0.0
     
     do i=1,MaxSC
       indx    = Nvar + i
       VarS(i) = indx
       write(key,'(''Scalar'',i3.3,''   '')') i 
       Variable(indx)  = key(1:12)
       Scalar(i)       = key(1:12)
     end do
     
     NScal = 0
     
   endif  
   
   Variable(VarU )   = 'U       '
   Variable(VarV )   = 'V       '
   Variable(VarW )   = 'W       '
   Variable(VarP )   = 'Pressure'

   Variable(VarTE)   = 'KE      '
   Variable(VarED)   = 'ED      '
   Variable(VarT )   = 'T       '
   Variable(VarSC)   = 'Scalar  '

   Variable(VarDen)  = 'Density '
   Variable(VarPP)   = 'PressCor'
   Variable(VarVis)  = 'VisEff  '
   Variable(VarLVis) = 'VisLam  '
   Variable(VarCP)   = 'Cp      '

   write(IOdbg,*)'Variable control arrays dimensioned'

   !
   ! particles
   !
   if( Npart > 0 )then
     allocate(Particle(NPart),stat=istat)
     call TrackMemory(istat,Npart,'Particle array allocated')
   endif

   if( Nsens > 0 )then
     allocate(Sensor(Nsens),stat=istat)
     call TrackMemory(istat,Nsens,'Sensor array allocated')
   endif

end subroutine ScanControlFile
subroutine ReadControlFile

   use constants
   use geometry
   use variables
   use scalars
   use particles
   use watches
   
   character(len=120) string

   character(len=32)  :: key, key2, key3, key4, blank

   integer, parameter :: MaxKeys = 20
   character(len=32)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer            :: ideltas(MaxKeys)
   real               :: rdeltas(MaxKeys)

   integer, save      :: MathMode = 0

   logical :: Messages = .false.
   logical :: Echo     = .false.
   
   logical :: LGenerator
   
  !call watch_enter('ReadControlFile')

   Noutlets = 0
   Split    = 0.0
   
   blank      = '                                '

   Small = tiny(small) * 1.e+3
   Large = huge(large) * 1.e-3
   
   SMALL = 1.e-12

   call SetUpVariables
   call ScanControlFile

   call SetUpVariables
   call openfile(IOinp,casename,'.din','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)
   
   do 
     read(IOinp,'(A)',end=20)    string   ! read string
     !write(IOdbg,'(1x,''>'',A)') string

     if( string(1:1) == '#' ) cycle    ! skip if comment line

     if( string(1:5) == 'title' .and.                     &
       ( string(6:6) == ',' .or. string(6:6) == ' ' .or.  &
         string(6:6) == '=') )then
       title = string(7:len(string))
       !write(*,'(1x,A,A)')'+ title ',title(1:lens(title))
       cycle
     endif

     indx = index(string,'#')
     if( indx >= 2 )then               ! strip trailing comments
       do i=indx,len(string)
         string(i:i) = ' '
       end do
     endif

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

     if( Nkeys == 0 ) cycle           ! skip blank line

     key = keys(1)
     call lowercase(key)

     !write(*,*)'using key:',key

     select case (key) 

       case ('debug')
         i = keyi(2)
         debug = i
         if( debug > 1 ) Echo = .true.
         if( Echo ) write(*,'(1x,A,i2)')'+ debug level = ',i

       case ('echo')

         call select_case_echo(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('steps')
         i = keyi(2)
         if( Echo ) write(*,*)'+ iterations/time steps = ',i
         Niter = i
  
         if( keyr(3) > 0.0 )then
           ResMax = keyr(3)
           if( Echo ) write(*,'(1x,A,1pe9.3)')'+ with target residual ',ResMax
         endif

       case ('pref')
         !
         ! nog chekken of in de juiste range/materiaal is enz.
         !
         i = keyi(2)
         if( i > Ncel )then
           write(*,*)'+++ Warning: Pressure reference cel out of range'
           i = Ncel 
         endif
         
         if( Echo ) write(*,*)'+ reference pressure in cell = ',i
         IPref = i

       case ('monitor')
         !
         ! nog chekken of in de juiste range/materiaal is enz.
         !
         i = keyi(2)
         if( i > Ncel )then
           write(*,*)'+++ Warning: Monitor cel out of range'
           i = Ncel 
         endif

         if( Echo ) write(*,*)'+ monitor cell = ',i
         IMoni = i

       case ('save')
         
         call select_case_save(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
         
       case ('output')
         
         call select_case_output(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
         
       case ('boundary')

         call select_case_boundary(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('ngradient')
         i = keyi(2)
         if( Echo ) write(*,'(1x,A,i3)')'+ gradient iterations = ',i
         Ngradient = i

       case ('thermal')
        
         call select_case_thermal(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('turbulence')
         if( Echo ) write(*,'(1x,A,i3)')'+ turbulence options'
         key = keys(2)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
         select case (key(1:3))
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ turbulence modelling mode switched off'
             SolveTurb       = .false.
             SolveVisc       = .false.
             SolveTurbEnergy = .false.
             SolveTurbDiss   = .false.
             TurbModel       =  TMnone
           case ('ke')
             if( Echo ) write(*,'(1x,A,i3)')'+ standard k-e turbulence model'
             SolveTurb       = .true.
             SolveVisc       = .true.
             SolveTurbEnergy = .true.
             SolveTurbDiss   = .true.
             TurbModel       =  TMkeps
             if( keyr(3) > Small )then
               TMLenSc = keyr(3)
             endif
           case ('rng')
             if( Echo ) write(*,'(1x,A,i3)')'+ rng k-e turbulence model'
             SolveTurb       = .true.
             SolveVisc       = .true.
             SolveTurbEnergy = .true.
             SolveTurbDiss   = .true.
             TurbModel       =  TMrng
         end select

       case ('scalars')
         !
         write(*,*)'*** Scalars: experimental feature ***'
         !
         if( Echo ) write(*,'(1x,A,i3)')'+ scalar options'
         key = keys(2)
         call lowercase(key)
         if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
         select case (key(1:3))
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ scalars switched off'
             SolveScalars = .false.
           case ('on ')
             if( Echo ) write(*,'(1x,A,i3)')'+ scalars switched on'
             SolveScalars = .true.
             UseScalars   = .true.
         end select

       case ('restart')

         call select_case_restart(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('gravity')
         Gravity(1) = keyr(2)
         Gravity(2) = keyr(3)
         Gravity(3) = keyr(4)
         if( Echo ) write(*,'(1x,A,3(1x,f6.3))')'+ gravity vector = ',Gravity(1:3)

       case ('beta')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ volumetric expansion coef = ',dummy
         beta = dummy
         
       case ('transient')
         dummy = keyr(2)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ transient, time step (dt) = ',dummy
         dt = dummy
         Transient = .true.
         Euler     = .true.
         QuadTime  = .false.
         if( keys(3) == 'quad' )then
           Euler    = .false.
           QuadTime = .true.
           GammaTime = keyr(4)
         endif
         !
         ! RESET RELAXATION FACTORS!
         !
         !URF(VarU ) = 1.0
         !URF(VarV ) = URF(VarU)
         !URF(VarW ) = URF(VarU)
         !URF(VarP ) = 0.5
         !URF(VarTE) = 1.0
         !URF(VarED) = URF(5)
         !URF(VarT ) = URF(5)
         !URF(VarSC) = URF(5)

       case ('vislam')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ laminar viscosity = ',dummy
         vislam = dummy

       case ('density')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ density = ',dummy
         DensRef = dummy

       case ('cp')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ specific heat Cp = ',dummy
         CpStd = dummy

       case ('cv')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ specific heat Cv = ',dummy
         CvStd = dummy

       case ('prandtl')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ Prandtl number = ',dummy
         Prandtl = dummy
         Lambda  = VisLam * CpStd / Prandtl

       case ('conductivity')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ conductivity = ',dummy
         Lambda = dummy
         Prandtl = VisLam * CpStd / Lambda

       case ('schmidt')
         dummy = keyr(2)
         call check_minmax(dummy,Small,Large)
         if( Echo ) write(*,'(1x,A,1pe9.3)')'+ Schmidt number = ',dummy
         Schmidt = dummy

       case ('relax*')
  
         call select_case_relax_star(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('relax')

         call select_case_relax(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
         
       case ('rtol') 

         call select_case_rtol(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
       
       case ('atol')

         call select_case_atol(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
       
       case ('init')
             
         call select_case_init(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
       
       case ('pcor')
         key2 = keys(2)
         call lowercase(key2)
         key3 = keys(3)
         call lowercase(key3)
         if( Echo ) &
           write(*,'(1x,A,A,A)')'+ pcor keyword with strings = ',key2,key3
         select case (key2(1:3))
           case ('max')
             MaxPCOR = keyi(3)
             if( Echo ) write(*,'(1x,A,i3)')'+ maximum corr. steps = ',MaxPCOR
           case ('fac')
             FactDPP = keyr(3)
             if( Echo ) write(*,'(1x,A,f6.3)')'+ target factor = ',FactDPP
           case default
             write(*,'(1x,A,i3)')'+ error in pcor command'
             messages = .true.
         end select

       case ('gamma*')

         call select_case_gamma_star(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
         
       case ('gamma')

         call select_case_gamma(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('double')
         !
         ! obsolete!!!
         !
         write(*,'(1x,A,A)')'+ double keyword obsolete! ignored'

       case ('solver')

         call select_case_solver(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
!
! choice of differencing schemes
!
       case ('scheme*')
         Scheme(VarU ) = keyi(2)                                      
         Scheme(VarV ) = keyi(3)                                      
         Scheme(VarW ) = keyi(4)                                      
        !Scheme(VarP ) = keyi(5)                                      
         Scheme(VarTE) = keyi(6)                                      
         Scheme(VarED) = keyi(7)                                        
         Scheme(VarT ) = keyi(8)                                      
         Scheme(VarSC) = keyi(9)                                     
         if( Echo ) write(*,'(1x,A,4(1x,i2))')'+ scheme (*) = ',Scheme(1:4)    
         if( Echo ) write(*,'(1x,A,4(1x,i2))')'               ',Scheme(5:8)    

       case ('scheme')

         call select_case_scheme(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

!
! set scalar limits, slope limiters and gradients
!
       case ('limit')

         call select_case_limit(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('slope')

         call select_case_slopelimiter(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('grad')

         call select_case_gradient(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

!
! special models, commands, switches
!
       case ('use')
         key2 = keys(2)
         call lowercase(key2)

         select case (key2(1:5))
           case ('fixab')
             i = keyi(3)
             if( i == 0 ) write(*,*)'WARNING ctid incorrect!'
             
             tmpu = keyr(4)
             tmpv = keyr(5)
             tmpw = keyr(6)
             tmpk = keyr(7)
             tmpe = keyr(8)
             
             write(*,'(1x,A,i3)')'Fixing ABL top layer, ctid =',i
             write(IOdbg,'(1x,A,i3)')'Fixing ABL top layer, ctid =',i
             
             UseFixABL = .true.
             IdFixABL  = i
             UFixABL   = tmpu
             VFixABL   = tmpv
             WFixABL   = tmpw
             TeFixABL  = tmpk
             EdFixABL  = tmpe
             
             write(*,*)'Using:',tmpu,tmpv,tmpw,tmpk,tmpe
             write(IOdbg,*)'Using:',tmpu,tmpv,tmpw,tmpk,tmpe
           case ('artif')
             key3 = keys(3)
             call lowercase(key3)
             select case (key3(1:3))
               case ('off')
                 if( Echo ) write(*,'(1x,A,i3)')'+ artificial compressibility off'
                 UseArtificialComp = .false.
               case ('on')     
                 if( Echo ) write(*,'(1x,A,i3)')'+ artificial compressibility on'
                 UseArtificialComp = .true.
               case ('   ')     
                 if( Echo ) write(*,'(1x,A,i3)')'+ artificial compressibility on'
                 UseArtificialComp = .true.
               case default
                 write(*,'(1x,A,A)')'+ use,artif: unknown key: ',key3(1:3)
             end select
           case ('vtk')
             UseVTK = .true.
             key3 = keys(3)
             call lowercase(key3)
            !if( Echo ) &
               write(*,'(1x,A,A)')'+ using VTK',key3
             select case (key3(1:5))
               case ('ascii')
                 if( Echo ) write(*,'(1x,A,i3)')'+ write out ascii file'
                 UseVTKbinary = .false.
               case ('binar')
                 if( Echo ) write(*,'(1x,A,i3)')'+ write out binary file'
                 UseVTKbinary = .true.
               case ('fluid')     
                !if( Echo )
                 write(*,'(1x,A,i3)')'+ write out only fluid cells'
                 UseVTKwalls = .false.
               case ('     ')     
                 if( Echo ) write(*,'(1x,A,i3)')'+ binary write out boundaries as cells'
                 UseVTKwalls  = .true.
                 UseVTKbinary = .true.
               case default
                !if( Echo ) 
                 write(*,'(1x,A,i3)')'+ default binary write out boundaries as cells'
                 UseVTKwalls  = .true.
                 UseVTKbinary = .true.
             end select
             
           case ('opend')
             if( Echo ) write(*,'(1x,A,i3)')'+ using OpenDX'
             UseOpenDX = .true.
           case ('tecpl')
             if( Echo ) write(*,'(1x,A,i3)')'+ using Tecplot'
             UseTECPLOT = .true.
           case ('fieldview')
             if( Echo ) write(*,'(1x,A,i3)')'+ using Fieldview'
             UseFIELDVIEW = .true.
           case ('ensight')
             if( Echo ) write(*,'(1x,A,i3)')'+ using Ensight'
             UseENSIGHT = .true.
           case ('gmv')
             if( Echo ) write(*,'(1x,A,i3)')'+ using GMV'
             UseGMV = .true.
           case ('gmsh')
             UseGMSH = .true.
             key3 = keys(3)
             call lowercase(key3)
            !if( Echo ) &
               write(*,'(1x,A,A)')'+ using GMSH ',key3
             select case (key3(1:5))
               case ('ascii')
                 if( Echo ) write(*,'(1x,A,i3)')'+ write out ascii file'
                !UseGMSHbinary = .false.
               case ('binar')
                 if( Echo ) write(*,'(1x,A,i3)')'+ write out binary file'
                !UseGMSHbinary = .true.
               case ('full ')
                 if( Echo ) write(*,'(1x,A,i3)')'+ write out boundaries as cells'
                 UseGMSHwalls    = .true.
                 UseGMSHdataonly = .false.
               case ('fluid')     
                !if( Echo )
                 write(*,'(1x,A,i3)')'+ write out only fluid cells'
                 UseGMSHwalls    = .false.
                 UseGMSHdataonly = .false.
               case ('     ')     
                 if( Echo ) write(*,'(1x,A,i3)')'+ binary write out boundaries as cells'
                 UseGMSHwalls    = .true.
                !UseGMSHbinary   = .true.
                 UseGMSHdataonly = .true.
               case default
                !if( Echo ) 
                 write(*,'(1x,A,i3)')'+ default binary write out boundaries as cells'
                 UseGMSHwalls    = .true.
                !UseGMSHbinary   = .true.
                 UseGMSHdataonly = .true.
             end select

           case ('patch')
             if( Echo ) write(*,'(1x,A,i3)')'+ using patches by. B. Tuinstra'
             UsePatches   = .true.  ! switch on
             
             if( keyi(3) > 0 )then
               UseScalars   = .true.  ! patches relies on scalars, so switch on
               SolveScalars = .true.  ! 
             endif
             
             ! from ScanControlfile:
             if( Echo ) write(*,'(1x,A,i3)')'+ max. scalar dimension set to',MaxSC
             write(IOdbg,'(1x,A,i3)')'+ scalars: MaxSC=',MaxSC              
             
           case ('lapac')
             if( Echo ) write(*,'(1x,A,i3)') &
                    '+ using LAPACK in GradientPhiLeastSquares'
             UseLapack = .true.
           case ('parti')
             if( Echo ) write(*,'(1x,A,i3)')'+ using particles'
             UseParticles = .true.
           case ('senso')
             if( Echo ) write(*,'(1x,A,i3)')'+ using sensors'
             UseSensors = .true.
           
           case ('least')
             if( Echo ) write(*,'(1x,A,i3)')'+ using least squares for gradient'
             GradAlg = GradLS
             GradVar(:) = GradLS
           
           case ('gauss')
             i = keyi(3)
             if( i == 0 ) i = 2
             Ngradient = i

             if( Echo ) write(*,'(1x,A,i3)')'+ using Gauss for gradient',&
                                            Ngradient
             GradAlg = GradGauss 
             GradVar(:)  = GradGauss 
             nGradVar(:) = Ngradient
             
           case default
             write(*,*)'+ use: unknown key: ',key2(1:5)
         end select
!
!===== postprocessors ==================================================
!
       case ('opendx')

         call select_case_opendx(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
       
       case ('post')

         call select_case_post(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

       case ('check')

         call select_case_check(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

!
!===== user ============================================================
!
       case('user')
       
         UserInt(1) = keyi(2)
         UserInt(2) = keyi(3)
         UserInt(3) = keyi(4)
         UserInt(4) = keyi(5)
 
         UserReal(1) = keyr(2)
         UserReal(2) = keyr(3)
         UserReal(3) = keyr(4)
         UserReal(4) = keyr(5)
          	       
         write(IOdef,*)'User data i:',UserInt(1:Nuser)
         write(IOdef,*)'User data r:',UserReal(1:Nuser)                      
!
!===== particles =======================================================
!
       case ('gene')
         write(*,'(1x,A,i3)')'+ generation ignored here'
         messages = .true.
       
       case ('part')
         if( .not. UseParticles )then
           write(*,*)'+ error use particles not set'                         
           messages = .true.
           cycle
         endif
         !
         ! two forms:
         !
         ! PARTicle,#,PROP,dens,diam,COUPled etc
         ! 
         ! PARTicle,#,INIT,icoor,x,y,z,u,v,w ; icoor = 1 (global cartesian)
         ! 
         ! followed by a generation command. example:
         !
         ! PART,   1,PROP,1000.0,100e-6    define 1 particle
         ! GENE 19 1,,       0.0, 10e-6    and generate 19 others 
         !
         
         ipart = keyi(2) 
         key2  = keys(3)         
         call lowercase(key2)
         
         select case (key2(1:4))

           case ('prop')
             r1 = keyr(4) ! dens     
             r2 = keyr(5) ! diam

             Particle(ipart)%dens0 = r1
             Particle(ipart)%diam0 = r2
             
             if( LGenerator() )then
               call GetGenerators(ngen,ideltas,rdeltas)

               if( (ipart+ngen*ideltas(1)) > Npart )then
                 write(*,'(1x,A,4(1x,i6))')'+ too many particles generated (1):', &
                                        ipart,ngen,ideltas(1),ipart+ngen*ideltas(1)
                 messages = .true.                
               else
                 do i=1,ngen
                   j = i*ideltas(1)
                   Particle(ipart+j)%dens0 = Particle(ipart)%dens0 + float(i)*rdeltas(3) 
                   Particle(ipart+j)%diam0 = Particle(ipart)%diam0 + float(i)*rdeltas(4) 
                 end do
               endif
             endif

           case ('init')
             i1 = keyi(4) ! coordinate system (still to be done)
             if( i1 <= 0 .or. i1 >= 3 )then
               write(*,'(1x,A,i3)')'+ coordinate system still 1 or missing'
               messages = .true.                
             endif

             if( i1 == 1 )then
               r1 = keyr( 5) ! x     
               r2 = keyr( 6) ! y  
               r3 = keyr( 7) ! z  
               r4 = keyr( 8) ! u  
               r5 = keyr( 9) ! v  
               r6 = keyr(10) ! w  
             else
               x0 = 0.5
               y0 = 0.5       ! simpel test hack
               z0 = 0.125
               u0 = 0.0
               v0 = 0.0       
               w0 = 0.0

               !r1 = x0 + keyr( 5)*cosd(keyr( 6))      
               !r2 = y0 + keyr( 5)*sind(keyr( 6))     
               !r3 = keyr( 7)
               !r4 =      keyr( 8)*cosd(keyr( 9))   
               !r5 =      keyr( 8)*sind(keyr( 9))   
               !r6 = keyr(10)  
               write(*,*)'C>',r1,r2,r3,r4,r5,r6             
             endif
             
             Particle(ipart)%x0(1) = r1
             Particle(ipart)%x0(2) = r2
             Particle(ipart)%x0(3) = r3

             Particle(ipart)%v0(1) = r4
             Particle(ipart)%v0(2) = r5
             Particle(ipart)%v0(3) = r6


             if( LGenerator() )then
               call GetGenerators(ngen,ideltas,rdeltas)

               if( (ipart+ngen*ideltas(1)) > Npart )then
                 write(*,'(1x,A,4(1x,i6))')'+ too many particles generated (1):', &
                                        ipart,ngen,ideltas(1),ipart+ngen*ideltas(1)
                 messages = .true.                
               else
                 if( i1 == 1 )then
                   do i=1,ngen
                     j = i*ideltas(1)
                     Particle(ipart+j)%x0(1) = Particle(ipart)%x0(1) + float(i)*rdeltas(4) 
                     Particle(ipart+j)%x0(2) = Particle(ipart)%x0(2) + float(i)*rdeltas(5) 
                     Particle(ipart+j)%x0(3) = Particle(ipart)%x0(3) + float(i)*rdeltas(6) 
                     Particle(ipart+j)%v0(1) = Particle(ipart)%v0(1) + float(i)*rdeltas(7) 
                     Particle(ipart+j)%v0(2) = Particle(ipart)%v0(2) + float(i)*rdeltas(8) 
                     Particle(ipart+j)%v0(3) = Particle(ipart)%v0(3) + float(i)*rdeltas(9) 
                   end do
                 else
                   do i=1,ngen
                     j = i*ideltas(1)
                     !Particle(ipart+j)%x0(1) = x0 + keyr(5)*cosd(keyr(6)+float(i)*rdeltas(5))
                     !Particle(ipart+j)%x0(2) = y0 + keyr(5)*sind(keyr(6)+float(i)*rdeltas(5))
                     !Particle(ipart+j)%x0(3) = keyr( 7) + float(i)*rdeltas(6) 
                     !Particle(ipart+j)%v0(1) =      keyr(7)*cosd(keyr(8)+float(i)*rdeltas(8))
                     !Particle(ipart+j)%v0(2) =      keyr(7)*sind(keyr(8)+float(i)*rdeltas(8))
                     !Particle(ipart+j)%v0(3) = keyr(10) + float(i)*rdeltas(9) 
                   end do
                 endif
                 
               endif
             endif
                      
           case default
             write(*,'(1x,A,i3)')'+ error in particle command'
             messages = .true.
         end select

       !
       ! sensors
       !       
       case ('sens')
         if( .not. UseSensors )then
           write(*,*)'+ error use sensors not set'                         
           messages = .true.
           cycle
         endif
         !
         ! only one form:
         !
         ! SENSor,#,icoor,x,y,z ; icoor = 1 (global cartesian)
         ! 
         ! followed by a generation command. example:
         !
         ! SENS,   1,1,0.1,0.1,0.1    define 1 sensor
         ! GENE 19 1, ,0.1,0.0,0.0    and generate 19 others 
         !
         isens = keyi(2) 
         
         i1 = keyi(3) ! coordinate system (still to be done)
         if( i1 <= 0 .or. i1 >= 3 )then
           write(*,'(1x,A,i3)')'+ coordinate system still 1 or missing'
           messages = .true.                
         endif

         if( i1 == 1 )then
           r1 = keyr(4) ! x     
           r2 = keyr(5) ! y  
           r3 = keyr(6) ! z  
         else
           x0 = 0.5
           y0 = 0.5       ! simpel test hack
           z0 = 0.125
           write(*,*)'C>',r1,r2,r3,r4,r5,r6             
         endif

         Sensor(isens)%x(1) = r1
         Sensor(isens)%x(2) = r2
         Sensor(isens)%x(3) = r3
         Sensor(isens)%defined = .true.

         if( LGenerator() )then
           ideltas =  0
	   rdeltas = 0.0
	   call GetGenerators(ngen,ideltas,rdeltas)

           if( (isens+ngen*ideltas(1)) > Nsens )then
             write(*,'(1x,A,4(1x,i6))')'+ too many sensors generated (1):', &
                                    isens,ngen,ideltas(1),isens+ngen*ideltas(1)
             messages = .true.                 
           else
             if( i1 == 1 )then
               do i=1,ngen
                 j = i*ideltas(1)
		!write(*,*)'s>',i,j,rdeltas(3:5)
                 Sensor(isens+j)%x(1) = Sensor(isens)%x(1) + float(i)*rdeltas(3) 
                 Sensor(isens+j)%x(2) = Sensor(isens)%x(2) + float(i)*rdeltas(4) 
                 Sensor(isens+j)%x(3) = Sensor(isens)%x(3) + float(i)*rdeltas(5) 
                 Sensor(isens+j)%defined = .true.
               end do
             else
               do i=1,ngen
                 j = i*ideltas(1)
                 ! todo other coor.
               end do
             endif

           endif
         endif
!
!===== special print options ===========================================
!
       case ('print')
   
         call select_case_print(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

!
!===== rarely used changed constants/switches ==========================
!
       case ('switch')
         key2 = keys(2)
         call lowercase(key2)
         !write(*,*)'switch> ',key2 
         select case (key2(1:4))
           case ('cour')
             r = keyr(3)
             if( r < 1e-6 .or. r > 2.0 )then
               write(*,*)'Invalid particle Courant number:',r   
               messages = .true.
             else    
               ParticleCourant = r                          
               write(IOdef,*)'Particle Courant number:',r
               write(IOdbg,*)'Switch: Particle Courant number:',r 
               write(IOrun,*)'Switch: Particle Courant number:',r 
             endif
           case ('maxo')
             i1 = keyi(3)
             if( i1 <= 0 )then
               write(*,*)'Invalid number of outer iterations (MaxOuter):',i1   
               messages = .true.
             else
               write(IOdef,*)'Number of outer iterations set to:',i1
               write(IOdbg,*)'Number of outer iterations set to:',i1
               write(IOrun,*)'Number of outer iterations set to:',i1
               MaxOuter = i1
             endif
           case ('pp  ')

             r = keyr(3)
             if( r < 1e-6 .or. r > 1.0 )then
               write(*,*)'Invalid pressure correction underrelaxation factor:',r   
               messages = .true.
             else    
               URF(VarPP) = r                          
               write(IOdef,*)'Pressure correction urf:',r
               write(IOdbg,*)'Switch: Pressure correction urf:',r 
               write(IOrun,*)'Switch: Pressure correction urf:',r 
             endif

           case ('cds ')
             key3 = keys(3)
             call lowercase(key3)
             if( key3(1:3) == 'cd1' )then 
	       DScdDefault = DScd1
	     elseif( key3(1:3) == 'cd2' )then 
	       DScdDefault = DScd2
	     elseif( key3(1:3) == 'cd3' )then 
	       DScdDefault = DScd3
	     else
	       write(*,*)'Unknown default CDS scheme ',key3
	     endif
             write(IOdef,*)'Default CDS:',DScdDefault
             write(IOdbg,*)'Switch: Default CDS:',DScdDefault
             write(IOrun,*)'Switch: Default CDS:',DScdDefault
	     
          !case ('cbc ')
          !  key3 = keys(3)
          !  call lowercase(key3)
          !  if( key3(1:3) == 'on ' )then 
	  !    UseCBC = .true.
	  !  elseif( key3(1:3) == 'off' )then 
	  !    UseCBC = .false.
	  !  else
	  !    write(*,*)'Switch CBC on or off ',key3
	  !  endif
          !  write(IOdef,*)'Convective Boundness Criterion (CBC):',key3(1:3)
          !  write(IOdbg,*)'Convective Boundness Criterion (CBC):',key3(1:3)
          !  write(IOrun,*)'Convective Boundness Criterion (CBC):',key3(1:3)

           case default
             write(*,'(1x,A)')'+ error in switch command'
             messages = .true.
         end select
!
!===== variables / library =============================================
!
       case ('set')
         key2 = keys(2)
         call lowercase(key2)

         r1  = keyr(3)
         r2  = keyr(4)
        
         call setvariable(key2(1:8),r1,r2)

       case ('math')
         key2 = keys(2)
         call lowercase(key2)

         if( key2(1:3) == 'deg' )then
           write(*,*)'Mathmode set to degrees'
           MathMode = 0
         else if( key2(1:3) == 'rad' )then
           write(*,*)'Mathmode set to radians'
           MathMode = 1         
         else
         
           key3 = keys(3)
           call lowercase(key3)

           r1  = keyr(4)
           r2  = keyr(5)
           r3  = keyr(6)
           r4  = keyr(7)

           call mathvariable(key3,MathMode,r1,r2,r3,r4)

           r2  = 0.0        
           call setvariable(key2(1:8),r1,r2)

         endif
!
!===== default =========================================================
!
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key(1:lens(key))
         write(*,*) '+ string: ',string
         messages = .true.

     end select
     cycle 

10   continue
     write(*,*)'+ error in: ',string(1:lens(string))
     messages = .true.
     cycle

11   continue
     write(*,*)'+ read error: something else expected'
     write(*,*)'+ trying to continue'
     messages = .true.

   end do
20 continue
   write(*,*)'End of input deck'
   close(IOinp)

   !
   ! some final checks
   !
   if( Noutlets == 0 )then
     write(*,*) 'No outlets'
   elseif( Noutlets == 1 )then
     write(*,*) 'One outlet'
     if( abs(Split-1.0) > Small )then
       write(*,*)'Error: Split not equal to 1.0',Split
       messages = .true.       
     endif     
   elseif( Noutlets > 1 )then
     write(*,*) 'Multiple outlets'
     if( abs(Split-1.0) > 1.e-6 )then
       write(*,*)'Error: Sum of splits greater than 1.0', Split
       write(*,*)' > > : ',abs(Split-1.0) 
       messages = .true.       
     else
       write(*,*)'Sum split ',Split     
     endif
   endif

   !
   ! fill in the Solve-array 
   !
   Solve(VarU)    = SolveU
   Solve(VarV)    = SolveV
   Solve(VarW)    = SolveW
   Solve(VarP)    = SolveP
   Solve(VarTE)   = SolveTurbEnergy 
   Solve(VarED)   = SolveTurbDiss
   Solve(VarT)    = SolveEnthalpy
   Solve(VarSC)   = SolveScalars
   Solve(VarDEN)  = .false.
   Solve(VarVis)  = SolveVisc
   Solve(VarLVis) = .false.
   Solve(VarCP)   = .false.

   if( SolveScalars )then
     do is=1,Nscal
       Solve(NVar+is) = .true.
     end do
   endif

   !
   ! nicer/cleaner
   !
   if( Scheme(VarU)  == DScds .and. Gamma(VarU)  <= 0.001 ) Scheme(VarU)  = DSud
   if( Scheme(VarV)  == DScds .and. Gamma(VarV)  <= 0.001 ) Scheme(VarV)  = DSud
   if( Scheme(VarW)  == DScds .and. Gamma(VarW)  <= 0.001 ) Scheme(VarW)  = DSud
   if( Scheme(VarTE) == DScds .and. Gamma(VarTE) <= 0.001 ) Scheme(VarTE) = DSud
   if( Scheme(VarED) == DScds .and. Gamma(VarED) <= 0.001 ) Scheme(VarED) = DSud
   if( Scheme(VarT)  == DScds .and. Gamma(VarT)  <= 0.001 ) Scheme(VarT)  = DSud
   if( Scheme(VarSC) == DScds .and. Gamma(VarSC) <= 0.001 ) Scheme(VarSC) = DSud

   if( Scheme(VarU)  == DSud ) Gamma(VarU)  = 0.0
   if( Scheme(VarV)  == DSud ) Gamma(VarV)  = 0.0
   if( Scheme(VarW)  == DSud ) Gamma(VarW)  = 0.0
   if( Scheme(VarTE) == DSud ) Gamma(VarTE) = 0.0
   if( Scheme(VarED) == DSud ) Gamma(VarED) = 0.0
   if( Scheme(VarT)  == DSud ) Gamma(VarT)  = 0.0
   if( Scheme(VarSC) == DSud ) Gamma(VarSC) = 0.0
   
  !call watch_leave('ReadControlFile')

   if( messages )then
     write(*,*)'***'
     write(*,*)'*** Errors and/or warnings in input deck.'
     write(*,*)'***'
     stop
   endif

end subroutine ReadControlFile
subroutine parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)

   character(len=*)   :: string
   character(len=120) :: tmpstring
   character(len=8)   :: name, value

   character(len=32)  :: keys(MaxKeys), blank
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   logical            :: math, found
   character(len=1)   :: s,s1 
 
   keyu(1:MaxKeys) = .false.
   math = .false.
   
   !
   ! test for blank line
   !
   if( lens(string) == 0 )then
     !write(*,*)'+ blank line'
     Nkeys = 0
     return
   endif

   !
   ! change first occurance of = into ,
   !
   indx = index(string,'=')
   if( indx > 0 ) string(indx:indx) = ' '
   !
   ! first remove extra spaces
   !  
   do i=1,len(string)
     s = string(i:i)
     if( s .eq. ' ' )then
       do j=i+1,len(string)
         s1 = string(j:j)
         if( s1 /= ' ' ) exit
       end do
       if( j /= (i+1) )then
         tmpstring = string(j:len(string)) 
         string(i+1:) = tmpstring
       endif
     endif
   end do

   !
   ! remove trailing komma's if any
   !
   do while( .true. )
     ilast = lens(string)
     if( string(ilast:ilast) == ',' )then
       string(ilast:ilast) = ' '
     else
       exit
     endif 
   end do

   !
   ! check for combination ', ' and replace by ',' 
   !
   do while( .true. )
     indx = index(string,', ')
     if( indx > 0 )then
       tmpstring = string(indx+2:len(string)) 
       string(indx:) = ','//tmpstring
     else
       exit
     endif 
   end do
   !
   ! test for leading space
   !
   if( string(1:1) == ' ' )then
     tmpstring = string(2:len(string))
     string(1:) = tmpstring
   endif

   !
   ! replace spaces by komma's
   !
   ilen = lens(string)
   do i=1,ilen
     if( string(i:i) == ' ' ) string(i:i) =','
   end do

   !
   ! the line is filled with komma's as separators
   !
   ! count keywords
   !
   ilen = lens(string)
   icnt = 1
   do i=1,ilen
     if( string(i:i) == ',' ) icnt = icnt + 1
   end do
   Nkeys = icnt
  !write(*,*)'p5: number of keys:',Nkeys
      
   !
   ! fill in the strings
   !
   blank = '                '
   do i=1,MaxKeys
     keys(i) = blank
     keyi(i) =  0
     keyr(i) = 0.0
   end do
   
   i0 = 1
   i1 = 1
   do i=1,Nkeys
     indx    = index(string(i0:ilen),',')
     i1      = indx - 1
     !write(*,*) 'ip:',i,i0,i1
     if( i1 > 0 .and. i < Nkeys )then
       keyu(i) = .true.
       keys(i) = string(i0:i0+i1-1)
       i0      = i0 + i1 + 1
     elseif( i1 < 0 .and. i == Nkeys )then
       keyu(i) = .true.
       keys(i) = string(i0:ilen)
     elseif( i1 == 0 )then
       ! placeholder
       i0 = i0 + 1
       !write(*,*)'px:',i,'- - -',i0,i1
       keyu(i) = .false.
     else
       write(*,*)'+ internal error parser'
     endif
     !write(*,*)'px:',i,keys(i),i0,i1
   end do
  
   do i=1,Nkeys
     if( keyu(i) )read(keys(i),*,err=1) keyi(i) 
 1   continue
     !write(*,*)'pi:',i,keyi(i)
   end do 

   do i=1,Nkeys
     if( keyu(i) )read(keys(i),*,err=2) keyr(i)
 2   continue
     !write(*,*)'pr:',i,keyr(i)
   end do 

   !
   ! check for variables and insert them from library
   ! (set by 'set' command)
   !
  !write(*,*)'>>>',string(1:lens(string))
   do i=1,Nkeys
     if( keys(i)(1:1) == '$' )then
       name    = keys(i)(2:9)
       call getvariable(name,r1,r2)
       keyu(i) = .true.
       keyr(i) = r1
       keyi(i) = int(r1)     
       !write(*,*)'ps:',i,keys(i),'>',r1
     endif
   end do 
    
   !
   ! check for simple math (*/+-)
   !
   icnt = 0
   do i=2,Nkeys
     if( keys(i)(1:4) == '-   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '+   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '/   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '*   ' ) icnt = icnt + 1
     if( keys(i)(1:4) == '**  ' ) icnt = icnt + 1
   end do
   
  !if( icnt > 0 ) write(*,*) 'Math mode',icnt
   if( icnt > 0 ) math = .true.
   
   do while( icnt > 0 )

     do i=2,Nkeys-1
       found = .false. 
       if( keys(i)(1:4) == '-   ' )then
         r1 = keyr(i-1) - keyr(i+1)
         i1 = keyi(i-1) - keyi(i+1)       
         found = .true.
       else if( keys(i)(1:4) == '+   ' )then
         r1 = keyr(i-1) + keyr(i+1)
         i1 = keyi(i-1) + keyi(i+1)
         found = .true.
       else if( keys(i)(1:4) == '*   ' )then
         r1 = keyr(i-1) * keyr(i+1)
         i1 = keyi(i-1) * keyi(i+1)
         found = .true.
       else if( keys(i)(1:4) == '/   ' )then
         if( keyr(i+1) /= 0.0 ) r1 = keyr(i-1) / keyr(i+1)
         if( keyi(i+1) /=  0  )i1 = keyi(i-1) / keyi(i+1)
         found = .true.
       else if( keys(i)(1:4) == '**  ' )then
         r1 = keyr(i-1) ** keyr(i+1)
         i1 = keyi(i-1) ** keyi(i+1)
         found = .true.
       endif

       if( found )then
           
         keyu(i-1) = .true.
         keyr(i-1) = r1
         keyi(i-1) = i1 
         keys(i-1) = '(var)'

         do j=i,Nkeys-2
           keyu(j) = keyu(j+2)
           keyr(j) = keyr(j+2)
           keyi(j) = keyi(j+2)
           keys(j) = keys(j+2)
         end do

         keyu(Nkeys-1:MaxKeys) = .false.
         keyr(Nkeys-1:MaxKeys) =  0.0
         keyi(Nkeys-1:MaxKeys) =   0
         keys(Nkeys-1:MaxKeys) =  blank

         Nkeys = Nkeys - 2
         icnt  = icnt - 1 
         exit  
       endif 
     end do
   end do 

   !if( math )then
   !  write(*,*)'String:',string(1:lens(string))
   !  write(*,*)'New Nkeys:',Nkeys
   !  write(*,*)'keys:',keys(1:Nkeys)(1:4)
   !  write(*,*)'ints:',keyi(1:Nkeys)
   !  write(*,*)'rals:',keyr(1:Nkeys)
   !  write(*,*)'used:',keyu(1:Nkeys)
   !endif


   !write(*,*)'=== parser'

   
end subroutine parser
subroutine getkeyword(string,key,ikey)

   character(len=*) :: string, key
   character(len=1) :: s
   
   do i=1,len(string)
     s = string(i:i)
     if( s .eq. ' ' .or. s .eq. ',' .or. s .eq. '=') exit
     key(i:i) = string(i:i)
   end do
   ikey = i
      
end subroutine getkeyword
subroutine check_minmax(variable,varmin,varmax)

   if( variable < varmin )then
     write(*,*) '+ error: invalid entry'
     variable = varmin
   endif

   if( variable > varmax )then
     write(*,*) '+ error: invalid entry'
     variable = varmax
   endif

end subroutine check_minmax
logical function LGenerator()

   use constants
   character(len=120) string

   read(IOinp,'(A)',end=1) string   

   if( string(1:4) == 'gene' )then   
     !write(*,*)'Generator'
     backspace(IOinp)
     LGenerator = .true.
   else
     !write(*,*)'No generator'
     backspace(IOinp)
     LGenerator = .false.
   endif
   
   return
   
 1 continue
   write(*,*)'+ read error: premature eof (1)'
   stop
   
end function LGenerator
subroutine GetGenerators(ngen,ideltas,rdeltas)

   use constants
   
   character(len=120) string

   integer, parameter :: MaxKeys = 20
   character(len=32)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer            :: ideltas(MaxKeys)
   real               :: rdeltas(MaxKeys)

   read(IOinp,'(A)',end=1) string    

   if( string(1:4) == 'gene' )then   
     !write(*,*)'Generator (2)'

     indx = index(string,'#')
     if( indx >= 2 )then               ! strip trailing comments
       do i=indx,len(string)
         string(i:i) = ' '
       end do
     endif

     call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
 
     ngen    = keyi(2)
     if( Nkeys > 2 )then
       ideltas(1:Nkeys-2) = keyi(3:Nkeys) 
       rdeltas(1:Nkeys-2) = keyr(3:Nkeys) 
     else
       ideltas =  0 
       rdeltas = 0.0
     endif
   else
     write(*,*)'Internal error'
   endif

   return
   
 1 continue
   write(*,*)'+ read error: premature eof (2)'
   stop
   
end subroutine GetGenerators
subroutine dinlibrary

  use constants, only: IOdbg

  integer, parameter :: MaxDinVar = 80
  integer, save      :: NDinVar   =  0

  character(len=8)   :: variable
  logical            :: found

  character(len=8),dimension(MaxDinVar), save :: name
  real, dimension(MaxDinVar), save            :: rval1, rval2

  !=====================================================================
  entry setupvariables
  
  NDinVar = 0

  name(1)  = 'pi      '
  rval1(1) = 3.14159265
  rval2(1) = 0.0

  name(2)  = 'e       '
  rval1(2) = 2.71828183
  rval2(2) = 0.0

  do i=3,MaxDinVar
    name(i)  = '        '
    rval1(i) = 0.0
    rval2(i) = 0.0
  end do

  write(IOdbg,*)'DinLibrary variables set up:',MaxDinVar

  return

  !=====================================================================
  entry setvariable(variable,r1,r2)

  write(IOdbg,*)'set ',variable,r1,r2

  call lowercase(variable)

  do i=1,NDinVar
    if( name(i) == variable )then
      rval1(i) = r1
      rval2(i) = r2
      
      write(IOdbg,*)'DinLibrary modified:',variable,r1,r2      
      return
    end if
  end do
  
  if( NDinVar == MaxDinVar )then
    write(*,*)'+ error maximum variables reached'
    return
  endif
 
  
  NDinVar = NDinVar + 1
  name(NDinVar) = variable
  rval1(NDinVar) = r1
  rval2(NDinVar) = r2

  write(IOdbg,*)'DinLibrary stored:',variable,r1,r2
   
  return
  
  !=====================================================================
  entry getvariable(variable,r1,r2)

  found = .false.

  call lowercase(variable)

  do i=1,MaxDinVar
    if( variable == name(i) )then
      r1 = rval1(i)
      r2 = rval2(i)
      found = .true.
      exit
    endif
  end do

  if( .not. found )then
    write(*,*) '+ error: variable not found/set: ',variable
    r1 = 0.0
    r2 = 0.0
    stop
  endif 

  return
  
end subroutine dinlibrary
subroutine mathvariable(funct,mode,arg1,arg2,arg3,arg4)

   character(len=*), intent(IN)    :: funct
   integer, intent(IN)             :: mode

   real, intent(INOUT)             :: arg1
   real, intent(IN)                :: arg2
   real, intent(IN)                :: arg3
   real, intent(IN)                :: arg4
   
   real, parameter                 :: PI = 3.14159265
   real, parameter                 :: E  = 2.71828183 

   select case (funct) 

     case ('cos')
       if( mode == 0 )then
         arg1 = cos(arg1 * PI/180.)
       else if( mode == 1 )then
         arg1 = cos(arg1)
       endif
     case ('sin')
       if( mode == 0 )then
         arg1 = sin(arg1 * PI/180.)
       else if( mode == 1 )then
         arg1 = sin(arg1)
       endif
     case ('tan')
       if( mode == 0 )then
         arg1 = tan(arg1 * PI/180.)
       else if( mode == 1 )then
         arg1 = tan(arg1)
       endif
     case ('abs')
       arg1 = abs(arg1)
     case ('exp')
       arg1 = exp(arg1)
     case ('ln')
       arg1 = log(arg1)
     case ('log10')
       arg1 = log10(arg1)
     case ('sqrt')
       arg1 = sqrt(arg1)

     case default

       write(*,*)'Unknown math subcommand:',funct,mode
       arg1 = 0.0

   end select

end subroutine mathvariable
subroutine PrettyPrint(IO)

   use constants
   use geometry
   use variables
   use scalars

   use PatchesModule

   logical            :: warning
   integer            :: datum(8)
   character (len=12) :: clock(3)

   call date_and_time(clock(1),clock(2),clock(3),datum)

   write(IO,'(/)')
   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A,A)')  'Case                 : ',Casename(1:lens(Casename))
   write(IO,'(1x,A,A)')  'Title                : ',title(1:lens(title))
   write(IO,'(1x,A,i2.2,A1,i2.2,A1,i4,A2,i2.2,A1,i2.2,A1)')   &
       'Date and time        : ',datum(3),'/',datum(2),'/',datum(1), &
       ' (',datum(5),':',datum(6),')'  
   write(IO,'(1x,A,i10)') 'Number of cells      : ',Ncel
   write(IO,'(1x,A,i10)') 'Number of vertices   : ',Nvrt
   write(IO,'(1x,A,i10)') 'Number of boundaries : ',Nbnd
   write(IO,'(1x,A,i10)') 'Number of faces      : ',Nfac
   write(IO,'(1x,A,i10)') 'Number of regions    : ',Nreg
   write(IO,'(1x,A,i10)') 'Number of scalars    : ',Nscal
   
   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A,i10)') 'Number of iterations : ',Niter
   if( Transient )then
     if( Euler )then
       write(IO,'(1x,A,1pe10.3,A)') 'Transient run, dt    : ',dt,' (Euler)'
     else
       write(IO,'(1x,A,1pe10.3,A)') 'Transient run, dt    : ',dt,' (Quadratic)'
       write(IO,'(1x,A,1pe10.3)')   'Euler blending factor: ',GammaTime
     endif
   else
     write(IO,'(1x,A,1pe10.3)')   'Steady state run'
   endif
   if( UsePatches )then
   write(IO,'(1x,A,1pe10.3,A)') '*** Using patches'
   endif
   write(IO,'(1x,A,1pe10.3)')   'Target residual tol. : ',ResMax
   write(IO,'(1x,A,i10)')       'Monitor cell         : ',Imoni
   write(IO,'(1x,A,i10)')       'Pressure ref. cell   : ',IPref
   write(IO,'(1x,A,1pe10.3,A)') 'Reference pressure   : ',Pref,' Pa'
   if( SolveEnthalpy )then
   write(IO,'(1x,A,i10)')       'Temperature ref. cell: ',IPref
   write(IO,'(1x,A,1pe10.3,A)') 'Reference Temperature: ',Tref,' K'
   endif
   write(IO,'(1x,A,1pe10.3,A)') 'Reference density    : ',DensRef,' kg/m3'
   write(IO,'(1x,A,1pe10.3,A)') 'Molecular viscosity  : ',VisLam,' Pa s'
   if( SolveEnthalpy )then
   write(IO,'(1x,A,1pe10.3,A)') 'Prandtl number       : ',Prandtl,' -'
   write(IO,'(1x,A,1pe10.3,A)') 'Specific heat coef.  : ',CpStd,' J/kgK'
   write(IO,'(1x,A,1pe10.3,A)') 'Conductivity         : ',Lambda,' W/mK'
   endif
   write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Gravity vector       : ',Gravity,' '

   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A,1X,A)')            'Turbulence Model     : ',TMnames(Turbmodel)

   write(IO,'(1x,A)')    '----------------------------------------------------------'
   if( GradAlg == GradGauss )then
     write(IO,'(1x,A,i2,A)') 'Gradients using Gauss with ',Ngradient,' passes'
   else
     write(IO,'(1x,A)') 'Gradients using Least Squares'
   endif
   if( DScdDefault == DScd1 )then
     write(IO,'(1x,A)') 'Default CDS scheme: distance weighted scheme (1)'
   elseif( DScdDefault == DScd2 )then
     write(IO,'(1x,A)') 'Default CDS scheme: gradient based scheme (2)'
   elseif( DScdDefault == DScd3 )then
     write(IO,'(1x,A)') 'Default CDS scheme: simple averaged scheme (3)'
   endif
   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A)') 'Var  SCV sch blend gr slp    urf        rtol       guess'
   write(IO,1) 'U  : ',&
     SolveU ,PostC(VarU) ,PostV(VarU) ,DSch(Scheme(VarU)),Gamma(VarU), &
     GradID(GradVar(VarU)),SlopeID(LimiterSlope(VarU)),                &
     URF(VarU) ,RTOL(VarU) ,Guess(VarU) 
   write(IO,1) 'V  : ',&
     SolveV ,PostC(VarV) ,PostV(VarV) ,DSch(Scheme(VarV)),Gamma(VarV), &
     GradID(GradVar(VarV)),SlopeID(LimiterSlope(VarV)),&
     URF(VarV) ,RTOL(VarV) ,Guess(VarV) 
   write(IO,1) 'W  : ',&
     SolveW ,PostC(VarW) ,PostV(VarW) ,DSch(Scheme(VarW)),Gamma(VarW), &
     GradID(GradVar(VarW)),SlopeID(LimiterSlope(VarW)),&
     URF(VarW) ,RTOL(VarW) ,Guess(VarW) 
   write(IO,1) 'P  : ',&
     SolveP ,PostC(VarP) ,PostV(VarP) ,'   ',Gamma(VarP),              &
     GradID(GradVar(VarP)),SlopeID(LimiterSlope(VarP)),&
     URF(VarP) ,RTOL(VarP) ,Guess(VarP) 
   write(IO,1) 'TE : ',&
     SolveTurbEnergy,PostC(VarTE),PostV(VarTE),DSch(Scheme(VarTE)),Gamma(VarTE), &
     GradID(GradVar(VarTE)),SlopeID(LimiterSlope(VarTE)),&
     URF(VarTE),RTOL(VarTE),Guess(VarTE)
   write(IO,1) 'ED : ',&
     SolveTurbDiss,PostC(VarED),PostV(VarED),DSch(Scheme(VarED)),Gamma(VarED),   &
     GradID(GradVar(VarED)),SlopeID(LimiterSlope(VarED)),&
     URF(VarED),RTOL(VarED),Guess(VarED)
   write(IO,1) 'T  : ',&
     SolveEnthalpy,PostC(VarT) ,PostV(VarT) ,DSch(Scheme(VarT)),Gamma(VarT),     &
     GradID(GradVar(VarT)),SlopeID(LimiterSlope(VarT)),&
     URF(VarT) ,RTOL(VarT) ,Guess(VarT) 
   write(IO,1) 'SC : ',&
     SolveScalars,PostC(VarSC),PostV(VarSC),DSch(Scheme(VarSC)),Gamma(VarSC),    &
     '--','---',& 
     URF(VarSC),RTOL(VarSC),Guess(VarSC) 
   if (SolveScalars) then
     do i = 1,Nscal
        write(IO,'(1x,A4,A,4L1,1x,4(1pe10.3,1x),A)') Variable(Nvar+i),' : ',     & ! format niet OK!!!
          Solve(Nvar+i),Store(Nvar+i),PostC(Nvar+i),PostV(Nvar+i),Gamma(Nvar+i), &
          GradID(GradVar(VarT)),SlopeID(LimiterSlope(VarT)), &
          URF(Nvar+i),RTOL(Nvar+i),Guess(Nvar+i)
     end do
   endif
 1 format(1x,A,3L1,1x,A3,1x,f5.3,1x,A2,1x,A3,1x,3(1pe10.3,1x),A)

   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A,3x, 1x,f5.3,A)') 'PP : ',URF(VarPP), &
            ' (pressure correction pp under-relaxation'
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'DEN: ',&
     SolveDensity,PostC(VarDEN),PostV(VarDEN)
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'VIS: ',&
     SolveVISC,PostC(VarVIS),PostV(VarVIS)
   write(IO,'(1x,A,3L1,1x,A3,1x,4(1pe10.3,1x),A)') 'CP : ',&
     SolveCP,PostC(VarCP),PostV(VarCP)

   if( Transient )then
     warning = .false.
     if( URF(VarU)  <= 0.99 ) warning = .true.
     if( URF(VarTE) <= 0.99 ) warning = .true.
     if( URF(VarT)  <= 0.99 ) warning = .true.
     if( URF(VarSC) <= 0.99 ) warning = .true.
     if( warning ) write(IO,'(/,1x,A)') &
        '  * Transient run with underrelaxtion factors < 1 '
   endif     

   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A)')    '--- boundary conditions ----------------------------------'
   write(IO,'(1x,A)')    '----------------------------------------------------------'
   write(IO,'(1x,A,i3)') 'Number of regions ',NReg+1
   do i=0,NReg
     
     if( Reg(i)%n == 0 ) cycle
     
     if( Reg(i)%typ > 0 )then
       write(IO,'(1x,A,i3,A)')      '-  r e g i o n ',i,' -'
       write(IO,'(1x,A,A,i3)')      'Name         :  ',Reg(i)%name
       write(IO,'(1x,A,A,i3)')      'Type         :  ',Region(reg(i)%typ)
       if( reg(i)%typ == RWall )then
         if( reg(i)%noslip )then
           if( reg(i)%std )then
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'U, V, W      : ',Reg(i)%uvw 
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'E log        : ',Reg(i)%elog
           else
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'U, V, W      : ',Reg(i)%uvw 
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Rough. z0, Cs: ',Reg(i)%z0,Reg(i)%cs
           endif             
           if( Reg(i)%adiab )then
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Adiabatic      ' 
           elseif( Reg(i)%flux )then
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'T flux, R    : ',Reg(i)%T,Reg(i)%R           
           else
             write(IO,'(1x,A,3(1pe10.3,1x),A)') 'T fixed, R   : ',Reg(i)%T,Reg(i)%R
           endif
         else
           write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Free slip wall'
         endif
       elseif( reg(i)%typ == RInlet )then
         write(IO,'(1x,A,3(1pe10.3,1x),L1)')'U, V, W      : ',Reg(i)%uvw,Reg(i)%user
         write(IO,'(1x,A,3(1pe10.3,1x),A)') 'Density, T   : ',Reg(i)%den,Reg(i)%T
         write(IO,'(1x,A,3(1pe10.3,1x),A)') 'k, e         : ',Reg(i)%k,Reg(i)%e
       elseif( reg(i)%typ == ROutlet )then
         write(IO,'(1x,A,3(1pe10.3,1x),L1)')'Split        : ',Reg(i)%splvl
       else
       endif
     endif
   end do

   if( UsePatches )then
     write(IO,'(1x,A)')    '----------------------------------------------------------'
     write(IO,'(1x,A)')    '--- patches ----------------------------------------------'
     write(IO,'(1x,A)')    '----------------------------------------------------------'
     write(IO,'(1x,A,i3)') 'Number of patches',NPatches
   
     do i=1,NPatches
       write(IO,'(1x,A,i3,A)')      '- p a t c h ',i,' -'
       write(IO,'(1x,A,A,A)')       'Name         : ',Patch(i)%sName
       write(IO,'(1x,A,i3)')        'Type         : ',Patch(i)%itype
       write(IO,'(1x,A,A,i3)')      'Variable     : ',Patch(i)%sVar,Patch(i)%iVar
       write(IO,'(1x,A,i3)')              'Cell ID      : ',Patch(i)%iCellID
       write(IO,'(1x,A,3(1pe10.3,1x),A)') 'x,y,z low    : ',Patch(i)%xyz1
       write(IO,'(1x,A,3(1pe10.3,1x),A)') 'x,y,z high   : ',Patch(i)%xyz2
       write(IO,'(1x,A,i3)')        'Data length  : ',Patch(i)%iDatLen
       do j=1,Patch(i)%iDatLen
         write(IO,'(1x,A,i3,A,1pe10.3)')  '   data ',j,'  : ',Patch(i)%dat(j)
       end do
     end do
   endif

   if( SolveScalars )then
     write(IO,'(1x,A)')    '----------------------------------------------------------'
     write(IO,'(1x,A)')    '--- scalars ----------------------------------------------'
     write(IO,'(1x,A)')    '----------------------------------------------------------'
     write(IO,'(1x,A,i3)') 'Number of scalars',NScal
     do i=1,NScal
       write(IO,'(1x,A,i3,A)')      '- s c a l a r ',i,' -'
       write(IO,'(1x,A,A,A)')       'Name         :  ',ScProp(i)%name,ScProp(i)%unit
       write(IO,'(1x,A,3L1)')       'Act/Sol/Sav  :  ',ScProp(i)%active,ScProp(i)%solve,ScProp(i)%save
       write(IO,'(1x,A,1pe10.3,A)') 'Lam. Prandtl : ',ScProp(i)%PrL,' -'
       write(IO,'(1x,A,1pe10.3,A)') 'Turb. Prandtl: ',ScProp(i)%PrT,' -'
     end do
     
     do i=1,NReg
       write(IO,'(1x,A,i3,A)')      '- r e g i o n ',i,' -'
       do j=1,NScal
         write(IO,'(i3,1x,A,3(1pe10.3,1x))')  j,ScProp(j)%name,ScReg(i,j)%fraction,ScReg(i,j)%value
       end do
     end do
   endif
   write(IO,'(1x,A)')    '----------------------------------------------------------'
   
end subroutine PrettyPrint
subroutine ReadCase

   use constants
   use geometry
   use variables

   logical :: Exists = .false.

   
   write(IOdef,'(1x,A,f5.3)') 'This is dolfyn version ',float(Version)*0.001
   write(IOdef,'(1x,A)') 'Copyright(C) 2002-2010 Cyclone Fluid Dynamics BV'
   write(IOdef,'(1x,A)') 'NL-5583 XM, Waalre, The Netherlands'
   write(IOdef,'(1x,A)') 'see http://www.cyclone.nl and http://www.dolfyn.net'
   write(IOdef,'(1x,A)') 'Using Sparsekit2 by Yousef Saad'
   write(IOdef,'(1x,A)') '(C) 2005, the Regents of the University of Minnesota'
   write(IOdef,'(1x,A)') 'Modules with patches (C) 2004-2010 by B. Tuinstra'
   write(IOdef,'(1x,A)') 'see http://www.home.zonnet.nl/bouke_1/dolfyn'
   write(IOdef,'(1x,A)') 'Tecplot interface (C) 2006 by S.B. Kuang'
   write(IOdef,'(1x,A)') 'VTK interface updated (C) 2007 by J. Jacobs'
    write(IOdef,'(/)')  
   

   inquire(file='dolfyn.cfg',exist=exists)
   if( exists )then
     call openfile(IOcfg,'dolfyn','.cfg','FORMATTED', &
                         'SEQUENTIAL','OLD',debug)
     write(IOdef,*) 'Reading case from dolfyn.cfg'
     read(IOcfg,'(A)') Casename 
     close(IOcfg)                         
   else
     write(*,*) 'Enter case:'
     read(*,'(A)') Casename 
   endif
   write(IOdef,*) 'Using case: ',Casename(1:lens(Casename))

   !
   ! check if geometry and din files exists
   !
   ierr = 0
   inquire(file=Casename(1:lens(Casename))//'.din',exist=Exists)
   if( .not. Exists )then
     write(*,*)'*** Error: Unknown File ',Casename(1:lens(Casename))//'.din'
     ierr = ierr + 1
   endif

   inquire(file=Casename(1:lens(Casename))//'.dge',exist=Exists)
   if( .not. Exists )then
     !
     ! check if old/pre-gmsh geometry file is used
     !
     write(*,*)'File ',Casename(1:lens(Casename))//'.dge',' unknown'
     write(*,*)'Checking for ',Casename(1:lens(Casename))//'.geo'
     inquire(file=Casename(1:lens(Casename))//'.geo',exist=Exists)
     if( .not. Exists )then
       write(*,*)'*** Error: File ',Casename(1:lens(Casename))//'.geo',' unknown too' 
       ierr = ierr + 1
     else
       Ext(IOgeo) = '.geo'       
     endif
   else
     Ext(IOgeo) = '.dge'  
   endif

   if( ierr > 0 ) stop
   !
   ! check if STOP-files are present; if so remove it or them
   ! (note: IOcfg used again)
   !
   inquire(file='STOP',exist=Exists)
   if( Exists )then
     call openfile(IOcfg,'STOP','','FORMATTED', &
                         'SEQUENTIAL','OLD',debug)
     write(*,*)'Removing STOP'
     close(IOcfg,status='DELETE') 
   endif

   inquire(file=casename(1:lens(casename))//'.STOP',exist=Exists)
   if( Exists )then
     call openfile(IOcfg,casename(1:lens(casename))//'.STOP','', &
                         'FORMATTED','SEQUENTIAL','OLD',debug)
     write(*,*)'Removing ',casename(1:lens(casename))//'.STOP'
     close(IOcfg,status='DELETE') 
   endif
   !
   ! echo important/relevant console output to a file as well
   !
   inquire(file=Casename(1:lens(Casename))//'.txt',exist=Exists)
   if( Exists )then
     write(*,*)'Appending echo to file ',Casename(1:lens(Casename))//'.txt'
   endif

   call openfile(IOrun,Casename,'.txt','FORMATTED', &
                       'APPEND','UNKNOWN',debug)

   write(IOrun,'(1x,A,f5.3)') 'This is dolfyn version ',float(Version)*0.001
   write(IOrun,'(1x,A)') 'Copyright(C) 2002-2010 Cyclone Fluid Dynamics BV'
   write(IOrun,'(1x,A)') 'NL-5583 XM, Waalre, The Netherlands'
   write(IOrun,'(1x,A)') 'see http://www.cyclone.nl and http://www.dolfyn.net'
   write(IOrun,'(1x,A)') 'Using Sparsekit2 by Yousef Saad'
   write(IOrun,'(1x,A)') '(C) 2005, the Regents of the University of Minnesota'
   write(IOrun,'(1x,A)') 'Modules with patches (C) 2004-2010 by B. Tuinstra'
   write(IOrun,'(1x,A)') 'see http://www.home.zonnet.nl/bouke_1/dolfyn'
   write(IOrun,'(1x,A)') 'Tecplot interface (C) 2006 by S.B. Kuang'
   write(IOrun,'(1x,A)') 'VTK interface updated (C) 2007 by J. Jacobs'
   write(IOrun,'(/)')  
   write(IOrun,*) 'Using case: ',Casename(1:lens(Casename))

end subroutine ReadCase


