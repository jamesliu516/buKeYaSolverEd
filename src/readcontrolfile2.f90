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
subroutine ReadControlFile2

   use constants
   use geometry
   use variables
   use scalars
   use particles
   
   character(len=120) string

   character(len=32)  :: key, key2, key3, key4, key5
   character(len=32)  :: blank = '                                '
   character(len=32)  :: str1, str2

   integer, parameter :: MaxKeys = 20
   character(len=32)  :: keys(MaxKeys) 
   integer            :: keyi(MaxKeys)
   real               :: keyr(MaxKeys)
   logical            :: keyu(MaxKeys)

   integer            :: ideltas(MaxKeys)
   real               :: rdeltas(MaxKeys)

   logical :: Messages 
   logical :: Echo     
   logical :: InputError    

   !
   !========================================================
   !
   entry select_case_echo(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
   
     key = keys(2)
     call lowercase(key)

     select case (key(1:3))
       case ('no')
         if( Echo ) write(*,'(1x,A,i3)')'+ echo input off'
         Echo = .false.
       case ('off')
         if( Echo ) write(*,'(1x,A,i3)')'+ echo input off'
         Echo = .false.
       case ('on')
         if( Echo ) write(*,'(1x,A,i3)')'+ echo input on'
         Echo = .true.
       case default
         if( Echo ) write(*,'(1x,A,i3)')'+ echo input on'
         Echo = .true.
     end select
   
   return

   !
   !========================================================
   !
   entry select_case_save(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
     !
     ! save restart data
     !
     key = keys(2)
     call lowercase(key)

     select case (key(1:4))
       case ('ever')
         i = keyi(3)
         if( i < 0 )then
           write(*,*)'+++ Error: invalid entry in "save,every,":',i
           messages = .true.
         else
           if( Echo ) write(*,'(1x,A,i3)')'+ saving every ',i,' steps'
           NSave = i
         endif
       case ('time')
         r = keyr(3)
         if( r < 0.0 )then
           write(*,*)'+++ Error: invalid entry in "save,time,":',r
           messages = .true.
         else
           if( Echo ) write(*,'(1x,A,1pe10.3)')'+ saving at time =',r
           TSave = r
         endif
       case ('iter')
         i = keyi(3)
         if( i < 0 )then
           write(*,*)'+++ Error: invalid entry in "save,iteration,":',i
           messages = .true.
         else
           if( Echo ) write(*,'(1x,A,i3)')'+ saving iteration ',i
           ISave = i
         endif
       case ('cpu')
         r = keyr(3)
         key2 = keys(4)
         call lowercase(key2)
         key3 = keys(5)
         call lowercase(key3)

         select case (key2(1:1))
           case ('s')
             fact =    1.0
           case ('m')
             fact =   60.0
           case ('h')
             fact = 3600.0
           case default
             fact =   1.0
         end select

         select case (key3(1:1))
           case ('s')
             CPUStop = .true.
           case default
             CPUStop = .false.
         end select             
         
         if( r < 0.0 )then
           write(*,*)'+++ Error: invalid entry in "save,cpu,":',r
           messages = .true.
         else
           r = fact * r
           if( Echo )  write(*,'(1x,A,f6.0)') &
                     '+ saving after elapsed CPU time seconds:',r
           TCPU = r
         endif
       case default
         write(*,*)'+++ Error: unknown option in save command.'
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_output(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
     !
     ! write output data...
     !
     key = keys(2)
     call lowercase(key)

     select case (key(1:4))
       case ('ever')
         i = keyi(3)
         if( i < 0 )then
           write(*,*)'+++ Error: invalid entry in "output,every,":',i
           messages = .true.
         else
           if( Echo ) write(*,'(1x,A,i3)')'+ output every ',i,' steps'
           NOutput = i
         endif
       case ('time')
         r = keyr(3)
         if( r < 0.0 )then
           write(*,*)'+++ Error: invalid entry in "output,time,":',r
           messages = .true.
         else
           if( Echo ) write(*,'(1x,A,1pe10.3)')'+ output at time =',r
           TOutput = r
         endif
       case ('iter')
         i = keyi(3)
         if( i < 0 )then
           write(*,*)'+++ Error: invalid entry in "output,iteration,":',i
           messages = .true.
         else
           if( Echo ) write(*,'(1x,A,i3)')'+ output iteration ',i
           IOutput = i
         endif
       case default
         write(*,*)'+++ Error: unknown option in output command.'
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_boundary(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
     
     InputError = .false.
     !
     ! boundaries section
     !
     ! first check if a name has been used instead of a number
     !
     ik = -1
     do i=0,Nreg
       str1 = Reg(i)%name
       call lowercase(str1)
       str2 = keys(2)
       call lowercase(str2)
       if( str1 == str2 )then
         ik = i
	 exit
       endif 
     end do

     if( ik /= -1 )then    !
       i = ik              ! use the name found
     else                  !
       i = keyi(2)         ! or directly the region number
     endif                 !
     
     if( Echo ) write(*,'(1x,A,i3)')'+ boundary category = ',i
     if( i < 0 .or. i > Nreg )then
       write(*,*)'Boundary out of range. Skipped.'
       Messages   = .true.
       InputError = .true.
     endif
     
     if( .not. InputError )then
       read(IOinp,'(A)',end=20) string
       key = blank
       call getkeyword(string,key,jkey)
       call lowercase(key)
       if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)

       select case (key(1:3))
         case ('sym')
           write(*,'(1x,A,i3,A)') 'Boundary ',i,' sym. conditions'
           Reg(i)%typ = RSymp
         case ('inl')
           !
           ! inlet
           !
           write(*,'(1x,A,i3,A)') 'Boundary ',i,' inlet conditions'
           Reg(i)%typ = RInlet
           if( keys(3) == 'user' ) Reg(i)%user  = .true. 
           if( keys(4) /= '    ' ) Reg(i)%name1 = keys(4)
           if( keys(5) /= '    ' ) Reg(i)%name2 = keys(5)
           if( Reg(i)%user ) write(*,*) '  *** User inlet ',&
                                        Reg(i)%name1,' ',Reg(i)%name2

           read(IOinp,'(A)',end=20) string    
           call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
           Uin = keyr(1)
           Vin = keyr(2)
           Win = keyr(3)
           Reg(i)%uvw(1) = Uin 
           Reg(i)%uvw(2) = Vin 
           Reg(i)%uvw(3) = Win 
           if( Echo ) write(*,*) '  UVW:',Reg(i)%uvw,' m/s'

           read(IOinp,'(A)',end=20) string    
           call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
           Reg(i)%den = keyr(1)
           if( Echo ) write(*,*) '  DEN:',Reg(i)%den,' kg/m3'

           read(IOinp,'(A)',end=20) string    
           call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
           Reg(i)%T= keyr(1)
           if( Echo ) write(*,*) '  T  :',Reg(i)%T,' K'

           Vmag2 = Uin**2 + Vin**2 + Win**2

           read(IOinp,'(A)',end=20) string
           key = blank
           call getkeyword(string,key,jkey)
           call lowercase(key)
           if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:4)
           if( key(1:4) == 'keps' )then
             if( Echo ) write(*,*) '  Using k and epsilon'
             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             TKin  = keyr(1)
             EPSin = keyr(2)

             Reg(i)%k = TKin
             Reg(i)%e = EPSin
             if( TKin <= 0.0 )then
               write(*,*)'Error: k should be positive! Reseting.'
               ti = 0.01
               Reg(i)%k = 3./2.*( ti**2 * vmag2 ) 
               messages = .true.
             endif
             if( EPSin <= 0.0 )then
               write(*,*)'Error: epsilon should be positive! Reseting.'
               tl = 1.0
               Reg(i)%e = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / tl
               messages = .true.
             endif
             TKin  = sqrt( 2./3.*Reg(i)%k / (Vmag2 + Small) )
             EPSin = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / (Reg(i)%e + Small)
             write(*,*) '  TE :',Reg(i)%k,'( i',TKin,')'
             write(*,*) '  ED :',Reg(i)%e,'( l',EPSin,')'
           elseif( key(1:4) == 'inle' )then
             if( Echo ) write(*,*) '  Using intensity and length scale'
             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             TKin  = keyr(1)
             EPSin = keyr(2)
             if( TKin <= 0.0 )then
               write(*,*)'Error: intensity should be positive! Reseting.'
               ti   = 0.01
               TKin = ti
               Reg(i)%k = 3./2.*( ti**2 * vmag2 ) 
               messages = .true.
             else
               Reg(i)%k = 3./2.*( TKin**2 * vmag2 ) 
             endif
             if( EPSin <= 0.0 )then
               write(*,*)'Error: length scale should be positive! Reseting.'
               tl    = 1.0
               EPSin = tl
               Reg(i)%e = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / tl
               messages = .true.
             else
               Reg(i)%e = (TMCmu)**(0.75) * Reg(i)%k**(1.5) / EPSin
             endif
             if( Echo ) write(*,*) '  TE :',TKin,'( k',Reg(i)%k,')'
             if( Echo ) write(*,*) '  ED :',EPSin,'( e',Reg(i)%e,')'
           else
             write(*,*)'Error: [keps*|inlen] expected'
             messages = .true.
           endif
         case ('out')
           ! 
           ! outlet 
           !
           write(*,'(1x,A,i3,A)') 'Boundary ',i,' outlet conditions'
           Reg(i)%typ    = ROutlet
           read(IOinp,*,err=11) string
           call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
           SplitIn = keyr(1)
           Reg(i)%splvl  = SplitIn  

           Noutlets = Noutlets + 1
           Split    = Split + SplitIn

           !write(*,*)'>> ',Noutlets,SplitIn,Split
           if( (Split-1.0) > 1.e-6 )then
             write(*,*)'Error: Sum of splits greater than 1.0 ',Split
             messages = .true.
           endif
         case ('wal')
           !
           ! wall
           !
           write(*,'(1x,A,i3,A)') 'Boundary ',i,' wall conditions'
           Reg(i)%typ = RWall
           read(IOinp,'(A)',end=20) string
           key = blank
          !call getkeyword(string,key,jkey)
           call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
           key = keys(1)
           call lowercase(key)
           if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:4)
           if( key(1:4) == 'slip' )then
             if( Echo ) write(*,*) '  Free slip wall'
             Reg(i)%noslip = .false.
             Reg(i)%adiab  = .true.
             Reg(i)%flux   = .false.
             InputError = .true.
           else if( key(1:4) == 'nosl' )then
             Reg(i)%noslip = .true.
           else if( key(1:4) == 'roug' )then
             Reg(i)%noslip = .true.
             Reg(i)%std    = .false.                    !<= veranderen in %rough?
             r = keyr(2)
             if( r < 0.0 )then
               write(*,*)'Error: roughness should be greater than zero'
               messages = .true.
               r = 0.03
             else if( r == 0.0 )then
               write(*,*)'Default roughness used'
               r = 0.03
             endif  
             c = keyr(3)
             if( c < 0.0 .or. c > 1.0 )then
               write(*,*)'Error: roughness factor out of range',c
               messages = .true.
               c = 0.3
             else if( c == 0.0 )then
               write(*,*)'Default roughness factor used'
               c = 0.3
             endif  
             Reg(i)%z0 = r
             Reg(i)%cs = c
             if( Echo ) write(*,*) ' Rough wall z0 [m]:',r,cs
           else
             write(*,*) '*** Error: [noslip*|slip|rough] expected'
             write(*,*) '*** Assuming noslip wall'
             Reg(i)%noslip = .false.
             messages = .true.
           endif
           if( .not. InputError )then
             !
             ! noslip wall
             !
             if( Echo ) write(*,*) '  No slip wall'
             read(IOinp,'(A)',end=20) string    
             call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
             Reg(i)%uvw(1) = keyr(1) 
             Reg(i)%uvw(2) = keyr(2) 
             Reg(i)%uvw(3) = keyr(3) 
             if( Echo ) write(*,*) '  UVW:',Reg(i)%uvw

             read(IOinp,'(A)',end=20) string
             key = blank
             call getkeyword(string,key,jkey)
             call lowercase(key)
             if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:4)
             if( Echo ) write(*,*)'+ note: resistance R added after T or q!'
             if( key(1:4) == 'adia' )then
               if( Echo ) write(*,*) '  Adiabatic wall'
               Reg(i)%adiab  = .true.
               Reg(i)%flux   = .false.
             elseif( key(1:4) == 'fixe' )then
               Reg(i)%adiab  = .false.
               Reg(i)%flux   = .false.
               read(IOinp,'(A)',end=20) string    
               call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
               Reg(i)%T = keyr(1) 
               Reg(i)%R = keyr(2) 

               if( Echo ) write(*,*) '  T,R:',Reg(i)%T,Reg(i)%R
             elseif( key(1:4) == 'flux' )then
               Reg(i)%adiab  = .false.
               Reg(i)%flux   = .true.
               read(IOinp,'(A)',end=20) string    
               call parser(MaxKeys,Nkeys,keys,keyi,keyr,keyu,string)
               Reg(i)%T = keyr(1) 
               Reg(i)%R = keyr(2) 

               if( Echo ) write(*,*) '  T  :',Reg(i)%T,Reg(i)%R
             else
               write(*,*)'Error: [adiabatic*|fixed|flux] expected'
               messages = .true.
             endif
           endif
       end select
     endif
   
   return
11   continue
     write(*,*)'+ read error: something else expected'
     write(*,*)'+ trying to continue'
     messages = .true.
20   continue
   return
   !
   !========================================================
   !
   entry select_case_thermal(Echo,Messages,KeyS,KeyI,KeyR,KeyU)
   
     if( Echo ) write(*,'(1x,A,i3)')'+ thermal options'
     key = keys(2)
     call lowercase(key)
     if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)

     select case (key(1:3))
       case ('off')
         if( Echo ) write(*,'(1x,A,i3)')'+ thermal mode switched off'
         SolveEnthalpy = .false.
       case ('on')
         if( Echo ) write(*,'(1x,A,i3)')'+ thermal mode switched on'
         SolveEnthalpy = .true.
       case default
         SolveEnthalpy = .true.
     end select
   
   return
   !
   !========================================================
   !
   entry select_case_restart(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key = keys(2)
     call lowercase(key)
     if( Echo ) write(*,'(1x,A,A)')'+ keyword = ',key(1:3)
     select case (key(1:3))
       case ('no')
         if( Echo ) write(*,'(1x,A,i3)')'+ no restart file used'
         Restart       = 0
       case ('off')
         if( Echo ) write(*,'(1x,A,i3)')'+ no restart file used'
         Restart       = 0
       case ('ini')
         !
         ! if ReadRestartField discovers that Nbnd and Nfaces
         ! have changed Restart is set there to 4 (massfluxes
         ! will not be read in and set to zero)
         !
         if( Echo ) write(*,'(1x,A,i3)')'+ use old results as initial geuss'
         Restart       = 2
         key2 = keys(3)
         call lowercase(key2)
         if( key2(1:3) == 'cel' )then
           !
           ! when Nbnd and Nfaces do not change but regions
           ! are added or changed (then the order of the bnd-list 
           ! will be changed and cannot be tracked -unless the
           ! bnd-list is written to the rst-file as well-)
           ! use flux key as an option 
           !
           key3 = keys(4)
           call lowercase(key2)
           if( key3(1:3) == 'flu' )then
             write(*,*)'Using only cell data with mass fluxes for restart '
             Restart     = 5
           else
             write(*,*)'Using only cell data for restart (ignoring mass fluxes)'
             Restart     = 4
           endif
         endif             
       case ('res')
         if( Echo ) write(*,'(1x,A,i3)')'+ just reset counters'
         Restart       = 3
       case default
         if( Echo ) write(*,'(1x,A,i3)')'+ run restarted'
         Restart       = 1
     end select

   return
   !
   !========================================================
   !
   entry select_case_relax_star(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     URF(VarU ) = keyr(2)
     URF(VarV ) = keyr(3)
     URF(VarW ) = keyr(4)
     URF(VarP ) = keyr(5)
     URF(VarTE) = keyr(6)
     URF(VarED) = keyr(7)
     URF(VarT ) = keyr(8)
     URF(VarSC) = keyr(9)

     if( URF(VarU ) == 0.0 ) URF(VarU ) = 0.7   ! set defaults
     if( URF(VarV ) == 0.0 ) URF(VarV ) = 0.7
     if( URF(VarW ) == 0.0 ) URF(VarW ) = 0.7
     if( URF(VarP ) == 0.0 ) URF(VarP ) = 0.35
     if( URF(VarTE) == 0.0 ) URF(VarTE) = 0.7
     if( URF(VarED) == 0.0 ) URF(VarED) = 0.7
     if( URF(VarT ) == 0.0 ) URF(VarT ) = 0.9
     if( URF(VarSC) == 0.0 ) URF(VarSC) = 0.9

     call check_minmax(URF(VarU ),Small,1.0)
     call check_minmax(URF(VarV ),Small,1.0)
     call check_minmax(URF(VarW ),Small,1.0)
     call check_minmax(URF(VarP ),Small,1.0)
     call check_minmax(URF(VarTE),Small,1.0)
     call check_minmax(URF(VarED),Small,1.0)
     call check_minmax(URF(VarT ),Small,1.0)
     call check_minmax(URF(VarSC),Small,1.0)
     if( Echo )then
       write(*,'(1x,A,3(1x,f6.3))')'+ under relaxation factors (*) => '
       write(*,*) '+ uvwp: ',(URF(iurf),iurf=1,4)
       write(*,*) '+ kets: ',(URF(iurf),iurf=5,8)
     endif
   
   return
   !
   !========================================================
   !
   entry select_case_relax(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     URF(VarU ) = keyr(2)
     URF(VarV ) = URF(VarU)
     URF(VarW ) = URF(VarU)
     URF(VarP ) = keyr(3)
     URF(VarTE) = keyr(4)
     URF(VarED) = URF(VarTE)
     URF(VarT ) = keyr(5)
     URF(VarSC) = keyr(6)

     if( URF(VarU ) == 0.0 ) URF(VarU ) = 0.7   ! set defaults
     if( URF(VarV ) == 0.0 ) URF(VarV ) = 0.7
     if( URF(VarW ) == 0.0 ) URF(VarW ) = 0.7
     if( URF(VarP ) == 0.0 ) URF(VarP ) = 0.35
     if( URF(VarTE) == 0.0 ) URF(VarTE) = 0.7
     if( URF(VarED) == 0.0 ) URF(VarED) = 0.7
     if( URF(VarT ) == 0.0 ) URF(VarT ) = 0.9
     if( URF(VarSC) == 0.0 ) URF(VarSC) = 0.9
     
     call check_minmax(URF(VarU) ,Small,1.0)
     call check_minmax(URF(VarP) ,Small,1.0)
     call check_minmax(URF(VarTE),Small,1.0)
     call check_minmax(URF(VarT) ,Small,1.0)
     call check_minmax(URF(VarSC),Small,1.0)
     if( Echo ) &
       write(*,'(1x,A,5(1x,f6.3))')'+ under relaxation factors  = ',&
       URF(VarU),URF(VarP),URF(VarTE),URF(VarT),URF(VarSC)

   return
   !
   !========================================================
   !
   entry select_case_rtol(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     RTOL(VarU ) = keyr(2)
     RTOL(VarV ) = RTOL(VarU)
     RTOL(VarW ) = RTOL(VarU)
     RTOL(VarP ) = keyr(3)
     RTOL(VarTE) = keyr(4)
     RTOL(VarED) = RTOL(VarTE)
     RTOL(VarT ) = keyr(5)
     RTOL(VarSC) = keyr(6)         

     if( RTOL(VarU ) == 0.0 ) RTOL(VarU ) = 0.1   ! set defaults
     if( RTOL(VarV ) == 0.0 ) RTOL(VarV ) = 0.1
     if( RTOL(VarW ) == 0.0 ) RTOL(VarW ) = 0.1
     if( RTOL(VarP ) == 0.0 ) RTOL(VarP ) = 0.05
     if( RTOL(VarTE) == 0.0 ) RTOL(VarTE) = 0.1
     if( RTOL(VarED) == 0.0 ) RTOL(VarED) = 0.1
     if( RTOL(VarT ) == 0.0 ) RTOL(VarT ) = 0.1
     if( RTOL(VarSC) == 0.0 ) RTOL(VarSC) = 0.1

     call check_minmax(RTOL(VarU) ,Small,1.0)
     call check_minmax(RTOL(VarP) ,Small,1.0)
     call check_minmax(RTOL(VarTE),Small,1.0)
     call check_minmax(RTOL(VarT) ,Small,1.0)
     call check_minmax(RTOL(VarSC),Small,1.0)
     if( Echo ) &
       write(*,'(1x,A,5(1x,f6.3))')'+ relative solver tolerance = ',&
       RTOL(VarU),RTOL(VarP),RTOL(VarTE),RTOL(VarT),RTOL(VarSC)
     
   return
   !
   !========================================================
   !
   entry select_case_atol(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     ATOL(VarU ) = keyr(2)
     ATOL(VarV ) = ATOL(VarU)
     ATOL(VarW ) = ATOL(VarU)
     ATOL(VarP ) = keyr(3)
     ATOL(VarTE) = keyr(4)
     ATOL(VarED) = ATOL(5)
     ATOL(VarT ) = keyr(5)
     ATOL(VarSC) = keyr(6)       
     call check_minmax(ATOL(1),-Small,1.0)
     call check_minmax(ATOL(4),-Small,1.0)
     call check_minmax(ATOL(5),-Small,1.0)
     if( Echo ) &
       write(*,'(1x,A,5(1x,1pe7.1))')'+ absolute solver tolerance = ',&
       ATOL(VarU),ATOL(VarP),ATOL(VarTE),ATOL(VarT),ATOL(VarSC)
   
   return
   !
   !========================================================
   !
   entry select_case_init(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)
     if( Echo ) &
       write(*,'(1x,A,A,A)')'+ init keyword with strings = ',key2,key3
     select case (key2(1:5))
       case ('field')
         Guess(VarU ) = keyr(3)
         Guess(VarV ) = keyr(4)
         Guess(VarW ) = keyr(5)
         Guess(VarP ) = keyr(6)
         Guess(VarTE) = keyr(7)
         Guess(VarED) = keyr(8)
         if( keyr(9) > Small ) Guess(VarT ) = keyr(9)
         Guess(VarSC) = keyr(10)
         if( Echo ) write(*,'(1x,A,4(1x,1pe8.1))')'+ initial fields = ',Guess(1:4)
         if( Echo ) write(*,'(1x,A,4(1x,1pe8.1))')'                   ',Guess(5:8)
       case ('user')
         if( Echo ) write(*,'(1x,A)') '+ use user subroutine'
         UserInitialisation = .true.
       case ('fact')
         VisInitialFactor = keyr(3)
         if( Echo ) write(*,'(1x,A,1pe8.1)')'+ initial viscosity factor  = ',&
                    VisInitialFactor
       case ('steps')
         InitialSteps = keyi(3)
         if( Echo ) write(*,'(1x,A,i3)')'+ initial steps = ',InitialSteps
       case default
         write(*,'(1x,A,i3)')'+ error in initial command'
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_gamma_star(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     Gamma(VarU ) = keyr(2)                                      
     Gamma(VarV ) = keyr(3)                                      
     Gamma(VarW ) = keyr(4)                                      
     Gamma(VarP ) = keyr(5)                                      
     Gamma(VarTE) = keyr(6)                                      
     Gamma(VarED) = keyr(7)                                        
     Gamma(VarT ) = keyr(8)                                      
     Gamma(VarSC) = keyr(9)                                     

     if( Gamma(VarU ) > 0.0 ) Scheme(VarU ) = DScds
     if( Gamma(VarV ) > 0.0 ) Scheme(VarV ) = DScds
     if( Gamma(VarW ) > 0.0 ) Scheme(VarW ) = DScds

     if( Gamma(VarTE) > 0.0 ) Scheme(VarTE) = DScds
     if( Gamma(VarED) > 0.0 ) Scheme(VarED) = DScds
     if( Gamma(VarT ) > 0.0 ) Scheme(VarT ) = DScds
     if( Gamma(VarSC) > 0.0 ) Scheme(VarSC) = DScds

     if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'+ blending (*) = ',Gamma(1:4)    
     if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'                 ',Gamma(5:8)    

   return
   !
   !========================================================
   !
   entry select_case_gamma(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     Gamma(VarU ) = keyr(2)                                      
     Gamma(VarV ) = Gamma(VarU )                                      
     Gamma(VarW ) = Gamma(VarU )                                      
     if( Gamma(VarU) > 0.0 )then
       Scheme(VarU) = DScds
       Scheme(VarV) = DScds
       Scheme(VarW) = DScds
     endif

     Gamma(VarP ) = keyr(3)         ! not used!                                  

     Gamma(VarTE) = keyr(4)                                      
     Gamma(VarED) = Gamma(VarTE)                                      
     if( Gamma(VarTE) > 0.0 )then
       Scheme(VarTE) = DScds
       Scheme(VarED) = DScds
     endif

     Gamma(VarT ) = keyr(5)                                      
     if( Gamma(VarT) > 0.0 ) Scheme(VarT) = DScds

     Gamma(VarSC) = keyr(6)                                     
     if( Gamma(VarSC) > 0.0 ) Scheme(VarSC) = DScds

     if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'+ blending = ',Gamma(1:4)    
     if( Echo ) write(*,'(1x,A,4(1x,1pe7.1))')'             ',Gamma(5:8)    

   return
   !
   !========================================================
   !
   entry select_case_solver(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)

     select case (key2)
       case ('u')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving U component'
             SolveU      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ U component not solved'
             SolveU      = .false.
           case ('sparse')
             if( Echo ) write(*,'(1x,A,i3)')'+ U component SparseKit2'
             solver(VarU) = SparseKit2 
           case ('bicg')
             if( Echo ) write(*,'(1x,A,i3)')'+ U component BiCGstabL'
             solver(VarU) = SolverBiCGstabL 
           case ('direct')
             if( Echo ) write(*,'(1x,A,i3)')'+ U component Direct'
             solver(VarU) = SolverDirect 
           case ('user')
             if( Echo ) write(*,'(1x,A,i3)')'+ U component User'
             solver(VarU) = SolverUser
           case default
             write(*,'(1x,A,i3)')'+ error in command solver'
             write(*,'(1x,A,i3)')'+ solving U component'
             SolveU      = .true.
             messages = .true.
         end select
       case ('v')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving V component'
             SolveV      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ V component not solved'
             SolveV      = .false.
           case ('sparse')
             if( Echo ) write(*,'(1x,A,i3)')'+ V component SparseKit2'
             solver(VarV) = SparseKit2 
           case ('bicg')
             if( Echo ) write(*,'(1x,A,i3)')'+ V component BiCGstabL'
             solver(VarV) = SolverBiCGstabL 
           case ('direct')
             if( Echo ) write(*,'(1x,A,i3)')'+ V component Direct'
             solver(VarV) = SolverDirect 
           case ('user')
             if( Echo ) write(*,'(1x,A,i3)')'+ V component User'
             solver(VarV) = SolverUser 
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
             if( Echo ) write(*,'(1x,A,i3)')'+ solving V component'
             SolveV      = .true.
             messages = .true.
         end select
       case ('w')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving W component'
             SolveW      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ W component not solved'
             SolveW      = .false.
           case ('sparse')
             if( Echo ) write(*,'(1x,A,i3)')'+ W component SparseKit2'
             solver(VarW) = SparseKit2 
           case ('bicg')
             if( Echo ) write(*,'(1x,A,i3)')'+ W component BiCGstabL'
             solver(VarW) = SolverBiCGstabL 
           case ('direct')
             if( Echo ) write(*,'(1x,A,i3)')'+ W component Direct'
             solver(VarW) = SolverDirect 
           case ('user')
             if( Echo ) write(*,'(1x,A,i3)')'+ W component User'
             solver(VarW) = SolverUser
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
             if( Echo ) write(*,'(1x,A,i3)')'+ solving W component'
             SolveW      = .true.
             messages = .true.
         end select
       case ('p')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving pressure'
             SolveP      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ pressure not solved'
             SolveP      = .false.
           case ('sparse')
             if( Echo ) write(*,'(1x,A,i3)')'+ pressure SparseKit2'
             solver(VarP) = SparseKit2 
           case ('bicg')
             if( Echo ) write(*,'(1x,A,i3)')'+ pressure BiCGstabL'
             solver(VarP) = SolverBiCGstabL 
           case ('direct')
             if( Echo ) write(*,'(1x,A,i3)')'+ pressure Direct'
             solver(VarP) = SolverDirect 
           case ('user')
             if( Echo ) write(*,'(1x,A,i3)')'+ pressure User'
             solver(VarP) = SolverUser 
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
             if( Echo ) write(*,'(1x,A,i3)')'+ solving pressure'
             SolveP      = .true.
             messages = .true.
         end select
       case ('t')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving temperature'
             SolveEnthalpy = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature not solved'
             SolveEnthalpy = .false.
           case ('sparse')
             if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature SparseKit2'
             solver(VarT) = SparseKit2 
           case ('bicg')
             if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature BiCGstabL'
             solver(VarT) = SolverBiCGstabL 
           case ('direct')
             if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature Direct'
             solver(VarT) = SolverDirect 
           case ('user')
             if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature User'
             solver(VarT) = SolverUser 
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ error in command solver'
             if( Echo ) write(*,'(1x,A,i3)')'+ enthalpy/temperature not solved'
             SolveEnthalpy = .false.
             messages = .true.
         end select
       case ('k')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving turb. kin. energy k'
             SolveTurbEnergy = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ turb. kin. energy k not solved'
             SolveTurbEnergy = .false.
           case default
             write(*,'(1x,A,i3)')'+ error in command solver'
             write(*,'(1x,A,i3)')'+ turb. kin. energy k not solved'
             SolveTurbEnergy = .false.
             messages = .true.
         end select
       case ('eps')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ solving turb. dissipation epsilon'
             SolveTurbDiss = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ turb. dissipation not solved'
             SolveTurbDiss = .false.
           case default
             write(*,'(1x,A,i3)')'+ error in command solver'
             write(*,'(1x,A,i3)')'+ turb. dissipation not solved'
             SolveTurbDiss = .false.
             messages = .true.
         end select
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_scheme(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)
     dum1 = keyr(4)
     dum2 = keyr(5)
     call check_minmax(dum1,0.0,1.0)
     call check_minmax(dum2,0.0,0.5)

     select case (key2)
       case ('uvw')
         select case (key3(1:3))
           case ('ud ')
             Scheme(VarU) = DSud
             Scheme(VarV) = DSud
             Scheme(VarW) = DSud
             Gamma(VarU)  = 0.0
             Gamma(VarV)  = 0.0
             Gamma(VarW)  = 0.0
             if( Echo ) write(*,*)'+ velocity components: UD'
           case ('cd1')
             Scheme(VarU) = DScd1
             Scheme(VarV) = DScd1
             Scheme(VarW) = DScd1
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
             if( Echo ) write(*,*)'+ velocity components: uncorr. CD, blend ',dum1
           case ('cd2')
             Scheme(VarU) = DScd2
             Scheme(VarV) = DScd2
             Scheme(VarW) = DScd2
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
           case ('cd3')
             Scheme(VarU) = DScd3
             Scheme(VarV) = DScd3
             Scheme(VarW) = DScd3
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
             if( Echo ) write(*,*)'+ velocity components: CD3, blend ',dum1
             if( Echo ) write(*,*)'+ velocity components: CD3, blend ',dum1
           case ('cd ')
             Scheme(VarU) = DScds
             Scheme(VarV) = DScds
             Scheme(VarW) = DScds
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
             if( Echo ) write(*,*)'+ velocity components: CD, blend ',dum1
           case ('lud')
             Scheme(VarU) = DSlud
             Scheme(VarV) = DSlud
             Scheme(VarW) = DSlud
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
             if( Echo ) write(*,*)'+ velocity components: LUD, blend ',dum1 
           case ('gam')
             Scheme(VarU) = DSgam
             Scheme(VarV) = DSgam
             Scheme(VarW) = DSgam
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
            !if( dum2 /= 0.0 ) Beta_m = dum2
             if( Echo ) write(*,*)'+ velocity components: Gamma, blend ',dum1!,dum2
           case ('min')
             Scheme(VarU) = DSmod
             Scheme(VarV) = DSmod
             Scheme(VarW) = DSmod
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
             if( Echo ) write(*,*)'+ velocity components: MinMod, blend ',dum1 
           case ('lux')
             Scheme(VarU) = DSlux
             Scheme(VarV) = DSlux
             Scheme(VarW) = DSlux
             Gamma(VarU)  = 1.0
             Gamma(VarV)  = 1.0
             Gamma(VarW)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarU) = dum1
               Gamma(VarV) = dum1
               Gamma(VarW) = dum1
             endif
             if( Echo ) write(*,*)'+ velocity components: LUX, blend ',dum1 
           case default
             Scheme(VarU) = DSud
             Scheme(VarV) = DSud
             Scheme(VarW) = DSud
             Gamma(VarU)  = 0.0
             Gamma(VarV)  = 0.0
             Gamma(VarW)  = 0.0
             if( Echo ) write(*,*)'+ velocity components: default (UD)'
         end select
       case ('keps')
         select case (key3(1:3))
           case ('ud ')
             Scheme(VarTE) = DSud
             Scheme(VarED) = DSud
             Gamma(VarTE)  = 0.0
             Gamma(VarED)  = 0.0
             if( Echo ) write(*,*)'+ turb. model: UD'
           case ('cd1')
             Scheme(VarTE) = DScd1
             Scheme(VarED) = DScd1
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: CD, blend ',dum1
           case ('cd2')
             Scheme(VarTE) = DScd2
             Scheme(VarED) = DScd2
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: CD2, blend ',dum1
           case ('cd3')
             Scheme(VarTE) = DScd3
             Scheme(VarED) = DScd3
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: CD3, blend ',dum1
           case ('cd ')
             Scheme(VarTE) = DScds
             Scheme(VarED) = DScds
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: CD, blend ',dum1
           case ('lud')
             Scheme(VarTE) = DSlud
             Scheme(VarED) = DSlud
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: LUD, blend ',dum1 
           case ('gam')
             Scheme(VarTE) = DSgam
             Scheme(VarED) = DSgam
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
            !if( dum2 /= 0.0 ) Beta_m = dum2
             if( Echo ) write(*,*)'+ turb. model: Gamma, blend ',dum1!,dum2
           case ('min')
             Scheme(VarTE) = DSmod
             Scheme(VarED) = DSmod
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: MinMod, blend ',dum1 
           case ('lux')
             Scheme(VarTE) = DSlux
             Scheme(VarED) = DSlux
             Gamma(VarTE)  = 1.0
             Gamma(VarED)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarTE) = dum1
               Gamma(VarED) = dum1
             endif
             if( Echo ) write(*,*)'+ turb. model: LUD uncor., blend ',dum1 
           case default
             Scheme(VarTE) = DSud
             Scheme(VarED) = DSud
             Gamma(VarTE)  = 0.0
             Gamma(VarED)  = 0.0
             if( Echo ) write(*,*)'+ turb. model: default (UD)'
         end select
       case ('t')
         select case (key3(1:3))
           case ('ud ')
             Scheme(VarT)  = DSud
             Scheme(VarSC) = DSud
             Gamma(VarT)   = 0.0
             Gamma(VarSC)  = 0.0
             if( Echo ) write(*,*)'+ temperature/scalars: UD'
           case ('cd1')
             Scheme(VarT)  = DScd1
             Scheme(VarSC) = DScd1
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: CD, blend ',dum1
           case ('cd2')
             Scheme(VarT)  = DScd2
             Scheme(VarSC) = DScd2
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: CD2, blend ',dum1
           case ('cd3')
             Scheme(VarT)  = DScd3
             Scheme(VarSC) = DScd3
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: CD3, blend ',dum1
           case ('cd ')
             Scheme(VarT)  = DScds
             Scheme(VarSC) = DScds
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: CD, blend ',dum1
           case ('lud')
             Scheme(VarT)  = DSlud
             Scheme(VarSC) = DSlud
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: LUD, blend ',dum1 
           case ('gam')
             Scheme(VarT)  = DSgam
             Scheme(VarSC) = DSgam
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
            !if( dum2 /= 0.0 ) Beta_m = dum2
             if( Echo ) write(*,*)'+ temperature/scalars: Gamma, blend ',dum1!,dum2
           case ('min')
             Scheme(VarT)  = DSmod
             Scheme(VarSC) = DSmod
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: MinMod, blend ',dum1 
           case ('lux')
             Scheme(VarT)  = DSlux
             Scheme(VarSC) = DSlux
             Gamma(VarT)   = 1.0
             Gamma(VarSC)  = 1.0
             if( dum1 >= 0.0 )then
               Gamma(VarT)  = dum1
               Gamma(VarSC) = dum1
             endif
             if( Echo ) write(*,*)'+ temperature/scalars: LUD uncor., blend ',dum1 
           case default
             Scheme(VarT)  = DSud
             Scheme(VarSC) = DSud
             Gamma(VarT)   = 0.0
             Gamma(VarSC)  = 0.0
             if( Echo ) write(*,*)'+ temperature/scalars: default (UD)'
         end select
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_limit(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)
     dum1 = keyr(4)
     key4 = keys(4)
     call lowercase(key4)

     if( Echo ) write(*,*)'+ limit: ',key2,key3,dum1,'|',key4
     select case (key2(1:3))

       case ('u  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowU = .false.             
               write(IOdef,*)'Lower limit U switched off'
               write(IOdbg,*)'Lower limit U switched off'
               write(IOrun,*)'Lower limit U switched off'
             else
               LimitLowU   = .true.
               LowerLimitU = dum1
               write(IOdef,*)'Lower limit U: ',dum1
               write(IOdbg,*)'Lower limit U: ',dum1
               write(IOrun,*)'Lower limit U: ',dum1
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpU = .false.             
               write(IOdef,*)'Upper limit U switched off'
               write(IOdbg,*)'Upper limit U switched off'
               write(IOrun,*)'Upper limit U switched off'
             else
               LimitUpU = .true.
               UpperLimitU = dum1
               write(IOdef,*)'Upper limit U: ',dum1
               write(IOdbg,*)'Upper limit U: ',dum1
               write(IOrun,*)'Upper limit U: ',dum1
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select
       case ('v  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowV = .false.             
               write(IOdef,*)'Lower limit V switched off'
               write(IOdbg,*)'Lower limit V switched off'
               write(IOrun,*)'Lower limit V switched off'
             else
               LimitLowV   = .true.
               LowerLimitV = dum1
               write(IOdef,*)'Lower limit V: ',dum1
               write(IOdbg,*)'Lower limit V: ',dum1
               write(IOrun,*)'Lower limit V: ',dum1
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpV = .false.             
               write(IOdef,*)'Upper limit V switched off'
               write(IOdbg,*)'Upper limit V switched off'
               write(IOrun,*)'Upper limit V switched off'
             else
               LimitUpV = .true.
               UpperLimitV = dum1
               write(IOdef,*)'Upper limit V: ',dum1
               write(IOdbg,*)'Upper limit V: ',dum1
               write(IOrun,*)'Upper limit V: ',dum1
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select
       case ('w  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowW = .false.             
               write(IOdef,*)'Lower limit W switched off'
               write(IOdbg,*)'Lower limit W switched off'
               write(IOrun,*)'Lower limit W switched off'
             else
               LimitLowW   = .true.
               LowerLimitW = dum1
               write(IOdef,*)'Lower limit W: ',dum1
               write(IOdbg,*)'Lower limit W: ',dum1
               write(IOrun,*)'Lower limit W: ',dum1
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpW = .false.             
               write(IOdef,*)'Upper limit W switched off'
               write(IOdbg,*)'Upper limit W switched off'
               write(IOrun,*)'Upper limit W switched off'
             else
               LimitUpW = .true.
               UpperLimitW = dum1
               write(IOdef,*)'Upper limit W: ',dum1
               write(IOdbg,*)'Upper limit W: ',dum1
               write(IOrun,*)'Upper limit W: ',dum1
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select
       case ('uvw')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowU = .false.             
               LimitLowV = .false.             
               LimitLowW = .false.             
               write(IOdef,*)'Lower limit UVW switched off'
               write(IOdbg,*)'Lower limit UVW switched off'
               write(IOrun,*)'Lower limit UVW switched off'
             else
               LimitLowU   = .true.
               LowerLimitU = dum1
               LimitLowV   = .true.
               LowerLimitV = dum1
               LimitLowW   = .true.
               LowerLimitW = dum1
               write(IOdef,*)'Lower limit UVW: ',dum1
               write(IOdbg,*)'Lower limit UVW: ',dum1
               write(IOrun,*)'Lower limit UVW: ',dum1
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpU = .false.             
               LimitUpV = .false.             
               LimitUpW = .false.             
               write(IOdef,*)'Upper limit UVW switched off'
               write(IOdbg,*)'Upper limit UVW switched off'
               write(IOrun,*)'Upper limit UVW switched off'
             else
               LimitUpU = .true.
               UpperLimitU = dum1
               LimitUpV = .true.
               UpperLimitV = dum1
               LimitUpW = .true.
               UpperLimitW = dum1
               write(IOdef,*)'Upper limit UVW: ',dum1
               write(IOdbg,*)'Upper limit UVW: ',dum1
               write(IOrun,*)'Upper limit UVW: ',dum1
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select

       case ('vma')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowVMag = .false.             
               write(IOdef,*)'Lower limit Vel. Magnitude switched off'
               write(IOdbg,*)'Lower limit Vel. Magnitude switched off'
               write(IOrun,*)'Lower limit Vel. Magnitude switched off'
             else
               LimitLowVMag   = .true.
               LowerLimitVMag = dum1
               write(IOdef,*)'Lower limit Vel. Magnitude: ',dum1
               write(IOdbg,*)'Lower limit Vel. Magnitude: ',dum1
               write(IOrun,*)'Lower limit Vel. Magnitude: ',dum1
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpVMag = .false.             
               write(IOdef,*)'Upper limit Vel. Magnitude switched off'
               write(IOdbg,*)'Upper limit Vel. Magnitude switched off'
               write(IOrun,*)'Upper limit Vel. Magnitude switched off'
             else
               LimitUpVMag = .true.
               UpperLimitVMag = dum1
               write(IOdef,*)'Upper limit Vel. Magnitude: ',dum1
               write(IOdbg,*)'Upper limit Vel. Magnitude: ',dum1
               write(IOrun,*)'Upper limit Vel. Magnitude: ',dum1
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select

       case ('p  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowP = .false.             
               write(IOdef,*)'Lower limit P switched off'
               write(IOdbg,*)'Lower limit P switched off'
               write(IOrun,*)'Lower limit P switched off'
             else
               LimitLowP   = .true.
               LowerLimitP = dum1
               write(IOdef,*)'Lower limit P: ',dum1
               write(IOdbg,*)'Lower limit P: ',dum1
               write(IOrun,*)'Lower limit P: ',dum1
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpP = .false.             
               write(IOdef,*)'Upper limit P switched off'
               write(IOdbg,*)'Upper limit P switched off'
               write(IOrun,*)'Upper limit P switched off'
             else
               LimitUpP = .true.
               UpperLimitP = dum1
               write(IOdef,*)'Upper limit P: ',dum1
               write(IOdbg,*)'Upper limit P: ',dum1
               write(IOrun,*)'Upper limit P: ',dum1
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select

       case ('k  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowTurbEnergy = .false.             
               write(IOdef,*)'Lower limit turb. kin. energy k switched off'
               write(IOdbg,*)'Lower limit turb. kin. energy k switched off'
               write(IOrun,*)'Lower limit turb. kin. energy k switched off'
             else
               if( dum1 > 0.0 )then
                 LimitLowTurbEnergy = .true.
                 LowerLimitTE = dum1
                 write(IOdef,*)'Lower limit turb. kin. energy k: ',dum1
                 write(IOdbg,*)'Lower limit turb. kin. energy k: ',dum1
                 write(IOrun,*)'Lower limit turb. kin. energy k: ',dum1
               else
                 write(*,*) '+ invalid limit: ',dum1
                 messages = .true.
               endif
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpTurbEnergy = .false.             
               write(IOdef,*)'Upper limit turb. kin. energy k switched off'
               write(IOdbg,*)'Upper limit turb. kin. energy k switched off'
               write(IOrun,*)'Upper limit turb. kin. energy k switched off'
             else
               if( dum1 > 0.0 )then
                 LimitUpTurbEnergy = .true.
                 UpperLimitTE = dum1
                 write(IOdef,*)'Upper limit turb. kin. energy k: ',dum1
                 write(IOdbg,*)'Upper limit turb. kin. energy k: ',dum1
                 write(IOrun,*)'Upper limit turb. kin. energy k: ',dum1
               else
                 write(*,*) '+ invalid limit: ',dum1
                 messages = .true.
               endif 
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select
       case ('e  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowTurbDiss = .false.             
               write(IOdef,*)'Lower limit turb. dissipation epsilon switched off'
               write(IOdbg,*)'Lower limit turb. dissipation epsilon switched off'
               write(IOrun,*)'Lower limit turb. dissipation epsilon switched off'
             else
               if( dum1 > 0.0 )then
                 LimitLowTurbDiss = .true.
                 LowerLimitED = dum1
                 write(IOdef,*)'Lower limit turb. dissipation epsilon: ',dum1
                 write(IOdbg,*)'Lower limit turb. dissipation epsilon: ',dum1
                 write(IOrun,*)'Lower limit turb. dissipation epsilon: ',dum1
               else
                 write(*,*) '+ invalid limit: ',dum1
                 messages = .true.
               endif
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpTurbDiss = .false.             
               write(IOdef,*)'Upper limit turb. dissipation epsilon switched off'
               write(IOdbg,*)'Upper limit turb. dissipation epsilon switched off'
               write(IOrun,*)'Upper limit turb. dissipation epsilon switched off'
             else
               if( dum1 > 0.0 )then
                 LimitUpTurbDiss = .true.
                 UpperLimitED= dum1
                 write(IOdef,*)'Upper limit turb. dissipation epsilon: ',dum1
                 write(IOdbg,*)'Upper limit turb. dissipation epsilon: ',dum1
                 write(IOrun,*)'Upper limit turb. dissipation epsilon: ',dum1
               else
                 write(*,*) '+ invalid limit: ',dum1
                 messages = .true.
               endif 
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select
       case ('t  ')
         select case (key3(1:3))
           case ('low')
             if( key4 == 'off' )then
               LimitLowEnthalpy = .false.             
               write(IOdef,*)'Lower limit enthalpy switched off'
               write(IOdbg,*)'Lower limit enthalpy switched off'
               write(IOrun,*)'Lower limit enthalpy switched off'
             else
               if( dum1 > 0.0 )then
                 LimitLowEnthalpy = .true.
                 LowerLimitT = dum1
                 write(IOdef,*)'Lower limit T: ',dum1
                 write(IOdbg,*)'Lower limit T: ',dum1
                 write(IOrun,*)'Lower limit T: ',dum1
               else
                 write(*,*) '+ invalid limit: ',dum1
                 messages = .true.
               endif
             endif
           case ('upp')
             if( key4 == 'off' )then
               LimitUpEnthalpy = .false.             
               write(IOdef,*)'Upper limit enthalpy switched off'
               write(IOdbg,*)'Upper limit enthalpy switched off'
               write(IOrun,*)'Upper limit enthalpy switched off'
             else
               if( dum1 > 0.0 )then
                 LimitUpEnthalpy = .true.
                 UpperLimitT= dum1
                 write(IOdef,*)'Upper limit T: ',dum1
                 write(IOdbg,*)'Upper limit T: ',dum1
                 write(IOrun,*)'Upper limit T: ',dum1
               else
                 write(*,*) '+ invalid limit: ',dum1
                 messages = .true.
               endif 
             endif
           case default
             write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
             messages = .true.
         end select
       case ('s  ')
         i1 = keyi(3)
         if( i1 <= 0  )then
           write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
           write(*,*) '+ or scalar outside of range',i1
           messages = .true.
         else
           key4 = keys(4)
           call lowercase(key4)
           dum1 = keyr(5)
           select case (key4(1:3))
             case ('low')
               LimitLowScalar(i1) = .true.
               LowerLimitSC(i1)   = dum1
               write(IOdef,*)'Lower limit scalar: ',i1,dum1
               write(IOdbg,*)'Lower limit scalar: ',i1,dum1
               write(IOrun,*)'Lower limit scalar: ',i1,dum1
             case ('upp')
               LimitUpScalar(i1) = .true.
               UpperLimitSC(i1)  = dum1
               write(IOdef,*)'Upper limit scalar: ',i1,dum1
               write(IOdbg,*)'Upper limit scalar: ',i1,dum1
               write(IOrun,*)'Upper limit scalar: ',i1,dum1
             case default
               write(*,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
               messages = .true.
           end select
         endif           
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
         messages = .true.
     end select
     
   return
   !
   !========================================================
   !
   entry select_case_slopelimiter(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)
     dum1 = keyr(4)

     if( Echo ) write(*,*)'+ slope limiter: ', &
                           key2(1:3),key3(1:3),dum1
     select case (key2(1:3))
       case ('u  ')
         i = VarU
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         
       case ('v  ')
         i = VarV 
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         
       case ('w  ')
         i = VarW
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         
       case ('uvw')
         i = VarU
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         i = VarV 
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         i = VarW
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)

       case ('p  ')
         i = VarP
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         i = VarPP
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)

       case ('k  ')
         i = VarTE
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
       case ('eps')
         i = VarED
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
       case ('kep')
         i = VarTE
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
         i = VarED
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
       case ('t  ')
         i = VarT
         call set_slopelimiter(i, UseSlopeLimiter(i), &
                                  URFSlopeLimiter(i), &
                                  LimiterSlope(i),    &
                                  Messages, key3, dum1)
       case ('sca')
         key3 = keys(3)
         call lowercase(key3)
         if( key3 == 'all' )then
           key4 = keys(4)
           call lowercase(key4)
           dum1 = keyr(5)           

           do j=MinSC,MaxVar
             call set_slopelimiter(j, UseSlopeLimiter(j), &
                                      URFSlopeLimiter(j), &
                                      LimiterSlope(j),    &
                                      Messages, key4, dum1)
           end do
         else       
           i = NVar + keyi(3)
           key4 = keys(4)
           call lowercase(key4)
           dum1 = keyr(5)
           
           call set_slopelimiter(i, UseSlopeLimiter(i), &
                                    URFSlopeLimiter(i), &
                                    LimiterSlope(i),    &
                                    Messages, key4, dum1)                                 
         endif
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_gradient(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)
     idum = keyi(4)

     if( Echo ) write(*,*)'+ set gradient: ', &
                           key2(1:3),key3(1:3),idum
     select case (key2(1:3))
       case ('u  ')
         i = VarU
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
       case ('v  ')
         i = VarV
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
       case ('w  ')
         i = VarW
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
       case ('uvw')
         i = VarU
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
         i = VarV
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
         i = VarW
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)

       case ('p  ')
         i = VarP
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
         i = VarPP
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)

       case ('k  ')
         i = VarTE
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
       case ('eps')
         i = VarED
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
       case ('kep')
         i = VarTE
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
         i = VarED
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)

       case ('t  ')
         i = VarT
         call set_gradient(i, GradVar(i), nGradVar(i), &
                              Messages, key3, idum)
                              
       case ('sca')
         key3 = keys(3)
         call lowercase(key3)
         if( key3 == 'all' )then
           key4  = keys(4)
           call lowercase(key4)
           idum1 = keyi(5)           

           do j=MinSC,MaxVar
             call set_gradient(j, GradVar(j), nGradVar(j), &
                                  Messages, key4, idum1)
           end do
         else       
           i = NVar + keyi(3)
           key4  = keys(4)
           call lowercase(key4)
           idum1 = keyi(5)
           
           call set_gradient(i, GradVar(i), nGradVar(i), &
                                Messages, key4, idum1)
         endif


       case default

         write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
         messages = .true.

     end select
                              
         
   return
   !
   !========================================================
   !
   entry select_case_opendx(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)

     if( Echo ) write(*,'(1x,A,i3)')'+ opendx...'
     if( Echo ) write(*,*)'+ opendx keyword with strings: >',&
                  key2(1:lens(key2)),'< and >',&
                  key3(1:lens(key3)),'<'

     select case (key2)
       case ('normals')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ writing DXnormals'
             DXnormals      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXnormals'
             DXnormals      = .false.
           case default
             write(*,'(1x,A,i3)')'+ error in command opendx'
             messages = .true.
         end select
       case ('centers')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ writing DXcenters'
             DXcenters      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXcenters'
             DXcenters      = .false.
           case default
             write(*,'(1x,A,i3)')'+ error in command opendx'
             messages = .true.
         end select
       case ('massflux')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ writing DXmassflux'
             DXmassflux      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXmassflux'
             DXmassflux      = .false.
           case default
             if( Echo ) write(*,'(1x,A,i3)')'+ error in command opendx'
             messages = .true.
         end select
       case ('debug')
         select case (key3(1:lens(key3)))
           case ('on')
             if( Echo ) write(*,'(1x,A,i3)')'+ writing DXdebugdata'
             DXdebugdata      = .true.
           case ('off')
             if( Echo ) write(*,'(1x,A,i3)')'+ not writing DXdebugdata'
             DXdebugdata      = .false.
           case default
             write(*,'(1x,A,i3)')'+ error in command opendx'
             messages = .true.
         end select
       case ('dump')
         DXdump = keyi(3)
         if( Echo ) write(*,'(1x,A,i3,A)')'+ call opendx every ',&
                    DXdump,' steps'
       case ('no')
         if( Echo ) write(*,'(1x,A,i3,A)')'+ no opendx output'
         UseOpenDX = .false.
       case ('off')
         if( Echo ) write(*,'(1x,A,i3,A)')'+ no opendx output'
         UseOpenDX = .false.
         
       case default
         write(*,*) '+ unknown/unexpected keyword: ',key2(1:lens(key2))
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_post(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)        ! variable
     call lowercase(key2)
     key3 = keys(3)        ! cell or vertex
     call lowercase(key3)
     key4 = keys(4)        ! yes or no
     call lowercase(key3)

     if( Echo ) &
       write(*,'(1x,A,A,A)')'+ post keyword with strings = ',key2,key3,key4
     select case (key2(1:3))
       case ('u  ')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarU) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarU) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarU) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarU) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('v  ')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarV) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarV) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarV) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarV) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('w  ')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarW) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarW) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarW) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarW) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('p  ')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarP) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarP) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarP) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarP) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('k  ')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarTE) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarTE) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarTE) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarTE) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('eps')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarED) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarED) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarED) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarED) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('t  ')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarT) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarT) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarT) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarT) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('sca')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarSC) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarSC) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarSC) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarSC) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('den')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarDen) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarDen) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarDen) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarDen) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('vis')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarVis) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarVis) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarVis) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarVis) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
       case ('lvi')
         if( key3(1:1) == 'c' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostC(VarLVis) = .true.
           else if( key4(1:1) == 'n' )then
             PostC(VarLVis) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif
         else if( key3(1:1) == 'v' )then
           if( key4(1:1) == 'y' .or. key4(1:1) == ' ' )then
             PostV(VarLVis) = .true.
           else if( key4(1:1) == 'n' )then
             PostV(VarLVis) = .false.
           else
             write(*,'(1x,A,i3)')'+ error in post command:',key2           
           endif               
         else
           write(*,'(1x,A,i3)')'+ error in post command:',key3
           messages = .true.
         endif
           
       case default
         write(*,'(1x,A,i3)')'+ error in post command'
         messages = .true.
     end select
 
   return
   !
   !========================================================
   !
   entry select_case_check(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     CheckOut = .true.

     key2 = keys(2)        ! variable
     call lowercase(key2)
     key3 = keys(3)        ! "range", lower, upper or 
     call lowercase(key3)  ! "average", average, band
     r1 = keyr(4)          !  
     r2 = keyr(5)          ! value
     key4 = keys(6)        ! "range", lower, upper or 
     call lowercase(key4)  ! report or do not report (silent)

     if( Echo ) &
       write(*,'(1x,A,A,A)')'+ check keyword with strings = ',key2,key3
     select case (key2(1:3))
       case ('u  ')
         CheckU = .true.
         if( key3(1:1) == 'r' )then
           CheckLowU = r1
           CheckUpU  = r2
           if( key4(1:1) == 'r' ) CheckReprtU = .true.
	 else if( key3(1:1) == 'a' )then
           CheckAverageU      =  r1
           CheckAverageBandU  =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtU = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('v  ')
         CheckV = .true.
         if( key3(1:1) == 'r' )then
           CheckLowV = r1
           CheckUpV  = r2
           if( key4(1:1) == 'r' ) CheckReprtV = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageV      =  r1
           CheckAverageBandV  =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtV = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('w  ')
         CheckW = .true.
         if( key3(1:1) == 'r' )then
           CheckLowW = r1
           CheckUpW  = r2
           if( key4(1:1) == 'r' ) CheckReprtW = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageW      =  r1
           CheckAverageBandW  =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtW = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('p  ')
         CheckP = .true.
         if( key3(1:1) == 'r' )then
           CheckLowP = r1
           CheckUpP  = r2
           if( key4(1:1) == 'r' ) CheckReprtP = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageP      =  r1
           CheckAverageBandP  =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtP = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('k  ')
         CheckTE = .true.
         if( key3(1:1) == 'r' )then
           CheckLowTE = r1
           CheckUpTE  = r2
           if( key4(1:1) == 'r' ) CheckReprtTE = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageTE      =  r1
           CheckAverageBandTE  =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtTE = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('eps')
         CheckED = .true.
         if( key3(1:1) == 'r' )then
           CheckLowED = r1
           CheckUpED  = r2
           if( key4(1:1) == 'r' ) CheckReprtED = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageED     =  r1
           CheckAverageBandED =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtED = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('t  ')
         CheckT = .true.
         if( key3(1:1) == 'r' )then
           CheckLowT = r1
           CheckUpT  = r2
           if( key4(1:1) == 'r' ) CheckReprtT = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageT     =  r1
           CheckAverageBandT =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtT = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('den')
         CheckDEN = .true.
         if( key3(1:1) == 'r' )then
           CheckLowDEN = r1
           CheckUpDEN  = r2
           if( key4(1:1) == 'r' ) CheckReprtDEN = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageDEN     =  r1
           CheckAverageBandDEN =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtDEN = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif
       case ('vis')
         CheckVIS = .true.
         if( key3(1:1) == 'r' )then
           CheckLowVIS = r1
           CheckUpVIS  = r2
           if( key4(1:1) == 'r' ) CheckReprtVIS = .true.
         else if( key3(1:1) == 'a' )then
           CheckAverageVIS     =  r1
           CheckAverageBandVIS =  abs(r2)
           if( key4(1:1) == 'r' ) CheckReprtVIS = .true.
         else
           write(*,'(1x,A,i3)')'+ error in check command:',key3
           messages = .true.
         endif

       case default
         write(*,'(1x,A,i3)')'+ error in check command'
         messages = .true.
     end select

   return
   !
   !========================================================
   !
   entry select_case_print(Echo,Messages,KeyS,KeyI,KeyR,KeyU)

     key2 = keys(2)
     call lowercase(key2)
     key3 = keys(3)
     call lowercase(key3)
     
     select case (key2(1:4))
       case ('cell')
         PrintCellVar     = .true.
         if( key3(1:4) == 'user' )then

           PrintCellVarUser = .true.             

         elseif( key3(1:3) == 'all' )then

           IPrintCellVarStart = 1
           IPrintCellVarEnd   = NCel
           IPrintCellVarInc   = 1

         else
         
           istart = keyi(3)
           iend   = keyi(4)
           iinc   = keyi(5)

           if( iend < istart )then
             write(*,'(1x,A)')'+ error in print cell command: end < start'
             messages = .true.
           elseif( iinc <= 0 )then
             write(*,'(1x,A)')'+ error in print cell command: invalid increment'
             messages = .true.
           elseif( istart <= 0 )then
             write(*,'(1x,A)')'+ error in print cell command: invalid start'
             messages = .true.
           elseif( iend <= 0 )then
             write(*,'(1x,A)')'+ error in print cell command: invalid end'
             messages = .true.
           endif

           IPrintCellVarStart = istart
           IPrintCellVarEnd   = iend
           IPrintCellVarInc   = iinc

           write(*,'(1x,A,3(i6))')'Printing cells ',istart,iend,iinc

         endif
       case ('wall')
         PrintWallVar     = .true.
         if( key3(1:4) == 'user' )then

           PrintWallVarUser = .true.             

         elseif( key3(1:3) == 'all' )then

           IPrintWallVarStart = 1
           IPrintWallVarEnd   = NBnd
           IPrintWallVarInc   = 1

         else
         
           istart = keyi(3)
           iend   = keyi(4)
           iinc   = keyi(5)

           if( iend < istart )then
             write(*,'(1x,A)')'+ error in print wall command: end < start'
             messages = .true.
           elseif( iinc <= 0 )then
             write(*,'(1x,A)')'+ error in print wall command: invalid increment'
             messages = .true.
           elseif( istart <= 0 )then
             write(*,'(1x,A)')'+ error in print wall command: invalid start'
             messages = .true.
           elseif( iend <= 0 )then
             write(*,'(1x,A)')'+ error in print wall command: invalid end'
             messages = .true.
           endif

           IPrintWallVarStart = istart
           IPrintWallVarEnd   = iend
           IPrintWallVarInc   = iinc

           write(*,'(1x,A,3(i6))')'Printing Walls ',istart,iend,iinc

         endif
       case ('file')
         PrintFile = keys(3)
         write(*,'(1x,A,A)')'Printing to file ',PrintFile
       case default
         write(*,'(1x,A)')'+ error in print command'
         messages = .true.
     end select
    
   return
    
end subroutine ReadControlFile2
!
!========================================================
!
subroutine set_slopelimiter(iVar, UseIt, UrfIt, LimitIt, Messages, &
                            key3, dummy)

   use constants
   use variables, only : Variable

   logical, intent(OUT)   :: UseIt
   real, intent(OUT)      :: UrfIt
   integer, intent(OUT)   :: LimitIt
   logical, intent(INOUT) :: Messages

   character(len=32), intent(IN) :: key3
   real, intent(INOUT)           :: dummy

   if( dummy >= 0.0 .and. dummy <= 1.0 )then     
     if( dummy <= Small ) dummy = 1.0
   else
     write(IOdef,*) '+ invalid value for slope limiter: ',dummy
     Messages = .true.
     dummy    = 1.0
   endif

   select case (key3(1:2))
   
     case ('of') 
       
       UseIt   = .false.
       UrfIt   = 0.0
       LimitIt = SlopeLimiterOFF

       write(IOdef,*)'Slope limiter switched off for ',Variable(iVar)
       write(IOdbg,*)'Slope limiter switched off for ',Variable(iVar)
       write(IOrun,*)'Slope limiter switched off for ',Variable(iVar)
 
     case ('bj')
       
       select case (key3(3:3))
       
         case ('c')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterBJc           
           write(IOdef,*)'Slope limiter BJ c ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter BJ c ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter BJ c ',Variable(iVar),' ',dummy

         case ('f')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterBJf           
           write(IOdef,*)'Slope limiter BJ f ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter BJ f ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter BJ f ',Variable(iVar),' ',dummy

         case ('n')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterBJn           
           write(IOdef,*)'Slope limiter BJ n ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter BJ n ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter BJ n ',Variable(iVar),' ',dummy

         case (' ')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterBJf           
           write(IOdef,*)'Slope limiter BJ ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter BJ ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter BJ ',Variable(iVar),' ',dummy

         case default
           
           write(IOdef,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
           Messages = .true.
  
       end select
       
     case ('vn')
       
       select case (key3(3:3))
       
         case ('c')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVNc           
           write(IOdef,*)'Slope limiter VN c ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VN c ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VN c ',Variable(iVar),' ',dummy

         case ('f')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVNf           
           write(IOdef,*)'Slope limiter VN f ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VN f ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VN f ',Variable(iVar),' ',dummy

         case ('n')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVNn           
           write(IOdef,*)'Slope limiter VN n ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VN n ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VN n ',Variable(iVar),' ',dummy

         case (' ')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVNf           
           write(IOdef,*)'Slope limiter VN ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VN ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VN ',Variable(iVar),' ',dummy

         case default
           
           write(IOdef,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
           Messages = .true.
  
       end select

     case ('va')
       
       select case (key3(3:3))
       
         case ('c')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVAc           
           write(IOdef,*)'Slope limiter VA c ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VA c ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VA c ',Variable(iVar),' ',dummy

         case ('f')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVAf           
           write(IOdef,*)'Slope limiter VA f ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VA f ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VA f ',Variable(iVar),' ',dummy

         case ('n')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVAn           
           write(IOdef,*)'Slope limiter VA n ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VA n ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VA n ',Variable(iVar),' ',dummy

         case (' ')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterVNf           
           write(IOdef,*)'Slope limiter VA ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter VA ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter VA ',Variable(iVar),' ',dummy

         case default
           
           write(IOdef,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
           Messages = .true.
  
       end select

     case ('p1')
       
       select case (key3(3:3))
       
         case ('c')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterP1c           
           write(IOdef,*)'Slope limiter P1 c ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter P1 c ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter P1 c ',Variable(iVar),' ',dummy

         case ('f')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterP1f           
           write(IOdef,*)'Slope limiter P1 f ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter P1 f ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter P1 f ',Variable(iVar),' ',dummy

         case ('n')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterP1n           
           write(IOdef,*)'Slope limiter P1 n ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter P1 n ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter P1 n ',Variable(iVar),' ',dummy

         case (' ')
         
           UseIt   = .true.             
           UrfIt   = dummy      
           LimitIt = SlopeLimiterP1f           
           write(IOdef,*)'Slope limiter P1 ',Variable(iVar),' ',dummy
           write(IOdbg,*)'Slope limiter P1 ',Variable(iVar),' ',dummy
           write(IOrun,*)'Slope limiter P1 ',Variable(iVar),' ',dummy

         case default
           
           write(IOdef,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
           Messages = .true.
  
       end select

     case default
           
       write(IOdef,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
       Messages = .true.
  
   end select
     
end subroutine set_slopelimiter
subroutine set_gradient(iVar, iChoice, nPasses, &
                        Messages, key3, idummy)

   use constants
   use variables, only : Variable

   integer, intent(OUT)   :: iChoice
   integer, intent(OUT)   :: nPasses
   logical, intent(INOUT) :: Messages

   character(len=32), intent(IN) :: key3
   integer, intent(INOUT)        :: idummy

   if( idummy == 0 )then     
     idummy = Ngradient
   else if( idummy < 0 )then  
     write(IOdef,*) '+ invalid number of passes: ',idummy
     Messages = .true.
     idummy    = Ngradient
   else if( idummy > 10 )then  
     write(IOdef,*) '+ are you sure of these gradient passes: ',idummy,'?'
   endif

   select case (key3(1:1))
   
     case ('g') 
       
       iChoice = GradGauss
       nPasses = idummy

       write(IOdef,*)'Using Gauss for ',Variable(iVar),' passes:',idummy
       write(IOdbg,*)'Using Gauss for ',Variable(iVar),' passes:',idummy
       write(IOrun,*)'Using Gauss for ',Variable(iVar),' passes:',idummy
 
     case ('l')

       iChoice = GradLS
       nPasses = idummy

       write(IOdef,*)'Using Least Squares for ',Variable(iVar),idummy
       write(IOdbg,*)'Using Least Squares for ',Variable(iVar),idummy
       write(IOrun,*)'Using Least Squares for ',Variable(iVar),idummy
   
     case default
           
       write(IOdef,*) '+ unknown/unexpected keyword: ',key3(1:lens(key3))
       Messages = .true.
  
   end select

end subroutine set_gradient
