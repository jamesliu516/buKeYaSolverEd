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
subroutine cross_product(a,b,c)
!
!  classic cross product: a x b = c
!
   real, dimension(3), intent(IN)  :: a, b
   real, dimension(3), intent(OUT) :: c
   
   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine cross_product
!========================================================================
subroutine normalise(a)
!
!  normalise a vector to length 1
!
   real, dimension(3), intent(INOUT) :: a
   real :: rlength
   
   rlength = 1.0/sqrt( dot_product(a,a) )
   
   a(1) = a(1)*rlength
   a(2) = a(2)*rlength
   a(3) = a(3)*rlength

end subroutine normalise
!========================================================================
real function vector_length(a)

   real, dimension(3), intent(IN) :: a

   vector_length = sqrt( dot_product(a,a) )
   
end function vector_length
!========================================================================
real function vector_cosangle(a,b)
!
!  the cosine of the angle between two vectors in the form of cos(angle)
!
   real, dimension(3), intent(IN) :: a, b

   vector_cosangle = dot_product(a,b)
   vector_cosangle = vector_cosangle/(vector_length(a)*vector_length(b))
   
end function vector_cosangle
!========================================================================
real function det33v(a,b,c) 
!
!  determinant of a 3 by 3 matrix
!
   real, dimension(3), intent(IN) :: a, b, c

   det33v = ( a(2)*b(3) - b(2)*a(3) )*c(1) &
          + ( b(1)*a(3) - a(1)*b(3) )*c(2) &
          + ( a(1)*b(2) - a(2)*b(1) )*c(3)

end function det33v
!========================================================================
real function det33m(a) 
!
!  same as the previous function, different form: classic determinant
!  of a 3 by 3 matrix
!
   real, dimension(3,3), intent(IN) :: a

   det33m = a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
          - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) &
          + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)
           
end function det33m
!========================================================================
real function det44m(a) 
!
!  determinant of a 4 by 4 matrix
!
   real, dimension(4,4), intent(IN) :: a
   real, dimension(3,3) :: b

   det44m = 0.0
   
   b(1,1) = a(2,2)
   b(1,2) = a(2,3)
   b(1,3) = a(2,4)
   b(2,1) = a(3,2)
   b(2,2) = a(3,3)
   b(2,3) = a(3,4)
   b(3,1) = a(4,2)
   b(3,2) = a(4,3)
   b(3,3) = a(4,4)
   
   det44m = det44m + a(1,1)*det33m(b)
   
   b(1,1) = a(2,1)
   b(1,2) = a(2,3)
   b(1,3) = a(2,4)
   b(2,1) = a(3,1)
   b(2,2) = a(3,3)
   b(2,3) = a(3,4)
   b(3,1) = a(4,1)
   b(3,2) = a(4,3)
   b(3,3) = a(4,4)
   
   det44m = det44m - a(1,2)*det33m(b)
   
   b(1,1) = a(2,1)
   b(1,2) = a(2,2)
   b(1,3) = a(2,4)
   b(2,1) = a(3,1)
   b(2,2) = a(3,2)
   b(2,3) = a(3,4)
   b(3,1) = a(4,1)
   b(3,2) = a(4,2)
   b(3,3) = a(4,4)
   
   det44m = det44m + a(1,3)*det33m(b)
   
   b(1,1) = a(2,1)
   b(1,2) = a(2,2)
   b(1,3) = a(2,3)
   b(2,1) = a(3,1)
   b(2,2) = a(3,2)
   b(2,3) = a(3,4)
   b(3,1) = a(4,1)
   b(3,2) = a(4,2)
   b(3,3) = a(4,3)
   
   det44m = det44m - a(1,4)*det33m(b)
      
end function det44m
!========================================================================
integer function lens(string)

   character(len=*) string
   
   do i=len(string),1,-1
     if( string(i:i) .ne. ' ')goto 10
   end do
   i = 0
10 continue

   lens = i
    
end function lens
!========================================================================
subroutine openfile(iunit,casename,extension,reqform,reqacc,status,idebug)

  character(len=*)   casename
  character(len=*)   extension
  character(len=*)   reqform
  character(len=*)   reqacc
  character(len=*)   status
  character(len=128) filename
  character(len=11)  form 

  logical exists

  filename = casename(1:lens(casename))//extension(1:lens(extension))
  length   = lens(filename)

  if( idebug > 2 )write(*,*) 'Opening ',filename(1:length),' l=',length
  
  if( status(1:3) == 'OLD' )then
    inquire(file=filename(1:length),exist=exists,form=form)
    if( .not. exists )then
      write(*,*) '*** Error: File ',filename(1:length),' does not exist'
      stop
    endif
  endif

  !write(*,*)'open:',iunit,filename(1:length),' ',reqacc

  if( reqacc == 'APPEND' )then
    open(iunit,file=filename(1:length),form=reqform, &
               access='SEQUENTIAL',position='APPEND',status=status)
  elseif( reqacc  == 'STREAM' .and.                  &
          reqform == 'UNFORMATTED_BIG_ENDIAN' )then
      open(iunit,file=filename(1:length),form='UNFORMATTED', &
                 convert='BIG_ENDIAN',                       &
		 access=reqacc,status='REPLACE')  
  else
    if( reqform == 'UNFORMATTED_BIG_ENDIAN' )then
      open(iunit,file=filename(1:length),form='UNFORMATTED', &
                 convert='BIG_ENDIAN',                       &
		 access=reqacc,status=status)
    else
      open(iunit,file=filename(1:length),form=reqform, &
                 access=reqacc,status=status)
    endif
  endif

  if( idebug >= 2 ) write(*,*) 'File ',filename(1:length),' opened'
  
end subroutine openfile
!========================================================================
subroutine feelfile(iunit,iform,casename,extension,reqform,reqacc,status,idebug)

   character(len=*)  casename
   character(len=*)  extension
   character(len=*)  reqform
   character(len=*)  reqacc
   character(len=*)  status
   character(len=48) filename

   character(len=12) string

   filename = casename(1:lens(casename))//extension(1:lens(extension))
   length   = lens(filename)

   if( idebug > 2 )write(*,*) 'Probing ',filename(1:length)

   call openfile(iunit,casename,extension(1:lens(extension)), &
                 'FORMATTED','SEQUENTIAL','OLD',idebug)

   read(iunit,'(a12)',err=2) string(1:12)
   
   if( string(1:12) == 'dolfyn asc g' )then
     write(*,*) filename(1:length),' ascii'
     iform = 1
   else 
     write(*,*) filename(1:length),' binary'
     iform = 0
   endif
   close(iunit)

 2 continue
   if( idebug >= 2 ) write(*,*) 'File ',filename(1:length),' probed'
  
end subroutine feelfile
!========================================================================
subroutine uppercase(string)

  character(len=*) :: string
  
  ial = ichar('a')
  izl = ichar('z') 
  iau = ichar('A')
  
  ioffset = iau - ial

  do i=1,len(string)
    if( ichar(string(i:i)) >= ial .and. ichar(string(i:i)) <= izl )then
      string(i:i) = char(ichar(string(i:i))+ioffset)
    endif      
  end do

end subroutine uppercase
!========================================================================
subroutine lowercase(string)

  character(len=*) :: string
  
  ial = ichar('a')
  izl = ichar('z')
  iau = ichar('A')
  izu = ichar('Z')
  
  ioffset = ial - iau

  do i=1,len(string)
    if( ichar(string(i:i)) >= iau .and. ichar(string(i:i)) <= izu )then
      string(i:i) = char(ichar(string(i:i))+ioffset)
    endif      
  end do

end subroutine lowercase
!========================================================================
real function random_valentclare()
!========================================================================
! random number generator 
! org. f77 code by: Dick Valent and Fred Clare addapted from
! "Random Number Generators: Good Ones are Hard to Find"
! by Steven Park, Keith Miller, Communications of the ACM, October, 1988.
!========================================================================
  
   integer, parameter :: MPLIER =  16807, MODLUS = 2147483647, &
                         MOBYMP = 127773, MOMDMP = 2836 

   integer :: hvlue, lvlue, testv, nextn = 0

   integer :: jseed = 123456789, ifirst = 0
   save       jseed, ifirst, nextnumber, nextn

   if( ifirst == 0 )then
     nextnumber = jseed
     ifirst = 1
   endif

   hvlue = nextnumber / mobymp
   lvlue = mod(nextn, mobymp)
   testv = mplier*lvlue - momdmp*hvlue
   if( testv > 0 )then
     nextn = testv
   else
     nextn = testv + modlus
   endif
   random_valentclare = real(nextn)/real(modlus)

   return

   entry RandomSeed(iseed)

   jseed  = iseed
   ifirst = 0  

   random_valentclare = 0.0

   return
   
end function random_valentclare   
subroutine CheckMassFlux

   use constants
   use geometry
   use variables

   write(*,*)'Checking mass fluxes'

   ierr = 0
   do ip=1,Ncel
     cmf = 0.0
     cmfmax = 1.e-12
     do j=1,NFaces(ip)
       k   = CFace(ip,j)
       tmp = MassFlux( k )
       if( tmp /= 0.0 )then
         cmfmax = max( cmfmax , abs(tmp) ) 
       endif
       if( ip == Face(k)%cell1 )then
         cmf = cmf + tmp
       else
         cmf = cmf - tmp
       endif
     end do

     if( abs(cmf) >= cmfmax )then
       ierr = ierr + 1
       if( ierr <= 10 )then
         write(*,*) '*** Check mass flux in cell =',ip,cmf,cmfmax
         write(*,*) 'volume:',Cell(ip)%vol
         cmf  = 0.0
         do j=1,NFaces(ip)
           k = CFace(ip,j)
           tmp = MassFlux( k )
           if( ip == Face(k)%cell1 )then
             cmf = cmf + tmp
             write(*,*) 'face ',k,' + ',tmp,' = ',cmf
           else
             cmf = cmf - tmp
             write(*,*) 'face ',k,' - ',tmp,' = ',cmf
           endif
         end do
       else if( ierr == 2 )then
         write(*,*)' (no more shown)'
       endif
     endif

   end do
   if( ierr > 0 )write(*,*)'Check :',ierr
    
end subroutine CheckMassFlux
subroutine TrackMemory(istat,i,string)
!========================================================================

   use constants, only: memory, debug, IOdbg

   integer, intent(in) :: i, istat
   character(len=*), intent(in) :: string 
   
   memory = memory + i
   
   change = 0.001*float(i)
   words  = 0.000001*float(memory)
   if( istat == 0 .and. debug > -1 )then

     write(IOdbg,'(1x,A,f8.3,A,f10.3,A,A)') &
       'MEMORY ',words,' MWords, change ',change,' kWords => ',string

   else
     if( istat /= 0 ) write(*,'(1x,A,A,A,i12,i12)') &
       '*** WARNING: ',string,', memory: ',words,' MW, status ',istat
   endif
   
end subroutine TrackMemory
subroutine CleanUp

   use constants
   use geometry
   use variables
   use scalars
   use particles
!   use searched

   if( allocated( U        ) ) deallocate ( U        )
   if( allocated( V        ) ) deallocate ( V        )
   if( allocated( W        ) ) deallocate ( W        )
   if( allocated( P        ) ) deallocate ( P        )
   if( allocated( PP       ) ) deallocate ( PP       )
   if( allocated( TE       ) ) deallocate ( TE       )
   if( allocated( ED       ) ) deallocate ( ED       )
   if( allocated( T        ) ) deallocate ( T        )

   if( allocated( Uold     ) ) deallocate ( Uold     )
   if( allocated( Vold     ) ) deallocate ( Vold     )
   if( allocated( Wold     ) ) deallocate ( Wold     )
   if( allocated( TEold    ) ) deallocate ( TEold    )
   if( allocated( EDold    ) ) deallocate ( EDold    )
   if( allocated( Told     ) ) deallocate ( Told     )

   if( allocated( Uold2    ) ) deallocate ( Uold2    )
   if( allocated( Vold2    ) ) deallocate ( Vold2    )
   if( allocated( Wold2    ) ) deallocate ( Wold2    )
   if( allocated( TEold2   ) ) deallocate ( TEold2   )
   if( allocated( EDold2   ) ) deallocate ( EDold2   )
   if( allocated( Told2    ) ) deallocate ( Told2    )

   if( allocated( Variable ) ) deallocate ( Variable )
   if( allocated( Solver   ) ) deallocate ( Solver   )
   if( allocated( Gamma    ) ) deallocate ( Gamma    )
   if( allocated( Scheme   ) ) deallocate ( Scheme   )
   if( allocated( PostC    ) ) deallocate ( PostC    )
   if( allocated( PostV    ) ) deallocate ( PostV    )
   if( allocated( Solve    ) ) deallocate ( Solve    )
   if( allocated( Store    ) ) deallocate ( Store    )

   if( allocated(UseSlopeLimiter)) deallocate ( UseSlopeLimiter )
   if( allocated(URFSlopeLimiter)) deallocate ( URFSlopeLimiter )
   if( allocated(LimiterSlope   )) deallocate ( LimiterSlope    )

   if( allocated( GradVar  ) ) deallocate ( GradVar  )
   if( allocated( nGradVar ) ) deallocate ( nGradVar )
 
   if( allocated( TurbP    ) ) deallocate ( TurbP    )
   if( allocated( DEN      ) ) deallocate ( DEN      )
   if( allocated( CP       ) ) deallocate ( CP       )
   if( allocated( VisEff   ) ) deallocate ( VisEff   )
   if( allocated( dUdX     ) ) deallocate ( dUdX     )
   if( allocated( dVdX     ) ) deallocate ( dVdX     )
   if( allocated( dWdX     ) ) deallocate ( dWdX     )
   if( allocated( dPdX     ) ) deallocate ( dPdX     )
   if( allocated( MassFlux ) ) deallocate ( MassFlux )
   if( allocated( dPhidXo  ) ) deallocate ( dPhidXo  )
   if( allocated( DXdebug  ) ) deallocate ( DXdebug  )
   if( allocated( DXgrad   ) ) deallocate ( DXgrad   )

   if( allocated( Cell     ) ) deallocate ( Cell     )
   if( allocated( NFaces   ) ) deallocate ( NFaces   )
   if( allocated( CFace    ) ) deallocate ( CFace    )
   if( allocated( Face     ) ) deallocate ( Face     )
   if( allocated( RFace    ) ) deallocate ( RFace    )
   if( allocated( Bnd      ) ) deallocate ( Bnd      )
   if( allocated( Reg      ) ) deallocate ( Reg      )
   if( allocated( Vert     ) ) deallocate ( Vert     )

   if( allocated( iFaces%i ) ) deallocate ( iFaces%i    )
   if( allocated( iFaces%list))deallocate ( iFaces%list )
   if( allocated( iNodes%i ) ) deallocate ( iNodes%i    )
   if( allocated( iNodes%list))deallocate ( iNodes%list )
   if( allocated( NNodes   ) ) deallocate ( NNodes      )
 
   if( allocated( Ar       ) ) deallocate ( Ar       )
   if( allocated( Au       ) ) deallocate ( Au       )
   if( allocated( Av       ) ) deallocate ( Av       )
   if( allocated( Aw       ) ) deallocate ( Aw       )
   if( allocated( Su       ) ) deallocate ( Su       )
   if( allocated( Sv       ) ) deallocate ( Sv       )
   if( allocated( Sw       ) ) deallocate ( Sw       )
   if( allocated( Res      ) ) deallocate ( Res      )
   if( allocated( Residual ) ) deallocate ( Residual )
   if( allocated( ResiNorm ) ) deallocate ( ResiNorm )
   if( allocated( Acoo     ) ) deallocate ( Acoo     )
   if( allocated( Arow     ) ) deallocate ( Arow     )
   if( allocated( Acol     ) ) deallocate ( Acol     )
   if( allocated( Acsr     ) ) deallocate ( Acsr     )
   if( allocated( Arwc     ) ) deallocate ( Arwc     )
   if( allocated( Aclc     ) ) deallocate ( Aclc     )

   if( allocated( Work     ) ) deallocate ( Work     )
   if( allocated( RHS      ) ) deallocate ( RHS      )
   if( allocated( SOL      ) ) deallocate ( SOL      )

   !
   ! scalars
   !
   if( allocated( Scalar   ) ) deallocate ( Scalar   )
   if( allocated( SC       ) ) deallocate ( SC       )
   if( allocated( SCold    ) ) deallocate ( SCold    )
   if( allocated( SCold2   ) ) deallocate ( SCold2   )
   if( allocated( ScReg    ) ) deallocate ( ScReg    )
   if( allocated( ScProp   ) ) deallocate ( ScProp   )
   if( allocated( VarS     ) ) deallocate ( VarS     )
   if( allocated( SolveScalar   ) ) deallocate ( SolveScalar )
   if( allocated( StoreScalar   ) ) deallocate ( StoreScalar )

   if( allocated( LimitLowScalar ) ) deallocate ( LimitLowScalar )
   if( allocated( LimitUpScalar  ) ) deallocate ( LimitUpScalar )
   if( allocated( LowerLimitSC   ) ) deallocate ( LowerLimitSC )
   if( allocated( UpperLimitSC   ) ) deallocate ( UpperLimitSC )

   !
   ! particles
   !
   if( UseParticles )then
     write(*,*)'Deallocate particle data'
     do ipart=1,Npart
       call ParticleClearTrack(ipart)
     end do
     do ipart=1,Npart
       if( associated( Tracks(ipart)%head ) ) deallocate ( Tracks(ipart)%head )
       if( associated( Tracks(ipart)%last ) ) deallocate ( Tracks(ipart)%last )
     end do

     if( allocated( Tracks   ) ) deallocate ( Tracks   )
     if( allocated( Particle ) ) deallocate ( Particle )

    !call CollectCells(-1)   ! <=== VISUAL FORTRAN WIL DIT NIET
     call CollectCells(-1,0,0,0)

   endif
   if( allocated( Sensor     ) ) deallocate ( Sensor )
   write(*,*)'Arrays cleaned up'

end subroutine CleanUp
subroutine A33xB3(A,RHS_A)
   
   real, dimension(3,3) :: A
   real, dimension(2,2) :: B
   real, dimension(3)   :: X
   real, dimension(3)   :: RHS_A
   real, dimension(2)   :: RHS_B

   logical :: UseX, UseY, UseZ
   
   UseX = .false.
   UseY = .false.
   UseZ = .false.
      
   X = 0.0

   irang = 0
   if( A(1,1) /= 0.0 )then  ! altijd >= 0
     UseX = .true.
     irang = irang + 1
   else
     X(1) = 0.0
   endif    
   if( A(2,2) /= 0.0 )then  ! altijd >= 0
     UseY = .true.
     irang = irang + 1
   else
     X(2) = 0.0
   endif  
   if( A(3,3) /= 0.0 )then  ! altijd >= 0
     UseZ = .true.
     irang = irang + 1
   else
     X(3) = 0.0
   endif  

   if( irang <= 1 )then
     !
     ! too simple, do nothing
     !

     X = 0.0

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
       write(*,1) 1
     endif
     if( B(2,2) /= 0.0 )then
       p2 = - B(1,2)/B(2,2)
       B(1,1) = B(1,1) + p2 * B(2,1)
       B(1,2) = B(1,2) + p2 * B(2,2)

       RHS_B(1) = RHS_B(1) + p2 * RHS_B(2)

     else
       write(*,1) 2
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
       write(*,1) 3
     endif
     if( B(2,2) /= 0.0 )then
       p4 = 1.0/B(2,2)
       B(2,1) = p4 * B(2,1)
       B(2,2) = p4 * B(2,2)

       RHS_B(2) = p4 * RHS_B(2)

     else
       write(*,1) 4
     endif
     !
     ! de gradient
     !
     if( UseX .and. UseY )then
       X(1) = RHS_B(1)
       X(2) = RHS_B(2)
       X(3) = 0.0
     else if( UseX .and. UseZ )then
       X(1) = RHS_B(1)
       X(2) = 0.0
       X(3) = RHS_B(2)
     else if( UseY .and. UseZ )then
       X(1) = 0.0
       X(2) = RHS_B(1)
       X(3) = RHS_B(2)
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
       write(*,1) 5
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
       write(*,1) 6
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
       write(*,1) 7
     endif
     !
     ! normeren
     !
     if( A(1,1) /= 0.0 )then
       p7 = 1.0/A(1,1)
       A(1,1) = p7 * A(1,1)
       A(1,2) = p7 * A(1,2)
       A(1,3) = p7 * A(1,3)
       RHS_A(1) = p7 * RHS_A(1)
     else
       write(*,1) 8
     endif
     if( A(2,2) /= 0.0 )then
       p8 = 1.0/A(2,2)
       A(2,1) = p8 * A(2,1)
       A(2,2) = p8 * A(2,2)
       A(2,3) = p8 * A(2,3)
       RHS_A(2) = p8 * RHS_A(2)
     else
       write(*,1) 8
     endif
     if( A(3,3) /= 0.0 )then
       p9 = 1.0/A(3,3)
       A(3,1) = p9 * A(3,1)
       A(3,2) = p9 * A(3,2)
       A(3,3) = p9 * A(3,3)
       RHS_A(3) = p9 * RHS_A(3)
     else
       write(*,1) 9
     endif
     !
     ! de gradient
     !
     X(1) = RHS_A(1)
     X(2) = RHS_A(2)
     X(3) = RHS_A(3)
   
   else
     write(*,*)'+ internal error A33xB3 catch 22:',irang
   endif
   
   RHS_A(1) = X(1)
   RHS_A(2) = X(2)
   RHS_A(3) = X(3)

 1 format('+ internal error: A33xB3 ',i2)
   
end subroutine A33xB3
