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
subroutine GradientPhi(ivar,Phi,dPhidX)
!========================================================================
!
!  Main calling subroutine to calculate the gradient of 'phi'.
!  Two methods are available:
!  1) using Gauss theorem, and
!  2) the least squares method
!
!  Once a gradient has been computed an optional check if the gradient
!  is reasonable is possible and if the gradient is too steep then
!  it is limited. Names associated to these slope limiters are:
!  Barth and Jespersen, Venkatakrishnan, Mavriplis, Aftosmis etc.
!  Slope limiters are especially useful for all-tetrahedral meshes
!  are used (only four relations for a 3D gradient, near walls the
!  the situation is worse), and for bad/awkward meshes.
!
   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd), intent(IN)    :: Phi
   real, dimension(Ncel+Nbnd,3), intent(OUT) :: dPhidX

   integer, intent(IN) :: ivar

  !call watch_enter('GradientPhi')

  !if( GradAlg == GradLS )then
  !  call GradientPhiLeastSquares(ivar,Phi,dPhidX)
  !else if( GradAlg == GradGauss )then
  !  call GradientPhiGauss(ivar,Phi,dPhidX)
  !else
  !  stop '*** internal error invalid gradient algorithm'
  !endif

   if( GradVar(iVar) == GradLS )then
     call GradientPhiLeastSquares(ivar,Phi,dPhidX)
   else if( GradVar(iVar) == GradGauss )then
     call GradientPhiGauss(ivar,Phi,dPhidX)
   else
     stop '*** internal error invalid gradient algorithm'
   endif
   
   if( UseSlopeLimiter(iVar) )then
     if( LimiterSlope(iVar) >= SlopeLimiterBJc .and. &
         LimiterSlope(iVar) <= SlopeLimiterBJn )then
       call GradientPhiLimiterBarthJespersen(ivar,LimiterSlope(iVar),Phi,dPhidX) 
     else if( LimiterSlope(iVar) >= SlopeLimiterVNc .and. &
         LimiterSlope(iVar) <= SlopeLimiterVNn )then
       call GradientPhiLimiterVenkatarishnan(ivar,LimiterSlope(iVar),Phi,dPhidX) 
     else if( LimiterSlope(iVar) >= SlopeLimiterVAc .and. &
         LimiterSlope(iVar) <= SlopeLimiterVAn )then
       call GradientPhiLimiterAlbada(ivar,LimiterSlope(iVar),Phi,dPhidX) 
     else if( LimiterSlope(iVar) >= SlopeLimiterP1c .and. &
         LimiterSlope(iVar) <= SlopeLimiterP1n )then
       call GradientPhiLimiterPolynomial(ivar,LimiterSlope(iVar),Phi,dPhidX) 
     endif
   endif
   
  !call watch_leave('GradientPhi')

end subroutine GradientPhi
subroutine GradientPhiGauss(ivar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables
   use watches
   
   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX

   real :: facn, facp
   real :: Xac(3), dPhidXac(3)

   real, dimension(3) :: Xp, Xn, Xf, Xnorm, Xpac, Xnac, delp, deln

   character(len=16) :: string

   if( Debug > 3 ) write(IOdef,*)'*** GradientPhiGauss: ',Variable(ivar)
   write(IOdbg,*)'*** GradientPhiGauss: ',Variable(ivar),Ngradient

   write(string,'(''Gauss'',A11)') Variable(ivar)
   call watch_enter(string)
!  call watch_enter('GradientPhiGauss')

   dPhidXo = 0.0
   !
   ! iterative calculation of gradients
   !
   gradientloop: do igrad=1,Ngradient
     dPhidX  = 0.0
     faceloop: do i=1,Nfac
       ip = Face(i)%cell1
       in = Face(i)%cell2

       if( in > 0 )then
         !
         ! internal cell face which points to two cells
         !
         facn = Face(i)%lambda
         facp = 1.0 - facn
         !
         ! for the next section see the discussion in the text
         !
         Xac      =     Cell(in)%x * facn +    Cell(ip)%x * facp
         dPhidXac =  dPhidXo(in,:) * facn + dPhidXo(ip,:) * facp
         !
         ! now the gradient at the shifted position is known
         ! correct the value at the cell face center
         !
         PhiFace = Phi(in) * facn + Phi(ip) * facp            ! the standard

         delta   = dot_product( dPhidXac , Face(i)%x - Xac )  ! correction
         PhiFace = PhiFace + delta

         !
         ! now only the value at the face center is known
         ! multiply it by the area-components
         ! this is basically Gauss' theorem
         !
         dPhidX(ip,:) = dPhidX(ip,:) + PhiFace * Face(i)%n    ! normal p -> n
         dPhidX(in,:) = dPhidX(in,:) - PhiFace * Face(i)%n    ! reversed

       else
         !
         ! boundary face
         !
         ib = Face(i)%bnd

         PhiFace = Phi(Ncel+ib)
         dPhidX(ip,:) = dPhidX(ip,:) + PhiFace * Face(i)%n

       endif

     end do faceloop

    !dPhidX(1:Ncel,1) = dPhidX(1:Ncel,1) / Cell(1:Ncel)%Vol
    !dPhidX(1:Ncel,2) = dPhidX(1:Ncel,2) / Cell(1:Ncel)%Vol
    !dPhidX(1:Ncel,3) = dPhidX(1:Ncel,3) / Cell(1:Ncel)%Vol

    !do i=1,Ncel
    !  fact = 1.0/Cell(i)%Vol              ! not faster
    !  dPhidX(i,:) =  fact*dPhidX(i,:)
    !end do
     
     do i=1,Ncel
       fact = 1.0/Cell(i)%Vol               
       dPhidX(i,1) =  fact*dPhidX(i,1)
       dPhidX(i,2) =  fact*dPhidX(i,2)
       dPhidX(i,3) =  fact*dPhidX(i,3)
     end do

     dPhidXo = dPhidX
     !dPhidXo = dPhidXo + 0.95*( dPhidX - dPhidXo) ! under relaxation

   end do gradientloop

   if( allocated(DXgrad) ) DXgrad = dPhidX

   if( Debug > 2 .and. iVar==VarT )then
     write(IOdef,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,1)),' <  dPhi/dX  < ',maxval(dPhidX(1:Ncel,1))
     write(IOdef,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,2)),' <  dPhi/dY  < ',maxval(dPhidX(1:Ncel,2))
     write(IOdef,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,3)),' <  dPhi/dZ  < ',maxval(dPhidX(1:Ncel,3))
   endif

   call watch_leave('GradientPhiGauss')

   if( Debug > 3 ) write(IOdef,*)'=== GradientPhi: ',Variable(ivar)

end subroutine GradientPhiGauss
subroutine GradientPhiLeastSquares(ivar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables
   use watches
   
   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX

   real, dimension(3)   :: GradPhi, Xf, Xac, Xp, dX, ds
   real, dimension(3,3) :: A
   real, dimension(3)   :: RHS_A

   integer, dimension(3):: IPIV    ! LAPACK

   character(len=16) :: string

   if( Debug > 3 ) write(IOdef,*)'*** GradientPhiLeastSquares ',Variable(ivar)
   write(IOdbg,*)'*** GradientPhiLeastSquares ',Variable(ivar)

   write(string,'(''Least'',A11)') Variable(ivar)
   call watch_enter(string)
  !call watch_enter('GradientPhiLeastSquares')

   dPhidX = 0.0

   do ip=1,Ncel

     Xp = Cell(ip)%x
     A     = 0.0
     RHS_A = 0.0

     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2

       if( ipn > 0 )then
         !
         ! original method Sum from P to N
         !
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

       else
         ! boundary face

         ib = Face(k)%bnd
         ir = Bnd(ib)%rid
         it = Reg(ir)%typ

         dPhi = Phi(Ncel+ib) - Phi(ipp)
         dX   = Face(k)%x - Cell(ipp)%x

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

     if( UseLapack )then

       call SGESV( 3, 1, A, 3, IPIV, RHS_A, 3, INFO )

       if( info /= 0 ) write(IOdef,*)'Lapack (sgesv) info(1):',INFO

       GradPhi(1) = RHS_A(1)
       GradPhi(2) = RHS_A(2)
       GradPhi(3) = RHS_A(3)

       dPhidX(ip,:) = GradPhi

     else

       call A33xB3(A,RHS_A)

       GradPhi(1) = RHS_A(1)
       GradPhi(2) = RHS_A(2)
       GradPhi(3) = RHS_A(3)

       dPhidX(ip,:) = GradPhi

     endif

   end do

  !if( allocated(DXdebug) .and. allocated(T) ) DXdebug =T(1:Ncel)
  !if( allocated(DXgrad) )                     DXgrad = dPhidX

   if( Debug > 2 .and. ivar == vart )then
     write(IOdef,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,1)),' <  dPhi/dX  < ',maxval(dPhidX(1:Ncel,1))
     write(IOdef,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,2)),' <  dPhi/dY  < ',maxval(dPhidX(1:Ncel,2))
     write(IOdef,'(1x,1pe10.3,A,e10.3)') &
       minval(dPhidX(1:Ncel,3)),' <  dPhi/dZ  < ',maxval(dPhidX(1:Ncel,3))
   endif

   call watch_leave('GradientPhiLeastSquares')

   if( Debug > 3 ) write(IOdef,*)'=== GradientPhiLeastSquares'

end subroutine GradientPhiLeastSquares
subroutine GradientPhiLeastSquaresN(ivar,Phi,dPhidX)
!========================================================================

   use constants
   use geometry
   use variables

   real, dimension(Ncel+Nbnd)   :: Phi
   real, dimension(Ncel+Nbnd,3) :: dPhidX

   real, dimension(3)   :: GradPhi, Xf, Xac, Xp, dX, ds
   real, dimension(3,3) :: A
   real, dimension(3)   :: RHS_A

   integer, dimension(3):: IPIV    ! LAPACK

   if( Debug > 1 ) write(IOdef,*)'*** GradientPhiLeastSquaresN ',Variable(ivar)

   dPhidXo = 0.0
   !
   ! iterative calculation of gradients
   !
   do igrad=1,Ngradient
     write(IOdef,*)'loop ',igrad
     dPhidX  = 0.0
     do ip=1,Ncel

       Xp = Cell(ip)%x
       A  = 0.0
       RHS_A = 0.0

       do j=1,Nfaces(ip)
         k  = CFace(ip,j)
         ipp = Face(k)%cell1
         ipn = Face(k)%cell2

         if( ipn > 0 )then

           facn = Face(k)%lambda
           facp = 1.0 - facn
           Xac  = Cell(ipn)%x * facn + Cell(ipp)%x * facp
           Xf   = Face(k)%x

           Phiac   = Phi(ipn) * facn + Phi(ipp) * facp
           GradPhi = dPhidXo(ipn,:) * facn + dPhidXo(ipp,:) * facp
           delta   = dot_product( GradPhi , Xf - Xac )  ! correction
           PhiFace = Phiac + delta

           !PhiFace = 0.5*( Phi(ipp) + Phi(ipn) + &
           !dot_product( dPhidXo(ipp,:) , Xf - Cell(ipp)%x ) + &
           !dot_product( dPhidXo(ipn,:) , Xf - Cell(ipn)%x ) )

           if( ipp == ip )then
             dPhi = PhiFace - Phi(ipp)
             dX   = Xf - Cell(ipp)%x
           else
             dPhi = Phi(ipp) - PhiFace
             dX   = Cell(ipp)%x - Xf
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

         else
           ! boundary face

           ib = Face(k)%bnd
           ir = Bnd(ib)%rid
           it = Reg(ir)%typ

           dPhi = Phi(Ncel+ib) - Phi(ipp)
           dX   = Face(k)%x - Cell(ipp)%x

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

       if( UseLapack )then
         call SGESV( 3, 1, A, 3, IPIV, RHS_A, 3, INFO )

         if( info /= 0 ) write(IOdef,*)'Lapack (sgesv) info(1):',INFO

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

         dPhidX(ip,:) = GradPhi

       else
         call A33xB3(A,RHS_A)

         GradPhi(1) = RHS_A(1)
         GradPhi(2) = RHS_A(2)
         GradPhi(3) = RHS_A(3)

         dPhidX(ip,:) = GradPhi

       endif

     end do ! cell loop

     dPhidXo = dPhidX
     !dPhidXo = dPhidXo + 0.75*( dPhidX - dPhidXo) ! under relaxation

     if( Debug > 2 )then
       write(IOdef,'(1x,1pe10.3,A,e10.3)') &
         minval(dPhidX(1:Ncel,1)),' <  dPhi/dX  < ',maxval(dPhidX(1:Ncel,1))
       write(IOdef,'(1x,1pe10.3,A,e10.3)') &
         minval(dPhidX(1:Ncel,2)),' <  dPhi/dY  < ',maxval(dPhidX(1:Ncel,2))
       write(IOdef,'(1x,1pe10.3,A,e10.3)') &
         minval(dPhidX(1:Ncel,3)),' <  dPhi/dZ  < ',maxval(dPhidX(1:Ncel,3))
     endif

   end do   ! gradient loop

   if( Debug > 3 ) write(IOdef,*)'=== GradientPhiLeastSquaresN'

end subroutine GradientPhiLeastSquaresN
!========================================================================
!========================================================================
! 
!   S L O P E  L I M I T E R S
!
!========================================================================
!========================================================================
subroutine GradientPhiLimiterBarthJespersen(ivar,limiter,Phi,dPhidX)
!========================================================================
!
!  AIAA-89-0366, The design and application of upwind schemes
!  on unstructured meshes, T.J.Barth, D.C.Jespersen, 1989
!
!  Basically: test the calculated gradient against the original
!  neighbor values.
!
!  Three different approaches are possible:
!  1) Test the gradient in the cell centre against the values 
!     in the neighboring centroids
!  2) Project the value of the neighbor on to the face center and 
!     use that value
!  3) Same as 2) but using the nodes instead of the face centres.
!
!  Approach 1) damps the gradient, 2) tends to overshoot and 
!  3) is somewere in the middle. Because 3) uses the nodes it 
!  is NOT suited for '2D' cases as the 'model thickness' enters 
!  the test.
!
   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd), intent(IN)    :: Phi
   real, dimension(Ncel+Nbnd,3), intent(OUT) :: dPhidX

   integer, intent(IN) :: ivar

   real, dimension(3) :: dS, xpn

   write(IOdbg,*) 'BarthJespersen: ',variable(ivar)
   
   call watch_enter('BarthJespersen')

   do ip=1,Ncel
     phimax = Phi(ip)
     phimin = Phi(ip)
     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2
       if( ipn > 0 )then
         if( ipp == ip )then
           phimax = max(phimax,Phi(ipn))
           phimin = min(phimin,Phi(ipn))  
         else
           phimax = max(phimax,Phi(ipp))
           phimin = min(phimin,Phi(ipp))  
         endif
       else
         ib = Face(k)%bnd
         phimax = max(phimax,Phi(Ncel+ib))
         phimin = min(phimin,Phi(Ncel+ib))     
       endif
     end do

     deltamax = phimax - Phi(ip)        
     deltamin = phimin - Phi(ip)

     alpha = 1.0
     ds    = 0.0
     xpn   = 0.0
     
     if( limiter == SlopeLimiterBJN )then
       !
       ! node based
       !
       do j=1,NItemsFromList(iNodes,ip)
         
         iv  = ItemFromList(iNodes,ip,j)

         ds  = Vert(iv,:) - Cell(ip)%x 
         delta_face = dot_product( dPhidX(ip,:) , ds )

         if( abs(delta_face) < 1.e-6 )then
           r = 1.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         alpha = min( alpha , r )

       end do     
     else
       !
       ! face or centroid based
       !
       do j=1,Nfaces(ip)
         k  = CFace(ip,j)

         if( limiter == SlopeLimiterBJF )then 

           ds = Face(k)%x - Cell(ip)%x 
           delta_face = dot_product( dPhidX(ip,:) , ds )

         else if( limiter == SlopeLimiterBJC )then 

           ipp = Face(k)%cell1
           ipn = Face(k)%cell2
           if( ipn > 0 )then
             if( ipp == ip )then
               xpn = Cell(ipn)%x - Cell(ipp)%x
             else
               xpn = Cell(ipp)%x - Cell(ipn)%x
             endif
           else
             xpn = Face(k)%x - Cell(ip)%x
           endif

           delta_face = dot_product( dPhidX(ip,:) , xpn )

         else

           write(IOdef,*)'+++ internal error: BJ limiter = ',limiter
           delta_face = 0.0
           
         endif

         if( abs(delta_face) < 1.e-6 )then
           r = 1.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         alpha = min( alpha , r )

       end do
     endif

     dPhidX(ip,:) = alpha*dPhidX(ip,:)

    !if( allocated(DXdebug) .and. iVar == VarT ) DXdebug(ip) = alpha

   end do

  !if( allocated(DXgrad) .and. iVar == VarT ) DXgrad = dPhidX
  !if( allocated(DXdebug) .and. iVar == VarT ) &
  !  write(*,*) minval(DXdebug(1:Ncel)),' < alpha < ',maxval(DXdebug(1:Ncel))

   call watch_leave('BarthJespersen')

end subroutine GradientPhiLimiterBarthJespersen
subroutine GradientPhiLimiterVenkatarishnan(ivar,limiter,Phi,dPhidX)
!========================================================================
!
!  AIAA-93-0880, On the accuracy of limiters and convergence
!  to steady state solutions, V.Venkatakrishnan, 1993
!
!  Basically: test the calculated gradient against the original
!  neighbor values as with Barth and Jespersen. However BJ uses
!  the min-function Venkatarishnan introduces instead of min(1,y):
!
!  phi(y) = (y**2+2*y)/(y**2+y+2)
!
!  Again three different approaches are possible:
!  1) Test the gradient in the cell centre against the values 
!     in the neighboring centroids
!  2) Project the value of the neighbor on to the face center and 
!     use that value
!  3) Same as 2) but using the nodes instead of the face centres.
!
!  Approach 1) damps the gradient, 2) tends to overshoot and 
!  3) is somewere in the middle. Because 3) uses the nodes it 
!  is NOT suited for '2D' cases as the 'model thickness' enters 
!  the test.
!
   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd), intent(IN)    :: Phi
   real, dimension(Ncel+Nbnd,3), intent(OUT) :: dPhidX

   integer, intent(IN) :: ivar

   real, dimension(3) :: dS, xpn

   write(IOdbg,*) 'Venkatarishnan: ',variable(ivar)
   
   call watch_enter('Venkatarishnan')

   alphamin = Large 
   do ip=1,Ncel
     phimax = Phi(ip)
     phimin = Phi(ip)
     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2
       if( ipn > 0 )then
         if( ipp == ip )then
           phimax = max(phimax,Phi(ipn))
           phimin = min(phimin,Phi(ipn))  
         else
           phimax = max(phimax,Phi(ipp))
           phimin = min(phimin,Phi(ipp))  
         endif
       else
         ib = Face(k)%bnd
         phimax = max(phimax,Phi(Ncel+ib))
         phimin = min(phimin,Phi(Ncel+ib))     
       endif
     end do

     deltamax = phimax - Phi(ip)        
     deltamin = phimin - Phi(ip)

     alpha = 1.0
     ds    = 0.0
     xpn   = 0.0
     
     if( limiter == SlopeLimiterVNN )then
       !
       ! node based
       !
       do j=1,NItemsFromList(iNodes,ip)
         
         iv  = ItemFromList(iNodes,ip,j)

         ds  = Vert(iv,:) - Cell(ip)%x 
         delta_face = dot_product( dPhidX(ip,:) , ds )

         if( abs(delta_face) < 1.e-6 )then
           r = 1000.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         alpha = min( alpha , (r**2+2.0*r)/(r**2+r+2.0) )

       end do
     
     else
       !
       ! face or centroid based
       !
       do j=1,Nfaces(ip)
         k  = CFace(ip,j)

         if( limiter == SlopeLimiterVNF )then 

           ds = Face(k)%x - Cell(ip)%x 
           delta_face = dot_product( dPhidX(ip,:) , ds )

         else if( limiter == SlopeLimiterVNC )then 

           ipp = Face(k)%cell1
           ipn = Face(k)%cell2
           if( ipn > 0 )then
             if( ipp == ip )then
               xpn = Cell(ipn)%x - Cell(ipp)%x
             else
               xpn = Cell(ipp)%x - Cell(ipn)%x
             endif
           else
             xpn = Face(k)%x - Cell(ip)%x
           endif

           delta_face = dot_product( dPhidX(ip,:) , xpn )

         else

           write(IOdef,*)'+++ internal error: VN limiter = ',limiter
           delta_face = 0.0
           
         endif

         if( abs(delta_face) < 1.e-6 )then
           r = 1000.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         alpha = min( alpha , (r**2+2.0*r)/(r**2+r+2.0) )

       end do
     endif

     dPhidX(ip,:) = alpha*dPhidX(ip,:)

     if( allocated(DXdebug) .and. iVar == VarT ) DXdebug(ip) = alpha

     alphamin = min(alphamin,alpha)

   end do

   if( allocated(DXgrad) .and. iVar == VarT ) DXgrad = dPhidX

  !if( allocated(DXdebug) .and. iVar == VarT ) &
  !  write(*,*) minval(DXdebug(1:Ncel)),' < alpha < ',maxval(DXdebug(1:Ncel))

   call watch_leave('Venkatarishnan')

end subroutine GradientPhiLimiterVenkatarishnan
subroutine GradientPhiLimiterAlbada(ivar,limiter,Phi,dPhidX)
!========================================================================
!
!  AIAA-93-0880, On the accuracy of limiters and convergence
!  to steady state solutions, V.Venkatakrishnan, 1993
!
!  Basically: test the calculated gradient against the original
!  neighbor values as with Barth and Jespersen. Instead of
!  BJ min(1,y) use van Albada's:
!
!  phi(y) = (y**2+y)/(y**2+1)
!
!  Again three different approaches are possible:
!  1) Test the gradient in the cell centre against the values 
!     in the neighboring centroids
!  2) Project the value of the neighbor on to the face center and 
!     use that value
!  3) Same as 2) but using the nodes instead of the face centres.
!
!  Approach 1) damps the gradient, 2) tends to overshoot and 
!  3) is somewere in the middle. Because 3) uses the nodes it 
!  is NOT suited for '2D' cases as the 'model thickness' enters 
!  the test.
!
   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd), intent(IN)    :: Phi
   real, dimension(Ncel+Nbnd,3), intent(OUT) :: dPhidX

   integer, intent(IN) :: ivar

   real, dimension(3) :: dS, xpn

   write(IOdbg,*) 'Albada: ',variable(ivar)
   
   call watch_enter('Albada')

   alphamin = Large 
   do ip=1,Ncel
     phimax = Phi(ip)
     phimin = Phi(ip)
     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2
       if( ipn > 0 )then
         if( ipp == ip )then
           phimax = max(phimax,Phi(ipn))
           phimin = min(phimin,Phi(ipn))  
         else
           phimax = max(phimax,Phi(ipp))
           phimin = min(phimin,Phi(ipp))  
         endif
       else
         ib = Face(k)%bnd
         phimax = max(phimax,Phi(Ncel+ib))
         phimin = min(phimin,Phi(Ncel+ib))     
       endif
     end do

     deltamax = phimax - Phi(ip)        
     deltamin = phimin - Phi(ip)

     alpha = 1.0
     ds    = 0.0
     xpn   = 0.0
     
     if( limiter == SlopeLimiterVAN )then
       !
       ! node based
       !
       do j=1,NItemsFromList(iNodes,ip)
         
         iv  = ItemFromList(iNodes,ip,j)

         ds  = Vert(iv,:) - Cell(ip)%x 
         delta_face = dot_product( dPhidX(ip,:) , ds )

         if( abs(delta_face) < 1.e-6 )then
           r = 1.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         alpha = min( alpha , (r**2+r)/(r**2+1.0) )

       end do
     
     else
       !
       ! face or centroid based
       !
       do j=1,Nfaces(ip)
         k  = CFace(ip,j)

         if( limiter == SlopeLimiterVAF )then 

           ds = Face(k)%x - Cell(ip)%x 
           delta_face = dot_product( dPhidX(ip,:) , ds )

         else if( limiter == SlopeLimiterVAC )then 

           ipp = Face(k)%cell1
           ipn = Face(k)%cell2
           if( ipn > 0 )then
             if( ipp == ip )then
               xpn = Cell(ipn)%x - Cell(ipp)%x
             else
               xpn = Cell(ipp)%x - Cell(ipn)%x
             endif
           else
             xpn = Face(k)%x - Cell(ip)%x
           endif

           delta_face = dot_product( dPhidX(ip,:) , xpn )

         else

           write(IOdef,*)'+++ internal error: vA limiter = ',limiter
           delta_face = 0.0
           
         endif

         if( abs(delta_face) < 1.e-6 )then
           r = 1.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         alpha = min( alpha , (r**2+r)/(r**2+1.0) )

       end do
     endif

     dPhidX(ip,:) = alpha*dPhidX(ip,:)

     if( allocated(DXdebug) .and. iVar == VarT ) DXdebug(ip) = alpha

     alphamin = min(alphamin,alpha)

   end do

   if( allocated(DXgrad) .and. iVar == VarT ) DXgrad = dPhidX

  !if( allocated(DXdebug) .and. iVar == VarT ) &
  !  write(*,*) minval(DXdebug(1:Ncel)),' < alpha < ',maxval(DXdebug(1:Ncel))

   call watch_leave('Albada')

end subroutine GradientPhiLimiterAlbada
subroutine GradientPhiLimiterPolynomial(ivar,limiter,Phi,dPhidX)
!========================================================================
!
!  K.Michalak, C. Ollivier-Gooch
!  Limiters for UnstructuredHigher-Order Accurate
!  Solutions of the Euler Equations
!
!  cubic polynomial: f(y)= A3 y**3 + A2 * y**2 + A1 * y + A0 
!
!  target: y=0, f(y)=0 => A0 = 0
!          df(y=0)/dy=1 => A1 = 1 
!
!  switch at point Yc where f(Yc)=1, df(y=Yc)/dy=0 
!
!  Again three different approaches are possible:
!  1) Test the gradient in the cell centre against the values 
!     in the neighboring centroids
!  2) Project the value of the neighbor on to the face center and 
!     use that value
!  3) Same as 2) but using the nodes instead of the face centres.
!
!  Approach 1) damps the gradient, 2) tends to overshoot and 
!  3) is somewere in the middle. Because 3) uses the nodes it 
!  is NOT suited for '2D' cases as the 'model thickness' enters 
!  the test.
!
   use constants
   use geometry
   use variables
   use watches

   real, dimension(Ncel+Nbnd), intent(IN)    :: Phi
   real, dimension(Ncel+Nbnd,3), intent(OUT) :: dPhidX

   integer, intent(IN) :: ivar

   real, dimension(3) :: dS, xpn

   write(IOdbg,*) 'Polynomial1: ',variable(ivar)
   
   call watch_enter('Polynomial1')

   Xc = 1.5 
   A3 = (Xc-2)/(Xc**3)
   A2 = -1.0*(1.0+3.0*A3*Xc**2)/(2.0*Xc) 

   alphamin = Large 
   do ip=1,Ncel
     phimax = Phi(ip)
     phimin = Phi(ip)
     do j=1,Nfaces(ip)
       k  = CFace(ip,j)
       ipp = Face(k)%cell1
       ipn = Face(k)%cell2
       if( ipn > 0 )then
         if( ipp == ip )then
           phimax = max(phimax,Phi(ipn))
           phimin = min(phimin,Phi(ipn))  
         else
           phimax = max(phimax,Phi(ipp))
           phimin = min(phimin,Phi(ipp))  
         endif
       else
         ib = Face(k)%bnd
         phimax = max(phimax,Phi(Ncel+ib))
         phimin = min(phimin,Phi(Ncel+ib))     
       endif
     end do

     deltamax = phimax - Phi(ip)        
     deltamin = phimin - Phi(ip)

     alpha = 1.0
     ds    = 0.0
     xpn   = 0.0
     
     if( limiter == SlopeLimiterP1N )then
       !
       ! node based
       !
       do j=1,NItemsFromList(iNodes,ip)
         
         iv  = ItemFromList(iNodes,ip,j)

         ds  = Vert(iv,:) - Cell(ip)%x 
         delta_face = dot_product( dPhidX(ip,:) , ds )

         if( abs(delta_face) < 1.e-6 )then
           r = 100.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         if( r >= Xc )then
	   alpha = min( alpha , 1.0 )
	 else
	   alpha = min( alpha , A3*r**3+A2*r**2+r ) 
         endif
	 
       end do
     
     else
       !
       ! face or centroid based
       !
       do j=1,Nfaces(ip)
         k  = CFace(ip,j)

         if( limiter == SlopeLimiterP1F )then 

           ds = Face(k)%x - Cell(ip)%x 
           delta_face = dot_product( dPhidX(ip,:) , ds )

         else if( limiter == SlopeLimiterP1C )then 

           ipp = Face(k)%cell1
           ipn = Face(k)%cell2
           if( ipn > 0 )then
             if( ipp == ip )then
               xpn = Cell(ipn)%x - Cell(ipp)%x
             else
               xpn = Cell(ipp)%x - Cell(ipn)%x
             endif
           else
             xpn = Face(k)%x - Cell(ip)%x
           endif

           delta_face = dot_product( dPhidX(ip,:) , xpn )

         else

           write(IOdef,*)'+++ internal error: P1 limiter = ',limiter
           delta_face = 0.0
           
         endif

         if( abs(delta_face) < 1.e-6 )then
           r = 100.0
         else if( delta_face > 0.0 )then
           r = deltamax/delta_face
         else
           r = deltamin/delta_face
         endif

         if( r >= Xc )then
	   alpha = min( alpha , 1.0 )
	 else
	   alpha = min( alpha , A3*r**3+A2*r**2+r ) 
         endif

       end do
     endif

     dPhidX(ip,:) = alpha*dPhidX(ip,:)

     if( allocated(DXdebug) .and. iVar == VarT ) DXdebug(ip) = alpha

     alphamin = min(alphamin,alpha)

   end do

   if( allocated(DXgrad) .and. iVar == VarT ) DXgrad = dPhidX

  !if( allocated(DXdebug) .and. iVar == VarT ) &
  !  write(*,*) minval(DXdebug(1:Ncel)),' < alpha < ',maxval(DXdebug(1:Ncel))

   call watch_leave('Polynomial1')

end subroutine GradientPhiLimiterPolynomial
