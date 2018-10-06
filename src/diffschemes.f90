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
subroutine SelectDiffSchemeScalar(iFace,iScheme,iP,iN, &
                                  Phi,dPhidX,PhiFace)
!========================================================================
!
!  selection of differencing schemes pulled out from
!  FluxScalar and FluxUVW. Two forms: a scalar version (this) and
!  a vector version 
!
!  all the trouble is about finding the value of Phi at the face
!
   use geometry
   use variables
   use constants

   real, intent(OUT)   :: PhiFace

   integer, intent(IN)                      :: iFace, iScheme, iP, iN   
   real, dimension(Ncel+Nbnd), intent(IN)   :: Phi
   real, dimension(Ncel+Nbnd,3), intent(IN) :: dPhidX

   real, dimension(3)  :: Xac, Xp, Xn, Xf, dPhiC, dPhiD, dPhidXac, &
                          d, d1, d2
      
   if( Debug > 3 ) write(IOdef,*)'*** SelectDiffSchemeScalar'

   Xp = Cell(ip)%x
   Xn = Cell(in)%x

   PhiP  = Phi(ip)
   PhiN  = Phi(in)
   !
   ! set the upwind C and downwind D values
   ! and PhiC-tilde
   !
   if( MassFlux(iFace) >= 0.0 )then
     iC    = iP

     PhiC  = PhiP
     PhiD  = PhiN
     dPhiC = dPhidX(ip,:)
     dPhiD = dPhidX(in,:)

     dPhi  = PhiD - PhiC
     d     = Xn - Xp
     ddPhi = 2.0 * dot_product( dPhiC , d )
     PhiCt = 1.0 - dphi/(ddphi+Small)
   else
     iC    = iN

     PhiC  = PhiN
     PhiD  = PhiP
     dPhiC = dPhidX(in,:)
     dPhiD = dPhidX(ip,:)

     dPhi  = PhiD - PhiC
     d     = Xp - Xn
     ddPhi = 2.0 * dot_product( dPhiC , d )
     PhiCt = 1.0 - dphi/(ddphi+Small)
   endif
   !
   ! now we know the up- and downwind values of Phi: PhiC resp. PhiD
   !
   ! most popular/used schemes are: 
   !
   !   UDS, CDS(1), Gamma, LUD, Minmod, CD1, CD2, CD3
   !
   ! so arrange them in the if-then-else tree accordingly
   !
   select case(iScheme)
     case( DSud ) 

       PhiFace = PhiC                                  ! standard upwind
       call DiffSchemeCounter(1,iScheme,1,0)

     case( DScds )

       PhiFace = SetFaceCDS(iFace,DScdDefault,iP,iN,Phi,dPhidX)
       call DiffSchemeCounter(1,DScds,1,0)

     case( DSgam )

       PhiCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,Phi,dPhidX)

       Beta_m = 0.1                                    ! Jasak's constant
       if( PhiCt <= 0.0 .or. PhiCt >= 1.0 )then        ! 0.1 <= Beta_m <= 0.5
                                                       !
         PhiFace = PhiC                                ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
                                                       !
       elseif( Beta_m <= PhiCt .and. PhiCt < 1.0 )then !
                                                       !
         PhiFace = PhiCDS                              ! select CDS
         call DiffSchemeCounter(1,DScds,1,0)           !
                                                       !
       elseif( 0.0 < PhiCt .and. PhiCt < Beta_m )then  !
                                                       !
         fg = PhiCt / Beta_m                           !
         PhiFace = fg*PhiCDS + (1.0-fg)*PhiC           ! select gamma
         call DiffSchemeCounter(1,DSgam,1,0)           !

       endif

     case( DSlud )                                     ! LUDS Scheme:
                                                       !  
       if( PhiCt <= 0.0 .or. PhiCt >= 1.0 )then        !
                                                       !
         PhiFace = PhiC                                ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !

       else                                         

         Xf      = Face(iFace)%x                    
         d       = Xf - Cell(iC)%x                  
         PhiFace = PhiC + dot_product( dPhiC , d )  
         call DiffSchemeCounter(1,DSlud,1,0)        

       endif

     case( DSmod )                                     ! MinMod Scheme:

       PhiCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,Phi,dPhidX)

       if( PhiCt <= 0.0 .or. PhiCt >= 1.0 )then        !
                                                       !
         PhiFace = PhiC                                ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
                                                       !
       elseif( PhiCt < 0.5 )then                       !
                                                       !
         Xf      = Face(iFace)%x                       !   
         d       = Xf - Cell(iC)%x                     !
         PhiFace = PhiC + dot_product( dPhiC , d )     ! select LUDS
         call DiffSchemeCounter(1,DSlud,1,0)           !
                                                       !
       else                                            !
                                                       !
         PhiFace = PhiCDS                              ! select CDS
         call DiffSchemeCounter(1,DScds,1,0)           !
                                                       !
       endif                                           !

     case( DScd1 )

       PhiFace = SetFaceCDS(iFace,DScd1,iP,iN,Phi,dPhidX)
       call DiffSchemeCounter(1,DScd1,1,0)

     case( DScd2 )

       PhiFace = SetFaceCDS(iFace,DScd2,iP,iN,Phi,dPhidX)
       call DiffSchemeCounter(1,DScd2,1,0)

     case( DScd3 )
                                                       !
       PhiFace  = 0.5*( Phi(in)+Phi(ip) )              ! very very simple average
       call DiffSchemeCounter(1,DScd3,1,0)             !
                                                       
     case( DSlux )                                     
                                                       
       Xf      = Face(iFace)%x                    
       d       = Xf - Cell(iC)%x                  
       PhiFace = PhiC + dot_product( dPhiC , d )  
       call DiffSchemeCounter(1,DSlux,1,0)        

     case default

       write(IOdef,*)'*** Error: Unknown differencing scheme'

   end select
   
   if( Debug > 3 ) write(IOdef,*)'=== SelectDiffSchemeScalar'

end subroutine SelectDiffSchemeScalar
real function SetFaceCDS(iFace,iScheme,iP,iN,Phi,dPhidX)

   use geometry
   use variables
   use constants

   integer, intent(IN)                      :: iFace, iScheme, iP, iN   
   real, dimension(Ncel+Nbnd), intent(IN)   :: Phi
   real, dimension(Ncel+Nbnd,3), intent(IN) :: dPhidX

   real, dimension(3)  :: Xac, Xp, Xn, Xf, dPhidXac, d1, d2
      
   if( Debug > 3 ) write(IOdef,*)'*** SetFacesCDS'

   Xp = Cell(ip)%x
   Xn = Cell(in)%x

   PhiP  = Phi(ip)
   PhiN  = Phi(in)

   if( iScheme == DScd1 )then
   
     facn  = Face(iFace)%lambda
     facp  = 1.0 - facn

     Xf    = Face(iFace)%x

     Xac   =      Xn * facn +      Xp * facp
     Phiac = Phi(in) * facn + Phi(ip) * facp 
     
     dPhidXac =  dPhidX(in,:) * facn + dPhidX(ip,:) * facp
     
     delta   = dot_product( dPhidXac , Xf - Xac )    ! non-orth. corr.
     PhiFace = Phiac + delta                         ! according to Ferziger

   else if( iScheme == DScd2 )then
     
     Xf       = Face(iFace)%x                        !
                                                     !
     d1       = Xf - Xp                              ! distance to cell face
     del1     = dot_product( dPhidX(ip,:) , d1 )     ! PhiP + grad_p*d1
                                                     !
     d2       = Xf - Xn                              ! distance to cell face
     del2     = dot_product( dPhidX(in,:) , d2 )     ! PhiN + grad_n*d2
                                                     !
     PhiFace  = 0.5*( PhiN + PhiP + del1 + del2 )    ! averaging

   elseif( iScheme == DScd3 )then
                                                     !
     PhiFace  = 0.5*( PhiN + PhiP )                  ! very very simple average
                                                     !
   else
       
     write(IOdef,*)'*** Error: In SetFaceCDS'
       
   endif
   
   SetFaceCDS = PhiFace

   if( Debug > 3 ) write(IOdef,*)'=== SetFacesCDS',PhiFace
   
end function SetFaceCDS
subroutine SelectDiffSchemeVector2(iFace,iScheme,iP,iN, &
                                  U,V,W,               &
				  dUdX,dVdX,dWdX,      &
				  UFace,VFace,Wface)
!========================================================================
!
!  selection of differencing schemes pulled out from
!  FluxScalar and FluxUVW. Two forms: a scalar version and
!  a vector version (this) (allowing to check the vector magnitude)
!
   use constants
   use geometry
   use variables, only : MassFlux

   real, intent(OUT)   :: UFace,VFace,WFace

   integer, intent(IN)                      :: iFace, iScheme, iP, iN   
   real, dimension(Ncel+Nbnd), intent(IN)   :: U,V,W
   real, dimension(Ncel+Nbnd,3), intent(IN) :: dUdX,dVdX,dWdX

   real, dimension(3)   :: Xac, Xp, Xn, Xf,          &
                           dUVW, ddUVW, d, d1, d2, ds

   real, dimension(3,3) :: TensorC 

   if( Debug > 3 ) write(IOdef,*)'*** SelectDiffSchemeVector'

   Xp = Cell(ip)%x
   Xn = Cell(in)%x

   !
   ! set the upwind and downwind values
   !
   if( MassFlux(iFace) >= 0.0 )then
     iC    = iP

     UC    = U(ip)
     VC    = V(ip)
     WC    = W(ip)

     UD    = U(in)
     VD    = V(in)
     WD    = W(in)

     TensorC(1:3,1) = dUdX(ip,:)
     TensorC(1:3,2) = dVdX(ip,:)
     TensorC(1:3,3) = dWdX(ip,:)

     dUVW(1) = UD - UC
     dUVW(2) = VD - VC
     dUVW(3) = WD - WC
     
     gradUVW = dot_product(dUVW,dUVW)
     
     d     = Xn - Xp

     ddUVW(1) = 2.0 * dot_product( TensorC(1:3,1) , d )
     ddUVW(2) = 2.0 * dot_product( TensorC(1:3,2) , d )
     ddUVW(3) = 2.0 * dot_product( TensorC(1:3,3) , d )

     gradCf   = dot_product(dUVW,ddUVW)
     
     UVWCt    = 1.0 - gradUVW/(gradCf+Small)

   else
     iC    = iN

     UC    = U(in)
     VC    = V(in)
     WC    = W(in)

     UD    = U(ip)
     VD    = V(ip)
     WD    = W(ip)

     TensorC(1:3,1) = dUdX(in,:)
     TensorC(1:3,2) = dVdX(in,:)
     TensorC(1:3,3) = dWdX(in,:)

     dUVW(1) = UD - UC
     dUVW(2) = VD - VC
     dUVW(3) = WD - WC
     
     gradUVW = dot_product(dUVW,dUVW)
     
     d     = Xn - Xp

     ddUVW(1) = 2.0 * dot_product( TensorC(1:3,1) , d )
     ddUVW(2) = 2.0 * dot_product( TensorC(1:3,2) , d )
     ddUVW(3) = 2.0 * dot_product( TensorC(1:3,3) , d )

     gradCf   = dot_product(dUVW,ddUVW)
     
     UVWCt    = 1.0 - gradUVW/(gradCf+Small)

   endif
   
   Beta_m   = 0.1                                    
   !
   ! most popular/used schemes are: 
   !
   !   UDS, CDS(1), Gamma, LUD, Minmod, CD2, CD3
   !
   ! so arrange them in the if-then-else tree accordingly
   !
   select case(iScheme)
     case( DSud ) 

       UFace = UC                                      !  
       VFace = VC                                      ! standard upwind
       WFace = WC                                      ! 
       call DiffSchemeCounter(1,iScheme,3,0)

     case( DScds ) 

       UFace = SetFaceCDS(iFace,DScdDefault,iP,iN,U,dUdX)
       VFace = SetFaceCDS(iFace,DScdDefault,iP,iN,V,dVdX)
       WFace = SetFaceCDS(iFace,DScdDefault,iP,iN,W,dWdX)
       call DiffSchemeCounter(1,DScds,3,0)

     case( DSgam )

       UCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,U,dUdX)   !
       VCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,V,dVdX)   ! later as 1 call
       WCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,W,dWdX)   !

       Beta_m = 0.1                                    ! Jasak's constant
       if( UVWCt <= 0.0 .or. UVWCt >= 1.0 )then        ! 0.1 <= Beta_m <= 0.5
         UFace = UC                                    !
         VFace = VC                                    ! select UDS
         WFace = WC                                    !
         call DiffSchemeCounter(1,DSud,3,0)             
       elseif( Beta_m <= UVWCt .and. UVWCt < 1.0  )then     
         UFace = UCDS                                  !
         VFace = VCDS                                  ! select CDS
         WFace = WCDS                                  !
         call DiffSchemeCounter(1,DScd1,3,0)            
       elseif(  0.0 < UVWCt .and. UVWCt < Beta_m )then 
         fg = UVWCt / Beta_m                           
         UFace = fg*UCDS + (1.0-fg)*UC                 !
         VFace = fg*VCDS + (1.0-fg)*VC                 ! select gamma
         WFace = fg*WCDS + (1.0-fg)*WC                 !
         call DiffSchemeCounter(1,DSgam,3,0)            
       endif

     case( DSlud )                                     ! LUDS Scheme:
                                                          
       if( UVWCt <= 0.0 .or. UVWCt >= 1.0  )then    
         UFace = UC                                    !
         VFace = VC                                    ! select UDS
         WFace = WC                                    !
         call DiffSchemeCounter(1,DSud,3,0)             
       else                                            
         Xf    = Face(iFace)%x                           
         d     = Xf - Cell(iC)%x                       
         UFace = UC + dot_product(TensorC(1:3,1),d)         
         VFace = VC + dot_product(TensorC(1:3,2),d)         
         WFace = WC + dot_product(TensorC(1:3,3),d)         
         call DiffSchemeCounter(1,DSlud,3,0)           
       endif

     case( DSmod )                                     ! MinMod Scheme:

       UCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,U,dUdX)
       VCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,V,dVdX)
       WCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,W,dWdX)

       if( UVWCt <= 0.0 .or. UVWCt >= 1.0  )then   
         UFace = UC                                    ! 
         VFace = VC                                    ! select UDS
         WFace = WC                                    !
         call DiffSchemeCounter(1,DSud,3,0)             
       elseif( UVWCt < 0.5 )then        
         Xf    = Face(iFace)%x                            
         d     = Xf - Cell(iC)%x                       
         UFace = UC + dot_product(TensorC(1:3,1),d)    !
         VFace = VC + dot_product(TensorC(1:3,2),d)    ! select LUDS
         WFace = WC + dot_product(TensorC(1:3,3),d)    !
         call DiffSchemeCounter(1,DSlud,3,0)           
       else                                            
         UFace = UCDS                                  !
         VFace = VCDS                                  ! select CDS
         WFace = WCDS                                  !
         call DiffSchemeCounter(1,DScd1,3,0)           
       endif                                           

     case( DScd1 )

       UFace = SetFaceCDS(iFace,DScd1,iP,iN,U,dUdX)
       VFace = SetFaceCDS(iFace,DScd1,iP,iN,V,dVdX)
       WFace = SetFaceCDS(iFace,DScd1,iP,iN,W,dWdX)
       call DiffSchemeCounter(1,DScd1,3,0)

     case( DScd2 )

       UFace = SetFaceCDS(iFace,DScd2,iP,iN,U,dUdX)
       VFace = SetFaceCDS(iFace,DScd2,iP,iN,V,dVdX)
       WFace = SetFaceCDS(iFace,DScd2,iP,iN,W,dWdX)
       call DiffSchemeCounter(1,DScd2,3,0)

     case( DScd3 )

       UFace = 0.5*( U(in)+U(ip) )                     
       VFace = 0.5*( V(in)+V(ip) )                     
       WFace = 0.5*( W(in)+W(ip) ) 
       call DiffSchemeCounter(1,DScd3,3,0)             

     case( DSlux )                                     ! LUDS Scheme:
                                                          
       Xf    = Face(iFace)%x                           
       d     = Xf - Cell(iC)%x                       
       UFace = UC + dot_product(TensorC(1:3,1),d)         
       VFace = VC + dot_product(TensorC(1:3,2),d)         
       WFace = WC + dot_product(TensorC(1:3,3),d)         
       call DiffSchemeCounter(1,DSlux,3,0)           

     case default

       write(IOdef,*)'*** Error: Unknown differencing scheme'

   end select

   if( Debug > 3 ) write(IOdef,*)'=== SelectDiffSchemeVector',UFace,VFace,WFace

end subroutine SelectDiffSchemeVector2
subroutine SelectDiffSchemeVector(iFace,iScheme,iP,iN, &
                                  U,V,W,               &
				  dUdX,dVdX,dWdX,      &
				  UFace,VFace,Wface)
!========================================================================
!
!  selection of differencing schemes pulled out from
!  FluxScalar and FluxUVW. Two forms: a scalar version and
!  a vector version (this) (allowing to check the vector magnitude)
!
   use constants
   use geometry
   use variables, only : MassFlux

   real, intent(OUT)   :: UFace,VFace,WFace

   integer, intent(IN)                      :: iFace, iScheme, iP, iN   
   real, dimension(Ncel+Nbnd), intent(IN)   :: U,V,W
   real, dimension(Ncel+Nbnd,3), intent(IN) :: dUdX,dVdX,dWdX

   real, dimension(3)  :: Xac, Xp, Xn, Xf,          &
			  dUdXac,dVdXac,dWdXac,     &
			  dUdXC,dVdXC,dWdXC,        &
			  dUdXD,dVdXD,dWdXD,        &
                          d, d1, d2, ds

   if( Debug > 3 ) write(IOdef,*)'*** SelectDiffSchemeVector'

   Xp = Cell(ip)%x
   Xn = Cell(in)%x

   UPs  = U(ip)
   VPs  = V(ip)
   WPs  = W(ip)
    
   UNs  = U(in)
   VNs  = V(in)
   WNs  = W(in)
   !
   ! set the upwind and downwind values
   !
   if( MassFlux(iFace) >= 0.0 )then
     iC    = iP

     UC    = UPs
     VC    = VPs
     WC    = WPs

     UD    = UNs
     VD    = VNs
     WD    = WNs

     dUdXC = dUdX(ip,:)
     dVdXC = dVdX(ip,:)
     dWdXC = dWdX(ip,:)

     dUdXD = dUdX(in,:)
     dVdXD = dVdX(in,:)
     dWdXD = dWdX(in,:)

     dU    = UD - UC
     dV    = VD - VC
     dW    = WD - WC
     
     d     = Xn - Xp
     ddU   = 2.0 * dot_product( dUdXC , d )
     UCt   = 1.0 - dU/(ddU+Small)
     ddV   = 2.0 * dot_product( dVdXC , d )
     VCt   = 1.0 - dV/(ddV+Small)
     ddW   = 2.0 * dot_product( dWdXC , d )
     WCt   = 1.0 - dW/(ddW+Small)

   else
     iC    = iN

     UC    = UNs
     VC    = VNs
     WC    = WNs

     UD    = UPs
     VD    = VPs
     WD    = WPs

     dUdXC = dUdX(in,:)
     dVdXC = dVdX(in,:)
     dWdXC = dWdX(in,:)

     dUdXD = dUdX(ip,:)
     dVdXD = dVdX(ip,:)
     dWdXD = dWdX(ip,:)

     dU    = UD - UC
     dV    = VD - VC
     dW    = WD - WC
     
     d     = Xn - Xp
     ddU   = 2.0 * dot_product( dUdXC , d )
     UCt   = 1.0 - dU/(ddU+Small)
     ddV   = 2.0 * dot_product( dVdXC , d )
     VCt   = 1.0 - dV/(ddV+Small)
     ddW   = 2.0 * dot_product( dWdXC , d )
     WCt   = 1.0 - dW/(ddW+Small)

   endif
   
   UVWCmag2 = UC*UC + VC*VC + WC*WC
   UVWDmag2 = UD*UD + VD*VD + WD*WD
   
   !
   ! most popular/used schemes are: 
   !
   !   UDS, CDS(1), Gamma, LUD, Minmod, CD2, CD3
   !
   ! so arrange them in the if-then-else tree accordingly
   !
   select case(iScheme)
     case( DSud ) 

       UFace = UC                                      !  
       VFace = VC                                      ! standard upwind
       WFace = WC                                      ! 
       call DiffSchemeCounter(1,iScheme,3,0)

     case( DScds ) 

       UFace = SetFaceCDS(iFace,DScdDefault,iP,iN,U,dUdX)
       VFace = SetFaceCDS(iFace,DScdDefault,iP,iN,V,dVdX)
       WFace = SetFaceCDS(iFace,DScdDefault,iP,iN,W,dWdX)

       call DiffSchemeCounter(1,DScds,3,0)

     case( DSgam )

       UCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,U,dUdX)   !
       VCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,V,dVdX)   ! later as 1 call
       WCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,W,dWdX)   !

       Beta_m = 0.1                                    ! Jasak's constant
       if( UCt <= 0.0 .or. UCt >= 1.0 )then            ! 0.1 <= Beta_m <= 0.5
         UFace = UC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       elseif( Beta_m <= UCt .and. UCt < 1.0 )then     !
         UFace = UCDS                                  ! select CDS
         call DiffSchemeCounter(1,DScd1,1,0)           !
       elseif( 0.0 < UCt .and. UCt < Beta_m )then      !
         fg = UCt / Beta_m                             !
         UFace = fg*UCDS + (1.0-fg)*UC                 ! select gamma
         call DiffSchemeCounter(1,DSgam,1,0)           !
       endif

       if( VCt <= 0.0 .or. VCt >= 1.0 )then            ! 0.1 <= Beta_m <= 0.5
         VFace = VC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       elseif( Beta_m <= VCt .and. VCt < 1.0 )then     !
         VFace = VCDS                                  ! select CDS
         call DiffSchemeCounter(1,DScd1,1,0)           !
       elseif( 0.0 < VCt .and. VCt < Beta_m )then      !
         fg = VCt / Beta_m                             !
         VFace = fg*VCDS + (1.0-fg)*VC                 ! select gamma
         call DiffSchemeCounter(1,DSgam,1,0)           !
       endif

       if( WCt <= 0.0 .or. WCt >= 1.0 )then            ! 0.1 <= Beta_m <= 0.5
         WFace = WC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       elseif( Beta_m <= WCt .and. WCt < 1.0 )then     !
         WFace = WCDS                                  ! select CDS
         call DiffSchemeCounter(1,DScd1,1,0)           !
       elseif( 0.0 < WCt .and. WCt < Beta_m )then      !
         fg = WCt / Beta_m                             !
         WFace = fg*WCDS + (1.0-fg)*WC                 ! select gamma
         call DiffSchemeCounter(1,DSgam,1,0)           !
       endif

     case( DSlud )                                     ! LUDS Scheme:
                                                       !  
       if( UCt <= 0.0 .or. UCt >= 1.0 )then            !
         UFace = UC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       else                                            !
         Xf    = Face(iFace)%x                         !   
         d     = Xf - Cell(iC)%x                       !
         UFace = UC + dot_product( dUdXC , d )         !
         call DiffSchemeCounter(1,DSlud,1,0)           !
       endif

       if( VCt <= 0.0 .or. VCt >= 1.0 )then            !
         VFace = VC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       else                                            !
         Xf    = Face(iFace)%x                         !   
         d     = Xf - Cell(iC)%x                       !
         VFace = VC + dot_product( dVdXC , d )         !
         call DiffSchemeCounter(1,DSlud,1,0)           !
       endif

       if( WCt <= 0.0 .or. WCt >= 1.0 )then            !
         WFace = WC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       else                                            !
         Xf    = Face(iFace)%x                         !   
         d     = Xf - Cell(iC)%x                       !
         WFace = WC + dot_product( dWdXC , d )         !
         call DiffSchemeCounter(1,DSlud,1,0)           !
       endif

     case( DSmod )                                     ! MinMod Scheme:

       UCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,U,dUdX)
       VCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,V,dVdX)
       WCDS = SetFaceCDS(iFace,DScdDefault,iP,iN,W,dWdX)

       if( UCt <= 0.0 .or. UCt >= 1.0 )then            !
         UFace = UC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       elseif( Uct < 0.5 )then                         !
         Xf    = Face(iFace)%x                         !   
         d     = Xf - Cell(iC)%x                       !
         UFace = UC + dot_product( dUdXC , d )         ! select LUDS
         call DiffSchemeCounter(1,DSlud,1,0)           !
       else                                            !
         UFace = UCDS                                  ! select CDS
         call DiffSchemeCounter(1,DScd1,1,0)           !
       endif                                           !

       if( VCt <= 0.0 .or. VCt >= 1.0 )then            !
         VFace = VC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       elseif( VCt < 0.5 )then                         !
         Xf    = Face(iFace)%x                         !   
         d     = Xf - Cell(iC)%x                       !
         VFace = VC + dot_product( dVdXC , d )         ! select LUDS
         call DiffSchemeCounter(1,DSlud,1,0)           !
       else                                            !
         VFace = VCDS                                  ! select CDS
         call DiffSchemeCounter(1,DScd1,1,0)           !
       endif                                           !

       if( WCt <= 0.0 .or. WCt >= 1.0 )then            !
         WFace = WC                                    ! select UDS
         call DiffSchemeCounter(1,DSud,1,0)            !
       elseif( WCt < 0.5 )then                         !
         Xf    = Face(iFace)%x                         !   
         d     = Xf - Cell(iC)%x                       !
         WFace = WC + dot_product( dWdXC , d )         ! select LUDS
         call DiffSchemeCounter(1,DSlud,1,0)           !
       else                                            !
         WFace = WCDS                                  ! select CDS
         call DiffSchemeCounter(1,DScd1,1,0)           !
       endif                                           !

     case( DScd1 )

       UFace = SetFaceCDS(iFace,DScd1,iP,iN,U,dUdX)
       VFace = SetFaceCDS(iFace,DScd1,iP,iN,V,dVdX)
       WFace = SetFaceCDS(iFace,DScd1,iP,iN,W,dWdX)
       call DiffSchemeCounter(1,DScd1,3,0)

     case( DScd2 )

       UFace = SetFaceCDS(iFace,DScd2,iP,iN,U,dUdX)
       VFace = SetFaceCDS(iFace,DScd2,iP,iN,V,dVdX)
       WFace = SetFaceCDS(iFace,DScd2,iP,iN,W,dWdX)
       call DiffSchemeCounter(1,DScd2,3,0)

     case( DScd3 )

       UFace = 0.5*( U(in)+U(ip) )                     
       VFace = 0.5*( V(in)+V(ip) )                     
       WFace = 0.5*( W(in)+W(ip) )                     
       call DiffSchemeCounter(1,DScd3,3,0)             

     case( DSlux )                                     
                                                    
       Xf    = Face(iFace)%x                        
       d     = Xf - Cell(iC)%x                      
       UFace = UC + dot_product( dUdXC , d )        
       VFace = VC + dot_product( dVdXC , d )        
       WFace = WC + dot_product( dWdXC , d )        
       call DiffSchemeCounter(1,DSlux,3,0)          

     case default

       write(IOdef,*)'*** Error: Unknown differencing scheme'

   end select

   if( Debug > 3 ) write(IOdef,*)'=== SelectDiffSchemeVector',UFace,VFace,WFace

end subroutine SelectDiffSchemeVector
subroutine DiffSchemeCounter(mode,iScheme,incr,iVar)
!========================================================================
!
!  counters in order to judge the relative number of scheme use
!
   use constants
   use variables
   
   integer, dimension(9), save :: icntrs
   real, dimension(9)          :: rcntrs

   if( Debug > 3 ) write(IOdef,*)'*** DiffSchemeCounter',mode,iScheme,iVar

   if( mode == 0 )then
    
     icntrs = 0
   
   elseif( mode == 1 )then
   
     icntrs(iScheme) = icntrs(iScheme) + 1
   
   elseif( mode == 2 )then
   
     isum = sum(icntrs)
     if( isum > 0 )then
       r = 1.0/float( sum(icntrs) )
     else
       r = 1.0
     endif
     
     rcntrs = 0.0
     rcntrs(DSud ) = float(icntrs(DSud ))*r
     rcntrs(DScds) = float(icntrs(DScds))*r
     rcntrs(DSgam) = float(icntrs(DSgam))*r
     rcntrs(DSlud) = float(icntrs(DSlud))*r
     rcntrs(DSmod) = float(icntrs(DSmod))*r
     rcntrs(DScd1) = float(icntrs(DScd1))*r
     rcntrs(DScd2) = float(icntrs(DScd2))*r
     rcntrs(DScd3) = float(icntrs(DScd3))*r
     rcntrs(DSlux) = float(icntrs(DSlux))*r
   
     if( iVar > 3 )then
       write(IOdbg,*)'Diff. scheme usage, variable ',Variable(iVar)
     else
       write(IOdbg,*)'Diff. scheme usage, velocity components'
     endif
     do i=1,9
       if( icntrs(i) > 0 ) write(IOdbg,'(1x, A,1x,i8,1x,f8.2,''%'')') &
                           DSnames(i),icntrs(i),rcntrs(i)*100.0
     end do
   else
   
     write(IOdef,*)'Foutje... counters...'
     
   endif

   if( Debug > 3 ) write(IOdef,*)'=== DiffSchemeCounter'

end subroutine DiffSchemeCounter
