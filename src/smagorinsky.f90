!
! Copyright 2009 Max Staufer
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
subroutine CalculateSmagorinskyViscosity

   use constants
   use geometry
   use variables
 
   implicit none

   real :: betaij(3,3), alpha

   real :: Cs = 0.1
   real :: C
   real :: FilterLength, mutsgs

   integer :: ip


   C = 2.5*Cs**2

   do ip=1,Ncel

     FilterLength = (Cell(ip)%vol)**(1./3.)
     FL2 = FilterLength*FilterLength
     
     betaij(1,1) = FL2*( dudx(ip,1)*dudx(ip,1) + dudx(ip,2)*dudx(ip,2) + dudx(ip,3)*dudx(ip,3) )
     betaij(1,2) = FL2*( dudx(ip,1)*dvdx(ip,1) + dudx(ip,2)*dvdx(ip,2) + dudx(ip,3)*dvdx(ip,3) )
     betaij(1,3) = FL2*( dudx(ip,1)*dwdx(ip,1) + dudx(ip,2)*dwdx(ip,2) + dudx(ip,3)*dwdx(ip,3) )

     betaij(2,1) = FL2*( dvdx(ip,1)*dudx(ip,1) + dvdx(ip,2)*dudx(ip,2) + dvdx(ip,3)*dudx(ip,3) )
     betaij(2,2) = FL2*( dvdx(ip,1)*dvdx(ip,1) + dvdx(ip,2)*dvdx(ip,2) + dvdx(ip,3)*dvdx(ip,3) )
     betaij(2,3) = FL2*( dvdx(ip,1)*dwdx(ip,1) + dvdx(ip,2)*dwdx(ip,2) + dvdx(ip,3)*dwdx(ip,3) )

     betaij(3,1) = FL2*( dwdx(ip,1)*dudx(ip,1) + dwdx(ip,2)*dudx(ip,2) + dwdx(ip,3)*dudx(ip,3) )
     betaij(3,2) = FL2*( dwdx(ip,1)*dvdx(ip,1) + dwdx(ip,2)*dvdx(ip,2) + dwdx(ip,3)*dvdx(ip,3) )
     betaij(3,3) = FL2*( dwdx(ip,1)*dwdx(ip,1) + dwdx(ip,2)*dwdx(ip,2) + dwdx(ip,3)*dwdx(ip,3) )

     alpha = dudx(ip,1)*dudx(ip,1) + dudx(ip,2)*dudx(ip,2) + dudx(ip,3)*dudx(ip,3) + &
             dvdx(ip,1)*dvdx(ip,1) + dvdx(ip,2)*dvdx(ip,2) + dvdx(ip,3)*dvdx(ip,3) + &
             dwdx(ip,1)*dwdx(ip,1) + dwdx(ip,2)*dwdx(ip,2) + dwdx(ip,3)*dwdx(ip,3)


     mutsgs = abs(&
                  betaij(1,1)*betaij(2,2) - betaij(1,2)*betaij(1,2) + &
                  betaij(1,1)*betaij(3,3) - betaij(1,3)*betaij(1,3) + &
                  betaij(2,2)*betaij(3,3) - betaij(2,3)*betaij(2,3)     )/alpha


     mutsgs =den(ip)*C*sqrt(mutsgs)

     if (alpha == 0.0 ) mutsgs = 0.0

     Viseff(ip) = Vislam + mutsgs

   end do

end subroutine CalculateSmagorinskyViscosity
