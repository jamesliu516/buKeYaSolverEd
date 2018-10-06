!
! Copyright 2009 Henk Krus, Cyclone Fluid Dynamics BV
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
! This is gmsh.f90, the dolfyn interface to gmsh 
!
! Gmsh is a preprocessor by Christophe Geuzaine and Jean-François Remacle. 
!
! See http://geuz.org/gmsh
!
subroutine dolfyn2gmsh

   use constants
   use geometry
   use variables
   use particles
   
   character(len=3)     :: string
   
   integer :: indx(8)
   
   real, allocatable    :: NodalData1(:)    ! array of data at vertices
   real, allocatable    :: NodalData2(:)    
   real, allocatable    :: NodalData3(:)    
   
   integer, allocatable :: NodalCounter(:)  

   if( .not. UseGMSHdataonly )then
     write(IOdef,*)  'Writing Full Gmsh file...'
     write(IOdbg,*)  'Writing Full Gmsh file...'
   else
 
     NCellOffset = count( Bnd(:)%rid > 0 )

     write(IOdef,*)  'Writing Gmsh data file...',NCellOffset
     write(IOdbg,*)  'Writing Gmsh data file...',NCellOffset  
   endif

   allocate( NodalData1(Nvrt),stat=istat)    
   allocate( NodalCounter(Nvrt),stat=istat)    
   call TrackMemory(istat,Ncel+Nbnd+2*Nvrt,'Work arrays for NodalData allocated')
      
   if( .not. UseGMSHdataonly )then
     !
     ! write everything into one single file
     !
     call openfile(IOpst,casename,'.msh','FORMATTED', &
                                         'SEQUENTIAL','UNKNOWN',debug)
   
     !
     ! note: in the binary format 'one-binary' is needed as well.
     !
     write(IOpst,'(A)')      '$MeshFormat'
     write(IOpst,'(A,A,A)')  '2.0 0 8' 
     write(IOpst,'(A)')      '$EndMeshFormat' 
     
     Nctid = maxval(cell(:)%ctid)
     NPhysicalNames = Nctid + (Nreg+1)
   
     write(IOpst,'(A)')      '$PhysicalNames' 
     write(IOpst,*) NPhysicalNames
     do i=1,Nctid
       write(string,'(i3.3)')   i
       write(IOpst,'(i4,1x,A)') i,'"Fluid'//string//'"'
     end do   
     write(IOpst,'(i4,1x,A)') Nctid+1,'"DefaultWall"'
     do i=1,Nreg
       write(IOpst,'(i4,1x,A)') &
         Nctid+1+i,'"'//Reg(i)%name(1:lens(Reg(i)%name))//'"' 
     end do
     write(IOpst,'(A)')      '$EndPhysicalNames' 

     !
     ! vertices (also called nodes)
     !
     write(IOpst,'(A)')  '$Nodes'
     write(IOpst,'(i8)') Nvrt 
     do i=1,Nvrt
       write(IOpst,'(i8,3(1x,1pe12.5))') i,Vert(i,:)
     end do
     write(IOpst,'(A)')  '$EndNodes'

     !
     ! the definitions are in the original *.cel file
     !
     call openfile(IOcel,casename,'.cel','FORMATTED', &
                                	 'SEQUENTIAL','OLD',debug)
     !
     ! gmsh elm-types 
     !
     !   2  3-node triangle
     !   3  4-node quadrangle
     !   4  4-node tet
     !   5  8-node hex
     !   6  6-node prism/wedge
     !   7  5-node pyramid
     !
     ! the types in *.cel
     !
     !   hexaeder:  1,2,3,4 5,6,7,8  
     !   pentaeder: 1,2,3,3 5,6,7,7  
     !   pyramid:   1,2,3,4 5,5,5,5  
     !   tetraeder: 1,2,3,3 4,4,4,4  
     !   quad:      1,2,3,4          
     !   triangle:  1,2,3,3          
     !
     ! 9.1 msh ascii file format:
     !
     !   $Elements
     !   number-of-elements
     !   elm-number elm-type number-of-tags < tag > ... node-number-list
     !   ....
     !   $EndElements
     !
     ! number-of-tags: gives the number of integer tags that follow for the
     ! n-th element. By default, the FIRST tag is the number of the physical
     ! entity to which the element belongs; the SECOND is the number of the
     ! elementary geometrical entity to which the element belongs; the THIRD
     ! is the number of a mesh partition to which the element belongs. All
     ! tags must be postive integers, or zero. A zero tag is equivalent to
     ! no tag.
     !
     !   Boundary region names: $PhysicalName
     !
     !   $PhysicalNames
     !   number-of-names
     !   physical-number "physical-name"
     !   ...
     !   $EndPhysicalNames
     !

     NHex  = 0
     NPen  = 0
     NPyr  = 0
     NTet  = 0
     NQuad = 0
     NTri  = 0

     write(IOpst,'(A)')  '$Elements'

     if( UseGMSHwalls )then
       write(IOpst,'(i8)') Ncel+Nbnd
     else
       write(IOpst,'(i8)') Ncel 
     endif

     do i=1,Ncel
       read(IOcel,*) idummy,(indx(j),j=1,8)
       i1 = indx(1)
       i2 = indx(2)
       i3 = indx(3)
       i4 = indx(4)
       i5 = indx(5)
       i6 = indx(6)
       i7 = indx(7)
       i8 = indx(8)

       if( i7 /= 0 )then
	 if( i7 /= i8 )then
           Nhex   = Nhex + 1
           write(IOpst,'(i8,2(1x,i1),i4,2(1x,i1),8(i8,1x))') &
	                             i,5,3,Cell(i)%ctid,0,1, &
                                     indx(1),indx(2),indx(3),indx(4), &                
                                     indx(5),indx(6),indx(7),indx(8)                  
	 else if( i5 /= i6 )then
           NPen = NPen + 1
	 else if( i3 /= i4 )then
           NPyr  = NPyr + 1
	 else if( i3 == i4 .and. i5 == i6 )then
           NTet   = NTet + 1
           !
           ! id type 1 id (rest) 1,2,3,3 4
           !
           write(IOpst,'(i8,2(1x,i1),i4,2(1x,i1),4(i8,1x))') &
	                             i,4,3,Cell(i)%ctid,0,1, &
                                     indx(1),indx(2),indx(3),indx(5)                  
	 else
           write(*,*)'Unknown gmsh shape, cell ',i,' : ',i1,i2,i3,i4,' ...'
	 endif
       endif
     end do

     if( UseGMSHwalls )then
       !
       ! boundaries
       !
       do ib=1,Nbnd
	 ir = bnd(ib)%rid 
	 indx(1:4) = bnd(ib)%vertices
	 if( indx(4) == -1 )then
	   Ntri  = Ntri + 1
	   write(IOpst,'(i8,2(1x,i1),i4,2(1x,i1),3(i8,1x))') &
                                     Ncel+ib,2,3,Nctid+1+ir,0,1, &
        			     indx(1),indx(2),indx(3)
	   icnt = icnt + 1
	 else
	   Nquad = Nquad + 1
	   write(IOpst,'(i8,2(1x,i1),i4,2(1x,i1),4(i8,1x))') &
        			     Ncel+ib,3,3,Nctid+1+ir,0,1, &
        			     indx(1),indx(2),indx(3),indx(4)
	   icnt = icnt + 1
	 endif
       end do
     endif
     
     write(IOpst,'(A)')  '$EndElements'

     if( NHex  > 0 )write(IOdbg,*)'Number of hexa''s     :',NHex 
     if( NPen  > 0 )write(IOdbg,*)'Number of prism''s    :',NPen 
     if( NPyr  > 0 )write(IOdbg,*)'Number of pyramids''s :',NPyr 
     if( NTet  > 0 )write(IOdbg,*)'Number of tetra''s    :',NTet 
     if( NQuad > 0 )write(IOdbg,*)'Number of quad''s     :',NQuad
     if( NTri  > 0 )write(IOdbg,*)'Number of triangles''s:',NTri 

     close(IOcel)
   
   endif
   
   if( PostC(VarU) .and. .not. UseGMSHdataonly )then

     write(IOdbg,*)'Writing Gmsh vectors cell'
     write(IOdef,*)'Writing Gmsh vectors cell'

     write(IOpst,'(A)')     '$ElementData'
     write(IOpst,'(A)')     '1'
     write(IOpst,'(A)')     '"Velocity vector cell"'
     write(IOpst,'(A)')     '1'
     write(IOpst,'(e10.3)') Time
     write(IOpst,'(A)')     '3'
     write(IOpst,'(i10)')    0     !  Iter
     write(IOpst,'(A)')     '3'
     if( UseGMSHwalls )then
       write(IOpst,'(i10)')   Ncel+Nbnd
     else
       write(IOpst,'(i10)')   Ncel  
     endif

     do i=1,Ncel
       write(IOpst,'(i8,3(1x,1pe12.5))') i, U(i),V(i),W(i) 
     end do
    
     if( UseGMSHwalls )then
       do i=1,Nbnd
	 write(IOpst,'(i8,3(1x,1pe12.5))') Ncel+i,U(Ncel+i),V(Ncel+i),W(Ncel+i)    
       end do
     endif 

     write(IOpst,'(A)')     '$EndElementData'

   elseif( PostC(VarU) .and. UseGMSHdataonly )then

     write(IOdbg,*)'Writing Gmsh data vectors cell'
     write(IOdef,*)'Writing Gmsh data vectors cell'

     call openfile(IOpst,casename,'_uvwc.msh','FORMATTED', &
                                              'SEQUENTIAL','UNKNOWN',debug)
   
     write(IOpst,'(A)')      '$MeshFormat'
     write(IOpst,'(A,A,A)')  '2.0 0 8' 
     write(IOpst,'(A)')      '$EndMeshFormat' 

     write(IOpst,'(A)')      '$ElementData'
     write(IOpst,'(A)')      '1'
     write(IOpst,'(A)')      '"Velocity vector cell"'
     write(IOpst,'(A)')      '1'
     write(IOpst,'(e10.3)')  Time
     write(IOpst,'(A)')      '3'
     write(IOpst,'(i10)')     0     !  Iter
     write(IOpst,'(A)')      '3'
     write(IOpst,'(i10)')    NCellOffset+Ncel

     icnt = 0
     do i=1,Nbnd
       if( bnd(i)%rid > 0 )then
         icnt = icnt + 1
         write(IOpst,'(i8,3(1x,1pe12.5))') icnt,U(Ncel+i),V(Ncel+i),W(Ncel+i)    
       endif
     end do

     do i=1,Ncel
       write(IOpst,'(i8,3(1x,1pe12.5))') NCellOffset+i, U(i),V(i),W(i) 
     end do

     write(IOpst,'(A)')     '$EndElementData'

     close(IOpst)

   endif

   !
   ! *** PRESSURE ***
   !
   if( PostC(VarP) .and. .not. UseGMSHdataonly .and. allocated(P) )then
     write(IOdbg,*)'Writing Gmsh pressure cell'
     write(IOdef,*)'Writing Gmsh pressure cell'
     if( UseGMSHbinary )then
     
     else
       write(IOpst,'(A)')     '$ElementData'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(A)')     '"Pressure cell"'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(e10.3)') Time
       write(IOpst,'(A)')     '3'
       write(IOpst,'(i10)')    0   !  Iter
       write(IOpst,'(A)')     '1'
       if( UseGMSHwalls )then
         write(IOpst,'(i10)')   Ncel+Nbnd
       else
         write(IOpst,'(i10)')   Ncel  
       endif

       do i=1,Ncel
         write(IOpst,'(i8,1x,1pe10.3)') i,P(i)
       end do
       if( UseGMSHwalls )then
         do i=1,Nbnd
           write(IOpst,'(i8,1x,1pe10.3)') Ncel+i,P(Ncel+i)
         end do
       endif 
     endif
     write(IOpst,'(A)')     '$EndElementData'
   
   elseif( PostC(VarP) .and. UseGMSHdataonly .and. allocated(P) )then   
     write(IOdbg,*)'Writing Gmsh data pressure cell'
     write(IOdef,*)'Writing Gmsh data pressure cell'

     if( UseGMSHbinary )then
     
     else
       call openfile(IOpst,casename,'_pc.msh','FORMATTED', &
                                              'SEQUENTIAL','UNKNOWN',debug)
       write(IOpst,'(A)')     '$MeshFormat'
       write(IOpst,'(A,A,A)') '2.0 0 8' 
       write(IOpst,'(A)')     '$EndMeshFormat' 

       write(IOpst,'(A)')     '$ElementData'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(A)')     '"Pressure cell"'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(e10.3)') Time
       write(IOpst,'(A)')     '3'
       write(IOpst,'(i10)')    0   !  Iter
       write(IOpst,'(A)')     '1'
       write(IOpst,'(i10)')    NCellOffset+Ncel

       icnt = 0
       do i=1,Nbnd
         if( bnd(i)%rid > 0 )then
	   icnt = icnt + 1
           write(IOpst,'(i8,1x,1pe10.3)') icnt,P(Ncel+i)
	 endif
       end do

       do i=1,Ncel
         write(IOpst,'(i8,1x,1pe10.3)') NCellOffset+i,P(i)
       end do

       write(IOpst,'(A)')     '$EndElementData'
       close(IOpst)
     endif  
   endif 

   if( PostV(VarP) .and. .not. UseGMSHdataonly .and. allocated(P) )then
     write(IOdbg,*)'Writing Gmsh pressure vertex cell'
     write(IOdef,*)'Writing Gmsh pressure vertex cell'

     call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)

     if( UseGMSHbinary )then
     
     else
       write(IOpst,'(A)')     '$NodeData'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(A)')     '"Pressure vertex"'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(e10.3)') Time
       write(IOpst,'(A)')     '3'
       write(IOpst,'(i10)')    0   !  Iter
       write(IOpst,'(A)')     '1'
       write(IOpst,'(i10)')   Nvrt

       do i=1,Nvrt
         write(IOpst,'(i8,1x,1pe10.3)') i,Nodaldata1(i)
       end do

       write(IOpst,'(A)')     '$EndNodeData'
     endif

   elseif( PostV(VarP) .and. UseGMSHdataonly .and. allocated(P) )then
     write(IOdbg,*)'Writing Gmsh data pressure vertex'
     write(IOdef,*)'Writing Gmsh data pressure vertex'

     call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)

     if( UseGMSHbinary )then
     
     else
       call openfile(IOpst,casename,'_pv.msh','FORMATTED', &
                                            'SEQUENTIAL','UNKNOWN',debug)
       write(IOpst,'(A)')     '$MeshFormat'
       write(IOpst,'(A,A,A)') '2.0 0 8' 
       write(IOpst,'(A)')     '$EndMeshFormat' 

       write(IOpst,'(A)')     '$NodeData'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(A)')     '"Pressure vertex"'
       write(IOpst,'(A)')     '1'
       write(IOpst,'(e10.3)') Time
       write(IOpst,'(A)')     '3'
       write(IOpst,'(i10)')    0   !  Iter
       write(IOpst,'(A)')     '1'
       write(IOpst,'(i10)')    Nvrt
 
       do i=1,Nvrt
         write(IOpst,'(i8,1x,1pe10.3)') i,Nodaldata1(i)
       end do
       write(IOpst,'(A)')     '$EndNodeData'
       close(IOpst)
     endif
   endif 
   
   !
   ! *** TURBULENCE ***
   !
   if( SolveTurb )then
     !
     ! *** TURBULENCE KINETIC ENERGY ***
     !
     if( PostC(VarTE) .and. .not. UseGMSHdataonly .and. allocated(TE) )then
       write(IOdbg,*)'Writing Gmsh turbulent kinetic energy cell'
       write(IOdef,*)'Writing Gmsh turbulent kinetic energy cell'
 
       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"TKE cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
         if( UseGMSHwalls )then
           write(IOpst,'(i10)')   Ncel+Nbnd
         else
           write(IOpst,'(i10)')   Ncel  
         endif

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,TE(i))
	 end do
	 if( UseGMSHwalls )then
	   do i=1,Nbnd
             write(IOpst,'(i8,1x,1pe10.3)') Ncel+i,max(0.0,TE(Ncel+i))
	   end do
	 endif 
       endif
       write(IOpst,'(A)')     '$EndElementData'
     
     elseif( PostC(VarTE) .and. UseGMSHdataonly .and. allocated(TE) )then
       write(IOdbg,*)'Writing Gmsh data turbulent kinetic energy cell'
       write(IOdef,*)'Writing Gmsh data turbulent kinetic energy cell'

       if( UseGMSHbinary )then

       else        
         call openfile(IOpst,casename,'_tkc.msh','FORMATTED', &
                                                 'SEQUENTIAL','UNKNOWN',debug) 
	 write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"TKE cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   NCellOffset+Ncel 

	 icnt = 0
	 do i=1,Nbnd
           if( bnd(i)%rid > 0 )then
	     icnt = icnt + 1
             write(IOpst,'(i8,1x,1pe10.3)') icnt,max(0.0,TE(Ncel+i))
	   endif
	 end do

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') NCellOffset+i,max(0.0,TE(i))
	 end do
         write(IOpst,'(A)')     '$EndElementData'
         close(IOpst) 
       endif
     endif

     if( PostV(VarTE) .and. .not. UseGMSHdataonly .and. allocated(TE) )then
       write(IOdbg,*)'Writing Gmsh turbulent kinetic energy vertex'
       write(IOdef,*)'Writing Gmsh turbulent kinetic energy vertex'
 
       call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"TKE vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,NodalData1(i))
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
       endif
     elseif( PostV(VarTE) .and. UseGMSHdataonly .and. allocated(TE) )then
       write(IOdbg,*)'Writing Gmsh data turbulent kinetic energy vertex'
       write(IOdef,*)'Writing Gmsh data turbulent kinetic energy vertex'
 
       call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_tkv.msh','FORMATTED', &
                                                 'SEQUENTIAL','UNKNOWN',debug)
         write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"TKE vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,NodalData1(i))
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
         close(IOpst) 
       endif
     endif

     !
     ! *** TURBULENCE DISSIPATION ***
     !
     if( PostC(VarED) .and. .not. UseGMSHdataonly .and. allocated(ED) )then
       write(IOdbg,*)'Writing Gmsh turbulent dissipation cell',Ncel+Nbnd,Ncel,Nbnd
       write(IOdef,*)'Writing Gmsh turbulent dissipation cell'

       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"ED cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
         if( UseGMSHwalls )then
           write(IOpst,'(i10)')   Ncel+Nbnd
         else
           write(IOpst,'(i10)')   Ncel  
         endif

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,ED(i))
	 end do
	 if( UseGMSHwalls )then
	   do i=1,Nbnd
             write(IOpst,'(i8,1x,1pe10.3)') Ncel+i,max(0.0,ED(Ncel+i))
	   end do
	 endif 
       endif
       write(IOpst,'(A)')     '$EndElementData'
     elseif( PostC(VarED) .and. UseGMSHdataonly .and. allocated(ED) )then
       write(IOdbg,*)'Writing Gmsh data turbulent dissipation cell'
       write(IOdef,*)'Writing Gmsh data turbulent dissipation cell'

       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_edc.msh','FORMATTED', &
                                                 'SEQUENTIAL','UNKNOWN',debug) 
         write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"ED cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   NCellOffset+Ncel 

	 icnt = 0
	 do i=1,Nbnd
           if( bnd(i)%rid > 0 )then
	     icnt = icnt + 1
             write(IOpst,'(i8,1x,1pe10.3)') icnt,max(0.0,ED(Ncel+i))
	   endif
	 end do

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,ED(i))
	 end do
         write(IOpst,'(A)')     '$EndElementData'
         close(IOpst)
       endif     
     endif

     if( PostV(VarED) .and. .not. UseGMSHdataonly .and. allocated(ED) )then
       write(IOdbg,*)'Writing Gmsh turbulent dissipation vertex'
       write(IOdef,*)'Writing Gmsh turbulent dissipation vertex'

       call InterpolateData(0,VarED,ED,dPdX,NodalData1,NodalCounter)

       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"ED vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,NodalData1(i))
	 end do
       endif
       write(IOpst,'(A)')     '$EndNodeData'
     elseif( PostV(VarED) .and. UseGMSHdataonly .and. allocated(ED) )then
       write(IOdbg,*)'Writing Gmsh data turbulent dissipation vertex'
       write(IOdef,*)'Writing Gmsh data turbulent dissipation vertex'

       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_edv.msh','FORMATTED', &
                                                 'SEQUENTIAL','UNKNOWN',debug) 
         write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"ED vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,NodalData1(i))
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
         close(IOpst)
       endif     
     endif

     !
     ! *** TURBULENT VISCOSITY ***
     !
     if( PostC(VarVis) .and. .not. UseGMSHdataonly .and. allocated(VisEff) )then
       write(IOdbg,*)'Writing Gmsh turbulent viscosity cell'
       write(IOdef,*)'Writing Gmsh turbulent viscosity cell'
       
       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"VisT cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
         if( UseGMSHwalls )then
           write(IOpst,'(i10)')   Ncel+Nbnd
         else
           write(IOpst,'(i10)')   Ncel  
         endif

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,VisEff(i) - VisLam)
	 end do
	 if( UseGMSHwalls )then
	   do i=1,Nbnd
             write(IOpst,'(i8,1x,1pe10.3)') Ncel+i,max(0.0,VisEff(Ncel+i) - VisLam)
	   end do
	 endif 
         write(IOpst,'(A)')     '$EndElementData'
       endif
     elseif( PostC(VarVis) .and. UseGMSHdataonly .and. allocated(VisEff) )then
       write(IOdbg,*)'Writing Gmsh data turbulent viscosity cell'
       write(IOdef,*)'Writing Gmsh data turbulent viscosity cell'
       
       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_vtc.msh','FORMATTED', &
                                                 'SEQUENTIAL','UNKNOWN',debug) 
	 write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"VisT cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   NCellOffset+Ncel 

	 icnt = 0
	 do i=1,Nbnd
           if( bnd(i)%rid > 0 )then
	     icnt = icnt + 1
             write(IOpst,'(i8,1x,1pe10.3)') icnt,max(0.0,VisEff(i) - VisLam)
	   endif
	 end do

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') NCellOffset+i,max(0.0,VisEff(i) - VisLam)
	 end do
         write(IOpst,'(A)')     '$EndElementData'
         close(IOpst) 
       endif
     endif

     if( PostV(VarVis) .and. .not. UseGMSHdataonly .and. allocated(VisEff) )then
       write(IOdbg,*)'Writing Gmsh turbulent viscosity vertex'
       write(IOdef,*)'Writing Gmsh turbulent viscosity vertex'
 
       call InterpolateData(0,VarVis,VisEff,dPdX,NodalData1,NodalCounter)
       
       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"VisT vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,NodalData1(i) - VisLam)
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
       endif
     elseif( PostV(VarVis) .and. UseGMSHdataonly .and. allocated(VisEff) )then
       write(IOdbg,*)'Writing Gmsh data turbulent viscosity vertex'
       write(IOdef,*)'Writing Gmsh data turbulent viscosity vertex'
 
       call InterpolateData(0,VarVis,VisEff,dPdX,NodalData1,NodalCounter)
       
       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_vtv.msh','FORMATTED', &
                                                 'SEQUENTIAL','UNKNOWN',debug) 
	 write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"VisT vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,max(0.0,NodalData1(i) - VisLam)
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
         close(IOpst) 
       endif
     endif

   endif

   !
   ! *** TEMPERATURES ***
   !
   if( SolveEnthalpy )then
     if( PostC(VarT) .and. .not. UseGMSHdataonly .and. allocated(T) )then
       write(IOdbg,*)'Writing GMSH temperature cell relative to',Tref
       write(IOdef,*)'Writing GMSH temperature cell relative to',Tref
       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"Temperature cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
         if( UseGMSHwalls )then
           write(IOpst,'(i10)')   Ncel+Nbnd
         else
           write(IOpst,'(i10)')   Ncel  
         endif

	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') i,T(i)-Tref
	 end do
	 if( UseGMSHwalls )then
	   do i=1,Nbnd
             write(IOpst,'(i8,1x,1pe10.3)') Ncel+i,T(Ncel+i)-Tref    
	   end do
	 endif 
       endif
       write(IOpst,'(A)')     '$EndElementData'
     elseif( PostC(VarT) .and. UseGMSHdataonly .and. allocated(T) )then
       write(IOdbg,*)'Writing GMSH data temperature cell relative to',Tref
       write(IOdef,*)'Writing GMSH data temperature cell relative to',Tref
       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_tc.msh','FORMATTED', &
                                                'SEQUENTIAL','UNKNOWN',debug) 
	 write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$ElementData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"Temperature cell"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   NCellOffset+Ncel 

	 icnt = 0
	 do i=1,Nbnd
           if( bnd(i)%rid > 0 )then
	     icnt = icnt + 1
             write(IOpst,'(i8,1x,1pe10.3)') icnt,T(Ncel+i)-Tref
	   endif
	 end do
	 do i=1,Ncel
           write(IOpst,'(i8,1x,1pe10.3)') NCellOffset+i,T(i)-Tref
	 end do
       endif
       write(IOpst,'(A)')     '$EndElementData'
       close(IOpst) 
     endif 

     if( PostV(VarT) .and. .not. UseGMSHdataonly .and. allocated(T) )then
       write(IOdbg,*)'Writing GMSH temperature vertex relative to',Tref
       write(IOdef,*)'Writing GMSH temperature vertex relative to',Tref

       call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

       if( UseGMSHbinary )then

       else
	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"Temperature vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,Nodaldata1(i)-Tref
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
       endif
     elseif( PostV(VarT) .and. UseGMSHdataonly .and. allocated(T) )then
       write(IOdbg,*)'Writing GMSH data temperature vertex relative to',Tref
       write(IOdef,*)'Writing GMSH data temperature vertex relative to',Tref

       call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

       if( UseGMSHbinary )then

       else
         call openfile(IOpst,casename,'_tv.msh','FORMATTED', &
                                                'SEQUENTIAL','UNKNOWN',debug) 
	 write(IOpst,'(A)')     '$MeshFormat'
         write(IOpst,'(A,A,A)') '2.0 0 8' 
         write(IOpst,'(A)')     '$EndMeshFormat' 

	 write(IOpst,'(A)')     '$NodeData'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(A)')     '"Temperature vertex"'
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(e10.3)') Time
	 write(IOpst,'(A)')     '3'
	 write(IOpst,'(i10)')    0   !  Iter
	 write(IOpst,'(A)')     '1'
	 write(IOpst,'(i10)')   Nvrt

	 do i=1,Nvrt
           write(IOpst,'(i8,1x,1pe10.3)') i,Nodaldata1(i)-Tref
	 end do
         write(IOpst,'(A)')     '$EndNodeData'
         close(IOpst) 
       endif
     endif 
   endif

   if( .not. UseGMSHdataonly ) close(IOpst)
   write(*,*)'Done'

end subroutine dolfyn2gmsh
