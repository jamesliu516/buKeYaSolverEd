!
! Copyright 2004-2009 Henk Krus, Cyclone Fluid Dynamics BV
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
! This is vtk.f90, the vtk interface to paraview, visit or mayavi
!
! Particles added by Johan Jacobs, Simex-Technology, summer 2007 
!
! Binary mode ideas and extensions by Max Staufer, Darmstadt, december 2008
! Relies on 'stream' and big endian IO. 
!
subroutine dolfyn2vtk

   use constants
   use geometry
   use variables
   use particles
   
   !
   ! binary extensions:
   !
  !logical              :: Binary = .true.   ! .true. or .false. (=ascii) 
   character(len=1), &                       ! end-character for binary-record
               parameter:: NWL    = char(10) ! finalize
   character(len=128)   :: sbuff             ! string buffer for binary IO
   !
   
   character(len=12)    :: string1, string2
   integer              :: indx(8)

   integer              :: datum(8)
   character (len=12)   :: clock(3)
   character (len=10)   :: datestring 
   character (len=13)   :: frmt

   character (len=96)   :: CaseNameTransient = ' '

   real, dimension(3)   :: Xp, X, ds 
   real, dimension(3)   :: Xpn, Normal 
   real, dimension(10)  :: tmp               ! temp array 
   
   real, allocatable    :: NodalData1(:)     ! array of data at vertices
   real, allocatable    :: NodalData2(:)    
   real, allocatable    :: NodalData3(:)    
   
   integer, allocatable :: NodalCounter(:)  

   integer, allocatable :: VTKcells(:)  

   integer, save        :: ICounter = 0      ! counter for transient data

   call date_and_time(clock(1),clock(2),clock(3),datum)
   write(datestring,'(i2.2,''/'',i2.2,''/'',i4)') datum(3),datum(2),datum(1)

   allocate( VTKcells(Ncel+Nbnd),stat=istat)    
   allocate( NodalData1(Nvrt),stat=istat)    
   allocate( NodalCounter(Nvrt),stat=istat)    
   call TrackMemory(istat,Ncel+Nbnd+2*Nvrt,'Work arrays for NodalData allocated')
   
   write(*,*)        'Writing VTK data file...'
   write(IOdbg,*)    'Writing VTK data file...'
   write(IOdbg,'(A,A)')  ' Case:  ',casename(1:lens(casename)) 
   write(IOdbg,'(A,A)')  ' Title: ',title(1:lens(title)) 
   write(IOdbg,'(A,A)')  ' Date : ',datestring
   
   if( Transient )then
     n = lens(CaseName)
     ICounter = ICounter + 1
     if( ICounter > 9999 )then 
       write(*,*) '+++ Warning: Transient VTK output data file counter reset'
       ICounter = 1
     endif
     write(string1,'(i4.4)') ICounter
     CaseNameTransient( 1 : n ) = CaseName(1:n)
     CaseNameTransient(n+1:n+1) = '_'
     CaseNameTransient(n+2:n+5) = string1

     write(IOdbg,'(A,A)')  ' File:  ',CaseNameTransient(1:n+5)//'.vtk'

     if( UseVTKbinary )then
       call openfile(IOpst,CaseNameTransient,'.vtk','UNFORMATTED_BIG_ENDIAN', &
                                             'STREAM','UNKNOWN',debug)
     else
       call openfile(IOpst,CaseNameTransient,'.vtk','FORMATTED', &
                                           'SEQUENTIAL','UNKNOWN',debug)
     endif
   else
     if( UseVTKbinary )then
       call openfile(IOpst,casename,'.vtk','UNFORMATTED_BIG_ENDIAN', &
                                           'STREAM','UNKNOWN',debug)
     else
       call openfile(IOpst,casename,'.vtk','FORMATTED', &
                                           'SEQUENTIAL','UNKNOWN',debug)
     endif
   endif

   if( UseVTKbinary )then
     write(IOpst)           '# vtk DataFile Version 3.0'//NWL
     write(sbuff,'(A,A,A)')  casename(1:lens(casename)),': ',title(1:lens(title)) 
     write(IOpst)            trim(sbuff)//NWL
     write(IOpst)           'BINARY'//NWL
   else
     write(IOpst,'(A)')     '# vtk DataFile Version 3.0'
     write(IOpst,'(A,A,A)')  casename(1:lens(casename)),': ',title(1:lens(title)) 
     write(IOpst,'(A)')     'ASCII'
   endif
   !
   ! algemene info
   !
   if( UseVTKbinary )then
     write(IOpst)           'DATASET UNSTRUCTURED_GRID'//NWL
   else
     write(IOpst,'(A)')     'DATASET UNSTRUCTURED_GRID'
   endif     
   !
   ! vertices
   !
   if( UseVTKbinary )then
     write(sbuff,'(A,i8,A)') 'POINTS ',Nvrt,' float'
     write(IOpst) trim(sbuff)//NWL
     
     write(IOpst) (Vert(i,1),Vert(i,2),Vert(i,3),i=1,Nvrt)

   else
     write(IOpst,'(A,i8,A)') 'POINTS ',Nvrt,' float'
     do i=1,Nvrt
       write(IOpst,'(3(1x,1pe12.5))') Vert(i,:)
     end do
   endif
   !
   ! the VTK file format needs to know in advance what is coming
   !
   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)

   !  first walk through the cells establish the type
   !
   !    6 hexaeder:  1,2,3,4 5,6,7,8  
   !    5 pentaeder: 1,2,3,3 5,6,7,7  
   !    5 pyramid:   1,2,3,4 5,5,5,5  
   !    4 tetraeder: 1,2,3,3 4,4,4,4  
   !    1 quad:      1,2,3,4          
   !    1 triangle:  1,2,3,3          
   !
   NHex  = 0
   NPen  = 0
   NPyr  = 0
   NTet  = 0
   NQuad = 0
   NTri  = 0
   
   do i=1,Ncel
     read(IOcel,*) idummy,(indx(j),j=1,8)
    !indx(1:8) = Original_Cells(i,1:8)
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
         VTKcells(i) = 12
       else if( i5 /= i6 )then
         NPen = NPen + 1
         VTKcells(i) = 13
       else if( i3 /= i4 )then
         NPyr  = NPyr + 1
         VTKcells(i) = 14
       else if( i3 == i4 .and. i5 == i6 )then
         NTet   = NTet + 1
         VTKcells(i) = 10
       else
         write(*,*)'Unknown VTK shape, cell ',i,' : ',i1,i2,i3,i4,' ...'
       endif
     else
       !
       ! Reserved for baffles 
       !
       write(*,*)'Error: Unsuppported feature'
              
       if( UseVTKwalls )then
         if( i4 /= -1 )then
           NQuad   = NQuad + 1
           VTKcells(i) = 9
         else 
           NTri    = NTri + 1
           VTKcells(i) = 5
         endif
       endif
     endif
   end do
   
   close(IOcel)

   if( UseVTKwalls )then
     do ib=1,Nbnd
       indx(1:4) = bnd(ib)%vertices
       if( indx(4) /= -1 )then
         NQuad   = NQuad + 1
         VTKcells(Ncel+ib) = 9    
       elseif( indx(4) == -1 )then
         NTri    = NTri + 1
         VTKcells(Ncel+ib) = 5    
       else
         write(*,*)'*** Error: internal VTK boundary cell format (1)'
       endif
     end do
   endif
   
   if( NHex  > 0 )write(IOdbg,*)'Number of VTK hexa''s     :',NHex 
   if( NPen  > 0 )write(IOdbg,*)'Number of VTK prism''s    :',NPen 
   if( NPyr  > 0 )write(IOdbg,*)'Number of VTK pyramids''s :',NPyr 
   if( NTet  > 0 )write(IOdbg,*)'Number of VTK tetra''s    :',NTet 
   if( NQuad > 0 )write(IOdbg,*)'Number of VTK quad''s     :',NQuad
   if( NTri  > 0 )write(IOdbg,*)'Number of VTK triangles''s:',NTri 

   NVTKc = NHex * 9  + NPen * 7 + NPyr * 6 + NTet * 5
   NVTKb = NQuad * 5 + NTri * 4
   !
   ! cells
   !
   if( UseVTKbinary )then
     if( UseVTKwalls )then
       write(sbuff,'(A,i12,i12)') 'CELLS ',Ncel+Nbnd,NVTKc + NVTKb
     else
       write(sbuff,'(A,i12,i12)') 'CELLS ',Ncel,NVTKc
     endif
     write(IOpst) trim(sbuff)//NWL
   else
     if( UseVTKwalls )then
       write(IOpst,'(A,i12,i12)') 'CELLS ',Ncel+Nbnd,NVTKc + NVTKb
     else
       write(IOpst,'(A,i12,i12)') 'CELLS ',Ncel,NVTKc
     endif
   endif
   
   call openfile(IOcel,casename,'.cel','FORMATTED', &
                                       'SEQUENTIAL','OLD',debug)
   rewind(IOcel)
   
!1 format(i2,8(1x,i4))
!2 format(i2,8(1x,i5))   
!3 format(i2,8(1x,i6))   
!4 format(i2,8(1x,i7))  
!5 format(i2,8(1x,i8))  
! 
!  if( Nvrt <= 99999999 ) assign 5 to ifmt 
!  if( Nvrt <=  9999999 ) assign 4 to ifmt 
!  if( Nvrt <=   999999 ) assign 3 to ifmt 
!  if( Nvrt <=    99999 ) assign 2 to ifmt 
!  if( Nvrt <=     9999 ) assign 1 to ifmt 
      
   frmt = '(i2,8(1x,i4))'
   if( Nvrt <=     9999 )then
     frmt = '(i2,8(1x,i4))'
   else if( Nvrt <=    99999 )then
     frmt = '(i2,8(1x,i5))'
   else if( Nvrt <=   999999 )then
     frmt = '(i2,8(1x,i6))'
   else if( Nvrt <=  9999999 )then
     frmt = '(i2,8(1x,i7))'
   else if( Nvrt <= 99999999 )then
     frmt = '(i2,8(1x,i8))'
   endif
   
   if( UseVTKbinary )then
     do i=1,Ncel
       read(IOcel,*) idummy,(indx(j),j=1,8)
      !indx(1:8) = Original_Cells(i,1:8)

       if( VTKcells(i) == 12 )then
         write(IOpst)      8, &
                                      indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1, &
                                      indx(5)-1,indx(6)-1,indx(7)-1,indx(8)-1
       elseif( VTKcells(i) == 13 )then
         write(IOpst)      6, &
                                      indx(1)-1,indx(2)-1,indx(3)-1, &
                                      indx(5)-1,indx(6)-1,indx(7)-1
       elseif( VTKcells(i) == 14 )then
         write(IOpst)      5, &
                                      indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1, &
                                      indx(5)-1 
       elseif( VTKcells(i) == 10 )then
         write(IOpst)      4, &
                                      indx(1)-1,indx(2)-1,indx(3)-1, &
                                      indx(5)-1 
       else
         write(*,*)'*** Error: internal VTK cell format (bin)'
       endif
     end do
   else
     do i=1,Ncel
       read(IOcel,*) idummy,(indx(j),j=1,8)
      !indx(1:8) = Original_Cells(i,1:8)

       if( VTKcells(i) == 12 )then
         write(IOpst,frmt) 8, &
                                      indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1, &
                                      indx(5)-1,indx(6)-1,indx(7)-1,indx(8)-1
       elseif( VTKcells(i) == 13 )then
         write(IOpst,frmt) 6, &
                                      indx(1)-1,indx(2)-1,indx(3)-1, &
                                      indx(5)-1,indx(6)-1,indx(7)-1
       elseif( VTKcells(i) == 14 )then
         write(IOpst,frmt) 5, &
                                      indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1, &
                                      indx(5)-1 
       elseif( VTKcells(i) == 10 )then
         write(IOpst,frmt) 4, &
                                      indx(1)-1,indx(2)-1,indx(3)-1, &
                                      indx(5)-1 
       else
         write(*,*)'*** Error: internal VTK cell format (asc)'
       endif
     end do
   endif
   !
   ! next loop over boundaries
   !
   if( UseVTKbinary )then
     if( UseVTKwalls )then
       do ib=1,Nbnd
         indx(1:4) = bnd(ib)%vertices
         if( indx(4) /= -1 )then
           write(IOpst)      4,indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1
         else 
           write(IOpst)      3,indx(1)-1,indx(2)-1,indx(3)-1
         endif
       end do
     endif
   else
     if( UseVTKwalls )then
       do ib=1,Nbnd
         indx(1:4) = bnd(ib)%vertices
         if( indx(4) /= -1 )then
           write(IOpst,frmt) 4,indx(1)-1,indx(2)-1,indx(3)-1,indx(4)-1
         else 
           write(IOpst,frmt) 3,indx(1)-1,indx(2)-1,indx(3)-1
         endif
       end do
     endif
   endif
   close(IOcel)
   if( Debug > 2 ) write(*,*)'Cells done'
   
   !
   ! cell types
   !
   if( UseVTKbinary )then
     if( UseVTKwalls )then
       write(sbuff,'(A,i8)') 'CELL_TYPES ',Ncel+Nbnd
     else
       write(sbuff,'(A,i8)') 'CELL_TYPES ',Ncel 
     endif
     write(IOpst) trim(sbuff)//NWL

     if( UseVTKwalls )then
       write(IOpst) (VTKcells(i),i=1,Ncel+Nbnd)
     else
       write(IOpst) (VTKcells(i),i=1,Ncel)
     endif

   else
     if( UseVTKwalls )then
       write(IOpst,'(A,i8)') 'CELL_TYPES ',Ncel+Nbnd
     else
       write(IOpst,'(A,i8)') 'CELL_TYPES ',Ncel 
     endif

     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(i2,9(1x,i2))') (VTKcells(j),j=i,k)
     end do
     if( UseVTKwalls )then
       do ib=1,Nbnd,10
         k = min(ib+9,Nbnd)
         write(IOpst,'(i2,9(1x,i2))') (VTKcells(Ncel+j),j=ib,k)
       end do
     endif

   endif
   !
   ! cell id's
   !
   write(IOdbg,*)'Writing VTK cell ids'
   
   if( UseVTKbinary )then
     if( UseVTKwalls )then
       write(sbuff,'(A,i8)')  'CELL_DATA ',Ncel+Nbnd
     else
       write(sbuff,'(A,i8)')  'CELL_DATA ',Ncel
     endif
     write(IOpst)	     trim(sbuff)//NWL
    
    !write(IOpst)	     'SCALARS materials int'//NWL
     write(IOpst)	     'SCALARS type_id int'//NWL
     write(IOpst)	     'LOOKUP_TABLE default'//NWL
   else
     if( UseVTKwalls )then
       write(IOpst,'(A,i8)')  'CELL_DATA ',Ncel+Nbnd
     else
       write(IOpst,'(A,i8)')  'CELL_DATA ',Ncel
     endif
    !write(IOpst,'(A)')      'SCALARS materials int'
     write(IOpst,'(A)')      'SCALARS type_id int'
     write(IOpst,'(A)')      'LOOKUP_TABLE default'
   endif
   
   if( UseVTKbinary )then

     if( UseVTKwalls )then
       IBoffset = maxval(cell(:)%ctid)
       write(IOpst) (Cell(i)%ctid,i=1,Ncel)
       write(IOpst) (Bnd(i)%rid+(IBoffset+1),i=1,Nbnd)
     else
       write(IOpst) (Cell(i)%ctid,i=1,Ncel)	
     endif
     
   else
   
     IBoffset = maxval(cell(:)%ctid)
     write(IOdbg,*) 'IBoffset: ',IBoffset

     if( IBoffset <= 9 )then
       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(i2))') (cell(j)%ctid,j=i,k)
       end do
     else if( IBoffset <= 99 )then
       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(i3))') (cell(j)%ctid,j=i,k)
       end do
     else
       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(i4))') (cell(j)%ctid,j=i,k)
       end do
     endif

     if( UseVTKwalls )then
       write(IOdbg,*) 'IBoffset+Nreg: ',IBoffset+Nreg
       if( IBoffset+Nreg+1 <= 9 )then
         do ib=1,Nbnd,10
           k = min(ib+9,Nbnd)
           write(IOpst,'(10(i2))') (bnd(j)%rid+(IBoffset+1),j=ib,k)
         end do
       else if( IBoffset+Nreg+1 <= 99 )then
         do ib=1,Nbnd,10
           k = min(ib+9,Nbnd)
           write(IOpst,'(10(i3))') (bnd(j)%rid+(IBoffset+1),j=ib,k)
         end do
       else
         do ib=1,Nbnd,10
           k = min(ib+9,Nbnd)
           write(IOpst,'(10(i4))') (bnd(j)%rid+(IBoffset+1),j=ib,k)
         end do
       endif
     endif
     
   endif
   !
   ! data
   !
   !
   ! VECTORS CELL DATA
   ! (in binary mode, the gradients as well)
   !
   write(IOdbg,*)'Writing VTK vectors'
   if( UseVTKbinary )then
     if( UseVTKwalls )then   
       write(IOpst) 'VECTORS velocity float'//NWL
       write(IOpst) (U(i),V(i),W(i),i=1,Ncel+Nbnd)
       
      !write(IOpst) 'VECTORS dUdX float'//NWL
      !write(IOpst) (dUdX(i,1),dUdX(i,2),dUdX(i,3),i=1,Ncel)
      !write(IOpst) (0.0,0.0,0.0,i=1,Nbnd)
       
      !write(IOpst) 'VECTORS dVdX float'//NWL
      !write(IOpst) (dVdX(i,1),dVdX(i,2),dVdX(i,3),i=1,Ncel)
      !write(IOpst) (0.0,0.0,0.0,i=1,Nbnd)
        
      !write(IOpst) 'VECTORS dWdX float'//NWL
      !write(IOpst) (dWdX(i,1),dWdX(i,2),dWdX(i,3),i=1,Ncel)
      !write(IOpst) (0.0,0.0,0.0,i=1,Nbnd)
        
      !call GradientPhi(VarP,P,dPdX)      

      !write(IOpst) 'VECTORS dPdX float'//NWL
      !write(IOpst) (dPdX(i,1),dPdX(i,2),dPdX(i,3),i=1,Ncel)
      !write(IOpst) (0.0,0.0,0.0,i=1,Nbnd)
 
      !write(IOpst) 'VECTORS dPPdx float'//NWL
      !write(IOpst) (dPPdX(i,1),dPPdX(i,2),dPPdX(i,3),i=1,Ncel)
      !write(IOpst) (0.0,0.0,0.0,i=1,Nbnd)
     else
       write(IOpst) 'VECTORS velocity float'//NWL
       write(IOpst) (U(i),V(i),W(i),i=1,Ncel)
       
      !write(IOpst) 'VECTORS dUdX float'//NWL
      !write(IOpst) (dUdX(i,1),dUdX(i,2),dUdX(i,3),i=1,Ncel)
      
      !write(IOpst) 'VECTORS dVdX float'//NWL
      !write(IOpst) (dVdX(i,1),dVdX(i,2),dVdX(i,3),i=1,Ncel)
       
      !write(IOpst) 'VECTORS dWdX float'//NWL
      !write(IOpst) (dWdX(i,1),dWdX(i,2),dWdX(i,3),i=1,Ncel)
       
      !call GradientPhi(VarP,P,dPdX)      

      !write(IOpst) 'VECTORS dPdX float'//NWL
      !write(IOpst) (dPdX(i,1),dPdX(i,2),dPdX(i,3),i=1,Ncel)

      !write(IOpst) 'VECTORS dPPdx float'//NWL
      !write(IOpst) (dPPdX(i,1),dPPdX(i,2),dPPdX(i,3),i=1,Ncel)
     endif     
   else
     write(IOpst,'(A)') 'VECTORS velocity float'

     do i=1,Ncel
       write(IOpst,'(3(1x,1pe10.3))') U(i),V(i),W(i)
     end do
     if( UseVTKwalls )then
       do i=1,Nbnd
         write(IOpst,'(3(1x,1pe10.3))') U(Ncel+i),V(Ncel+i),W(Ncel+i)
       end do
     endif 
   endif  

   !
   ! SCALARS CELL DATA
   !
   write(IOdbg,*)'Writing VTK magnitude'
   if( UseVTKbinary )then
     write(IOpst)	'SCALARS velmag float 1'//NWL
     write(IOpst)	'LOOKUP_TABLE default'//NWL

     if( UseVTKwalls )then
       write(IOpst) &
	 (sqrt( U(i)**2 + V(i)**2 + W(i)**2 ),i=1,Ncel+Nbnd)
     else
       write(IOpst) &
	 (sqrt( U(i)**2 + V(i)**2 + W(i)**2 ),i=1,Ncel)
     endif
     
   else
     write(IOpst,'(A)') 'SCALARS velmag float 1'
     write(IOpst,'(A)') 'LOOKUP_TABLE default'
   
     do i=1,Ncel,10
       k = min(i+9,Ncel)
       write(IOpst,'(10(1x,1pe10.3))') &
         (sqrt( U(j)**2 + V(j)**2 + W(j)**2 ),j=i,k)
     end do

     if( UseVTKwalls )then
       do i=1,Nbnd,10
         k = min(i+9,Nbnd)
         write(IOpst,'(10(1x,1pe10.3))') &
           (sqrt( U(Ncel+j)**2 + V(Ncel+j)**2 + W(Ncel+j)**2 ),j=i,k)
       end do
     endif
     
   endif  
   !
   ! pressure cell data
   !
   if( PostC(VarP) )then
     write(IOdbg,*)'Writing VTK pressure'
     if( UseVTKbinary )then
       write(IOpst)	  'SCALARS pressure float 1'//NWL
       write(IOpst)	  'LOOKUP_TABLE default'//NWL

       if( UseVTKwalls )then
	 write(IOpst) (P(i),i=1,Ncel+Nbnd)
       else
	 write(IOpst) (P(i),i=1,Ncel)
       endif	 
     else
       write(IOpst,'(A)') 'SCALARS pressure float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (P(j),j=i,k)
       end do

       if( UseVTKwalls )then
         do i=1,Nbnd,10
           k = min(i+9,Nbnd)
           write(IOpst,'(10(1x,1pe10.3))') (P(Ncel+j),j=i,k)
         end do
       endif
     endif
   endif

   !
   ! tke cell data
   !
   if( PostC(VarTE) .and. allocated(TE) )then
     write(IOdbg,*)'Writing VTK turbulent kinetic energy'
     if( UseVTKbinary )then
       write(IOpst)       'SCALARS k float 1'//NWL
       write(IOpst)       'LOOKUP_TABLE default'//NWL

       if( UseVTKwalls )then
         write(IOpst) (max(0.0,TE(i)),i=1,Ncel+Nbnd)
       else
         write(IOpst) (max(0.0,TE(i)),i=1,Ncel)
       endif       
     else
       write(IOpst,'(A)') 'SCALARS k float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (max(0.0,TE(j)),j=i,k)
       end do

       if( UseVTKwalls )then
         do i=1,Nbnd,10
           k = min(i+9,Nbnd)
           write(IOpst,'(10(1x,1pe10.3))') (max(0.0,TE(Ncel+j)),j=i,k)
         end do
       endif
     endif
   endif
   !
   ! turbulent dissipation cell data
   !
   if( PostC(VarED) .and. allocated(ED) )then
     write(IOdbg,*)'Writing VTK turbulent dissipation'
     if( UseVTKbinary )then
       write(IOpst)       'SCALARS eps float 1'//NWL
       write(IOpst)       'LOOKUP_TABLE default'//NWL

       if( UseVTKwalls )then
         write(IOpst) (max(0.0,ED(i)),i=1,Ncel+Nbnd)
       else
         write(IOpst) (max(0.0,ED(i)),i=1,Ncel)
       endif
     
     else
       write(IOpst,'(A)') 'SCALARS eps float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (max(0.0,ED(j)),j=i,k)
       end do

       if( UseVTKwalls )then
         do i=1,Nbnd,10
           k = min(i+9,Nbnd)
           write(IOpst,'(10(1x,1pe10.3))') (max(0.0,ED(Ncel+j)),j=i,k)
         end do
       endif
     endif
   endif
   !
   ! temperature cell data
   !
   if( PostC(VarT) .and. allocated(T) )then
     write(IOdbg,*)'Writing VTK temperature relative to',Tref
     if( UseVTKbinary )then
       write(IOpst)       'SCALARS temperature float 1'//NWL
       write(IOpst)       'LOOKUP_TABLE default'//NWL

       if( UseVTKwalls )then
         write(IOpst) (T(i)-Tref,i=1,Ncel+Nbnd)
       else
         write(IOpst) (T(i)-Tref,i=1,Ncel)
       endif
     
     else
       write(IOpst,'(A)') 'SCALARS temperature float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (T(j)-Tref,j=i,k)
       end do

       if( UseVTKwalls )then
         do i=1,Nbnd,10
           k = min(i+9,Nbnd)
           write(IOpst,'(10(1x,1pe10.3))') (T(Ncel+j)-Tref,j=i,k)
         end do
       endif
     endif
   endif 
   !
   ! effective viscosity cell data
   !
   if( PostC(VarVIS) )then
     write(IOdbg,*)'Writing VTK effective viscosity'
     if( UseVTKbinary )then
       write(IOpst)       'SCALARS viscosity float 1'//NWL
       write(IOpst)       'LOOKUP_TABLE default'//NWL

       if( UseVTKwalls )then
         write(IOpst) (VisEff(i),i=1,Ncel+Nbnd)
       else
         write(IOpst) (VisEff(i),i=1,Ncel)
       endif
 
     else
       write(IOpst,'(A)') 'SCALARS viscosity float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (VisEff(j),j=i,k)
       end do

       if( UseVTKwalls )then
         do i=1,Nbnd,10
           k = min(i+9,Nbnd)
           write(IOpst,'(10(1x,1pe10.3))') (VisEff(Ncel+j),j=i,k)
         end do
       endif
     endif
   endif 
   !
   ! variable density cell data
   !
   if( PostC(VarDEN) )then
     write(IOdbg,*)'Writing VTK density'
     if( UseVTKbinary )then
       write(IOpst)       'SCALARS density float 1'//NWL
       write(IOpst)       'LOOKUP_TABLE default'//NWL

       if( UseVTKwalls )then
         write(IOpst) (Den(i),i=1,Ncel+Nbnd)
       else
         write(IOpst) (Den(i),i=1,Ncel)
       endif
     else
       write(IOpst,'(A)') 'SCALARS density float 1'
       write(IOpst,'(A)') 'LOOKUP_TABLE default'

       do i=1,Ncel,10
         k = min(i+9,Ncel)
         write(IOpst,'(10(1x,1pe10.3))') (Den(j),j=i,k)
       end do

       if( UseVTKwalls )then
         do i=1,Nbnd,10
           k = min(i+9,Nbnd)
           write(IOpst,'(10(1x,1pe10.3))') (Den(Ncel+j),j=i,k)
         end do
       endif
     endif
   endif   
   !
   ! TODO: Binary 
   !
   if( SolveScalars )then
     if( PostC(VarSC) )then
       do is=1,NScal
         write(IOdbg,*)'Writing VTK scalar(s) ',Variable(NVar+is)
         write(IOpst,'(A,A,A)') 'SCALARS ',&
           Variable(NVar+is)(1:lens(Variable(NVar+is))),' float 1'
         write(IOpst,'(A)') 'LOOKUP_TABLE default'

         do i=1,Ncel,10
           k = min(i+9,Ncel)
           write(IOpst,'(10(1x,1pe10.3))') (max(0.0,SC(j,is)),j=i,k)
         end do

         if( UseVTKwalls )then
           do i=1,Nbnd,10
             k = min(i+9,Nbnd)
             write(IOpst,'(10(1x,1pe10.3))') (max(0.0,SC(Ncel+j,is)),j=i,k)
           end do
         endif
       end do
     endif
   endif

   !
   ! SCALARS VERTEX DATA
   !
   if( count( PostV(:)) > 0 )then
     if( UseVTKbinary )then
       write(sbuff,'(A,i8)') 'POINT_DATA ',Nvrt
       write(IOpst) trim(sbuff)//NWL
     else
       write(IOpst,'(A,i8)') 'POINT_DATA ',Nvrt       
     endif
   endif

   if( PostV(VarU) .or. PostV(VarV) .or. PostV(VarW) )then

     allocate( NodalData2(Nvrt),stat=istat)    
     allocate( NodalData3(Nvrt),stat=istat)    

     write(IOdbg,*)'Writing VTK data on nodes'
     call TrackMemory(istat,2*Nvrt,'Work arrays for extra NodalData allocated')

     call InterpolateData(0,VarU,U,dUdX,NodalData1,NodalCounter)
     call InterpolateData(0,VarV,V,dVdX,NodalData2,NodalCounter)
     call InterpolateData(0,VarW,W,dWdX,NodalData3,NodalCounter)
 
     if( UseVTKbinary )then
       !
       ! dump vectors on vertices
       !
       write(IOdbg,*)'Writing VTK vectors on nodes ',Nvrt
       write(IOpst)  'VECTORS nvelocity float'//NWL

       write(IOpst) (Nodaldata1(i),NodalData2(i),NodalData3(i),i=1,Nvrt)

       !
       ! dump velmagnitude on vertices
       !
       write(IOdbg,*)'Writing VTK velocity magnitude on vertices'
       write(IOpst)  'SCALARS nvelmag float 1'//NWL
       write(IOpst)  'LOOKUP_TABLE default'//NWL

       do i=1,Nvrt
         write(IOpst) sqrt( Nodaldata1(i)**2 + NodalData2(i)**2 + &
                            NodalData3(i)**2 )
       end do
     
     else
       !
       ! dump vectors on vertices
       !
       write(IOdbg,*)'Writing VTK vectors on nodes ',Nvrt
       write(IOpst,'(A)')    'VECTORS nvelocity float'

       do i=1,Nvrt
         if( abs(Nodaldata1(i)) < Small ) Nodaldata1(i) = 0.0
         if( abs(Nodaldata2(i)) < Small ) Nodaldata2(i) = 0.0
         if( abs(Nodaldata3(i)) < Small ) Nodaldata3(i) = 0.0

         write(IOpst,'(3(1x,1pe10.3))') Nodaldata1(i),NodalData2(i),NodalData3(i)
       end do

       !
       ! dump velmagnitude on vertices
       !
       write(IOdbg,*)'Writing VTK velocity magnitude on vertices'
       write(IOpst,'(A)')    'SCALARS nvelmag float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       do i=1,Nvrt
         if( abs(Nodaldata1(i)) < Small ) Nodaldata1(i) = 0.0
         if( abs(Nodaldata2(i)) < Small ) Nodaldata2(i) = 0.0
         if( abs(Nodaldata3(i)) < Small ) Nodaldata3(i) = 0.0

         tmpmag = sqrt( Nodaldata1(i)**2 + NodalData2(i)**2 + NodalData3(i)**2 )

         write(IOpst,'(1x,1pe10.3)') tmpmag
       end do
     endif
   endif

   if( PostV(VarP) )then   
     write(IOdbg,*)'Writing VTK pressure on vertices'
     if( UseVTKbinary )then
       write(IOpst)          'SCALARS npressure float 1'//NWL
       write(IOpst)          'LOOKUP_TABLE default'//NWL

       call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)

       write(IOpst) (Nodaldata1(i),i=1,Nvrt)
     
     else
       write(IOpst,'(A)')    'SCALARS npressure float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       call InterpolateData(0,VarP,P,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt,10
         k = min(i+9,Nvrt)
         write(IOpst,'(10(1x,1pe10.3))') (Nodaldata1(j),j=i,k)
       end do
     endif
   endif

   if( PostV(VarT) .and. allocated(T) )then
     write(IOdbg,*)'Writing VTK temperature on vertices relative to ',Tref
     if( UseVTKbinary )then
       write(IOpst)          'SCALARS ntemperature float 1'//NWL
       write(IOpst)          'LOOKUP_TABLE default'//NWL

       call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

       write(IOpst) (Nodaldata1(i),i=1,Nvrt)

     else
       write(IOpst,'(A)')    'SCALARS ntemperature float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       call InterpolateData(0,VarT,T,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt,10
         k = min(i+9,Nvrt)
         write(IOpst,'(10(1x,1pe11.4))')  (Nodaldata1(j)-Tref,j=i,k)
       end do
     endif
   endif

   if( PostV(VarTE) .and. allocated(TE) )then
     write(IOdbg,*)'Writing VTK turbulent kinetic energy on vertices'
     if( UseVTKbinary )Then
       write(IOpst)          'SCALARS nk float 1'//NWL
       write(IOpst)          'LOOKUP_TABLE default'//NWL

       call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

       write(IOpst) (Nodaldata1(i),i=1,Nvrt)

     else
       write(IOpst,'(A)')    'SCALARS nk float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       call InterpolateData(0,VarTE,TE,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt,10
         k = min(i+9,Nvrt)
         write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
       end do
     endif
   endif

   if( PostV(VarED) .and. allocated(ED) )then
     write(IOdbg,*)'Writing VTK turbulent dissipation on vertices'
     if( UseVTKbinary )then
       write(IOpst)          'SCALARS neps float 1'//NWL
       write(IOpst)          'LOOKUP_TABLE default'//NWL

       call InterpolateData(0,VarED,ED,dPdX,NodalData1,NodalCounter)

       write(IOpst) (Nodaldata1(i),i=1,Nvrt)
 
     else
       write(IOpst,'(A)')    'SCALARS neps float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       call InterpolateData(0,VarED,ED,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt,10
         k = min(i+9,Nvrt)
         write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
       end do
     endif
   endif

   if( PostV(VarVIS) )then
     write(IOdbg,*)'Writing VTK effective viscosity on vertices'
     if( UseVTKbinary )then
       write(IOpst)          'SCALARS neffvis float 1'//NWL
       write(IOpst)          'LOOKUP_TABLE default'//NWL

       call InterpolateData(0,VarVIS,VisEff,dPdX,NodalData1,NodalCounter)

       write(IOpst) (Nodaldata1(i),i=1,Nvrt)

     else
       write(IOpst,'(A)')    'SCALARS neffvis float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       call InterpolateData(0,VarVIS,VisEff,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt,10
         k = min(i+9,Nvrt)
         write(IOpst,'(10(1x,1pe11.4))') (max(0.0,Nodaldata1(j)),j=i,k)
       end do
     endif
   endif

   if( PostV(VarDEN) .and. allocated(DEN) )then
     write(IOdbg,*)'Writing VTK density on vertices'
     if( UseVTKbinary )then
       write(IOpst)          'SCALARS ndensity float 1'//NWL
       write(IOpst)          'LOOKUP_TABLE default'//NWL

       call InterpolateData(0,VarDEN,DEN,dPdX,NodalData1,NodalCounter)

       write(IOpst) (Nodaldata1(i),i=1,Nvrt)

     else
       write(IOpst,'(A)')    'SCALARS ndensity float 1'
       write(IOpst,'(A)')    'LOOKUP_TABLE default'

       call InterpolateData(0,VarDEN,DEN,dPdX,NodalData1,NodalCounter)

       do i=1,Nvrt,10
         k = min(i+9,Nvrt)
         write(IOpst,'(10(1x,1pe11.4))')  (max(0.0,Nodaldata1(j)),j=i,k)
       end do
     endif
   endif

   !
   ! TODO: BINARY
   !
   if( SolveScalars )then
     if( PostV(VarSC) )then
         do is=1,NScal
         
           call InterpolateData(0,VarS(is),SC(:,is),dPdX,NodalData1,NodalCounter)

           write(IOdbg,*)'Writing VTK scalar(s) ',&
                                     Variable(NVar+is),' on vertices'
           write(IOpst,'(A,A,A)') 'SCALARS ',&
             Variable(NVar+is)(1:lens(Variable(NVar+is))),' float 1'
           write(IOpst,'(A)') 'LOOKUP_TABLE default'

           do i=1,NVrt,10
             k = min(i+9,Nvrt)
             write(IOpst,'(10(1x,1pe10.3))') (max(0.0,Nodaldata1(j)),j=i,k)
           end do
        
         end do
     endif
   endif

   if( allocated( NodalData1   ) ) deallocate ( NodalData1   )
   if( allocated( NodalData2   ) ) deallocate ( NodalData2   )
   if( allocated( NodalData3   ) ) deallocate ( NodalData3   )
   if( allocated( NodalCounter ) ) deallocate ( NodalCounter )
   if( allocated( VTKcells     ) ) deallocate ( VTKcells     )

   !
   ! do not forget to close, the unit number IOpst might be reused
   !
   close(IOpst)
   
   !
   ! VTK particle track plot, added by Johan Jacobs, Simex-Technology
   !
   if( UseParticles )then
     if( Npart >= 1000 )then 
       write(*,*) '+++ Warning: Number of particles must be lower then 1000'
       write(*,*) 'File tracks.vtk not written'
     else
       call openfile(IOpst,'tracks','.vtk','FORMATTED', &
                                      'SEQUENTIAL','UNKNOWN',debug)
       write(IOpst,'(A)')      '# vtk DataFile Version 3.0'
       write(IOpst,'(A,A,A)')   casename(1:lens(casename)),': ',title(1:lens(title)) 
       write(IOpst,'(A)')      'ASCII'
       write(IOpst,'(A)')      'DATASET POLYDATA'

       n_part_points = 0
       do i=1,Npart
         if( .not. associated( Tracks(i)%head ) )then
           write(*,*)'ParticlePrint: error head not associated',i
           cycle
         endif
         n_part_points = n_part_points + Tracks(i)%n      
        !write(*,*)'>>>',i,Tracks(i)%n, n_part_points 
      end do

       write(IOpst,'(A,i8,A)') 'POINTS ',n_part_points,' float'
       izz = 0
       do i=1,Npart
         Track => Tracks(i)%head
         ipn = 0
         do while( associated(Track) )
           ipn = ipn + 1
           if( ipn <= Tracks(i)%n ) write(IOpst,'(3(1x,1pe12.5))') Track%x
           izz = izz + 1
           Track => Track%next
         end do
        !write(*,*)'n>>',i,ipn,Tracks(i)%n
       end do
       write(*,*)'Punten geschreven:',izz
       
       write(IOpst,'(A,i8,i8)') 'LINES ',Npart,n_part_points+Npart
       i_part_point=0
       do i=1,Npart
         write(IOpst,'(i8)') Tracks(i)%n
         do j=1,Tracks(i)%n
           write(IOpst,'(i8)') i_part_point+j-1
         end do
         i_part_point=i_part_point+Tracks(i)%n   
       end do
       write(IOpst,'(A,i8)') 'POINT_DATA', n_part_points
       write(IOpst,'(A,i8)') 'SCALARS part_inr int'
       write(IOpst,'(A,i8)') 'LOOKUP_TABLE default'
       do i=1,Npart
         do j=1,Tracks(i)%n 
           write(IOpst,'(i8)') i
         end do
       end do
       close(IOpst)
     endif
   endif
   !
   ! end addition Johan Jacobs, Simex-Technology
   !
   
   if( Debug > 2 ) write(*,*)'Done'

end subroutine dolfyn2vtk
