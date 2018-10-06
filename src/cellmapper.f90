!
! Copyright 2008-2009 Henk Krus, Cyclone Fluid Dynamics BV
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
! old cell: center and faces with vertices (x,y,z)
! new cell: center only needed
!
! map the center of new cell to old cell; keep this map 
!
subroutine mapper


   logical :: BinaryGeometry = .false.

   character(len=32) :: NameOldGeometryFile
   character(len=32) :: NameOldRestartFile
   

   NameOldGeometryFile = 'mand.old.geo'
   NameOldRestartFile  = 'mand.old.rst'
   
   NameOldGeometryFile = 'mand.old'
   NameOldRestartFile  = 'mand.old.rst'
   !
   ! inquire for the type of dolfyn geometry file
   !
   call feelfile(IOoldgeo,IFORM,NameOldGeometryFile,&
                         '.geo','FORMATTED','SEQUENTIAL','OLD',debug)

   if( IFORM == 0 )then
     BinaryGeometry = .true.
   else
     BinaryGeometry = .false.
   endif

   !
   ! read dolfyn old geomtery file
   !
   write(IOdef,*) 'Opening old geometry file'
   write(IOdbg,*) 'Opening old geometry file'
   write(IOrun,*) 'Opening old geometry file'

   if( BinaryGeometry )then
     call openfile(IOoldgeo,NameOldGeometryFile, &
                           '.geo','UNFORMATTED','SEQUENTIAL','OLD',debug)
   else
     call openfile(IOoldgeo,NameOldGeometryFile, &
                           '.geo','FORMATTED','SEQUENTIAL','OLD',debug)
   endif

   if( BinaryGeometry ) then
     !
     ! binary
     !
     read(IOoldgeo) string(1:12)
     if( string(1:12) /= 'dolfyn bin g' )then
       write(IOdef,*)'*** Error: Not a dolfyn geometry file'
       stop
     endif
     read(IOoldgeo) geoversion
     if( geoversion /= version )then
       write(IOdef,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOdbg,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOrun,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
     endif
     read(IOoldgeo) ScaleFactor
     write(IOdef,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOdbg,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOrun,*) 'Using ScaleFactor = ',ScaleFactor

     read(IOoldgeo) string(1:4)
     if( string(1:3) /= 'ce:' )then
       write(IOdef,*)'*** Error: Cell definitions expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     read(IOoldgeo) icel
     write(IOdef,*)'Reading cells ',icel,' BINARY'
     write(IOdbg,*)'Reading cells ',icel,' BINARY'
     write(IOrun,*)'Reading cells ',icel,' BINARY'

     Noldcel = icel

     do i=1,Noldcel
       read(IOoldgeo) j
       if( j /= i ) write(IOdef,*)'*** Warning: cell array corrupted'
     end do
   else
     !
     ! ascii
     !
     read(IOoldgeo,'(a12)') string(1:12)
     if( string(1:12) /= 'dolfyn asc g' )then
       write(IOdef,*)'*** Error: Not a dolfyn geometry file'
       stop
     endif
     read(IOoldgeo,'(8x,i6)') geoversion
     if( geoversion /= version )then
       write(IOdef,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOdbg,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
       write(IOrun,'(1x,a,f6.3)') &
         '*** Warning: Geometry file of version: ',geoversion*.001
     endif
     read(IOoldgeo,'(6x,e12.5)') ScaleFactor
     write(IOdef,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOdbg,*) 'Using ScaleFactor = ',ScaleFactor
     write(IOrun,*) 'Using ScaleFactor = ',ScaleFactor

     read(IOoldgeo,'(a3)') string
     if( string(1:3) /= 'ce:' )then
       write(IOdef,*)'*** Error: Cell definitions expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     read(IOoldgeo,*) icel
     write(IOdef,*)'Reading cells ',icel
     write(IOdbg,*)'Reading cells ',icel
     write(IOrun,*)'Reading cells ',icel

     Noldcel = icel

     do i=1,Noldcel
       read(IOoldgeo,'(2(i8,1x),4(1pe16.9,1x))') j
       if( j /= i ) write(IOdef,*)'*** Warning: cell array corrupted'
     end do
   endif

   !
   ! list of cell faces
   !
   if( BinaryGeometry )then
     read(IOoldgeo) string(1:4)
     if( string(1:4) /= 'cf: ' )then
       write(IOdef,*)'*** Error: List of cell faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading list of cell faces'

     read(IOoldgeo) maxfaces
     if( maxfaces > 6 )then
       write(IOdef,*)'*** Warning: Unsupported feature'
     endif
     allocate(NOldFaces(Noldcel),stat=istat)
     call TrackMemory(istat,Ncel,'NOldFaces array allocated')

     allocate(COldFace(Ncel,12),stat=istat)
     call TrackMemory(istat,Ncel*12,'COldFace array allocated')

     do i=1,Noldcel
       read(IOoldgeo) j,NOldFaces(j)
       do k=1,Nfaces(j)
         read(IOgeo) COldFace(j,k)
       end do
       if( j /= i ) write(IOdef,*)'*** Warning: cell face array corrupted'
     end do

   else
     !
     ! ascii
     !
     read(IOoldgeo,'(a3)') string
     if( string(1:3) /= 'cf:' )then
       write(IOdef,*)'*** Error: List of cell faces expected'
       write(IOdef,*)'    Found: ',string(1:3)
     endif
     write(IOdbg,*)'Reading list of cell faces'

     read(IOoldgeo,*) maxfaces
     if( maxfaces > 6 )then
       write(IOdef,*)'*** Warning: Unsupported feature'
     endif
     allocate(NOldFaces(Ncel),stat=istat)
     call TrackMemory(istat,Noldcel,'NFaces array allocated')

     allocate(COldFace(Ncel,12),stat=istat)
     call TrackMemory(istat,Noldcel*12,'COldFace array allocated')

     do i=1,Ncel
       read(IOoldgeo,'(12(i8,1x))') j,NOldFaces(j),(COldFace(j,k),k=1,NOldFaces(j))
       if( j /= i ) write(IOdef,*)'*** Warning: cell face array corrupted'
     end do
   endif



end

