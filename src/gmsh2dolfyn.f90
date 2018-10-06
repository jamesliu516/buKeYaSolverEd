! ==================== 
!  README Gmsh2Dolfyn 
! ==================== 
! 
! Original version by Navier, July 2007
! Adapted to Gmsh 2.3.1, March 2009
! Adapted to Gmsh 2.4.2, October 2009
!
! Purpose 
! -------
! This Fortran95 program translates a mesh file from Gmsh (.msh) format 
! to Dolfyn (.vrt, .cel, .bnd) format. 
! 
! Input and Output 
! ---------------- 
! Input : A Gmsh .msh file (version 2.0, ascii format). 
! Output: Dolfyn .vrt, .cel and .bnd files. 
! 
! Running the Program 
!-------------------- 
! First compile it using a Fortran95 compiler (eg g95 or gfortran).  
! Run it from the command line. 
! The program prompts for the name of the input file. 
! 
! Bug reports 
! ----------- 
! Please report bugs to http://www.dolfyn.net/dolfyn/forum/
! 
! Important note about the Gmsh msh format and Physical Groups. 
! ------------------------------------------------------------- 
! In order to define boundary conditions, the Gmsh geometry-builder allows a 
! group of faces to be assigned a common 'physical group' label.  The mesh 
! inherits this label, and the label is used in the Dolfyn .din file. 
! 
! When saving the mesh, the default is to save only mesh elements with a 
! physical group label.  This means that some mesh elements will be missing, 
! unless every mesh element belongs to a physical group.  
!
! For example in the adapted gmsh tutorial t2.geo enter:
!
! Physical Volume ("Fluid") = {119,120};
! Physical Surface("Inlet") = {111};
! Physical Surface("Outlet") = {132};
!
! Mesh 3D and save it as t2.msh
!
!======================================================================== 
!======================================================================== 
!======================================================================== 
integer function lens(string) 
 
   character(len=*) string 
     
   do i=len(string),0,-1 
     if( string(i:i) .ne. ' ') goto 10 
   end do 
   i = 0 
10 continue 
 
   lens = i 
     
end function lens 
!====================================================================== 
subroutine openfile(iunit,casename,extension,reqform,status,idebug) 
 
   character(len=*)  casename 
   character(len=*)  extension 
   character(len=*)  reqform 
   character(len=*)  status 
   character(len=48) filename 
   character(len=11) form 
 
   logical exists 
 
   filename = casename(1:lens(casename))//extension(1:lens(extension)) 
   length   = lens(filename) 
 
   if( idebug > 2 )write(*,*) 'Opening ',filename(1:length) 
   
   if( status(1:3) == 'OLD' )then 
     inquire(file=filename(1:length),exist=exists,form=form) 
     if( .not. exists )then 
       write(*,*) '*** Error: File ',filename(1:length),' does not exist' 
       stop 
     endif 
   endif 
 
   open(iunit,file=filename(1:length),form=reqform,status=status) 
 
   if( idebug >= 2 ) write(*,*) 'File ',filename(1:length),' opened' 
   
end subroutine openfile 
!============================================================================== 
program gmsh2dolfyn 

   implicit none 

   integer, parameter :: IOinp = 13, IOcel = 14, IOvrt = 15, IObnd = 23  ! I/O file numbers 
   integer, parameter :: IOdbg = 63, IOcfg = 12    
   integer, parameter :: IOgmsh= 24                                      !Gmsh mesh file 
 
   integer :: debug = 0 
 
   integer, parameter :: version = 0530 

   character(len=128) :: line 

   integer, parameter :: MaxNames = 100
   character(len=64), dimension(MaxNames) :: Names 
   character(len=64), dimension(MaxNames) :: Regions  
   integer, dimension(MaxNames)           :: ICTID    = -1
   integer, dimension(MaxNames)           :: Partition=  1

   logical, dimension(MaxNames)           :: Fluid    = .false.
   logical, dimension(MaxNames)           :: Boundary = .false.

   character(len=64)  :: casename = 'dolfyn' 
   character(len=72)  :: c_input1, c_input2, c_input3 

   integer i, j, k, ie, icel, ibnd, iloop 

   !
   ! nodes/vertices 
   !
   integer n_nodes,inode 
   real  node(3) 

   !
   ! there are 19 gmsh element types: \
   !
   !    1 : 2-node line 
   !    2 : 3-node triangle (face) 
   !    3 : 4-node quadrangle (face) 
   !    4 : 4-node tetrahedron 
   !    5 : 8-node hexahedron (eg cube) 
   !    6 : 6-node triangular-prism 
   !    7 : 5-node pyramid 
   !
   !  8-14: 'second-order' elements.  Ref Gmsh manual. 
   !   15 : 1-node point
   ! 16-19: more second-order FEM elements
   !
   ! the nodes/vertices for each element are read into the 
   ! v array. 
   !
   ! each element can have several tags.   
   ! the first tag gives the physical_group number.  
   ! all elements on the same boundary have the same physical_group number. 
   !
   integer, parameter :: element_type(19) = & 
                      (/ 2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13 /) 
			
   integer :: n_elements, ielement, ielement_type, n_tags, n_names, lens
   integer :: tags(64), v(27) 
   integer :: i3, i4q, i4, i5, i6, i8
   integer :: ivs = 0
   
   if( size(v) /= maxval(element_type) )then      
     stop'bug: error in dimensions of array v' 
   endif
   
   !
   ! read the gmsh filename, then open the .msh file 
   !
   write(*,*) 'Gmsh2Dolfyn: Converts a Gmsh mesh file to Dolfyn format.' 
   write(*,*) '(Input must be in Gmsh version 2.0 ascii format.' 
   write(*,*) ' Output is in Dolfyn version Jan 2009 format.)' 
   write(*,*) ' ' 
   write(*,*) 'Input Gmsh filename, excluding the .msh suffix' 
   read(*,'(A)') casename 

   write(*,*) 'Opening the Gmsh file' 
   call openfile(IOgmsh,casename,'.msh','FORMATTED','OLD',debug) 

   !
   ! read the Gmsh file header 
   !
   write(*,*)'Reading MeshFormat'
   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$MeshFormat') 
   
   read(IOgmsh,*) c_input1,c_input2,c_input3  
   if( c_input1 == '2.1' )then
     ivs = 21
   else if( c_input1 == '2' )then
     ivs = 20
   else
     write(*,*) '*** WARNING: unknown Gmsh version'
     write(*,*) '*** Unexpected results might happen'
     ivs = 21
   endif
   
   if( ivs == 20 )then
     call check_input_character(c_input1,'2') 
   else if( ivs == 21 )then
     call check_input_character(c_input1,'2.1') 
   else
     write(*,*) '*** Version found ',c_input1 
   endif
   
   call check_input_character(c_input2,'0') 
   call check_input_character(c_input3,'8') 
   
   write(*,*) 'MeshFormat: ', c_input1(1:lens(c_input1)),' ', &
                              c_input2(1:lens(c_input2)),' ', &
                              c_input3(1:lens(c_input3)) 

   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$EndMeshFormat') 

   !
   ! read the Gmsh PhysicalNames 
   !
   write(*,*)'Reading PhysicalNames'
   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$PhysicalNames') 

   read(IOgmsh,*) n_names 
   if( n_names <= 0 )then
     write(*,*) 'error: number of names must be a positive number' 
     stop
   endif

   if( ivs == 20 )then
     do i=1,n_names 
      !read(IOgmsh,*) k,j,c_input1
       read(IOgmsh,*)   j,c_input1
       write(*,*) 'Name ',j,'-> ', c_input1      
       Names(j) = c_input1  
     end do
   else if( ivs == 21 )then
     do i=1,n_names 
       read(IOgmsh,*) k,j,c_input1
       write(*,*) 'Name ',j,'-> ', c_input1      
       Names(j) = c_input1  
     end do
   else
     do i=1,n_names 
       read(IOgmsh,*) k,j,c_input1
       write(*,*) 'Name ',j,'-> ', c_input1      
       Names(j) = c_input1  
     end do   
   endif
   
   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$EndPhysicalNames') 
   !
   ! read the nodes from the .msh file and write them 
   ! to the .vrt file. 
   !
   write(*,*)'Reading Nodes'
   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$Nodes') 
   
   read(IOgmsh,*) n_nodes 
   if( n_nodes <= 0 )then
     write(*,*) 'error: number of nodes must be a positive number' 
     stop
   endif

   !
   ! open the dolfyn .vrt file 
   !
   write(*,*) 'Creating the dolfyn .vrt file' 
   call openfile(IOvrt,casename,'.vrt','FORMATTED','UNKNOWN',debug) 

   nodes: do iloop=1,n_nodes 
     read(IOgmsh,*) inode,(node(i), i=1,3) 
     write(IOvrt,'(i9,6x,3g16.9)') inode,(node(i),i=1,3) 
   enddo nodes 

   write(*,*) 'Nodes written ',n_nodes
   !
   ! close the dolfyn .vrt file 
   !
   close(IOvrt) 

   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$EndNodes') 

   !
   ! read the elements from the .msh file and write them 
   ! to the .cel and .bnd files. 
   !
   read(IOgmsh,*) c_input1  
   call check_input_character(c_input1,'$Elements') 
   
   read(IOgmsh,*) n_elements 
   if( n_elements <= 0 )then
     write(*,*) 'error: number of elements must be a positive number' 
     stop
   endif

   write(*,*) 'Total Gmsh elements to be read in:',n_elements 
   !
   ! open the dolfyn .cel and .bnd files 
   !
   write(*,*) 'Creating the dolfyn .cel and .bnd files' 
   call openfile(IOcel,casename,'.cel','FORMATTED','UNKNOWN',debug) 
   call openfile(IObnd,casename,'.bnd','FORMATTED','UNKNOWN',debug) 

   !
   ! note in Gmsh fluid cells and boudaries can be mixed
   ! we just keep track on them both
   ! remind default region is not assigned 
   !
   icel = 0
   ibnd = 0
   
   i3   = 0
   i4q  = 0
   i4   = 0
   i5   = 0
   i6   = 0
   i8   = 0

   do ie=1,n_elements 
   
     read(IOgmsh,*) ielement, ielement_type, n_tags      
     if( n_tags /= 3 ) write(*,*) 'tag error n_tags /= 3:',ielement,n_tags
     call check_element_type(ielement_type,element_type) 
     call check_n_tags(n_tags,tags) 
     backspace(IOgmsh) 

     !
     ! we need to circumvent backspace but
     ! advance='no' requires fixed format
     ! just keep it for now.
     !

     !
     ! now we know what to to expect to find on the line
     !
     read(IOgmsh,*) ielement, ielement_type, &
     		    n_tags, (tags(i),i=1,n_tags),& 
     		   (v(i),i=1,element_type(ielement_type)) 
      
     if( 4 <= ielement_type .and. ielement_type <= 7 )then
       
       icel = icel + 1

       if( .not. Fluid(tags(1)) ) Fluid(tags(1)) = .true.
       if( Boundary(tags(1)) )then
         write(*,*) 'Inconsistent data: Physical names ids overlap 1' 
       endif 

     1 format(i8,8(1x,i8),2(1x,i4))
       select case(ielement_type)
       
         case(4) ! 4-node tet

           write(IOcel,1) icel,& 
     	     	v(1),v(2),v(3),v(3),  v(4),v(4),v(4),v(4),& 
     	     	tags(1),tags(3)
           i4 = i4 + 1

         case(5) ! 8-node hex

           write(IOcel,1) icel,& 
     	     	v(1),v(2),v(3),v(4),  v(5),v(6),v(7),v(8),& 
     	     	tags(1),tags(3)
           i8 = i8 + 1

         case(6) ! 6-node prism or wedge

           write(IOcel,1) icel,& 
     	     	v(1),v(2),v(3),v(3),  v(4),v(5),v(6),v(6),& 
     	     	tags(1),tags(3)
           i6 = i6 + 1

         case(7) ! 5-node pyramid

           write(IOcel,1) icel,& 
     	     	v(1),v(2),v(3),v(3),  v(5),v(5),v(5),v(5),& 
     	     	tags(1),tags(3)
           i5 = i5 + 1
       
         case default 
	 
	   write(*,*)'internal error 1'

       end select 
       
     elseif( ielement_type == 2 .or. ielement_type == 3 )then
     
       ibnd = ibnd + 1
       
       if( .not. Boundary(tags(1)) ) Boundary(tags(1)) = .true.
       if( Fluid(tags(1)) )then
         write(*,*) 'Inconsistent data: Physical names ids overlap 2' 
       endif 

       select case(ielement_type)
       
         case(2) ! 3-node tri

           write(IObnd,1) ibnd, v(1),v(2),v(3),v(3), tags(1) 
           i3 = i3 + 1

         case(3) ! 4-node quad

           write(IObnd,1) ibnd, v(1),v(2),v(3),v(4), tags(1) 
           i4q = i4q + 1

         case default 
	 
	   write(*,*)'internal error 2'

       end select 

     else

       write(*,*)'internal error 3'
     
     endif
     

   end do

   close(IOcel) 
   close(IObnd) 

   read(IOgmsh,*) c_input1 
   call check_input_character(c_input1,'$EndElements') 

   if( i3  > 0 ) write(*,*) 'Triangle boundaries: ',i3
   if( i4q > 0 ) write(*,*) 'Quad boundaries:     ',i4q
   if( i4  > 0 ) write(*,*) 'Tetrahedral cells:   ',i4
   if( i5  > 0 ) write(*,*) 'Pyramid cells:       ',i5
   if( i6  > 0 ) write(*,*) 'Prism cells:         ',i6
   if( i8  > 0 ) write(*,*) 'Hexahedral cells:    ',i8

   !
   ! finally write out the boundary names
   !
   write(*,*) 'Writing the .inp file' 
   call openfile(IOinp,casename,'.inp','FORMATTED','UNKNOWN',debug) 

   do i=1,MaxNames
     if( Boundary(i) )then
       if( i <= 9 )then
         write(IOinp,'(''rname,'',i1,'','',A64)') i,Names(i)
       elseif( i <= 99 )then
         write(IOinp,'(''rname,'',i2,'','',A64)') i,Names(i)
       elseif( i <= 999 )then
         write(IOinp,'(''rname,'',i3,'','',A64)') i,Names(i)
       elseif( i <= 9999 )then
         write(IOinp,'(''rname,'',i4,'','',A64)') i,Names(i)
       elseif( i <= 99999 )then
         write(IOinp,'(''rname,'',i5,'','',A64)') i,Names(i)
       else
         write(IOinp,'(''rname,'',i6,'','',A64)') i,Names(i)      
       endif
     endif
   end do
   
   close(IOinp)
   
   write(*,*) 'Done gmsh2dolfyn'

   contains      
   !------------------------------------------------------------------------------------ 
   subroutine check_input_character(c1,c2) 

     implicit none 

     character (len=*) :: c1, c2 

     if( c1(1:len(c2)) /= c2 )then 
       write(*,*)  'error reading Gmsh input file: ',& 
                   'the following two characters should be the ',&
		   'same but differ ',c1(1:len(c2)),c2 
       stop 
     endif 
    
   end subroutine 
   subroutine check_element_type(ielement_type,element_type) 

     implicit none 
     integer ielement_type 
     integer element_type(:) 

     if( ielement_type < 0 )then 
       write(*,*) 'error reading Gmsh file: element type must be positive' 
       write(*,*) 'element type = ',ielement_type 
       stop 
     endif  

     if( ielement_type > size(element_type) )then 
       write(*,*) 'error reading Gmsh file: unrecognised element type' 
       write(*,*) 'element type ',ielement_type 
       write(*,*) 'max recognised element type ',size(element_type) 
       stop
     endif 
    
   end subroutine 
   subroutine check_n_tags(ntags,itags) 

     implicit none 

     integer ntags 
     integer itags(:) 

     if( ntags > size(itags) )then 
       write(*,*) 'error: The Gmsh file contains ',ntags,' tags per element' 
       write(*,*) 'Gmsh2Dolfyn is hard-wired for a maximum of ',size(itags),& 
        	  'tags.  The dimension of this array needs to be increased.' 
       stop         
     endif 
      
   end subroutine 
end 
 
