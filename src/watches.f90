module watches

   public  :: watch_setup, &
              watch_enter, &
              watch_leave, &
              watch_print, &
              watch_destroy

   private
   
   logical :: UseWatches = .true.
   logical :: WatchesSet = .false.
   
   integer :: CurrentLevel   = 0
   integer :: DefinedWatches = 0
   integer :: MaxWatch       = 0
   integer :: RunningWatch   = 0

   integer, parameter :: NameLength  =   64
   integer, parameter :: AbsMaxWatch = 9999
   
   type watchtype
     character(len=NameLength) :: Name      =  ''
     integer                   :: Level     =  -1
     integer                   :: Parent    =  -1
     real                      :: T         = 0.0
     real                      :: T0        = 0.0
     real                      :: T1        = 0.0
     logical                   :: Running   = .false.
     logical                   :: Active    = .false.
   end type

   type(watchtype), save :: MainWatch

   type(watchtype), allocatable, dimension(:) :: Watch
   
   contains
   
   !========================================================
   real function GetCpuTime()
        
     call cpu_time( GetCpuTime )
   
   end function GetCpuTime
   !========================================================
   subroutine watch_setup(N)

     integer, optional, intent(IN) :: N
     
     if( .not. UseWatches ) return

     if( present(N) .and. N <= 0 )then
       UseWatches = .false.
       return
     endif

     if( .not. present(N) )then
       MaxWatch = 20
     else
       MaxWatch =  N 
     endif

     allocate(Watch(N),stat=istat)   
     WatchesSet     = .true.
     CurrentLevel   =    1
     DefinedWatches = DefinedWatches + 1
     
     MainWatch%Name    = '- - main - -'
     MainWatch%T0      =  GetCpuTime()
     MainWatch%T1      =  GetCpuTime()
     MainWatch%T       =   0.0

     Watch(1)%Name    = 'Master'
     Watch(1)%Level   =    1
     Watch(1)%Parent  =    0
     Watch(1)%T0      =  GetCpuTime()
     Watch(1)%T1      =  GetCpuTime()
     Watch(1)%T       =   0.0
     Watch(1)%Active  =  .true.
     Watch(1)%Running =  .true.

     RunningWatch     = 1
     
   end subroutine watch_setup
   !========================================================
   subroutine watch_destroy

     if( .not. UseWatches ) return

     if( allocated(Watch) ) deallocate(Watch)

     MaxWatch       =  0
     WatchesSet     = .false.
     DefinedWatches =  0
     CurrentLevel   = -1
     RunningWatch   = -1
     
   end subroutine watch_destroy
   !========================================================
   subroutine watch_enter(NameIn)

     character(len=*), optional, intent(IN) :: NameIn
     
     character(len=NameLength) :: Name, String
     logical :: NameFound, WatchFound
     
     if( .not. UseWatches ) return
     if( DefinedWatches+1 > MaxWatch ) return
 
     if( .not. present(NameIn) )then
       if( DefinedWatches+1 <= MaxWatch )then
         write(String,'(''Watch_'',i4.4)') DefinedWatches+1
         Name = String
       else
         return
       endif
     else
       Name = NameIn
     endif

     NameFound  = .false.
     WatchFound = .false.
     do i=1,DefinedWatches
       if( Watch(i)%Name == Name )then
         NameFound = .true.
         do j=1,DefinedWatches
           if( Watch(i)%Level == CurrentLevel+1 )then
             ! already known watch
             WatchFound = .true.
             iWatch = i
             exit
           endif
         end do
         if( WatchFound ) exit
       endif  
     end do

     !
     ! halt running watch
     !      
     Watch(RunningWatch)%T1 = GetCpuTime()
     Watch(RunningWatch)%T  = Watch(RunningWatch)%T  + &
                            ( Watch(RunningWatch)%T1 - &
                              Watch(RunningWatch)%T0 )
     Watch(RunningWatch)%Running = .false.
     
!write(*,*)'E stop',RunningWatch,Watch(RunningWatch)%Name
     !
     ! next level
     !
     CurrentLevel = CurrentLevel + 1

     !
     ! start watch
     !
     if( WatchFound )then

       Watch(iWatch)%Name   = Name
       Watch(iWatch)%Level  = CurrentLevel
       Watch(iWatch)%Parent = RunningWatch
       Watch(iWatch)%T0     = GetCpuTime()
       Watch(iWatch)%Active = .true.
      
       RunningWatch = iWatch
       
     else
       
       DefinedWatches = DefinedWatches+1
       Watch(DefinedWatches)%Name     = Name
       Watch(DefinedWatches)%Level    = CurrentLevel
       Watch(DefinedWatches)%Parent   = RunningWatch
       Watch(DefinedWatches)%T0       = GetCpuTime()
       Watch(DefinedWatches)%Active = .true.

       RunningWatch = DefinedWatches

     endif

!write(*,*)'E star',RunningWatch,Watch(RunningWatch)%Name

   end subroutine watch_enter
   !========================================================
   subroutine watch_leave(Name)

     character(len=*), optional, intent(IN) :: Name
     logical :: WatchFound

     if( .not. UseWatches ) return

     !
     ! halt running watch
     !
     Watch(RunningWatch)%T1 = GetCpuTime()
     Watch(RunningWatch)%T  = Watch(RunningWatch)%T  + &
                            ( Watch(RunningWatch)%T1 - &
                              Watch(RunningWatch)%T0 )

     Watch(RunningWatch)%Active  = .false.
     Watch(RunningWatch)%Running = .false.

!write(*,*)'L stop',RunningWatch,Watch(RunningWatch)%Name
     !
     ! restart previous watch
     !    
     iWatch = Watch(RunningWatch)%Parent
     if( iWatch > 0 )then
       CurrentLevel = CurrentLevel - 1

       Watch(iWatch)%T0      = GetCpuTime()
       Watch(iWatch)%Active  = .true.
       Watch(iWatch)%Running = .true.

       RunningWatch = iWatch
     endif
!write(*,*)'L star',iWatch,Watch(iWatch)%Name
     
   end subroutine watch_leave
   !========================================================
   subroutine watch_print(IOin) 

     integer, optional, intent(IN) :: IOin
     
     integer, parameter :: MaxStar = 7
     character(len=MaxStar)   ::  Star 
     character(len=MaxStar+1) ::  Blank
     logical, save :: Done = .false.
     
     integer, allocatable, dimension(:) :: list
     logical, allocatable, dimension(:) :: mask
     
     if( .not. Done )then
       do i=1,MaxStar
         Star(i:i)  = '*'
         Blank(i:i) = ' '
       end do
       Blank(MaxStar+1:MaxStar+1) = ' '
       Done = .true.
     endif

     if( .not. UseWatches ) return

    !if( .not. present(IOin) )then
    !  IO = 6
    !else
    !  if( IO <= 0 )then
    !    IO = 6
    !  else
         IO = IOin
    !  endif
    !endif

     allocate(list(1:DefinedWatches),stat=istat)
     allocate(mask(1:DefinedWatches),stat=istat)
     list = -1
     mask = .true.

     do i=1,DefinedWatches
       j = maxloc(Watch(1:DefinedWatches)%T,1,MASK=mask)
       list(i) = j
       mask(j) = .false.
     end do

     MainWatch%T1 = GetCpuTime()
     MainWatch%T  = MainWatch%T1-MainWatch%T0

     fact = 1.0
     if( MainWatch%T > 0.0 ) fact = 1.0/MainWatch%T*100.0
     
     sum1 = 0.0
     sum2 = 0.0
     
     write(IO,*)'Watches     : ',DefinedWatches,'>',IO
     write(IO,*)'Elapsed time: ',MainWatch%T
     write(IO,*) &
     ' W L RA         :                             cpu       %'
     do iw=1,DefinedWatches
     
       i    = list(iw)
       
       sum1 = sum1 + Watch(i)%T
       sum2 = sum2 + fact*Watch(i)%T
       
       j   = min(Watch(i)%Level,MaxStar)
       
       write(IO,1) i,Watch(i)%Level, &
                   Watch(i)%Running,Watch(i)%Active, &
                   Star(1:j),Blank(j+1:MaxStar+1),   &
                   Watch(i)%Name(1:24),Watch(i)%T,fact*Watch(i)%T 
     
     end do
           !      lvl   rn,ac    stars    naam     tijd
   1 format(1x,i2,i2,1x,L1,L1,1x,A,A,': ',A24,1x,1pe9.3,1x,0pf6.2)
   
     write(IO,'(44x,''---------  -----'')')  
     write(IO,'(44x,1pe9.3,1x,0pf6.2)') sum1, sum2
   
     if( allocated(list) ) deallocate(list)
     if( allocated(mask) ) deallocate(mask)

   end subroutine watch_print
   
end module
