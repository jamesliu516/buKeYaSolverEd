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
! With contributions by B. Tuinstra
! (www.home.zonnet.nl/bouke_1/dolfyn)
!
module constants
   !
   ! global constants
   !
   integer, parameter :: version = 0527

   integer, parameter :: Solid = 1, Fluid = 2, Ghost = 3, Baffle = 4, Shell = 5

   integer, parameter :: Hexa = 1, Prism = 2, Pyramid  = 3, Tetra = 5,     &
                         Poly = 6, Quad  = 7, Triangle = 8, Polygon = 9

   integer, parameter :: NFhex  = 6, NFprism = 5, NFpyrmd = 5, NFtet = 4,  &
                         NFpoly = 8, NFquad  = 2, NFtri = 2

   integer, parameter :: IOgeo =  8, IOcel = 14, IOvrt = 15, IObnd = 23,   &
                         IOdbg = 30, IOinp = 10, IOpst = 13, IOrst =  9,   &
                         IOmtx = 11, IOres =  7, IOdef =  6, IOcfg = 12,   &
                         IOext = 17, IOmon = 12, IOprt = 18, IOrun = 19,   &
			 IOchk = 20,                                       &
			 IOoldgeo = 21, IOoldrst = 22

   character(len=4)   :: Ext(1:30)
   data Ext(IOcel), Ext(IOvrt), Ext(IObnd) / '.cel', '.vrt', '.bnd' /
   data Ext(IOinp), Ext(IOrst), Ext(IOgeo) / '.din', '.rst', '.dge' /  !<== .geo or .dge
   !
   ! boundary regions
   !
   integer, parameter :: RWall = 1, RInlet = 2,  ROutlet = 3, RPressure = 4,  &
                         RSymp = 5, RThin = 6,   RCyclic = 7, RDomain = 8

   integer, parameter :: NRegionTypes = 7
   character(len=12)  :: Region(1:NRegionTypes)

                 !123456789 12'
   data Region / 'Wall        ', &
                 'Inlet       ', &
		 'Outlet      ', &
		 'Pressure    ', &
                 'Symm. plane ', &
                 'Cyclic      ', &		
                 'Dom.Boundary'  /

   !
   ! turbulence models
   !
   integer, parameter :: TMnone = 0,                        &  
                         TMkeps = 1, TMrng = 2, TMcube = 3, & ! k-e models
                         TMrsm  = 4 ,                       & ! RSM models
			 TMles  = 5, TMlesdyn = 6             ! LES models >= 5

   integer            :: TurbModel = TMnone                   ! set default

   character(len=16)  :: TMnames(0:6)
                  !123456789 123456
   data TMnames / '-               ', &
                  'k-epsilon std   ', &
                  'k-epsilon rng   ', &
                  'k-epsilon cubic ', &
                  'Reynolds stress ', &
                  'LES Smag. const.', &
                  'LES Smag. dynam.'  /
   !
   ! standard k-epsilon turbulence model constants
   !
   real      :: Kappa   = 0.419
   real      :: TMCmu   = 0.09
   real      :: TMCeps1 = 1.44
   real      :: TMCeps2 = 1.92
   real      :: TMLenSc = 0.1
   !
   ! gradient methods
   !
   integer, parameter   :: GradLS    = 1
   integer, parameter   :: GradGauss = 2
   
   character(len=2)     :: GradID(1:2)
   data GradID / 'LS', 'GS' /
   
   integer              :: GradAlg    = GradGauss ! Gauss default again
   integer              :: Ngradient  =   2       ! default passes for Gauss' gradient

   integer, allocatable :: GradVar(:)             ! holds gradient type per variable 
   integer, allocatable :: nGradVar(:)            ! passes for Gauss 

   !
   ! some constants
   !
   real                 :: Small = 1.e-12
   real                 :: Large = 1.e+24

   real, parameter      :: Zero = 0.0
   real, parameter      :: One  = 1.0
   real, parameter      :: Two  = 2.0

   !
   ! differencing schemes
   !
   integer, parameter :: DSud  = 1, DScds = 2,            &
                         DSgam = 3, DSlud = 4, DSmod = 5, &
                         DScd1 = 6, DScd2 = 7, DScd3 = 8, &
                         DSlux = 9 
			  
   character(len=3)   :: DSch(9)
   data DSch / 'UD ','CD ','GAM','LUD','MMD','CD1','CD2','CD3','LUX' /

   integer            :: DScdDefault = DScd1
   
   character(len=16)  :: DSnames(9)
                  !123456789 123456
   data DSnames / 'Standard Upwind ', &  ! classic UD
                  'Central Differ. ', &  ! based on lambda
                  'Gamma Differenc.', &  ! gamma scheme
                  'Linear Upwind   ', &  ! lin. upwind based on CBC
                  'MinMod Scheme   ', &  ! 
                  'CDS uncorrected ', &  ! CDS
                  'CDS Gradients   ', &  ! based on gradients
                  'CDS Simple Aver.', &  ! (too) simple average
                  'LUDS uncorrected'  /  ! lin. upwind straight 2nd order
   !
   ! character array with konsole output flags
   !
   integer, parameter    :: NFlags = 5
   character*1, dimension(NFlags) :: Flags

   integer, parameter    :: IFlagP = 1, IFlagInflow = 2, IFlagMass = 3, &
                            IFlagTE = 4, IFlagVis = 5
   !
   ! keywords for variable density, viscosity, Cp
   !
   integer, parameter :: Constant = 1, Ideal = 2, Isobaric = 3, Inviscid = 4, &
                         Multicomponent = 4, Nonnewtonian = 5, User = 6,      &
                         Polynomial = 7
   !
   ! materials
   !
   integer            :: Density = Constant
   integer            :: LaminarViscosity = Constant

   !
   ! symbolic constants for number of variables
   !
   integer, parameter :: NVar  = 13               ! number of standard variables
   integer, parameter :: MinSC = NVar + 1         ! start of scalars
   integer            :: MaxVar, MaxSC, MaxMat    ! no parameter, will be set

   integer, parameter :: VarU   =  1, VarV   =  2, VarW   =  3, VarP   =  4,  &
                         VarTE  =  5, VarED  =  6, VarT   =  7, VarSC  =  8,  &
                         VarDen =  9, VarPP  = 10, VarVis = 11, VarLvis= 12,  &
                         VarCP  = 13

   integer, allocatable :: VarS(:)

   character(len=64) :: casename = 'dolfyn'
   character(len=96) :: title    = 'dolfyn'

   integer   :: memory     =   0
   integer   :: Debug      =   0

   integer   :: NSave      = 500                 ! save every
   integer   :: ISave      =  -1                 ! save just iter i
   real      :: TSave      =  -1.                ! save only at time
   logical   :: TSaveDone  = .false.
   real      :: TCPU       =  -1.                ! save every elapsed  cpu time
   logical   :: CPUStop    = .false.             ! stop or continue

   integer   :: Nouter     =   1                 ! outer step counter
   integer   :: Niter      = 100                 ! number of steps/iterations

   integer   :: MaxOuter   =  40                 ! maximum number of outer steps
   integer   :: NOutput    = 500                 ! write every
   integer   :: IOutput    =  -1                 ! write just iter i
   real      :: TOutput    =  -1.                ! write only at time
   logical   :: TOutputDone= .false.

   logical   :: StoreMesh  = .false.             ! write mesh to rst file 

   logical   :: Transient  = .false.
   real      :: dt         = 0.001
   logical   :: Euler      = .true.
   logical   :: QuadTime   = .false.
   real      :: GammaTime  = 1.0

   integer   :: Restart    = 0

   integer, parameter :: RestartNone = 0
   integer, parameter :: RestartFull = 1
   integer, parameter :: RestartInitial = 2
   integer, parameter :: RestartCounterReset = 3
   

   logical   :: MovingGrid      = .false.

   logical   :: SolveU          = .true.
   logical   :: SolveV          = .true.
   logical   :: SolveW          = .true.
   logical   :: SolveUVW        = .true.

   logical   :: SolveP          = .true.
   logical   :: SolveEnthalpy   = .false.
   logical   :: SolveTurb       = .false.
   logical   :: SolveVisc       = .false.
   logical   :: SolveScalars    = .false.
   logical   :: SolveTurbEnergy = .false.
   logical   :: SolveTurbDiss   = .false.

   logical   :: SolveDensity    = .false.
   logical   :: SolveCP         = .false.

   logical   :: Initialisation     = .false.
   logical   :: UserInitialisation = .false.
   integer   :: InitialSteps       =  10
   real      :: VisInitialFactor   =   1.0

   !
   ! additional limiters
   !
   logical   :: LimitLowU          = .false.
   logical   :: LimitUpU           = .false.
   real      :: LowerLimitU        =  -1.e9
   real      :: UpperLimitU        =   1.e9
   logical   :: LimitLowV          = .false.
   logical   :: LimitUpV           = .false.
   real      :: LowerLimitV        =  -1.e9
   real      :: UpperLimitV        =   1.e9
   logical   :: LimitLowW          = .false.
   logical   :: LimitUpW           = .false.
   real      :: LowerLimitW        =  -1.e9
   real      :: UpperLimitW        =   1.e9

   logical   :: LimitLowVMag       = .false.
   logical   :: LimitUpVMag        = .false.
   real      :: LowerLimitVMag     =    0.0
   real      :: UpperLimitVmag     =   1.e9

   logical   :: LimitLowP          = .false.
   logical   :: LimitUpP           = .false.
   real      :: LowerLimitP        =  -1.e9
   real      :: UpperLimitP        =   1.e9

   logical   :: LimitLowTurbEnergy = .true.
   logical   :: LimitUpTurbEnergy  = .false.
   real      :: LowerLimitTE       =  1.e-9
   real      :: UpperLimitTE       =  1.e-9
   logical   :: LimitLowTurbDiss   = .true.
   logical   :: LimitUpTurbDiss    = .false.
   real      :: LowerLimitED       =  1.e-12
   real      :: UpperLimitED       =  1.e-12

   logical   :: LimitLowEnthalpy   = .false.
   logical   :: LimitUpEnthalpy    = .false.
   real      :: LowerLimitT     =      0.0
   real      :: UpperLimitT     =  10000.0

   logical, allocatable :: LimitLowScalar(:)
   logical, allocatable :: LimitUpScalar(:)
   real, allocatable    :: LowerLimitSC(:)
   real, allocatable    :: UpperLimitSC(:)


   !
   ! slope limiters
   !
   integer, parameter   :: SlopeLimiterOFF =  0 ! off
   integer, parameter   :: SlopeLimiterBJc =  1 ! Barth,Jespersen, cell centroids
   integer, parameter   :: SlopeLimiterBJf =  2 ! Barth,Jespersen, face based
   integer, parameter   :: SlopeLimiterBJn =  3 ! Barth,Jespersen, node based
   integer, parameter   :: SlopeLimiterVNc =  4 ! Venkatarishnan, cell centroid
   integer, parameter   :: SlopeLimiterVNf =  5 ! Venkatarishnan, face based
   integer, parameter   :: SlopeLimiterVNn =  6 ! Venkatarishnan, node based
   integer, parameter   :: SlopeLimiterVAc =  7 ! van Albada, cell centroid
   integer, parameter   :: SlopeLimiterVAf =  8 ! .. face based
   integer, parameter   :: SlopeLimiterVAn =  9 ! .. node based
   integer, parameter   :: SlopeLimiterP1c = 10 ! cell centroid, polynomial vs. 1, see 
   integer, parameter   :: SlopeLimiterP1f = 11 ! Michalak,Ollivier-Gooch, face based
   integer, parameter   :: SlopeLimiterP1n = 12 ! and node based
   
   logical, allocatable :: UseSlopeLimiter(:)            
   real, allocatable    :: URFSlopeLimiter(:)            
   integer, allocatable :: LimiterSlope(:)         

   character(len=3)     :: SlopeID(0:12)
   data SlopeID / 'off', 'BJc','BJf','BJn', 'VNc','VNf','VNn' , &
                         'vAc','vAf','vAn', 'P1c','P1f','P1n'   /

   !
   ! checkout testing stuff (ombouwen tot zoals de slopelimiters)
   !
   logical   :: CheckOut     = .false.

   logical   :: CheckU       = .false.
   logical   :: CheckV       = .false.
   logical   :: CheckW       = .false.
   logical   :: CheckP       = .false.
   logical   :: CheckTE      = .false.
   logical   :: CheckED      = .false.
   logical   :: CheckT       = .false.
   logical   :: CheckDEN     = .false.
   logical   :: CheckVIS     = .false.

   logical   :: CheckReprtU  = .false.
   logical   :: CheckReprtV  = .false.
   logical   :: CheckReprtW  = .false.
   logical   :: CheckReprtP  = .false.
   logical   :: CheckReprtTE = .false.
   logical   :: CheckReprtED = .false.
   logical   :: CheckReprtT  = .false.
   logical   :: CheckReprtDEN= .false.
   logical   :: CheckReprtVIS= .false.

   real      :: CheckLowU    = -1.e24
   real      :: CheckLowV    = -1.e24
   real      :: CheckLowW    = -1.e24
   real      :: CheckLowP    = -1.e24
   real      :: CheckLowTE   = -1.e24
   real      :: CheckLowED   = -1.e24
   real      :: CheckLowT    = -1.e24
   real      :: CheckLowDEN  = -1.e24
   real      :: CheckLowVIS  = -1.e24

   real      :: CheckUpU     =  1.e24
   real      :: CheckUpV     =  1.e24
   real      :: CheckUpW     =  1.e24
   real      :: CheckUpP     =  1.e24
   real      :: CheckUpTE    =  1.e24
   real      :: CheckUpED    =  1.e24
   real      :: CheckUpT     =  1.e24
   real      :: CheckUpDEN   =  1.e24
   real      :: CheckUpVIS   =  1.e24

   real      :: CheckAverageU	    =  0.0
   real      :: CheckAverageV	    =  0.0
   real      :: CheckAverageW	    =  0.0
   real      :: CheckAverageP	    =  0.0
   real      :: CheckAverageTE      =  0.0
   real      :: CheckAverageED      =  0.0
   real      :: CheckAverageT	    =  0.0
   real      :: CheckAverageDEN     =  0.0
   real      :: CheckAverageVIS     =  0.0

   real      :: CheckAverageBandU   = -1.0
   real      :: CheckAverageBandV   = -1.0
   real      :: CheckAverageBandW   = -1.0
   real      :: CheckAverageBandP   = -1.0
   real      :: CheckAverageBandTE  = -1.0
   real      :: CheckAverageBandED  = -1.0
   real      :: CheckAverageBandT   = -1.0
   real      :: CheckAverageBandDEN = -1.0
   real      :: CheckAverageBandVIS = -1.0

   !
   ! data print options
   !
   character(len=64) :: PrintFile  = '-'

   logical   :: PrintCellVar       = .false.
   logical   :: PrintCellVarUser   = .false.

   integer   :: IPrintCellVarStart = 0
   integer   :: IPrintCellVarEnd   = 0
   integer   :: IPrintCellVarInc   = 1

   logical   :: PrintWallVar       = .false.
   logical   :: PrintWallVarUser   = .false.

   integer   :: IPrintWallVarStart = 0
   integer   :: IPrintWallVarEnd   = 0
   integer   :: IPrintWallVarInc   = 1
   !
   ! general options
   !
   logical   :: UsePatches         = .false.
   logical   :: UseScalars         = .false.
   logical   :: UseParticles       = .false.
   logical   :: UseSensors         = .false.

   logical   :: UseArtificialComp  = .false.
   logical   :: UseGMV             = .false.

   logical   :: UseVTK             = .false.
   logical   :: UseVTKwalls        = .true.
   logical   :: UseVTKbinary       = .true.

   logical   :: UseGMSH            = .false.
   logical   :: UseGMSHwalls       = .true.
   logical   :: UseGMSHbinary      = .false.
   logical   :: UseGMSHdataonly    = .true.

   logical   :: UseTECPLOT         = .false.
   logical   :: UseFIELDVIEW       = .false.
   logical   :: UseENSIGHT         = .false.
   logical   :: UseLapack          = .false.
   logical   :: UseOpenDX          = .true.
   !
   ! special ABL options
   !
   logical   :: UseFixABL          = .false.
   integer   :: IdFixABL           = 0
   real      :: UFixABL            = 0.0
   real      :: VFixABL            = 0.0
   real      :: WFixABL            = 0.0
   real      :: TeFixABL           = 0.0
   real      :: EdFixABL           = 0.0

  !integer, parameter :: SolverSparseKit2=  1  ! Sparsekit2 BiCGstab by Youcef Saad
   integer, parameter :: SparseKit2      =  1  ! Sparsekit2 BiCGstab by Youcef Saad
   integer, parameter :: SolverHypre     =  2  ! Hypre
   integer, parameter :: SolverBiCGstabL =  3  ! van der Vorst BiCGstab Ell
   integer, parameter :: SolverGMRESR    =  4  ! van der Vorst GMRESr
   integer, parameter :: SolverICCG      =  5  !
   integer, parameter :: SolverDirect    =  6  ! Harwell HSL MA27 (testing only)
   integer, parameter :: SolverSuperLU   =  7  ! superLU
   integer, parameter :: SolverPARMS     =  8  ! pARMS
   integer, parameter :: SolverMUMPS     =  9  ! the MUMPS suite
   integer, parameter :: SolverUser      = 10  ! easy way to allow for your own solver

   !
   ! hypre options (Thomas' work comes back in)
   !
   integer, parameter :: HypreBiCGstab   =  1  ! Hypre solver BiCGstab
   integer, parameter :: HypreGMRES      =  2  ! Hypre solver GMRES
   integer, parameter :: HyprePCG        =  3  ! Hypre CG solver 
   
   integer, parameter :: HypreBoomerAMG  =  4  ! BoomerAMG preconditioner
   integer, parameter :: HyprePilut      =  5  ! Pilut preconditioner
   integer, parameter :: HypreParasails  =  6  ! Parasails preconditioner

   integer            :: HypreSolver = HypreGMRES
   integer            :: HyprePreCon = HyprePilut

   !
   ! rarely changed switches/constants
   !
   real      :: ParticleCourant = 0.35         ! default Courant number
			
end module constants
!========================================================================
module geometry
   !
   ! the basic geometry related features
   ! the rest of the program should rely on this geometry related
   ! module only for its data
   !
   use constants

   public
   !
   ! general stuff
   !
   integer :: Nvrt = 0, Ncel = 0, Nbnd = 0, Nreg  = 0, &
              Nfac = 0, Nint = 0, NNZ = 0
   !
   ! vertex section
   !
   real, allocatable :: Vert(:,:)   ! array of vertices

   !
   ! scalefactor
   !
   real :: ScaleFactor = 1.0

   !
   ! face section
   !
   integer, parameter :: MaxFaceNodes = 4

   type FaceData
     integer bnd                 ! internal (0), boundary (#) ...
     integer cell1               ! first cell
     integer cell2               ! second cell
     integer vertices(4)         ! still max 4 vertices; extension needed
     real    area                ! area
     real    n(3)                ! normal / surface area components
     real    x(3)                ! center coordinates
     real    lambda              ! interpolation factor cell1 -> cell2
     real    rlencos             ! impl. coeff.: area/|Xpn|/vect_cosangle(n,Xpn)
     real    xnac(3)             ! auxiliary vectors
     real    xpac(3)             ! 
   end type

   type(FaceData), allocatable   :: Face(:)

   real, allocatable :: FaceNormal(:,:) ! Face normals

   !
   !--- for each computational cell there will be a list of faces
   !
   type StoredList
     integer, allocatable, dimension(:) :: i    ! index of item i (items + 1)
     integer, allocatable, dimension(:) :: list ! list for item i
   end type

   integer, allocatable :: NFaces(:)     ! number of faces for a computational cell
   integer, allocatable :: NNodes(:)     ! number of nodes for a computational cell
   
   integer, allocatable :: CFace(:,:)    ! list of faces (old style)

   type(StoredList)     :: iFaces        ! list of faces (new compressed style)
   type(StoredList)     :: iNodes        ! list of nodes (new compressed style)
   
   real, allocatable    :: RFace(:,:)    ! Ae, Aw coef. list (face -> cell)

   real, allocatable    :: Ar(:)         ! reciprocal coef. of A (pres. cor.)

   real, allocatable    :: Au(:)         ! coef. for the cell (Au)
   real, allocatable    :: Su(:)         ! cell source term (Su)

   real, allocatable    :: Av(:)         ! coef. for the cell (Av)
   real, allocatable    :: Sv(:)         ! cell source term (Sv)

   real, allocatable    :: Aw(:)         ! coef. for the cell (Aw)
   real, allocatable    :: Sw(:)         ! cell source term (Sw)

   !
   ! SparsKit2 / Hypre solver work arrays
   !
   double precision, allocatable :: Work(:,:)

   double precision, allocatable :: Acoo(:)   ! A-matrix in COO-format
   integer, allocatable          :: Arow(:)   ! A-matrix row entries
   integer, allocatable          :: Acol(:)   ! A-matrix column entries

   double precision, allocatable :: Acsr(:)   ! A-matrix in CSR-format
   integer, allocatable          :: Arwc(:)   ! A-matrix row entries
   integer, allocatable          :: Aclc(:)   ! A-matrix column entries

   double precision, allocatable :: Acsc(:)   ! A-matrix in CSC-format
   integer, allocatable          :: Arwi(:)   ! A-matrix row entries
   integer, allocatable          :: Acle(:)   ! A-matrix column entries

   double precision, allocatable :: RHS(:)    ! righthand side
   double precision, allocatable :: SOL(:)    ! solution
   !
   ! cell section
   !
   type :: CellData
     real    :: x(3)                      ! center coordinates
     real    :: vol                       ! volume
     integer :: ctid = 1                  ! fluid type id as set in fluid table
   end type

   type(CellData), allocatable :: Cell(:)! array of computational cells

   !
   ! boundary section (both set by input and default boundaries)
   !
   type :: BoundaryData
     integer :: face                     ! belongs to face...
     integer, dimension(4) :: vertices   ! the 4 vertices, to be done allocatable
     integer :: rid                      ! region id as set in rtable
     real    :: distance                 ! normal distance from cell face
                                         ! center to cell center
     real    :: yplus                    ! y+
     real    :: uplus                    ! u+
     real    :: shear(3) = 0.0           ! shearstress components
    !real    :: normal(3) = 0.0          ! normalstress components
     real    :: h                        ! local heattransfer coef.
     real    :: q                        ! local heat flux (in W/m2)
     real    :: T                        ! local wall temperature
   end type

   type(BoundaryData), allocatable :: Bnd(:) ! array of boundaries

   real      :: Split    = 0.0  ! moet anders in de calls...
   integer   :: NOutlets =  0
   !
   ! boundary region definitions
   !
   type :: RegionData
     character(len=32) :: name = '-'     ! region name
     integer :: typ =       RWall        ! type
     real    :: uvw(3) =     0.0         ! place for velocities
     real    :: den =        1.205       ! place for density
     real    :: T =        293.0         ! place for temperature
     real    :: P =          0.0         ! place for pressure
     real    :: R =          0.0         ! place for resistance
     real    :: k =          1.e-6       ! place for turb. kin. energy
     real    :: e =          1.0         ! place for turb. diss.
     logical :: noslip =   .true.        ! slip/noslip switch
     logical :: std  =     .true.        ! standard/roughness switch
     real    :: elog =       9.0         ! E-parameter
     real    :: ylog =      11.0         ! y+ value where u+=y+ matches wall law
     real    :: z0   =       0.03        ! roughness parameter
     real    :: cs   =       0.3         ! roughness parameter factor
     logical :: split =    .true.        ! split/fixed flow rate switch
     real    :: splvl =      1.0
     logical :: adiab =    .true.        ! adiab./fixed temperature switch
     logical :: flux  =    .false.       ! fixed flux temperature switch
     real    :: area  =      0.0         ! total region area (set in ...)
     integer :: n =           0          ! number of boundaries using this region

     real    :: massflow =   0.0         ! total mass flow of region
     real    :: flowfact =   1.0         ! massflux corr. factor
     real    :: mffixed =    0.0         ! prescribed mass flow
     logical :: fluxmass = .false.       ! true for presc. mass flow
     logical :: table =    .false.       ! use tabelised profiles (still in userinlet)

     real    :: ForceNormal(3) = 0.0
     real    :: ForceTangnt(3) = 0.0 

     logical :: perfect =  .true.        ! perfect matching (default) or search 
     real    :: cnt(3) = (/0.0,0.0,0.0/) ! central point
     real    :: dir(3) = (/1.0,0.0,0.0/) ! direction vector
     real    :: trn(3) = (/0.0,0.0,0.0/) ! translation vector
     real    :: angle  =     0.0         ! rotation angle (=0.0 => translation)

     logical :: user  =    .false.       ! call special user subroutine?
     character(len=4) :: name1 = '    '  ! IDs of special user test subroutines
     character(len=4) :: name2 = '    '  ! IDs of special user test subroutines
   end type

   type(RegionData), allocatable :: Reg(:) ! array of regions

end module geometry
module particles

   integer :: Nsens = 0

   type SensorProp
     real    :: x(3)    =  0.0                        ! coordinates
     real    :: v(3)    =  0.0                        ! velocity
     real    :: p       =  0.0                        ! variable
     real    :: te      =  0.0                        ! variable
     real    :: ed      =  0.0                        ! variable
     real    :: T       =  0.0                        ! variable
     integer :: cell    =   0                         ! cell the sensor is in
     integer :: face    =   0                         ! face the sensor is on
     logical :: wall    = .false.
     logical :: defined = .false.
   end type

   type(SensorProp), allocatable    :: Sensor(:)

   integer :: Npart = 0

   type ParticleProp
     real    :: x0(3) = 0.0                        ! starting coordinates
     real    :: v0(3) = 0.0                        ! starting velocity
     real    :: x1(3) = 0.0                        ! current coordinates
     real    :: v1(3) = 0.0                        ! current velocity
     real    :: a1(3) = 0.0                        ! current acceleration
     integer :: cell0 =  0                         ! starting cell the particle is in
     integer :: cell1 =  0                         ! current cell the particle is in
     integer :: face  =  0                         ! current face the particle is on
     real    :: dens0 = -1.0  ! <= undefined flag  ! current density
     real    :: diam0 =  1.0
     real    :: mass0 =  1.0
     logical :: wall  = .false.
   end type

   type ParticleTrackData
     real    x(3)                                       ! coordinates
     real    v(3)                                       ! velocity
     integer cell                                       ! cell the particle is in
     real    time                                       ! travelling time since start
     type(ParticleTrackData), pointer :: next => NULL() ! list members
   end type

   type ParticleList
     type(ParticleTrackData), pointer :: head => NULL() ! start of list i
     type(ParticleTrackData), pointer :: last => NULL() ! last entry
     integer                          :: n              ! number of entries
   end type ParticleList

   type(ParticleList), allocatable    :: Tracks(:)      ! all the tracks
   
   type(ParticleTrackData), pointer   :: Track, TrackEntry

   type(ParticleProp), allocatable    :: Particle(:)

   type(ParticleProp), allocatable    :: ParticleTmp(:) ! temporary array to store

end module particles
module scalars
   !
   ! scalar boundary region definitions
   !
   type :: ScalarRegionData
     integer :: typ      =  -1            ! type
     real    :: fraction =  0.0           ! boundary mass fraction
     logical :: adiab    = .true.         ! adiab./fixed temperature switch
     logical :: flux     = .false.        ! fixed flux temperature switch
     real    :: value    =  0.0           ! value
     logical :: user     = .false.        ! call special user subroutine?
   end type

   type(ScalarRegionData), allocatable :: ScReg(:,:)  ! array of scalar regions

   real, allocatable :: ScFluxIn(:,:)     ! store influx for each scalar
   real, allocatable :: ScFluxOut(:,:)    ! store outflux for each scalar

   !
   ! molecular scalar properties
   !
   type :: ScalarData
     character(len=12) name              ! name
     character(len=12) unit              ! unit

     logical :: active                   ! passive (simple tranport) or active
     logical :: solve                    ! use it or not
     logical :: save                     ! if true, variable field will be written on restart file

     real    :: PrL = 1.0                ! laminar Prandtl number
     real    :: PrT = 1.0                ! turbulent Prandtl number

     real    :: den                      ! boundary mass fraction (kg/m3)
     real    :: mol                      ! molecular weight (kg/kmol)
     real    :: beta                     ! thermal expansion coefficient (1/K)
     real    :: vis                      ! molecular viscosity (kg/ms)
     real    :: cp                       ! specific heat (J/kgK)
     real    :: lambda                   ! conductivity (W/mK)
     real    :: heat                     ! heat of formation (J/kg)
     real    :: temp                     ! temp. of formation (K)
   end type

   type(ScalarData), allocatable :: ScProp(:)  ! array of scalar properties

end module scalars
!========================================================================
module artificial_compressibility

    !double precision, save :: delta_t !the pseudo timestep
    double precision, dimension(3) :: sum_VS          !factor for pseudo time step
    double precision, dimension(3) :: sum_area        !the sum of the areas in the direction of the flow
    double precision, dimension(3) :: U_P, V_P, W_P   !local velocity vectors
    double precision :: cfl = 1d50                     !courant criterium
    double precision, allocatable :: art_snd_speed(:) !artificial sound speed
    double precision :: ac_beta = 1d50                 !art snd speed prefactor
    double precision, allocatable :: delta_t(:)       !the pseudo timestep
    double precision, allocatable :: density_change(:)!artificial change in density
    integer :: iter

end module artificial_compressibility
!========================================================================
module variables
   !
   ! main variables defined at the cell centers
   !
   integer :: Nscal = 0, Nmat = 0


   character(len=12), allocatable  :: Variable(:)
   character(len=12), allocatable  :: Scalar(:)

   real, allocatable :: U(:),V(:),W(:), TE(:), ED(:), &
                        VisEff(:), TurbP(:), T(:), Res(:)

   real, allocatable :: Cp(:),DEN(:), DENold(:), DENdp(:)

   real, allocatable :: Uold(:), Vold(:), Wold(:),  TEold(:),  EDold(:),  Told(:)
   real, allocatable :: Uold2(:),Vold2(:),Wold2(:), TEold2(:), EDold2(:), Told2(:)

   real, allocatable :: SC(:,:), SCold(:,:), SCold2(:,:)

   real, allocatable :: P(:), PP(:)
   !double precision, allocatable  :: P(:) <=moet op een gegeven moment

   real, allocatable    :: Residual(:)
   real, allocatable    :: ResiNorm(:)

   integer, allocatable :: Solver(:)           ! blending factor UDS en HigherOrder DS.
   real, allocatable    :: Gamma(:)            ! blending factor UDS en HigherOrder DS.
   integer, allocatable :: Scheme(:)           ! choice of DS

   logical, allocatable :: Solve(:)            ! simple array filled in after readcontrolfile
   logical, allocatable :: Store(:)            ! simple array filled in after readcontrolfile

   logical, allocatable :: SolveScalar(:)      ! solve extra scalars
   logical, allocatable :: StoreScalar(:)      ! solve extra scalars

   logical, allocatable :: PostC(:)            ! postprocessing arrays, cells
   logical, allocatable :: PostV(:)            ! vertices

   !
   ! variables defined at the cell faces only
   !
   real, allocatable :: MassFlux(:)
   real, allocatable :: MassFluxDebug(:,:)

   real, allocatable :: DXdebug(:), DXgrad(:,:)
   !
   ! gradients
   !
   real, allocatable :: dudx(:,:), dvdx(:,:), dwdx(:,:), &
                        dpdx(:,:)

   real, allocatable :: dPhidXo(:,:)
   !
   ! materials section
   !
   !real, allocatable :: Tref(:), Pref(:), DensRef(:), VisMol(:)

   integer :: IPref  =   1        ! ref. cell for pressure (material 1)
   integer :: IMoni  =   1        ! monitor cell for pressure (material 1)
   integer :: Iter   =   0        ! iteration counter
   real    :: Time   =   0.0      ! current time (if applicable)

   real :: Tref      = 273.       ! ref. cell for temperature (material 1)
   real :: Pref      =   0.0      ! ref. cell for pressure (material 1)

   real :: DensRef   =   1.2      ! <= later lucht op 20C
   real :: VisLam    =   0.001    ! <= later lucht op 20C

   real :: Prandtl   =   0.6905   ! lucht op 20C
   real :: Schmidt   =   0.9

   real :: Sigma_T   =   0.9      !
   real :: Sigma_k   =   1.0      ! turbulence diff. coef. factors
   real :: Sigma_e   =   1.219    ! ie. turbulent Prandtl numbers
   real :: Sigma_s   =   0.9      !

   real :: Gravity(3) =     0.0   ! <= gravity vector
   real :: Beta       =     0.001 ! expansie coef.
   real :: CpStd      =  1006.0   !
   real :: CvStd      =  1006.0   !
   real :: Lambda     =     0.02637 ! warmtegeleiding lucht

   real :: Qtransfer  = 0.0
   !                   u    v    w    p     k    e      T    Sc    Den  PP
   real :: URF(10) =(/0.5, 0.5, 0.5, 0.2,  0.5, 0.5,   0.95, 0.95, 1.0, 0.8 /)
   real :: RTOL(8) =(/0.1, 0.1, 0.1, 0.05, 0.1, 0.1,   0.1 , 0.1            /)
   real :: ATOL(8) =(/0.0, 0.0, 0.0, 0.0,  0.0, 0.0,   0.0 , 0.0            /)
   real :: Guess(8)=(/0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 293.0 , 0.0            /) !<==!!!!

   real    :: ResMax   = 1.e-4
   integer :: MaxPCOR  =    4
   real    :: FactDPP  = 0.25

   !
   ! common block section for user data related stuff
   !
   integer, parameter :: Nuser = 4
   integer, dimension(Nuser) :: UserInt
   real, dimension(Nuser)    :: UserReal
   !
   ! special post processing section
   !
   logical :: DXnormals     = .false.
   logical :: DXcenters     = .false.
   logical :: DXmassflux    = .false.
   logical :: DXtemperature = .false.
   logical :: DXdebugdata   = .false.
   integer :: DXdump        =    100

end module variables
