!==================================================================================================
MODULE Setup_fine !Defines global variables, reads input from file, allocates arrays
!==================================================================================================
USE SetPrecision
USE Setup

IMPLICIT NONE

!**************************** Global Simulation Quanitities ***************************************

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LBM Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! NumDistDirs - same as Setup

REAL(dbl),		ALLOCATABLE :: f_fine(:,:,:,:)							! distribution function
REAL(dbl), 		ALLOCATABLE :: fplus_fine(:,:,:,:)						! post-collision distribution function
REAL(dbl), 		ALLOCATABLE :: u_fine(:,:,:),v_fine(:,:,:),w_fine(:,:,:)		! x,y, and z components of the fluid velocity vector
REAL(dbl), 		ALLOCATABLE :: rho_fine(:,:,:)							! density
INTEGER(lng), 	ALLOCATABLE :: node_fine(:,:,:)    					! node flags (FLUID/SOLID)
! ex(:),ey(:),ez(:) - same as Setup
! bb_fine(:), sym_fine(:,:) - same as Setup
REAL(dbl), 		ALLOCATABLE :: wt_fine(:)    							! weighting coefficients for the equilibrium distribution functions

! den, denL - same as Setup 
! nu, nuL - same as Setup
! cs - same as Setup
REAL(dbl)		:: tau_fine            ! relaxation parameters of coarse and fine blocks
REAL(dbl)		:: oneOVERtau_fine     ! reciprical of tau
INTEGER(lng)	:: nx_fine,ny_fine,nz_fine     ! number of global nodes in the x, y, and z directions respectively
INTEGER(lng)	:: nxSub_fine,nySub_fine,nzSub_fine   ! number of local nodes in the each direction
INTEGER(lng)	:: iter0_fine,iter_fine,nt_fine	      ! initial time step, timestep index, total number of timesteps
INTEGER         :: subIter                            ! Sub-iteration counter
! domaintype - same as Setup
! FLUID - same as Setup
! SOLID - same as Setup
INTEGER(lng), PARAMETER :: COARSEMESH		= -1_lng					! coarseMesh

! restart - same as Setup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Scalar Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REAL(dbl), ALLOCATABLE :: phi_fine(:,:,:)			! passive scalar
REAL(dbl), ALLOCATABLE :: delphi_particle_fine(:,:,:)	! passive scalar contribution from particles
REAL(dbl), ALLOCATABLE :: phiTemp_fine(:,:,:)		! temporary storage of passive scalar
! Sc - same as Setup
! Dm,Dmcf - same as Setup
! Delta_fine - same as Setup
! phiIC, phiWall - same as Setup
REAL(dbl) :: phiAbsorbed_fine							! total amount of scalar absorbed up to current time
REAL(dbl) :: phiAbsorbedS_fine							! total amount of scalar absorbed up to current time - through the macroscopic surface
REAL(dbl) :: phiAbsorbedV_fine							! total amount of scalar absorbed up to current time - through the villi
REAL(dbl) :: phiInOut_fine							! total amount of scalar leaving/entering the domain
REAL(dbl) :: phiTotal_fine							! total intial amount of scalar in the domain
! sigma - same as Setup
! phiPer - same as Setup
! phiStart - same as Setup
! sclrIC - same as Setup

! ee - same as Setup

! BLOB - same as Setup
! LINE - same as Setup
! INLET - same as Setup
! UNIFORM - same as Setup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parallel (MPI) Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Number of communication directions (3D LBM)
! NumCommDirs           =  26_lng - same as Setup
! MaxDistFns		= 5_lng	- same as Setup
! NumFs_face		= 5_lng - same as Setup
! NumFs_side		= 2_lng - same as Setup
! NumFs_corner	        = 1_lng - same as Setup

! MPI Arrays (arranged by descending size for storage efficiency)
! f_Comps(:,:)  - same as Setup
INTEGER(lng), ALLOCATABLE :: Corner_SendIndex_fine(:,:)						! i, j, and k indices for each corner
INTEGER(lng), ALLOCATABLE :: Corner_RecvIndex_fine(:,:)						! i, j, and k indices for each corner (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: Z_SendIndex_fine(:,:)						! i and j indices for each Z side 
INTEGER(lng), ALLOCATABLE :: Z_RecvIndex_fine(:,:)						! i and j indices for each Z side (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: X_SendIndex_fine(:,:)						! j and k indices for each X side 
INTEGER(lng), ALLOCATABLE :: X_RecvIndex_fine(:,:)						! j and k indices for each X side (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: Y_SendIndex_fine(:,:)						! i and k indices for each Y side 
INTEGER(lng), ALLOCATABLE :: Y_RecvIndex_fine(:,:)						! i and k indices for each Y side (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: YZ_SendIndex_fine(:)										! i index for each YZ face 
INTEGER(lng), ALLOCATABLE :: YZ_RecvIndex_fine(:)						! i index for each YZ face (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: ZX_SendIndex_fine(:)						! j index for each ZX face 
INTEGER(lng), ALLOCATABLE :: ZX_RecvIndex_fine(:)						! j index for each ZX face (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: XY_SendIndex_fine(:)										! k index for each XY face 
INTEGER(lng), ALLOCATABLE :: XY_RecvIndex_fine(:)						! k index for each XY face (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: SubID_fine(:)							! id number of neighboring subdomains (same as rank of processing unit working on domain)
! OppCommDir_fine - Same as coarse mesh
INTEGER(lng), ALLOCATABLE :: CommDataStart_f_fine(:)						! array of starting indices in the send arrays for the distribution functions from each communication direction 
INTEGER(lng), ALLOCATABLE :: CommDataStart_rho_fine(:)						! array of starting indices in the send arrays for the density from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_phi_fine(:)						! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_u_fine(:)						! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_v_fine(:)						! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_w_fine(:)						! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: fSize_fine(:)							! array of the number of elements sent for each communication direction (distribution functions)
INTEGER(lng), ALLOCATABLE :: dsSize_fine(:)							! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: uvwSize_fine(:)							! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: msgSize_fine(:)							! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: req_fine(:)							! array of MPI send/receive requests 
INTEGER(lng), ALLOCATABLE :: waitStat_fine(:,:)							! array of MPI_WAITALL status objects
REAL(dbl), ALLOCATABLE :: msgSend_fine(:)							! array of ALL of the sent information (total)
REAL(dbl), ALLOCATABLE :: msgRecv_fine(:)							! array of ALL of the received information (total)

! MPI Variables
! master = 0_lng - same as Setup
INTEGER(lng) :: numprocs_fine, myid_fine, mySub_fine			! number of processing units, rank of current processing unit, subdomain of current processing unit

REAL(dbl) :: CommTime_f0_fine, CommTime_fEnd_fine, CommTime_f_fine	! communication time - distribution functions: start time, end time, current time
REAL(dbl) :: CommTime_ds0_fine, CommTime_dsEnd_fine, CommTime_ds_fine   ! communication time - distribution functions: start time, end time, current time

! Starting/Ending indices for each subdomain
INTEGER(lng) :: iMin_fine				! starting local i index
INTEGER(lng) :: iMax_fine				! ending local i index
INTEGER(lng) :: jMin_fine				! starting local j index
INTEGER(lng) :: jMax_fine				! ending local j index
INTEGER(lng) :: kMin_fine				! starting local k index
INTEGER(lng) :: kMax_fine				! ending local k index

! extension for output files
! sub	- same as Setup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Geometry Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REAL(dbl), ALLOCATABLE	:: villiLoc_fine(:,:)		    ! location of each villous
REAL(dbl), ALLOCATABLE 	:: x_fine(:),y_fine(:),z_fine(:)    ! physical coordinate arrays
REAL(dbl), ALLOCATABLE	:: xx_fine(:),yy_fine(:),zz_fine(:) ! x,y,z arrays (global)
REAL(dbl), ALLOCATABLE 	:: ub_fine(:),vb_fine(:),wb_fine(:)		! x,y, and z components of the solid boundary velocity vector
REAL(dbl), ALLOCATABLE 	:: rDom0_fine(:),rDom_fine(:),r_fine(:)        ! initial, and current radius at each z-location (global), radius at each location (local)
REAL(dbl), ALLOCATABLE	:: velDom_fine(:),vel_fine(:)		! global and local wall velocities 
REAL(dbl), ALLOCATABLE	:: rnd_fine(:)			! array of random numbers for random villi phase angles
INTEGER(lng), ALLOCATABLE :: villiGroup_fine(:)		! array of which groups the villi are in 
REAL(dbl)		:: xcf_fine, ycf_fine, zcf_fine	! length conversion factors
REAL(dbl)		:: dcf_fine, vcf_fine, pcf_fine	! density, velocity, pressure conversion factors
REAL(dbl)		:: tcf_fine			! time conversion factor
REAL(dbl)		:: nPers_fine			! number of time periods simulated
! Lv - same as Setup
! Rv - same as Setup
! villiAngle - same as Setup
INTEGER(lng)	:: iLv_fine	! length of the villi in lattice units
INTEGER(lng)	:: Ci_fine,Cj_fine,Ck_fine ! center node location (global)

! randORord								! flag to determine if the villous motion is random or ordered
! RANDOM=1  - same as Setup
! ORDERED=2 - same as Setup

! PI = 3.1415926535897932384626433832 - same as Setup
! D, L - same as Setup
REAL(dbl) :: fractionDfine						! fraction of the diameter used for the fine mesh
INTEGER ::   gridRatio						! Ratio of coarse/fine grid size 
! a1, a2 - same as Setup
! eps1, eps2 - same as Setup
! amp1, amp2 - same as Setup
! epsOVERa1, epsOVERa2 - same as Setup
! aOVERlam1,	aOVERlam2 - same as Setup
! lambda1, lambda2 - same as Setup
! kw1 - same as Setup												! wave number (peristalsis)
! s1, s2 - same as Setup
! Ts, Tp - same as Setup
! wc1, wc2 - same as Setup
! Re1, Re2 - same as Setup
! shift2 - same as Setup
! Tmix - same as Setup
! period - same as Setup
! freqRatioT - same as Setup
! freqRatioZ - same as Setup
! vFreqT - same as Setup
! vFreqZ - same as Setup
! activeVflagT - same as Setup
! activeVflagZ - same as Setup
! nlambda2 - same as Setup
! numw1, numw2 - same as Setup
! nzSL, nzSR - same as Setup
! segment - same as Setup
! seg1L, seg1R - same as Setup
! seg2L, seg2R - same as Setup
! numVilliZ, numVilliTheta - same as Setup
! numVilli, numVilliActual - same as Setup
! numVilliGroups - same as Setup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~Multi-grid algorithm variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
real(dbl), allocatable, dimension(:,:,:,:) :: fCtoF_topXZ
real(dbl), allocatable, dimension(:,:,:,:) :: fCtoF_bottomXZ
real(dbl), allocatable, dimension(:,:,:,:) :: fCtoF_frontYZ
real(dbl), allocatable, dimension(:,:,:,:) :: fCtoF_backYZ                                

real(dbl), allocatable, dimension(:,:,:) :: feqFC_topXZ
real(dbl), allocatable, dimension(:,:,:) :: feqFC_bottomXZ
real(dbl), allocatable, dimension(:,:,:) :: feqFC_frontYZ
real(dbl), allocatable, dimension(:,:,:) :: feqFC_backYZ

real(dbl), allocatable, dimension(:,:,:) :: fFtoC_topXZ
real(dbl), allocatable, dimension(:,:,:) :: fFtoC_bottomXZ
real(dbl), allocatable, dimension(:,:,:) :: fFtoC_frontYZ
real(dbl), allocatable, dimension(:,:,:) :: fFtoC_backYZ             

real(dbl), allocatable, dimension(:,:,:) :: feqFF_topXZ
real(dbl), allocatable, dimension(:,:,:) :: feqFF_bottomXZ
real(dbl), allocatable, dimension(:,:,:) :: feqFF_frontYZ
real(dbl), allocatable, dimension(:,:,:) :: feqFF_backYZ
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REAL(dbl), ALLOCATABLE	   :: radius_fine(:,:)  ! radius stored during output iterations
INTEGER(lng), ALLOCATABLE  :: filenum_fine(:)   ! array of output file numbers
INTEGER(lng)               :: numOuts_fine      ! number of output files
INTEGER(lng)	           :: fileCount_fine    ! current output file number (out of total number of output files)
INTEGER(lng)	           :: outFlag_fine      ! specifies whether to output in readable format (1), binaries (2), or both (3)
INTEGER(lng)               :: radcount_fine     ! counts the number of output iterations for storing the radius

! System Clock Variables (for PrintStatus)
! start, current, final, rate - same as Setup
! tStart,tEnd,tTotal,tRecv,tSum,tAvg - same as Setup

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Particle Tracking Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Removing all particle tracking variables
!************************************************

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Global_Setup_Fine		! sets up simulation
! The new code will now call this routine instead of the older one in Setup.f90. Since this module 'uses' the older Setup module, it will call the subroutines of the older module instead.    
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

CALL ReadInput_fine			! read input from file
CALL SubDomainSetup_fine		! set up the MPI subdomains
CALL AllocateArrays_fine		! allocate global variable arrays

!------------------------------------------------
END SUBROUTINE Global_Setup_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ReadInput_Fine			! read the input file
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Read input from input file
OPEN(10,FILE='input.txt')
READ(10,*) domaintype	 				! a flag to denote domain type - 0 for 1/4th cylinder and 1 for full cylinder
READ(10,*) nx	 				! number of nodes in the x-direction
READ(10,*) ny					! number of nodes in the y-direction
READ(10,*) nz					! number of nodes in the z-direction

READ(10,*) NumSubsX			! number of subdomains in the X direction
READ(10,*) NumSubsY			! number of subdomains in the Y direction
READ(10,*) NumSubsZ			! number of subdomains in the Z direction

READ(10,*) L					! length
READ(10,*) D					! diameter
READ(10,*) fractionDfine			! fraction of the diameter used for the fine mesh
READ(10,*) gridRatio			        ! Ratio of coarse/fine grid size

IF (gridRatio .gt. 1) THEN
   nx_fine = ( fractionDfine * (nx - 1) + 2 ) * gridRatio + 1
   ny_fine = ( fractionDfine * (ny - 1) + 2 ) * gridRatio + 1
   nz_fine = nz * gridRatio
END IF


READ(10,*) epsOVERa1			! peristaltic occlusion ratio (distance of occlusion/mean half-width)
READ(10,*) s1					! peristaltic wave speed
READ(10,*) numw1				! number of peristaltic waves
READ(10,*) wc1					! peristaltic weighting coefficient

READ(10,*) epsOVERa2			! segmental occlusion ratio (distance of occlusion/mean half-width)
READ(10,*) Ts					! segmental contraction period
READ(10,*) numw2				! number of segmental waves
READ(10,*) wc2					! segmental weighting coefficient

READ(10,*) Tmix				! period of mixed mode simulation

READ(10,*) numVilliZ			! number of rows of villi in the Z-direction
READ(10,*) numVilliTheta	! number of villi in each row (theta-direction)
READ(10,*) numVilliGroups	! number of villi in each row (theta-direction)

READ(10,*) Lv					! length of the villi (micrometers)
READ(10,*) Rv					! radius of the villi (micrometers)

READ(10,*) freqRatioT		! villous frequency to macroscopic contraction frequency (azimuthal, theta direction)
READ(10,*) freqRatioZ		! villous frequency to macroscopic contraction frequency (axial, z direction)

READ(10,*) randORord			! flag to determine if the villous motion is random or ordered
READ(10,*) villiAngle		! maximum angle of active villous travel (degrees) 

READ(10,*) den					! density
READ(10,*) nu					! kinematic viscosity
READ(10,*) tau					! relaxation parameter
READ(10,*) tau_fine				! relaxation parameter - Fine mesh

READ(10,*) Sc					! Schmidt number
READ(10,*) sclrIC				! initial/maintained scalar distribution (1=BLOB,2=LINE,3=INLET,4=UNIFORM)
READ(10,*) phiPer				! period at which to start the scalar
READ(10,*) phiIC				! maximum scalar concentration

READ(10,*) nPers				! total number of periods to run
READ(10,*) numOuts			! number of output files (roughly)
READ(10,*) restart			! use restart file? (0 if no, 1 if yes)
READ(10,*) ParticleTrack		! A flag to indicate if particle is on or off (0 if off, 1 if on)
CLOSE(10)

tau=1.0_dbl

! Check to make sure the number of processors (numprocs) and the number of subdomains are equal
NumSubsTotal = NumSubsX*NumSubsY*NumSubsZ
IF(NumSubsTotal .NE. numprocs) THEN
  OPEN(1000,FILE="error.dat")
  WRITE(1000,*) 'NumSubsTotal is .NE. to numprocs'
  WRITE(1000,*) 'NumSubsTotal', NumSubsTotal
  WRITE(1000,*) 'NumSubsX', NumSubsX
  WRITE(1000,*) 'NumSubsY', NumSubsY
  WRITE(1000,*) 'NumSubsZ', NumSubsZ
  WRITE(1000,*) 'numprocs', numprocs
  WRITE(1000,*) 'check input file and PBS queueing file...'
  CLOSE(1000)
  STOP
END IF

! Check to make sure the number of villi in the axial direction is a multiple of the number of villi groups
IF(MOD(numVilliZ,numVilliGroups) .NE. 0) THEN
  OPEN(1000,FILE="error.dat")
  WRITE(1000,*) 'numVilliZ is not a multiple of numVilliGroups'
  WRITE(1000,*) 'check input file...'
  CLOSE(1000)
  STOP
END IF

!------------------------------------------------
END SUBROUTINE ReadInput_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SubDomainSetup_Fine	! generates the information (ID number, starting/ending indices) of each neighboring subdomain (using the subroutine SetSubIDBC in this module)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: CDx(NumCommDirs), CDy(NumCommDirs), CDz(NumCommDirs) ! communication direction vectors in the x, y, and z directions respectively
INTEGER(lng) :: thisSub	 ! ID of the current subdomain
INTEGER(lng) :: iComm,iSub,jSub,kSub,iiSub,jjSub,kkSub ! index variables
INTEGER(lng) :: quotientX, quotientY, quotientZ	 ! variables for determining the local subdomain bounds

ALLOCATE(SubID_fine(NumCommDirs)) ! id number of neighboring subdomains (same as rank of processing unit working on domain)

! fill out communication direction vectors
CDx(1) =   1_lng
CDy(1) =   0_lng
CDz(1) =   0_lng

CDx(2) =  -1_lng
CDy(2) =   0_lng
CDz(2) =   0_lng

CDx(3) =   0_lng
CDy(3) =   1_lng
CDz(3) =   0_lng

CDx(4) =   0_lng
CDy(4) =  -1_lng
CDz(4) =   0_lng

CDx(5) =   0_lng
CDy(5) =   0_lng
CDz(5) =   1_lng

CDx(6) =   0_lng
CDy(6) =   0_lng
CDz(6) =  -1_lng

CDx(7) =   1_lng
CDy(7) =   1_lng
CDz(7) =   0_lng

CDx(8) =  -1_lng
CDy(8) =  -1_lng
CDz(8) =   0_lng

CDx(9) =   1_lng
CDy(9) =  -1_lng
CDz(9) =   0_lng

CDx(10) = -1_lng
CDy(10) =  1_lng
CDz(10) =  0_lng

CDx(11) =  0_lng
CDy(11) =  1_lng
CDz(11) =  1_lng

CDx(12) =  0_lng
CDy(12) = -1_lng
CDz(12) = -1_lng

CDx(13) =  0_lng
CDy(13) =  1_lng
CDz(13) = -1_lng

CDx(14) =  0_lng
CDy(14) = -1_lng
CDz(14) =  1_lng

CDx(15) =  1_lng
CDy(15) =  0_lng
CDz(15) =  1_lng

CDx(16) = -1_lng
CDy(16) =  0_lng
CDz(16) = -1_lng

CDx(17) = -1_lng
CDy(17) =  0_lng
CDz(17) =  1_lng

CDx(18) =  1_lng
CDy(18) =  0_lng
CDz(18) = -1_lng

CDx(19) =  1_lng
CDy(19) =  1_lng
CDz(19) =  1_lng

CDx(20) = -1_lng
CDy(20) = -1_lng
CDz(20) = -1_lng

CDx(21) =  1_lng
CDy(21) =  1_lng
CDz(21) = -1_lng

CDx(22) = -1_lng
CDy(22) = -1_lng
CDz(22) =  1_lng

CDx(23) = -1_lng
CDy(23) =  1_lng
CDz(23) =  1_lng

CDx(24) =  1_lng
CDy(24) = -1_lng
CDz(24) = -1_lng

CDx(25) =  1_lng
CDy(25) = -1_lng
CDz(25) =  1_lng

CDx(26) = -1_lng
CDy(26) =  1_lng
CDz(26) = -1_lng

! Number of the current subdomain
mySub = myid + 1_lng							! subdomain number
!WRITE(sub(1:2),'(I2.2)') mySub			! write subdomain number to 'sub' for output file exentsions
WRITE(sub(1:5),'(I5.5)') mySub			! write subdomain number to 'sub' for output file exentsions

! Loop through the subdomains
DO kSub=1,NumSubsZ
  DO jSub=1,NumSubsY
    DO iSub=1,NumSubsX

      thisSub = iSub + (jSub-1)*NumSubsX + (kSub-1)*NumSubsX*NumSubsY	! get the ID of the current Subdomain

      IF(mySub .EQ. thisSub) THEN													! fill out the SubID array of the current subdomain is the 
     
        ! Loop through the communication directions for the current subdomain
        DO iComm=1,NumCommDirs

          iiSub = iSub + CDx(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
          jjSub = jSub + CDy(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
          kkSub = kSub + CDz(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
      
          CALL SetSubID_fine(iComm,iiSub,jjSub,kkSub)								! identify the neighboring subdomains (SubID)

        END DO

      END IF

    END DO
  END DO
END DO

! Define the local computational domain bounds (iMin:iMax,jMin:jMax,kMin:kMax)
quotientX	= CEILING(REAL(nx_fine)/NumSubsX)	! divide the number of nodes by the number of subdomains (round up)
quotientY	= CEILING(REAL(ny_fine)/NumSubsY)	! divide the number of nodes by the number of subdomains (round up)
quotientZ	= CEILING(REAL(nz_fine)/NumSubsZ)	! divide the number of nodes by the number of subdomains (round up)

iMin_fine = MOD(myid,NumSubsX)*quotientX + 1				! starting local i index 
iMax_fine = iMin_fine + (quotientX - 1)					! ending local i index

jMin_fine = MOD((myid/NumSubsX),NumSubsY)*quotientY + 1	                ! starting local j index
jMax_fine = jMin_fine + (quotientY - 1)					! ending local j index

kMin_fine = (myid/(NumSubsX*NumSubsY))*quotientZ + 1		        ! starting local k index 
kMax_fine = kMin_fine + (quotientZ - 1)					! ending local k index

! Check the bounds
IF(iMax_fine .GT. nx_fine) THEN
  iMax_fine = nx_fine																! if iMax is greater than nx, correct it
END IF

IF(jMax_fine .GT. ny_fine) THEN
  jMax_fine = ny_fine																! if jMax is greater than ny, correct it
END IF

IF(kMax_fine .GT. nz_fine) THEN
  kMax_fine = nz_fine																! if kMax is greater than nz, correct it
END IF

! Determine the number of nodes in each direction
nxSub_fine = (iMax_fine - iMin_fine) + 1_lng
nySub_fine = (jMax_fine - jMin_fine) + 1_lng
nzSub_fine = (kMax_fine - kMin_fine) + 1_lng

! Write the local bounds to a file [TEST]
!OPEN(171,FILE='localBounds-'//sub//'.dat')
!WRITE(171,*) 'iMin =', iMin, 'iMax=', iMax
!WRITE(171,*) 'jMin =', jMin, 'jMax=', jMax 
!WRITE(171,*) 'kMin =', kMin, 'kMax=', kMax
!WRITE(171,*) 
!WRITE(171,*) 'nx =', nx, 'ny=', ny, 'nz=', nz
!WRITE(171,*) 'nxSub_fine =', nxSub, 'nySub=', nySub, 'nzSub=', nzSub
!WRITE(171,*) 
!WRITE(171,*) 'quotientX =', quotientX
!WRITE(171,*) 'quotientY =', quotientY
!WRITE(171,*) 'quotientZ =', quotientZ
!CLOSE(171)

! Write the subID to a file [TEST]
!OPEN(172,FILE='sub-'//sub//'.dat')
!DO iComm=1,NumCommDirs
!  WRITE(172,*)'iComm=',iComm,'SubID(iComm)=', SubID(iComm)
!END DO
!CLOSE(172)
!STOP

!------------------------------------------------
END SUBROUTINE SubDomainSetup_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SetSubID_Fine(iComm,iiSub,jjSub,kkSub)									! sets SubID based on neighboring subdomains
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng), INTENT(IN) :: iComm, iiSub, jjSub, kkSub 					! index variables
INTEGER(lng) :: nSub, kkSub2														! neighboring subdomain ID, kkSub (reset for periodicity)

IF(((jjSub .LT. 1) .OR. (jjSub .GT. NumSubsY))	.OR.	&
!   ((kkSub .LT. 1) .OR. (kkSub .GT. NumSubsZ))	.OR. 	& 					! comment out for periodic BCs in the k-direction
   ((iiSub .LT. 1) .OR. (iiSub .GT. NumSubsX)))	THEN 

  SubID_fine(iComm) = 0_lng																! no neighbor

ELSE IF((kkSub .LT. 1)) THEN

  kkSub2 = NumSubsZ																	! reset kkSub for periodicity in the z-direction
  nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
  SubID_fine(iComm) = nSub																! set SubID(iComm) to neighboring sudomain ID

ELSE IF((kkSub .GT. NumSubsZ)) THEN

  kkSub2 = 1_lng																		! reset kkSub for periodicity in the z-direction
  nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
  SubID_fine(iComm) = nSub																! set SubID(iComm) to neighboring sudomain ID
    
ELSE
  
  nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub-1)*NumSubsX*NumSubsY		! neighboring subdomain ID
  SubID_fine(iComm) = nSub																! set SubID(iComm) to neighboring sudomain ID

END IF

!------------------------------------------------
END SUBROUTINE SetSubID_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AllocateArrays_Fine	! allocates array space
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Distribution Functions
ALLOCATE(f_fine(0:NumDistDirs,0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1),			&
         fplus_fine(0:NumDistDirs,0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1))
! Velocity, Density
ALLOCATE(u_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1),							&
         v_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1),							&
         w_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1))
ALLOCATE(rho_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1))

! Scalar
ALLOCATE(phi_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1), 						&
         phiTemp_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1))
ALLOCATE(delphi_particle_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1))

! Node Flags
ALLOCATE(node_fine(0:nxSub_fine+1,0:nySub_fine+1,0:nzSub_fine+1))

! MPI Communication Arrays
ALLOCATE(Corner_SendIndex_fine(19:26,3))							! i, j, and k indices for each corner
ALLOCATE(Corner_RecvIndex_fine(19:26,3))							! i, j, and k indices for each corner (phantom node for recieving data)
ALLOCATE(Z_SendIndex_fine(7:10,2))									! i and j indices for each Z side 
ALLOCATE(Z_RecvIndex_fine(7:10,2))									! i and j indices for each Z side (phantom node for recieving data)
ALLOCATE(X_SendIndex_fine(11:14,2))								! j and k indices for each X side 
ALLOCATE(X_RecvIndex_fine(11:14,2))								! j and k indices for each X side (phantom node for recieving data)
ALLOCATE(Y_SendIndex_fine(15:18,2))								! i and k indices for each Y side 
ALLOCATE(Y_RecvIndex_fine(15:18,2))								! i and k indices for each Y side (phantom node for recieving data)
ALLOCATE(YZ_SendIndex_fine(1:2))									! i index for each YZ face 
ALLOCATE(YZ_RecvIndex_fine(1:2))									! i index for each YZ face (phantom node for recieving data)
ALLOCATE(ZX_SendIndex_fine(3:4))									! j index for each ZX face 
ALLOCATE(ZX_RecvIndex_fine(3:4))									! j index for each ZX face (phantom node for recieving data)
ALLOCATE(XY_SendIndex_fine(5:6))									! k index for each XY face 
ALLOCATE(XY_RecvIndex_fine(5:6))									! k index for each XY face (phantom node for recieving data)
ALLOCATE(CommDataStart_f_fine(NumCommDirs))						! array of starting indices in the send arrays for the distribution functions from each communication direction 
ALLOCATE(CommDataStart_rho_fine(NumCommDirs))					! array of starting indices in the send arrays for the density from each communication direction
ALLOCATE(CommDataStart_phi_fine(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_u_fine(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_v_fine(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_w_fine(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(fSize_fine(NumCommDirs))									! array of the number of elements sent for each communication direction (distribution functions)
ALLOCATE(dsSize_fine(NumCommDirs))									! array of the number of elements sent for each communication direction (density and scalar)
ALLOCATE(uvwSize_fine(NumCommDirs))									! array of the number of elements sent for each communication direction (density and scalar)
ALLOCATE(msgSize_fine(NumCommDirs))								! array of the number of elements sent for each communication direction (total)
ALLOCATE(req_fine(2*NumCommDirs))									! allocate the MPI send request array

! Geometry Arrays
ALLOCATE(rDom0_fine(0:nz_fine+1),rDom_fine(0:nz_fine+1),r_fine(0:nzSub_fine+1))		! intial and current radius (global), current radius (local)
ALLOCATE(velDom_fine(0:nz_fine+1),vel_fine(0:nzSub_fine+1))					! global and local wall velocities
ALLOCATE(x_fine(0:nxSub_fine+1),y_fine(0:nySub_fine+1),z_fine(-gridRatio+1:nzSub_fine+gridRatio))		! x, y, z, physical coordinate arrays (local)
ALLOCATE(xx_fine(0:nx_fine+1),yy_fine(0:ny_fine+1),zz_fine(0:nz_fine+1))				! x, y, z, physical coordinate arrays (global)


!Multigrid algorithm variables

!Post collision density distribution from the fine mesh for the coarse mesh
allocate(fFtoC_topXZ(1:14,46:56,nzSub))     !Includes the ends - Indices are directionalDensity, x Index, z Index
allocate(fFtoC_bottomXZ(1:14,46:56,nzSub))  !Includes the ends - Indices are directionalDensity, x Index, z Index
allocate(fFtoC_frontYZ(1:14,47:55,nzSub))     !Does not include the ends - Indices are directionalDensity, y Index, z Index
allocate(fFtoC_backYZ(1:14,47:55,nzSub))     !Does not include the ends - Indices are directionalDensity, y Index, z Index

!Equilibrium density distribution from the fine mesh for the coarse mesh
allocate(feqFC_topXZ(1:14,46:56,nzSub))     !Includes the ends - Indices are directionalDensity, x Index, z Index
allocate(feqFC_bottomXZ(1:14,46:56,nzSub))  !Includes the ends - Indices are directionalDensity, x Index, z Index
allocate(feqFC_frontYZ(1:14,47:55,nzSub))     !Does not include the ends - Indices are directionalDensity, y Index, z Index
allocate(feqFC_backYZ(1:14,47:55,nzSub))     !Does not include the ends - Indices are directionalDensity, y Index, z Index


!Post collision density distribution from the coarse mesh for the fine mesh
allocate(fCtoF_topXZ(1:14,3,nxSub_fine,-gridRatio+1:nzSub_fine+gridRatio))     !Includes the ends - Indices are directionalDensity, timeLevel, x Index, z Index
allocate(fCtoF_bottomXZ(1:14,3,nxSub_fine,-gridRatio+1:nzSub_fine+gridRatio))  !Includes the ends - Indices are directionalDensity, timeLevel, x Index, z Index
allocate(fCtoF_frontYZ(1:14,3,2:nySub_fine-1,-gridRatio+1:nzSub_fine+gridRatio))   !Does not include the ends - Indices are directionalDensity, timeLevel, y Index, z Index
allocate(fCtoF_backYZ(1:14,3,2:nySub_fine-1,-gridRatio+1:nzSub_fine+gridRatio))    !Does not include the ends - Indices are directionalDensity, timeLevel, y Index, z Index

!Equilibrium density distribution from the coarse mesh for the fine mesh
allocate(feqFF_topXZ(1:14,46:56,nzSub))     !Includes the ends - Indices are directionalDensity, x Index, z Index
allocate(feqFF_bottomXZ(1:14,46:56,nzSub))  !Includes the ends - Indices are directionalDensity, x Index, z Index
allocate(feqFF_frontYZ(1:14,47:55,nzSub))     !Does not include the ends - Indices are directionalDensity, y Index, z Index
allocate(feqFF_backYZ(1:14,47:55,nzSub))     !Does not include the ends - Indices are directionalDensity, y Index, z Index


!------------------------------------------------
END SUBROUTINE AllocateArrays_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE DEAllocateArrays_Fine	! allocates array space
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Distribution Functions
DEALLOCATE(f_fine,fplus_fine)

! Velocity, Density
DEALLOCATE(u_fine,v_fine,w_fine,rho_fine)

! Scalar
DEALLOCATE(phi_fine,phiTemp_fine,delphi_particle_fine)

! Node Flags
DEALLOCATE(node_fine)

! MPI Communication Arrays
DEALLOCATE(Corner_SendIndex_fine)		! i, j, and k indices for each corner
DEALLOCATE(Corner_RecvIndex_fine)		! i, j, and k indices for each corner (phantom node for recieving data)
DEALLOCATE(Z_SendIndex_fine)				! i and j indices for each Z side 
DEALLOCATE(Z_RecvIndex_fine)				! i and j indices for each Z side (phantom node for recieving data)
DEALLOCATE(X_SendIndex_fine)				! j and k indices for each X side 
DEALLOCATE(X_RecvIndex_fine)				! j and k indices for each X side (phantom node for recieving data)
DEALLOCATE(Y_SendIndex_fine)				! i and k indices for each Y side 
DEALLOCATE(Y_RecvIndex_fine)				! i and k indices for each Y side (phantom node for recieving data)
DEALLOCATE(YZ_SendIndex_fine)			! i index for each YZ face 
DEALLOCATE(YZ_RecvIndex_fine)			! i index for each YZ face (phantom node for recieving data)
DEALLOCATE(ZX_SendIndex_fine)			! j index for each ZX face 
DEALLOCATE(ZX_RecvIndex_fine)			! j index for each ZX face (phantom node for recieving data)
DEALLOCATE(XY_SendIndex_fine)			! k index for each XY face 
DEALLOCATE(XY_RecvIndex_fine)			! k index for each XY face (phantom node for recieving data)
DEALLOCATE(CommDataStart_f_fine)		! array of starting indices in the send arrays for the distribution functions from each communication direction 
DEALLOCATE(CommDataStart_rho_fine)		! array of starting indices in the send arrays for the density from each communication direction
DEALLOCATE(CommDataStart_phi_fine)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_u_fine)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_v_fine)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_w_fine)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(fSize_fine)						! array of the number of elements sent for each communication direction (distribution functions)
DEALLOCATE(dsSize_fine)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(uvwSize_fine)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(msgSize_fine)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(req_fine)						! array of MPI send/receive requests
DEALLOCATE(waitStat_fine)					! array of MPI_WAITALL requests

! Geometry Arrays
DEALLOCATE(rDom0_fine,rDom_fine,r_fine)			! intial and current radius (global), current radius (local)
DEALLOCATE(velDom_fine,vel_fine)				! global and local wall velocities
DEALLOCATE(x_fine,y_fine,z_fine)						! x, y, z, physical coordinate arrays (local)
!------------------------------------------------
END SUBROUTINE DEAllocateArrays_Fine
!------------------------------------------------

!================================================
END MODULE Setup_fine
!================================================
