!==================================================================================================
MODULE Geometry_fine	! Defines the geometry for the fine mesh in the simulation
						! Subroutines (NodeFlags, BoundaryVelocity)
!==================================================================================================
USE SetPrecision      
USE Setup
USE Setup_fine
USE LBM_fine
USE ICBC_fine
USE MPI

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Geometry_Setup_fine					! sets up the geometry
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER :: isize,idate(8)									! size of seed array for the random number genreator, array for output of DATE_AND_TIME
INTEGER,ALLOCATABLE  :: iseed(:)							! seeds for random number generator
INTEGER(lng) :: i,j,k,kk,iCon,it,iPer,nPers_INT		! index variables
INTEGER(lng) :: nvz,nvt,n,g								! index variables
INTEGER(lng) :: mpierr										! MPI standard error variable 
REAL(dbl) :: macroFreq										! macroscopic contraction frequency
INTEGER(lng) :: xaxis,yaxis								! axes index variables

! Define the lattice <=> physical conversion factors
IF(domaintype .EQ. 0) THEN
        xcf_fine	= (0.5_lng*D)/(nx-1_lng)		! length conversion factor: x-direction
        ycf_fine	= (0.5_lng*D)/(ny-1_lng)		! length conversion factor: y-direction
ELSE
        ! begin Balaji added
        xcf_fine	= (1.0_lng*D*fractionDfine)/(nx_fine-1_lng - 2*gridRatio)	! length conversion factor: x-direction
        ycf_fine	= (1.0_lng*D*fractionDfine)/(ny_fine-1_lng - 2*gridRatio)	! length conversion factor: y-direction
        ! end Balaji added
ENDIF

zcf_fine 		= L/nz_fine					! length conversion factor: z-direction
tcf_fine 		= nuL_fine*((xcf_fine*xcf_fine)/nu)				! time conversion factor
dcf_fine 		= den/denL							! density conversion factor
vcf_fine 		= xcf_fine/tcf_fine							! velocity conversion factor
pcf_fine 		= cs*cs*vcf_fine*vcf_fine					! pressure conversion factor

! Initialize arrays
node_fine	= -99_lng							! node flag array
rDom_fine	= 0.0_dbl							! radius at each z-location
r_fine		= 0.0_dbl							! temporary radius array for entire computational domain
velDom_fine	= 0.0_dbl							! wall velocity at each z-location (global)
vel_fine	= 0.0_dbl							! wall velocity at each z-location (local)

! Check to ensure xcf_fine=ycf_fine=zcf_fine (LBM grid must be cubic)
IF((ABS(xcf_fine-ycf_fine) .GE. 1E-8) .OR. (ABS(xcf_fine-zcf_fine) .GE. 1E-8) .OR. (ABS(ycf_fine-zcf_fine) .GE. 1E-8)) THEN
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "Conversion factors not equal... Geometry_Setup.f90: Line 93."
  WRITE(1000,*) "xcf_fine=", xcf_fine, "ycf_fine=", ycf_fine, "zcf_fine=", zcf_fine
  WRITE(1000,*) "L=", L, "D=", D, " Fraction of dimeter for fine mesh = ", fractionDfine
  WRITE(1000,*) "nx_fine=", nx_fine, "ny_fine=", ny_fine, "nz_fine=", nz_fine
  CLOSE(1000)
  STOP
END IF


! IF CONDITION TO CHECK IF THE DOMAIN TO BE MODELLED IS FULL CYLINDER OR JUST A QUARTER OF A CYLINDER
IF(domaintype .EQ. 0) THEN 
      ! Fill out x,y,z arrays (local)
      DO i=0,nxSub_fine+1
        x_fine(i) = ((iMin - 1_lng) + (i-1_lng))*xcf_fine
      END DO
      
      DO j=0,nySub_fine+1
        y_fine(j) = ((jMin - 1_lng) + (j-1_lng))*ycf_fine
      END DO
      
      DO k=0,nzSub_fine+1
        z_fine(k) = (((kMin - 1_lng) + k) - 0.5_dbl)*zcf_fine
      END DO
      
      ! Fill out xx,yy,zz arrays (global)
      DO i=0,nx_fine+1
        xx_fine(i) = (i-1_lng)*xcf_fine
      END DO
      
      DO j=0,ny_fine+1
        yy_fine(j) = (j-1_lng)*ycf_fine
      END DO
      
      DO k=0,nz_fine+1
        zz_fine(k) = (k - 0.5_dbl)*zcf_fine
      END DO
      
      ! Center node locations
      Ci = 1	
      Cj = 1
      Ck = ANINT(0.5_dbl*nz_fine)

ELSE
      ! begin Balaji added 
      !INTEGER(lng) :: xaxis,yaxis								! axes index variables
      xaxis=ANINT(0.5_dbl*(nx_fine+1))
      yaxis=ANINT(0.5_dbl*(ny_fine+1))
      
      ! Fill out x,y,z arrays (local)
      DO i=0,nxSub_fine+1
        x_fine(i) = ((iMin_fine - 1_lng - (xaxis-1_lng)) + (i-1_lng))*xcf_fine
      END DO
      
      DO j=0,nySub_fine+1
        y_fine(j) = ((jMin_fine - 1_lng - (yaxis-1_lng)) + (j-1_lng))*ycf_fine
      END DO
      
      DO k=-gridRatio+1,nzSub_fine+gridRatio
!        z_fine(k) = (((kMin_fine - 1_lng) + k) - 0.5_dbl)*zcf_fine
        z_fine(k) = z(1) + (k-1)*zcf_fine
      END DO

      ! Fill out xx,yy,zz arrays (global)
      DO i=0,nx_fine+1
        xx_fine(i) = (i-1_lng-(xaxis-1_lng))*xcf_fine
      END DO
      
      DO j=0,ny_fine+1
        yy_fine(j) = (j-1_lng-(yaxis-1_lng))*ycf_fine
      END DO
      
      DO k=0,nz_fine+1
        zz_fine(k) = (k - 0.5_dbl)*zcf_fine
      END DO
      
      ! Center node locations
      ! Ci = xaxis
      ! Cj = yaxis
      ! Ck = ANINT(0.5_dbl*nz)
      ! end Balaji added 
ENDIF

IF(restart .EQV. .FALSE.) THEN
  ! Initialize the Geometry
  CALL AdvanceGeometry_fine
END IF

!------------------------------------------------
END SUBROUTINE Geometry_Setup_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AdvanceGeometry_fine												! advances the geometry in time
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

! Calculate the radius at the current time step
CALL BoundaryPosition_fine

! Calculate the velocity at boundary point
CALL BoundaryVelocity_fine

! Flag the fluid/solid nodes based on the new geometry
CALL SetNodes_fine

!------------------------------------------------
END SUBROUTINE AdvanceGeometry_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BoundaryPosition_fine		! Calculates the position of the wall at the current time step
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl) :: h1(0:nz_fine+1)				! Mode 1 (peristalsis)
REAL(dbl) :: h2(0:nz_fine+1)				! Mode 2	(segmental)
REAL(dbl) :: Ac, lambdaC, shiftC	! temporary variables for the cos slopes
REAL(dbl) :: time						! time
INTEGER(lng) :: i,j,ii,k			! indices

! Initialize Variables
time 	= 0.0_dbl						! time					
!rDom	= 0.5_dbl*D						! summed height
!h1 	= 0.5_dbl*D						! mode 1 height
!h2 	= 0.5_dbl*D						! mode 2 height
h1 	= 0.0_dbl						! mode 1 height
h2 	= 0.0_dbl						! mode 2 height
rDom_fine	= 0.0_dbl						! summed height

! Current Physical Time
time	= iter*tcf + subIter*tcf_fine

!------------------------- Mode 1 - peristalsis -----------------------------
DO i=0,nz_fine-1

  h1(i) 	= amp1*(COS(kw1*(zz_fine(i) - (s1*time)))) + (0.5_dbl*D - amp1)

END DO

! since PI cannot be stored exactly, the wavelength(s) does/do not EXACTLY span the domain...
! set h1(nz) to h1(0) and h1(nz+1) to h(1) to ensure periodicity
h1(nz_fine) 	= h1(0)
h1(nz_fine+1)= h1(1)

!------------------- Mode 2 - segmental contractions ------------------------

! Calculate the geometry for the first wave
! First Straight Piece
DO i=0,seg1L

  h2(i) = amp2*(COS(((2.0_dbl*PI)/Ts)*time)) + shift2
  
END DO

! Second Straight Piece
DO i=seg1R,seg2L

  h2(i) = amp2*(COS(((2.0_dbl*PI)/Ts)*(time-(Ts/2.0_dbl)))) + shift2
  
END DO

! Third Straight Piece
DO i=seg2R,nlambda2+1

  h2(i) = amp2*(COS(((2.0_dbl*PI)/Ts)*time)) + shift2
  
END DO

! First Cos Piece
Ac	= 0.5_dbl*(h2(seg1L)-h2(seg1R))
lambdaC	= 2.0_dbl*(zz(seg1L)-zz(seg1R))
shiftC	= 0.5_dbl*(h2(seg1L)+h2(seg1R))
DO i=seg1L+1,seg1R-1

  h2(i) = Ac*COS((2.0_dbl*PI/lambdaC)*(zz_fine(i)-zz_fine(seg1L))) + shiftC
  
END DO

! Second Cos Piece
Ac			= 0.5_dbl*(h2(seg2L)-h2(seg2R))
lambdaC	= 2.0_dbl*(zz_fine(seg2L)-zz_fine(seg2R))
shiftC	= 0.5_dbl*(h2(seg2L)+h2(seg2R))
DO i=seg2L+1,seg2R-1

  h2(i) = Ac*COS((2.0_dbl*PI/lambdaC)*(zz_fine(i)-zz_fine(seg2L))) + shiftC
  
END DO

! Repeat for the rest of the waves
DO j=1,(numw2-1)
  DO i=0,nlambda2+1

    ii = i + j*nlambda2
    h2(ii) = h2(i)

  END DO
END DO

! "fudging" to make sure that the whole domain is filled (and periodic) - more logic (and computational expense would be
! necessary to do this correctly: ideally, one would determine if an even or odd number of waves was specified
! and then work from either end, and meet in the middle to ensure a symetric domain...
!h2(nz-1:nz+1) = h2(1)

!----------------------------------------------------------------------------

!-------------------------------- Mode Sum  ---------------------------------

! Sum the modes in a weighted linear combination
DO i=0,nz_fine+1
  rDom_fine(i) = wc1*h1(i) + wc2*h2(i)
END DO

!----------------------------------------------------------------------------

! Fill out the local radius array
r_fine(0:nzSub_fine+1) = rDom_fine(kMin_fine-1:kMax_fine+1)


!------------------------------------------------
END SUBROUTINE BoundaryPosition_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BoundaryVelocity_fine	! defines the velocity of the solid boundaries (fills "ub", "vb", and "wb" arrays)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

REAL(dbl) :: v1(0:nz_fine+1), v2(0:nz_fine+1)	! velocity arrays for each mode
REAL(dbl) :: lambdaC						! wavelength of the cos segments (mode 2)
REAL(dbl) :: time							! time
INTEGER(lng) :: i,j,ii					! indices

! Initialize Variables
time		= 0.0_dbl						! time
velDom_fine	= 0.0_dbl						! summed velocity
v1		= 0.0_dbl						! mode 1 velocity
v2 		= 0.0_dbl						! mode 2 velocity				

! Current Physical Time
time = iter*tcf

!------------------------- Mode 1 - peristalsis -----------------------------
!DO i=1,nz
DO i=0,nz_fine-1 ! Balaji added to ensure periodicity just like in h1. 

  v1(i)	= kw1*s1*amp1*(SIN(kw1*(zz_fine(i) - (s1*time))))

END DO

! Balaji added
v1(nz_fine)=v1(0)
v1(nz_fine+1)=v1(1)
!----------------------------------------------------------------------------

!------------------- Mode 2 - segmental contractions  -----------------------

! Calculate the wall velocity for the first wave
! First Straight Piece
DO i=0,seg1L

  v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*time))*((2.0_dbl*PI)/Ts)
  
END DO

! Second Straight Piece
DO i=seg1R,seg2L

  v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*(time-(Ts/2.0_dbl))))*((2.0_dbl*PI)/Ts)
  
END DO

! Third Straight Piece
DO i=seg2R,nlambda2

  v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*time))*((2.0_dbl*PI)/Ts)
  
END DO

! First Cos Piece
lambdaC	= 2.0_dbl*(zz_fine(seg1L)-zz_fine(seg1R))
DO i=seg1L+1,seg1R-1

  v2(i) = (0.5_dbl*(v2(seg1L)-v2(seg1R)))*COS((2.0_dbl*PI/lambdaC)*(zz_fine(i)-zz_fine(seg1L))) &
        + (0.5_dbl*(v2(seg1L)+v2(seg1R)))
    
END DO

! Second Cos Piece
lambdaC	= 2.0_dbl*(zz_fine(seg2L)-zz(seg2R))
DO i=seg2L+1,seg2R-1

  v2(i) = (0.5_dbl*(v2(seg2L)-v2(seg2R)))*COS((2.0_dbl*PI/lambdaC)*(zz_fine(i)-zz_fine(seg2L))) &
        + (0.5_dbl*(v2(seg2L)+v2(seg2R)))

END DO

! Repeat for the rest of the waves
DO j=1,(numw2-1)
  DO i=1,nlambda2+1

    ii = i + j*nlambda2
    v2(ii) = v2(i)

  END DO
END DO

! "fudging" to make sure that the whole domain is filled (and periodic) - more logic (and computational expense would be
! necessary to do this correctly: ideally, one would determine if an even or odd number of waves was specified
! and then work from either end, and meet in the middle to ensure a symetric domain...
v2(nz_fine-1:nz_fine+1) = v2(1)

!----------------------------------------------------------------------------

!-------------------------------- Mode Sum  ---------------------------------

! Sum the modes in a weighted linear combination
DO i=0,nz_fine+1
  velDom_fine(i) = wc1*v1(i) + wc2*v2(i)
END DO

!----------------------------------------------------------------------------

! Fill out the local velocity array
vel_fine(0:nzSub_fine+1) = velDom_fine(kMin_fine-1:kMax_fine+1)/vcf_fine

!------------------------------------------------
END SUBROUTINE BoundaryVelocity_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SetNodes_fine					! defines the geometry via the "node" array of flags
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

INTEGER(lng)	:: i,j,k,m,iComm	! index variables
REAL(dbl)		:: rijk				! radius of the current node
REAL(dbl)      :: ubx,uby,ubz		! boundary velocity
INTEGER(lng) :: mpierr										! MPI standard error variable 

! Flag the interior nodes and give values to nodes that just came in
DO k=1,nzSub_fine
  DO j=1,nySub_fine
    DO i=1,nxSub_fine

      rijk = SQRT(x_fine(i)*x_fine(i) + y_fine(j)*y_fine(j))

      IF(rijk .LT. r_fine(k)) THEN

         IF( ((i .eq. 1) .or. (i .eq. nx_fine)) .or. ((j .eq. 1) .or. (j .eq. ny_fine)) ) THEN !Trying to find the outermost node on the fine mesh, set that as COARSEMESH
            node_fine(i,j,k) = COARSEMESH !No computations to be carried out in these nodes

         ELSE

            IF(node_fine(i,j,k) .EQ. SOLID) THEN   ! just came into the domain
               ! calculate the wall velocity (boundary)
               
               ubx = vel_fine(k)*(x_fine(i)/rijk)
               uby = vel_fine(k)*(y_fine(j)/rijk)
               ubz = 0.0_dbl
               
               CALL SetProperties_fine(i,j,k,ubx,uby,ubz)

            END IF
	
            node_fine(i,j,k) = FLUID ! reset the SOLID node that just came in to FLUID
            
         END IF

      ELSE

        node_fine(i,j,k) = SOLID        ! if rijk is GT r(k) then it's a SOLID node

      END IF

    END DO
  END DO
END DO

! Loop through the phantom nodes, and set the entity, but do not give values
! YZ Faces
DO iComm=1,2

  i = YZ_RecvIndex_fine(OppCommDir(iComm))															! i index of the phantom nodes
 	
  DO j=0,nySub_fine+1_lng

    rijk = SQRT(x_fine(i)*x_fine(i) + y_fine(j)*y_fine(j))

    DO k=0,nzSub_fine+1_lng

      IF(rijk .LT. r_fine(k)) THEN

         IF(rijk .GT. (0.5*fractionDfine*D - 0.1*ycf_fine) ) THEN !Trying to find the outermost node on the fine mesh, set that as COARSEMESH

            node_fine(i,j,k) = COARSEMESH !No computations to be carried out in these nodes

         ELSE

            node_fine(i,j,k) = FLUID																		! set the SOLID node that just came in to FLUID
         END IF
      ELSE
        node_fine(i,j,k) = SOLID																		! if rijk is GT r_fine(k) then it's a SOLID node

      END IF
        
    END DO
  END DO

END DO

! ZX Faces
DO iComm=3,4

  j = ZX_RecvIndex_fine(OppCommDir(iComm))															! j index of the phantom nodes

  DO i=0,nxSub_fine+1_lng

    rijk = SQRT(x_fine(i)*x_fine(i) + y_fine(j)*y_fine(j))

    DO k=0,nzSub_fine+1_lng

      IF(rijk .LT. r_fine(k)) THEN

         IF(rijk .GT. (0.5*fractionDfine*D - 0.1*ycf_fine) ) THEN !Trying to find the outermost node on the fine mesh, set that as COARSEMESH

            node_fine(i,j,k) = COARSEMESH !No computations to be carried out in these nodes

         ELSE

            node_fine(i,j,k) = FLUID																		! set the SOLID node that just came in to FLUID

         END IF
      ELSE
        node_fine(i,j,k) = SOLID																		! if rijk is GT r_fine(k) then it's a SOLID node

      END IF
        
    END DO

  END DO

END DO

! XY Faces
DO iComm=5,6

  k = XY_RecvIndex_fine(OppCommDir(iComm))															! k index of the phantom nodes

  DO j=0,nySub_fine+1_lng
    DO i=0,nxSub_fine+1_lng

      rijk = SQRT(x_fine(i)*x_fine(i) + y_fine(j)*y_fine(j))

      IF(rijk .LT. r_fine(k)) THEN

         IF(rijk .GT. (0.5*fractionDfine*D - 0.1*ycf_fine) ) THEN !Trying to find the outermost node on the fine mesh, set that as COARSEMESH
            
            node_fine(i,j,k) = COARSEMESH !No computations to be carried out in these nodes

         ELSE
         
            node_fine(i,j,k) = FLUID																		! set the SOLID node that just came in to FLUID
         END IF
         
      ELSE
        node_fine(i,j,k) = SOLID																		! if rijk is GT r_fine(k) then it's a SOLID node
      END IF

    END DO
  END DO

END DO

!------------------------------------------------
END SUBROUTINE SetNodes_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SetProperties_fine(i,j,k,ubx,uby,ubz)	! give properties to nodes that just came into the fluid domain (uncovered)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

INTEGER(lng), INTENT(IN) :: i,j,k				! current node location
REAL(dbl), INTENT(IN) :: ubx,uby,ubz			! velocity of the boundary
INTEGER(lng)	:: m,ii,jj,kk						! index variables
INTEGER(lng)	:: numFLUIDs						! number of fluid nodes
REAL(dbl)		:: rhoSum, rhoTemp				! sum of the densities of the neighboring fluid nodes, pre-set density
REAL(dbl)		:: feq								! equilibrium distribution function
CHARACTER(7)	:: iter_char						! iteration stored as a character

! initialize the sum of surrounding densities
rhoSum = 0.0_dbl
numFLUIDs = 0_lng

! calculate the average density of the current node's neighbors
DO m=1,NumDistDirs

  ii = i + ex(m)
  jj = j + ey(m)
  kk = k + ez(m)

  IF(((ii .GE. 0) .AND. (ii .LE. nxSub_fine+1_lng)) .AND.	&
     ((jj .GE. 0) .AND. (jj .LE. nySub_fine+1_lng)) .AND.	&
     ((kk .GE. 0) .AND. (kk .LE. nzSub_fine+1_lng))) THEN

    IF(node_fine(ii,jj,kk) .EQ. FLUID) THEN
      rhoSum = rhoSum + rho_fine(ii,jj,kk)
      numFLUIDs = numFLUIDs + 1_lng     
    END IF       

  END IF

END DO

! This should rarely happen...
IF(numFLUIDs .NE. 0_lng) THEN

  rho_fine(i,j,k) = rhoSum/numFLUIDs

ELSE

  rho_fine(i,j,k) = denL

END IF

! velocity and scalar (use boundary conditions)
u_fine(i,j,k) 	= ubx									! wall velocity			
v_fine(i,j,k) 	= uby														
w_fine(i,j,k) 	= ubz
phi_fine(i,j,k)	= phiWall							! scalar			

! distribution functions (set to equilibrium)
DO m=0,NumDistDirs
  CALL Equilibrium_LOCAL_fine(m,rho_fine(i,j,k),u_fine(i,j,k),v_fine(i,j,k),w_fine(i,j,k),feq)	! distribution functions
  f_fine(m,i,j,k) = feq
END DO

!------------------------------------------------
END SUBROUTINE SetProperties_Fine
!------------------------------------------------

!================================================
END MODULE Geometry_Fine
!================================================
