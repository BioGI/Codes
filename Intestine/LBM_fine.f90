!==================================================================================================
MODULE LBM_fine				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE Setup_fine
USE ICBC_fine

IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE LBM_Setup_fine	! setup the LBM simulation
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Initialize variables and arrays
f_fine		= 0.0_dbl		! distribution functions
fplus_fine	= 0.0_dbl		! post-collision distribution functions
u_fine  	= 0.0_dbl		! x-velocity
v_fine     = 0.0_dbl		! y-velocity
w_fine     = 0.0_dbl		! z-velocity
rho_fine   = 0.0_lng		! density

! Define other simulation parameters
nuL_fine   		= (2.0_dbl*tau_fine - 1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
oneOVERtau_fine 	= 1.0_dbl/tau_fine			! reciprical of tau

! Initialize timestep
iter_fine = 0_lng												! intialize the starting timestep to 0 - will get reset in 'ICs' in ICBCM.f90

! Calculate feq for initial condition
CALL Equilibrium_fine

!------------------------------------------------
END SUBROUTINE LBM_Setup_Fine
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE Equilibrium_fine		! calculate the equilibrium distribution function and set f to feq (initial condition)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,m					! index variables
REAL(dbl)		:: uu,ue,ve,we,Usum		! precalculated quantities for use in the feq equation
REAL(dbl)		:: feq						! equilibrium distribution function

! Balaji modified to change indices form 0 to nzSub+1
DO k=1,nzSub_fine+0
  DO j=1,nySub_fine+0
    DO i=1,nxSub_fine+0
!DO k=0,nzSub+1
!  DO j=0,nySub+1
!    DO i=0,nxSub+1

      IF( (node_fine(i,j,k) .EQ. FLUID) ) THEN
      
        uu = u_fine(i,j,k)*u_fine(i,j,k) + v_fine(i,j,k)*v_fine(i,j,k) + w_fine(i,j,k)*w_fine(i,j,k)						! u . u
      
        DO m=0,NumDistDirs
        
          ue	= u_fine(i,j,k)*ex(m)																			! u . e
          ve	= v_fine(i,j,k)*ey(m)																			! v . e
          we	= w_fine(i,j,k)*ez(m)																			! w . e

          Usum	= ue + ve + we																				! U . e
        
          feq = (wt(m)*rho_fine(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function

          f_fine(m,i,j,k) = feq    

        END DO

      END IF
      
    END DO
  END DO
END DO

!------------------------------------------------
END SUBROUTINE Equilibrium_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Collision_Fine		! calculates equilibrium distribution function AND collision step for each node
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,m					! index variables
REAL(dbl)		:: UU,ue,ve,we,Usum		! precalculated quantities for use in the feq equation
REAL(dbl)		:: feq						! equilibrium distribution function

! Balaji modified to change indices form 0 to nzSub+1
DO k=1,nzSub_fine+0
  DO j=1,nySub_fine+0
    DO i=1,nxSub_fine+0
!DO k=0,nzSub+1
!  DO j=0,nySub+1
!    DO i=0,nxSub+1

      IF(node_fine(i,j,k) .EQ. FLUID) THEN

        UU = u_fine(i,j,k)*u_fine(i,j,k) + v_fine(i,j,k)*v_fine(i,j,k) + w_fine(i,j,k)*w_fine(i,j,k)						! U . U
      
        DO m=0,NumDistDirs
        
          ue	= u_fine(i,j,k)*ex(m)																			! u . e
          ve	= v_fine(i,j,k)*ey(m)																			! v . e
          we	= w_fine(i,j,k)*ez(m)																			! w . e

          Usum	= ue + ve + we																				! U . e
        
          feq	= (wt(m)*rho_fine(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5*Usum*Usum - 1.5*uu)	! equilibrium distribution function
          f_fine(m,i,j,k)		= f_fine(m,i,j,k) - oneOVERtau_fine*(f_fine(m,i,j,k) - feq)			! collision
        
        END DO 

      END IF
      
    END DO
  END DO
END DO


!------------------------------------------------
END SUBROUTINE Collision_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Stream_fine	! stream the distribution functions between neighboring nodes (stream - using Lallemand 2nd order moving BB)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: fbb								! bounced back distribution function

fplus_fine = f_fine										! store the post-collision distribution function

! interior nodes (away from other subdomains)
DO k=2,nzSub_fine-1
  DO j=2,nySub_fine-1
    DO i=2,nxSub_fine-1

      IF(node_fine(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)
    
          IF(node_fine(im1,jm1,km1) .EQ. FLUID) THEN
             f_fine(m,i,j,k) = fplus_fine(m,im1,jm1,km1)             
          ELSE IF(node_fine(im1,jm1,km1) .EQ. COARSEMESH) THEN
             f_fine(m,i,j,k) = fplus_fine(m,im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN														 ! macro- boundary
            ! CALL BounceBack2(m,i,j,k,im1,jm1,km1,fbb)					  										 ! implement the bounceback BCs [MODULE: ICBC]
	    ! Balaji added after commenting out the earlier method
            CALL BounceBack2New_fine(m,i,j,k,im1,jm1,km1,fbb)					  										! implement the bounceback BCs [MODULE: ICBC]
            f_fine(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 160: node_fine(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x_fine(i)=",x_fine(i),"y_fine(j)=",y_fine(j),"z_fine(k)=",z_fine(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x_fine(im1)=",x_fine(im1),"y_fine(jm1)=",y_fine(jm1),"z_fine(km1)=",z_fine(km1)
            WRITE(1000,*) "node_fine(i,j,k)=",node_fine(i,j,k)
            WRITE(1000,*) "node_fine(im1,jm1,km1)=",node_fine(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF
 
    END DO
  END DO
END DO

! XY faces
DO k=1,nzSub_fine,(nzSub_fine-1)
  DO j=1,nySub_fine
    DO i=1,nxSub_fine

      IF(node_fine(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

          !IF(km1.eq.0) km1=nzSub ! Balaji added
          !IF(km1.eq.nzSub+1) km1=1 ! Balaji added

          IF(node_fine(im1,jm1,km1) .EQ. FLUID) THEN 
            f_fine(m,i,j,k) = fplus_fine(m,im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. COARSEMESH) THEN
            f_fine(m,i,j,k) = fplus_fine(m,im1,jm1,km1)            
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN
            ! macro- boundary
            CALL BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)			  												! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f_fine(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM_fine.f90 at Line 208: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x_fine(i)=",x_fine(i),"y_fine(j)=",y_fine(j),"z_fine(k)=",z_fine(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x_fine(im1)=",x_fine(im1),"y_fine(jm1)=",y_fine(jm1),"z_fine(km1)=",z_fine(km1)
            WRITE(1000,*) "node_fine(i,j,k)=",node_fine(i,j,k)
            WRITE(1000,*) "node_fine(im1,jm1,km1)=",node_fine(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF

    END DO
  END DO
END DO

! XZ faces
DO j=1,nySub_fine,(nySub_fine-1)
  DO k=1,nzSub_fine
    DO i=1,nxSub_fine

      IF(node_fine(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node_fine
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

        
          IF(node_fine(im1,jm1,km1) .EQ. FLUID) THEN
            f_fine(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. COARSEMESH) THEN 
            f_fine(m,i,j,k) = fplus_fine(m,im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
            CALL BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f_fine(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM_fine.f90 at Line 253: node_fine(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x_fine(i)=",x_fine(i),"y_fine(j)=",y_fine(j),"z_fine(k)=",z_fine(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x_fine(im1)=",x_fine(im1),"y_fine(jm1)=",y_fine(jm1),"z_fine(km1)=",z_fine(km1)
            WRITE(1000,*) "node_fine(i,j,k)=",node_fine(i,j,k)
            WRITE(1000,*) "node_fine(im1,jm1,km1)=",node_fine(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF

    END DO
  END DO
END DO

! YZ faces
DO i=1,nxSub_fine,(nxSub_fine-1)
  DO k=1,nzSub_fine
    DO j=1,nySub_fine

      IF(node_fine(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node_fine
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)
        
          IF(node_fine(im1,jm1,km1) .EQ. FLUID) THEN 
            f_fine(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. COARSEMESH) THEN 
            f_fine(m,i,j,k) = fplus_fine(m,im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
            CALL BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f_fine(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM_finef90 at Line 297: node_fine(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x_fine(i)=",x_fine(i),"y_fine(j)=",y_fine(j),"z_fine(k)=",z_fine(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x_fine(im1)=",x_fine(im1),"y_fine(jm1)=",y_fine(jm1),"z_fine(km1)=",z_fine(km1)
            WRITE(1000,*) "node_fine(i,j,k)=",node_fine(i,j,k)
            WRITE(1000,*) "node_fine(im1,jm1,km1)=",node_fine(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF

    END DO
  END DO
END DO

!------------------------------------------------
END SUBROUTINE Stream_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Macro_fine	! calculate the macroscopic quantities
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m						! index variables

INTEGER(lng) :: ii,jj,kk

! Balaji modified to include 0 to nzSub+1
DO k=1,nzSub_fine
  DO j=1,nySub_fine
    DO i=1,nxSub_fine
!DO k=0,nzSub+1
!  DO j=0,nySub+1
!    DO i=0,nxSub+1
      
      IF(node_fine(i,j,k) .EQ. FLUID) THEN

        ! initialize arrays
        rho_fine(i,j,k)		= 0.0_dbl								! density
        u_fine(i,j,k)		= 0.0_dbl								! x-velocity
        v_fine(i,j,k)		= 0.0_dbl								! y-velocity
        w_fine(i,j,k)		= 0.0_dbl								! z-velocity     

        DO m=0,NumDistDirs  
          rho_fine(i,j,k)	= rho_fine(i,j,k) + f_fine(m,i,j,k)			! density
          u_fine(i,j,k)	= u_fine(i,j,k)   + f_fine(m,i,j,k)*ex(m)	! x-velocity
          v_fine(i,j,k)	= v_fine(i,j,k)   + f_fine(m,i,j,k)*ey(m)	! y-velocity
          w_fine(i,j,k)	= w_fine(i,j,k)   + f_fine(m,i,j,k)*ez(m)	! z-velocity
        END DO

        IF(rho_fine(i,j,k) .NE. 0) THEN
          u_fine(i,j,k) = u_fine(i,j,k)/rho_fine(i,j,k)					! x-velocity
          v_fine(i,j,k) = v_fine(i,j,k)/rho_fine(i,j,k)					! y-velocity
          w_fine(i,j,k) = w_fine(i,j,k)/rho_fine(i,j,k)					! z-velocity
        ELSE          

          OPEN(6678,FILE='error.'//sub//'.txt')
          WRITE(6678,*) 'rho_fine(i,j,k) = 0: Line 362 in Macro in LBM_fine.f90'
          WRITE(6678,*) 'iter_fine', iter_fine
          WRITE(6678,*) 'i,j,k:', i,j,k
          WRITE(6678,*) 'node_fine(i,j,k)', node_fine(i,j,k)
          WRITE(6678,*) 'rho_fine(i,j,k)', rho_fine(i,j,k)
          WRITE(6678,*)
          WRITE(6678,*)
          DO m=1,NumDistDirs
            ii = i + ex(m)
            jj = j + ey(m)
            kk = k + ez(m)       
            WRITE(6678,*) 'ii,jj,kk:', ii,jj,kk
            WRITE(6678,*) 'node(ii,jj,kk)', node_fine(ii,jj,kk)
            WRITE(6678,*) 'rho(ii,jj,kk)', rho_fine(ii,jj,kk)
            WRITE(6678,*)
          END DO
          CLOSE(6678)

          OPEN(1001,FILE='rhoMacro.'//sub//'.dat')
          WRITE(1001,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
          WRITE(1001,'(8E15.5,I6)') x_fine(i), y_fine(j), z_fine(k), u_fine(i,j,k), v_fine(i,j,k), w_fine(i,j,k), (rho_fine(i,j,k)-denL)*dcf_fine*pcf_fine, phi_fine(i,j,k), node_fine(i,j,k)
          CLOSE(1001)

          CALL PrintFieldsTEST_fine										! output the velocity, density, and scalar fields [MODULE: Output]

          STOP

        END IF

      ELSE IF (node_fine(i,j,k) .EQ. SOLID) THEN
      
        rho_fine(i,j,k)	= denL									! density (zero gauge pressure)
        u_fine(i,j,k)		= 0.0_dbl								! x-velocity
        v_fine(i,j,k)		= 0.0_dbl								! y-velocity
        w_fine(i,j,k)		= 0.0_dbl								! z-velocity
        phi_fine(i,j,k)	= phiWall								! scalar
        
	

      END IF

    END DO
  END DO
END DO   

!------------------------------------------------
END SUBROUTINE Macro_Fine
!------------------------------------------------

FUNCTION lowerCoarseXindex(xf)
  !Returns the lower coarse X index
  INTEGER :: xf
  INTEGER :: lowerCoarseXindex
  lowerCoarseXindex = 45 + FLOOR((xf - 1)/dble(gridRatio))
  RETURN 
END FUNCTION lowerCoarseXindex

FUNCTION lowerCoarseYindex(yf)
  !Returns the lower coarse Y index
  INTEGER :: yf
  INTEGER :: lowerCoarseYindex
  lowerCoarseYindex = 45 + FLOOR((yf - 1)/dble(gridRatio))
  RETURN 
END FUNCTION lowerCoarseYindex

FUNCTION closestCoarseZindex(zf)
  !Returns the lower coarse Z index
  REAL(dbl) :: zf
  INTEGER :: closestCoarseZindex
  closestCoarseZindex = ANINT( ((zf-z(1))/zcf) ) + 1  
  RETURN 
END FUNCTION closestCoarseZindex

FUNCTION closestFineIindex(x)
  !Returns the closest fine mesh I index
  REAL(dbl) :: x
  INTEGER :: closestFineIindex
  closestFineIindex = ANINT( ( (x-x_fine(1))/xcf_fine) ) + 1
  RETURN 
END FUNCTION closestFineIindex

FUNCTION closestFineJindex(y)
  !Returns the closest fine mesh I index
  REAL(dbl) :: y
  INTEGER :: closestFineJindex
  closestFineJindex = ANINT( ( (y-y_fine(1))/ycf_fine) ) + 1
  RETURN 
END FUNCTION closestFineJindex

FUNCTION closestFineKindex(z)
  !Returns the closest fine mesh K index
  REAL(dbl) :: z
  INTEGER :: closestFineKindex
  closestFineKindex = ANINT( ((z - z_fine(1))/zcf_fine) ) + 1
  RETURN 
END FUNCTION closestFineKindex


SUBROUTINE ComputeEquilibriumForFineGrid
  !!! Compute the equilibrium distribution function at the coarse grid interface for the fine grid 

  INTEGER   :: i,j,k,m
  REAL(dbl) :: uu,ue,ve,we,Usum		! precalculated quantities for use in the feq equation
  REAL(dbl) :: feq			! equilibrium distribution function

  !Do the bottom and top x-z planes first
  do k = 1, nzSub
     do i = 44, 59
        IF(node(i,45,k) .EQ. FLUID) THEN
           
           uu = u(i,45,k)*u(i,45,k) + v(i,45,k)*v(i,45,k) + w(i,45,k)*w(i,45,k)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u(i,45,k)*ex(m)		! u . e
              ve	= v(i,45,k)*ey(m)		! v . e
              we	= w(i,45,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(i,45,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function

              feqFF_bottomXZ(m,i,k) = feq

           END DO
           
        END IF

        IF(node(i,57,k) .EQ. FLUID) THEN
           
           uu = u(i,57,k)*u(i,57,k) + v(i,57,k)*v(i,57,k) + w(i,57,k)*w(i,57,k)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u(i,57,k)*ex(m)		! u . e
              ve	= v(i,57,k)*ey(m)		! v . e
              we	= w(i,57,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(i,57,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFF_topXZ(m,i,k) = feq
              
           END DO
           
        END IF
        
     end do
  end do

  !Fill in the remaining points on the front and back planes
  do k = 1, nzSub
     do j = 44, 58
        IF(node(45,j,k) .EQ. FLUID) THEN
           
           uu = u(45,j,k)*u(45,j,k) + v(45,j,k)*v(45,j,k) + w(45,j,k)*w(45,j,k)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u(45,j,k)*ex(m)		! u . e
              ve	= v(45,j,k)*ey(m)		! v . e
              we	= w(45,j,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(45,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFF_frontYZ(m,j,k) = feq
              
           END DO
           
        END IF

        IF(node(57,j,k) .EQ. FLUID) THEN
           
           uu = u(57,j,k)*u(57,j,k) + v(57,j,k)*v(57,j,k) + w(57,j,k)*w(57,j,k)	! u . u
           DO m=0,NumDistDirs
              
              ue	= u(57,j,k)*ex(m)		! u . e
              ve	= v(57,j,k)*ey(m)		! v . e
              we	= w(57,j,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(57,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFF_backYZ(m,j,k) = feq
              
           END DO
           
        END IF
        
     end do
  end do
  
END SUBROUTINE ComputeEquilibriumForFineGrid

SUBROUTINE XYSpatialInterpolateBufferToFineGrid    ! Interpolate required variables on the two buffer rows to fine grid

  IMPLICIT NONE
  INTEGER :: i,j,k,m
  REAL(dbl) :: xInterp, yInterp, zInterp
  INTEGER :: lCxIndex, lCzIndex, lCyIndex, lFzIndex
  REAL(dbl) :: f1,f2,f3,f4              ! Temporary variables for interpolation.
  
  !Do the bottom and top x-z planes first
  !Only x  - interpolation 
  do k=1, gridRatio+1, gridRatio
     do i=1,nxSub_fine

        lCxIndex = lowerCoarseXindex(i)  ! Lower Coarse x Index
        
        xInterp = dble( MODULO(i-1, gridRatio) ) / dble(gridRatio)
        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        do m=0,NumDistDirs	  

           lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
           fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step           
           f1 =  feqFF_bottomXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex-1,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_bottomXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_bottomXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+1,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_bottomXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+2,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step

           lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index           
           fCtoF_bottomXZ(m,1,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_bottomXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_bottomXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_bottomXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_bottomXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex-1,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_bottomXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_bottomXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+1,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_bottomXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+2,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_bottomXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step

           lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
           fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_topXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex-1,57,lCzIndex) - feqFF_topXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_topXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex,57,lCzIndex) - feqFF_topXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_topXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+1,57,lCzIndex) - feqFF_topXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_topXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+2,57,lCzIndex) - feqFF_topXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_topXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step

           lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index           
           fCtoF_topXZ(m,1,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_topXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_topXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_topXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_topXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex-1,57,lCzIndex) - feqFF_topXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_topXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex,57,lCzIndex) - feqFF_topXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_topXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+1,57,lCzIndex) - feqFF_topXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_topXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+2,57,lCzIndex) - feqFF_topXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_topXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step

        end do

        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        dsCtoF_bottomXZ(:,1,i,k) = dsCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
        dsCtoF_bottomXZ(:,2,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
        f1 =  rho(lCxIndex-1,45,lCzIndex)
        f2 =  rho(lCxIndex,45,lCzIndex) 
        f3 =  rho(lCxIndex+1,45,lCzIndex)
        f4 =  rho(lCxIndex+2,45,lCzIndex)
        dsCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(lCxIndex-1,45,lCzIndex)
        f2 =  phi(lCxIndex,45,lCzIndex) 
        f3 =  phi(lCxIndex+1,45,lCzIndex)
        f4 =  phi(lCxIndex+2,45,lCzIndex)
        dsCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step

        lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index           
        dsCtoF_bottomXZ(:,1,i,nzSub_fine-gridRatio+1-k+1) = dsCtoF_bottomXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
        dsCtoF_bottomXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) = dsCtoF_bottomXZ(:,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
        f1 =  rho(lCxIndex-1,45,lCzIndex)
        f2 =  rho(lCxIndex,45,lCzIndex)
        f3 =  rho(lCxIndex+1,45,lCzIndex)
        f4 =  rho(lCxIndex+2,45,lCzIndex)
        dsCtoF_bottomXZ(1,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(lCxIndex-1,45,lCzIndex)
        f2 =  phi(lCxIndex,45,lCzIndex)
        f3 =  phi(lCxIndex+1,45,lCzIndex)
        f4 =  phi(lCxIndex+2,45,lCzIndex)
        dsCtoF_bottomXZ(2,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time ste        
        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        dsCtoF_topXZ(:,1,i,k) = dsCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
        dsCtoF_topXZ(:,2,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
        f1 =  rho(lCxIndex-1,57,lCzIndex) 
        f2 =  rho(lCxIndex,57,lCzIndex) 
        f3 =  rho(lCxIndex+1,57,lCzIndex) 
        f4 =  rho(lCxIndex+2,57,lCzIndex) 
        dsCtoF_topXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(lCxIndex-1,57,lCzIndex) 
        f2 =  phi(lCxIndex,57,lCzIndex) 
        f3 =  phi(lCxIndex+1,57,lCzIndex) 
        f4 =  phi(lCxIndex+2,57,lCzIndex) 
        dsCtoF_topXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
        
        lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index           
        dsCtoF_topXZ(:,1,i,nzSub_fine-gridRatio+1-k+1) = dsCtoF_topXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
        dsCtoF_topXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) = dsCtoF_topXZ(:,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
        f1 =  rho(lCxIndex-1,57,lCzIndex) 
        f2 =  rho(lCxIndex,57,lCzIndex) 
        f3 =  rho(lCxIndex+1,57,lCzIndex) 
        f4 =  rho(lCxIndex+2,57,lCzIndex) 
        dsCtoF_topXZ(1,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(lCxIndex-1,57,lCzIndex) 
        f2 =  phi(lCxIndex,57,lCzIndex) 
        f3 =  phi(lCxIndex+1,57,lCzIndex) 
        f4 =  phi(lCxIndex+2,57,lCzIndex) 
        dsCtoF_topXZ(2,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
           
     end do
  end do
 
  !Fill out the remaining points on the front and back y-z planes
  !Only y-interpolation
  do k=1, gridRatio+1, gridRatio
     do j=2,nySub_fine-1

        lCyIndex = lowerCoarseYindex(j)  ! Lower Coarse x Index
        
        yInterp = dble( MODULO(j-1, gridRatio) ) / dble(gridRatio)
        
        do m=0,NumDistDirs	  
           
           lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
           fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_frontYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex-1,lCzIndex) - feqFF_frontYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_frontYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex,lCzIndex) - feqFF_frontYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_frontYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex+1,lCzIndex) - feqFF_frontYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_frontYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex+2,lCzIndex) - feqFF_frontYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step

           lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index - No interpolation in z
           fCtoF_frontYZ(m,1,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_frontYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_frontYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_frontYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_frontYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex-1,lCzIndex) - feqFF_frontYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_frontYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex,lCzIndex) - feqFF_frontYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_frontYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex+1,lCzIndex) - feqFF_frontYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_frontYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex+2,lCzIndex) - feqFF_frontYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_frontYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step

           lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
           fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_backYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(57,lCyIndex-1,lCzIndex),node(57,lCyIndex,lCzIndex),node(57,lCyIndex+1,lCzIndex),node(57,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step

           lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index - No interpolation in z
           fCtoF_backYZ(m,1,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_backYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_backYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_backYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_backYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(57,lCyIndex-1,lCzIndex),node(57,lCyIndex,lCzIndex),node(57,lCyIndex+1,lCzIndex),node(57,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
           
        end do

        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
        dsCtoF_frontYZ(:,1,j,k) = dsCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
        dsCtoF_frontYZ(:,2,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
        f1 =  rho(45,lCyIndex-1,lCzIndex) 
        f2 =  rho(45,lCyIndex,lCzIndex) 
        f3 =  rho(45,lCyIndex+1,lCzIndex) 
        f4 =  rho(45,lCyIndex+2,lCzIndex) 
        dsCtoF_frontYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(45,lCyIndex-1,lCzIndex) 
        f2 =  phi(45,lCyIndex,lCzIndex) 
        f3 =  phi(45,lCyIndex+1,lCzIndex) 
        f4 =  phi(45,lCyIndex+2,lCzIndex) 
        dsCtoF_frontYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        
        
        lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index - No interpolation in z
        dsCtoF_frontYZ(:,1,j,nzSub_fine-gridRatio+1-k+1) = dsCtoF_frontYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
        dsCtoF_frontYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) = dsCtoF_frontYZ(:,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
        f1 =  rho(45,lCyIndex-1,lCzIndex) 
        f2 =  rho(45,lCyIndex,lCzIndex) 
        f3 =  rho(45,lCyIndex+1,lCzIndex) 
        f4 =  rho(45,lCyIndex+2,lCzIndex) 
        dsCtoF_frontYZ(1,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(45,lCyIndex-1,lCzIndex) 
        f2 =  phi(45,lCyIndex,lCzIndex) 
        f3 =  phi(45,lCyIndex+1,lCzIndex) 
        f4 =  phi(45,lCyIndex+2,lCzIndex) 
        dsCtoF_frontYZ(2,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        
        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
        dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
        dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
        f1 =  rho(57,lCyIndex-1,lCzIndex) 
        f2 =  rho(57,lCyIndex,lCzIndex) 
        f3 =  rho(57,lCyIndex+1,lCzIndex) 
        f4 =  rho(57,lCyIndex+2,lCzIndex) 
        dsCtoF_backYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(57,lCyIndex-1,lCzIndex) 
        f2 =  phi(57,lCyIndex,lCzIndex) 
        f3 =  phi(57,lCyIndex+1,lCzIndex) 
        f4 =  phi(57,lCyIndex+2,lCzIndex) 
        dsCtoF_backYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        
        lCzIndex = closestCoarseZindex(z_fine(nzSub_fine-gridRatio+1-k+1))  ! Lower Coarse z Index - No interpolation in z
        dsCtoF_backYZ(:,1,j,nzSub_fine-gridRatio+1-k+1) = dsCtoF_backYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
        dsCtoF_backYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) = dsCtoF_backYZ(:,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
        f1 =  rho(57,lCyIndex-1,lCzIndex) 
        f2 =  rho(57,lCyIndex,lCzIndex) 
        f3 =  rho(57,lCyIndex+1,lCzIndex) 
        f4 =  rho(57,lCyIndex+2,lCzIndex) 
        dsCtoF_backYZ(1,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(57,lCyIndex-1,lCzIndex) 
        f2 =  phi(57,lCyIndex,lCzIndex) 
        f3 =  phi(57,lCyIndex+1,lCzIndex) 
        f4 =  phi(57,lCyIndex+2,lCzIndex) 
        dsCtoF_backYZ(2,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
     end do
  end do
  
END SUBROUTINE XYSpatialInterpolateBufferToFineGrid

SUBROUTINE XYSpatialInterpolateInternalNodesToFineGrid    ! Interpolate required variables on the two buffer rows to fine grid

  IMPLICIT NONE
  INTEGER :: i,j,k,m
  REAL(dbl) :: xInterp, yInterp, zInterp
  INTEGER :: lCxIndex, lCzIndex, lCyIndex, lFzIndex
  REAL(dbl) :: f1,f2,f3,f4              ! Temporary variables for interpolation.
  
  !Do the bottom and top x-z planes first
  !Only x  - interpolation 
  do k=2*gridRatio+1,nzSub_fine-3*gridRatio+1, gridRatio
     lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
     do i=1,nxSub_fine

        lCxIndex = lowerCoarseXindex(i)  ! Lower Coarse x Index
        xInterp = dble( MODULO(i-1, gridRatio) ) / dble(gridRatio)        

        do m=0,NumDistDirs	  
           fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_bottomXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex-1,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_bottomXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_bottomXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+1,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_bottomXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+2,45,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+2,lCzIndex))
           fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step

           fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_topXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex-1,57,lCzIndex) - feqFF_topXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_topXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex,57,lCzIndex) - feqFF_topXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_topXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+1,57,lCzIndex) - feqFF_topXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_topXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,lCxIndex+2,57,lCzIndex) - feqFF_topXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_topXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
        end do

        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        dsCtoF_bottomXZ(:,1,i,k) = dsCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
        dsCtoF_bottomXZ(:,2,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
        f1 =  rho(lCxIndex-1,45,lCzIndex)
        f2 =  rho(lCxIndex,45,lCzIndex) 
        f3 =  rho(lCxIndex+1,45,lCzIndex)
        f4 =  rho(lCxIndex+2,45,lCzIndex)
        dsCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(lCxIndex-1,45,lCzIndex)
        f2 =  phi(lCxIndex,45,lCzIndex) 
        f3 =  phi(lCxIndex+1,45,lCzIndex)
        f4 =  phi(lCxIndex+2,45,lCzIndex)
        dsCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
        
        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        dsCtoF_topXZ(:,1,i,k) = dsCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
        dsCtoF_topXZ(:,2,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
        f1 =  rho(lCxIndex-1,57,lCzIndex) 
        f2 =  rho(lCxIndex,57,lCzIndex) 
        f3 =  rho(lCxIndex+1,57,lCzIndex) 
        f4 =  rho(lCxIndex+2,57,lCzIndex) 
        dsCtoF_topXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(lCxIndex-1,57,lCzIndex) 
        f2 =  phi(lCxIndex,57,lCzIndex) 
        f3 =  phi(lCxIndex+1,57,lCzIndex) 
        f4 =  phi(lCxIndex+2,57,lCzIndex) 
        dsCtoF_topXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step

     end do
  end do
 
  !Fill out the remaining points on the front and back y-z planes
  !Only y-interpolation
  do k=2*gridRatio+1,nzSub_fine-3*gridRatio+1, gridRatio
     do j=2,nySub_fine-1

        lCyIndex = lowerCoarseYindex(j)  ! Lower Coarse x Index
        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
        
        yInterp = dble( MODULO(j-1, gridRatio) ) / dble(gridRatio)
        
        do m=0,NumDistDirs	  
           fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_frontYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex-1,lCzIndex) - feqFF_frontYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_frontYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex,lCzIndex) - feqFF_frontYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_frontYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex+1,lCzIndex) - feqFF_frontYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_frontYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,45,lCyIndex+2,lCzIndex) - feqFF_frontYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
           
           fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_backYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(57,lCyIndex-1,lCzIndex),node(57,lCyIndex,lCzIndex),node(57,lCyIndex+1,lCzIndex),node(57,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step

        end do


        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
        dsCtoF_frontYZ(:,1,j,k) = dsCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
        dsCtoF_frontYZ(:,2,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
        f1 =  rho(45,lCyIndex-1,lCzIndex) 
        f2 =  rho(45,lCyIndex,lCzIndex) 
        f3 =  rho(45,lCyIndex+1,lCzIndex) 
        f4 =  rho(45,lCyIndex+2,lCzIndex) 
        dsCtoF_frontYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(45,lCyIndex-1,lCzIndex) 
        f2 =  phi(45,lCyIndex,lCzIndex) 
        f3 =  phi(45,lCyIndex+1,lCzIndex) 
        f4 =  phi(45,lCyIndex+2,lCzIndex) 
        dsCtoF_frontYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step

        lCzIndex = closestCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
        dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
        dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
        f1 =  rho(57,lCyIndex-1,lCzIndex) 
        f2 =  rho(57,lCyIndex,lCzIndex) 
        f3 =  rho(57,lCyIndex+1,lCzIndex) 
        f4 =  rho(57,lCyIndex+2,lCzIndex) 
        dsCtoF_backYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        f1 =  phi(57,lCyIndex-1,lCzIndex) 
        f2 =  phi(57,lCyIndex,lCzIndex) 
        f3 =  phi(57,lCyIndex+1,lCzIndex) 
        f4 =  phi(57,lCyIndex+2,lCzIndex) 
        dsCtoF_backYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
        
     end do
  end do

END SUBROUTINE XYSpatialInterpolateInternalNodesToFineGrid

SUBROUTINE ZSpatialInterpolateToFineGrid    ! Interpolate required variable to fine grid

  IMPLICIT NONE
  INTEGER :: i,j,k,m
  REAL(dbl) :: xInterp, yInterp, zInterp
  INTEGER :: lCxIndex, lCzIndex, lCyIndex, lFzIndex
  REAL(dbl) :: f1,f2,f3,f4              ! Temporary variables for interpolation.
  
  !Do the bottom and top x-z planes first
  !z - interpolation
  do k=1,nzSub_fine
     IF ( MODULO(k-1, gridRatio) .gt. 0 ) THEN
        do i=1,nxSub_fine
           lFzIndex = k - MODULO(k-1, gridRatio)  ! Lower Fine z Index 
           
           zInterp = dble(MODULO(k-1, gridRatio)) / dble(gridRatio)
           do m=0,NumDistDirs
              fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,2,i,k) !Cycle the second time step to the first time step
              fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
              fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(fCtoF_bottomXZ(m,3,i,lFzIndex-gridRatio),fCtoF_bottomXZ(m,3,i,lFzIndex),fCtoF_bottomXZ(m,3,i,lFzIndex+gridRatio),fCtoF_bottomXZ(m,3,i,lFzIndex+2*gridRatio), node_fine(i,1,lFzIndex-gridRatio), node_fine(i,1,lFzIndex), node_fine(i,1,lFzIndex+gridRatio), node_fine(i,1,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step

              fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
              fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
              fCtoF_topXZ(m,3,i,k) = spatialInterpolate(fCtoF_topXZ(m,3,i,lFzIndex-gridRatio),fCtoF_topXZ(m,3,i,lFzIndex),fCtoF_topXZ(m,3,i,lFzIndex+gridRatio),fCtoF_topXZ(m,3,i,lFzIndex+2*gridRatio), node_fine(i,ny_fine,lFzIndex-gridRatio), node_fine(i,ny_fine,lFzIndex), node_fine(i,ny_fine,lFzIndex+gridRatio), node_fine(i,ny_fine,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
              
           end do

           dsCtoF_bottomXZ(:,1,i,k) = dsCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
           dsCtoF_bottomXZ(:,2,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
           dsCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(dsCtoF_bottomXZ(1,3,i,lFzIndex-gridRatio),dsCtoF_bottomXZ(1,3,i,lFzIndex),dsCtoF_bottomXZ(1,3,i,lFzIndex+gridRatio),dsCtoF_bottomXZ(1,3,i,lFzIndex+2*gridRatio),node_fine(i,1,lFzIndex-gridRatio), node_fine(i,1,lFzIndex), node_fine(i,1,lFzIndex+gridRatio), node_fine(i,1,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
           dsCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(dsCtoF_bottomXZ(2,3,i,lFzIndex-gridRatio),dsCtoF_bottomXZ(2,3,i,lFzIndex),dsCtoF_bottomXZ(2,3,i,lFzIndex+gridRatio),dsCtoF_bottomXZ(2,3,i,lFzIndex+2*gridRatio),node_fine(i,1,lFzIndex-gridRatio), node_fine(i,1,lFzIndex), node_fine(i,1,lFzIndex+gridRatio), node_fine(i,1,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
           
           dsCtoF_topXZ(:,1,i,k) = dsCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
           dsCtoF_topXZ(:,2,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
           dsCtoF_topXZ(1,3,i,k) = spatialInterpolate(dsCtoF_topXZ(1,3,i,lFzIndex-gridRatio),dsCtoF_topXZ(1,3,i,lFzIndex),dsCtoF_topXZ(1,3,i,lFzIndex+gridRatio),dsCtoF_topXZ(1,3,i,lFzIndex+2*gridRatio), node_fine(i,ny_fine,lFzIndex-gridRatio), node_fine(i,ny_fine,lFzIndex), node_fine(i,ny_fine,lFzIndex+gridRatio), node_fine(i,ny_fine,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
           dsCtoF_topXZ(2,3,i,k) = spatialInterpolate(dsCtoF_topXZ(2,3,i,lFzIndex-gridRatio),dsCtoF_topXZ(2,3,i,lFzIndex),dsCtoF_topXZ(2,3,i,lFzIndex+gridRatio),dsCtoF_topXZ(2,3,i,lFzIndex+2*gridRatio), node_fine(i,ny_fine,lFzIndex-gridRatio), node_fine(i,ny_fine,lFzIndex), node_fine(i,ny_fine,lFzIndex+gridRatio), node_fine(i,ny_fine,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
           
        end do
     END IF
  end do
  
  !Fill out the remaining points on the front and back y-z planes
  !z - interpolation
  do k=1,nzSub_fine
     IF ( MODULO(k-1, gridRatio) .gt. 0) THEN
        do j=2,nySub_fine-1
           lFzIndex = k - MODULO(k-1, gridRatio)  ! Lower Fine z Index 
           
           zInterp = dble( MODULO(k-1, gridRatio) ) / dble(gridRatio)
           do m=0,NumDistDirs              
              fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,2,j,k) !Cycle the second time step to the first time step
              fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
              fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(fCtoF_frontYZ(m,3,j,lFzIndex-gridRatio),fCtoF_frontYZ(m,3,j,lFzIndex),fCtoF_frontYZ(m,3,j,lFzIndex+gridRatio),fCtoF_frontYZ(m,3,j,lFzIndex+2*gridRatio), node_fine(1,j,lFzIndex-gridRatio), node_fine(1,j,lFzIndex), node_fine(1,j,lFzIndex+gridRatio), node_fine(1,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
              
              fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
              fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
              fCtoF_backYZ(m,3,j,k) = spatialInterpolate(fCtoF_backYZ(m,3,j,lFzIndex-gridRatio),fCtoF_backYZ(m,3,j,lFzIndex),fCtoF_backYZ(m,3,j,lFzIndex+gridRatio),fCtoF_backYZ(m,3,j,lFzIndex+2*gridRatio), node_fine(nx_fine,j,lFzIndex-gridRatio), node_fine(nx_fine,j,lFzIndex), node_fine(nx_fine,j,lFzIndex+gridRatio), node_fine(nx_fine,j,lFzIndex+2*gridRatio), zInterp)              
             
           end do

           dsCtoF_frontYZ(:,1,j,k) = dsCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
           dsCtoF_frontYZ(:,2,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
           dsCtoF_frontYZ(1,3,j,k) = spatialInterpolate(dsCtoF_frontYZ(1,3,j,lFzIndex-gridRatio),dsCtoF_frontYZ(1,3,j,lFzIndex),dsCtoF_frontYZ(1,3,j,lFzIndex+gridRatio),dsCtoF_frontYZ(1,3,j,lFzIndex+2*gridRatio), node_fine(1,j,lFzIndex-gridRatio), node_fine(1,j,lFzIndex), node_fine(1,j,lFzIndex+gridRatio), node_fine(1,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
           dsCtoF_frontYZ(2,3,j,k) = spatialInterpolate(dsCtoF_frontYZ(2,3,j,lFzIndex-gridRatio),dsCtoF_frontYZ(2,3,j,lFzIndex),dsCtoF_frontYZ(2,3,j,lFzIndex+gridRatio),dsCtoF_frontYZ(2,3,j,lFzIndex+2*gridRatio), node_fine(1,j,lFzIndex-gridRatio), node_fine(1,j,lFzIndex), node_fine(1,j,lFzIndex+gridRatio), node_fine(1,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
           
           dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
           dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
           dsCtoF_backYZ(1,3,j,k) = spatialInterpolate(dsCtoF_backYZ(1,3,j,lFzIndex-gridRatio),dsCtoF_backYZ(1,3,j,lFzIndex),dsCtoF_backYZ(1,3,j,lFzIndex+gridRatio),dsCtoF_backYZ(1,3,j,lFzIndex+2*gridRatio),  node_fine(nx_fine,j,lFzIndex-gridRatio), node_fine(nx_fine,j,lFzIndex), node_fine(nx_fine,j,lFzIndex+gridRatio), node_fine(nx_fine,j,lFzIndex+2*gridRatio), zInterp) 
           dsCtoF_backYZ(2,3,j,k) = spatialInterpolate(dsCtoF_backYZ(2,3,j,lFzIndex-gridRatio),dsCtoF_backYZ(2,3,j,lFzIndex),dsCtoF_backYZ(2,3,j,lFzIndex+gridRatio),dsCtoF_backYZ(2,3,j,lFzIndex+2*gridRatio),  node_fine(nx_fine,j,lFzIndex-gridRatio), node_fine(nx_fine,j,lFzIndex), node_fine(nx_fine,j,lFzIndex+gridRatio), node_fine(nx_fine,j,lFzIndex+2*gridRatio), zInterp) 
           
        end do
     END IF
  end do
  
END SUBROUTINE ZSpatialInterpolateToFineGrid

SUBROUTINE TemporalInterpolateToFineGrid

  INTEGER :: i,j,k,m
  REAL(dbl) :: tInterp !The time to which the temporal interpolation has to be done - Non-dimensionalized by the coarse mesh time step.
  REAL(dbl) :: tmp

  write(*,*) 'subIter in Temporal Interpolate = ', subIter  
  tInterp = dble(subIter-1)/dble(gridRatio)
  
  !Do the bottom and top x-z planes first
  do k=1,nzSub_fine
     do i=1,nxSub_fine
        do m=0,NumDistDirs
           f_fine(m,i,1,k) = temporalInterpolate(fCtoF_bottomXZ(m,1,i,k),fCtoF_bottomXZ(m,2,i,k),fCtoF_bottomXZ(m,3,i,k), node_fine_bottomXZ(1,i,k), node_fine_bottomXZ(2,i,k), node_fine_bottomXZ(3,i,k), tInterp)
           f_fine(m,i,ny_fine,k) = temporalInterpolate(fCtoF_topXZ(m,1,i,k),fCtoF_topXZ(m,2,i,k),fCtoF_topXZ(m,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
        end do
        rho_fine(i,1,k) = temporalInterpolate(dsCtoF_bottomXZ(1,1,i,k),dsCtoF_bottomXZ(1,2,i,k),dsCtoF_bottomXZ(1,3,i,k),node_fine_bottomXZ(1,i,k), node_fine_bottomXZ(2,i,k), node_fine_bottomXZ(3,i,k), tInterp)
        phi_fine(i,1,k) = temporalInterpolate(dsCtoF_bottomXZ(2,1,i,k),dsCtoF_bottomXZ(2,2,i,k),dsCtoF_bottomXZ(2,3,i,k), node_fine_bottomXZ(1,i,k), node_fine_bottomXZ(2,i,k), node_fine_bottomXZ(3,i,k), tInterp)
        rho_fine(i,ny_fine,k) = temporalInterpolate(dsCtoF_topXZ(1,1,i,k),dsCtoF_topXZ(1,2,i,k),dsCtoF_topXZ(1,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
        phi_fine(i,ny_fine,k) = temporalInterpolate(dsCtoF_topXZ(2,1,i,k),dsCtoF_topXZ(2,2,i,k),dsCtoF_topXZ(2,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)

     end do
  end do

  !Fill out the remaining points on the front and back y-z planes
  do k=1,nzSub_fine
     do j=2,nySub_fine-1
        do m=0,NumDistDirs
           f_fine(m,1,j,k) = temporalInterpolate(fCtoF_frontYZ(m,1,j,k),fCtoF_frontYZ(m,2,j,k),fCtoF_frontYZ(m,3,j,k), node_fine_frontYZ(1,j,k), node_fine_frontYZ(2,j,k), node_fine_frontYZ(3,j,k), tInterp)
           f_fine(m,nx_fine,j,k) = temporalInterpolate(fCtoF_backYZ(m,1,j,k),fCtoF_backYZ(m,2,j,k),fCtoF_backYZ(m,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)
        end do
        rho_fine(1,j,k) = temporalInterpolate(dsCtoF_frontYZ(1,1,j,k),dsCtoF_frontYZ(1,2,j,k),dsCtoF_frontYZ(1,3,j,k), node_fine_frontYZ(1,j,k), node_fine_frontYZ(2,j,k), node_fine_frontYZ(3,j,k), tInterp)
        phi_fine(1,j,k) = temporalInterpolate(dsCtoF_frontYZ(2,1,j,k),dsCtoF_frontYZ(2,2,j,k),dsCtoF_frontYZ(2,3,j,k), node_fine_frontYZ(1,j,k), node_fine_frontYZ(2,j,k), node_fine_frontYZ(3,j,k), tInterp)
        rho_fine(nx_fine,j,k) = temporalInterpolate(dsCtoF_backYZ(1,1,j,k),dsCtoF_backYZ(1,2,j,k),dsCtoF_backYZ(1,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)
        phi_fine(nx_fine,j,k) = temporalInterpolate(dsCtoF_backYZ(2,1,j,k),dsCtoF_backYZ(2,2,j,k),dsCtoF_backYZ(2,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)        
     end do
  end do

END SUBROUTINE TemporalInterpolateToFineGrid

SUBROUTINE ComputeEquilibriumForCoarseGrid

  !!! Compute the equilibrium distribution function at the fine grid interface for the coarse grid 
  INTEGER   :: i,j,k,m,iFine,jFine,kFine
  REAL(dbl) :: uu,ue,ve,we,Usum		! precalculated quantities for use in the feq equation
  REAL(dbl) :: feq			! equilibrium distribution function
  REAL(dbl) :: f1,f2,f3,f4              ! Temporary variables for interpolation.
  
  !Do the bottom and top x-z planes first
  do k = 1, nzSub
     kFine = 1 + (k-1)*gridRatio
     do i = 46, 56
        iFine = 1 + (i-45)*gridRatio
        IF(node_fine(iFine,1+gridRatio,kFine) .EQ. FLUID) THEN

           uu = u_fine(iFine,1+gridRatio,kFine)*u_fine(iFine,1+gridRatio,kFine) + v_fine(iFine,1+gridRatio,kFine)*v_fine(iFine,1+gridRatio,kFine) + w_fine(iFine,1+gridRatio,kFine)*w_fine(iFine,1+gridRatio,kFine)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u_fine(iFine,1+gridRatio,kFine)*ex(m)		! u . e
              ve	= v_fine(iFine,1+gridRatio,kFine)*ey(m)		! v . e
              we	= w_fine(iFine,1+gridRatio,kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(iFine,1+gridRatio,kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_bottomXZ(m,i,k) = feq
              
           END DO
           
        END IF

        IF(node_fine(iFine,nySub_fine-gridRatio,kFine) .EQ. FLUID) THEN
           
           uu = u_fine(iFine,nySub_fine-gridRatio,kFine)*u_fine(iFine,nySub_fine-gridRatio,kFine) + v_fine(iFine,nySub_fine-gridRatio,kFine)*v_fine(iFine,nySub_fine-gridRatio,kFine) + w_fine(iFine,nySub_fine-gridRatio,kFine)*w_fine(iFine,nySub_fine-gridRatio,kFine)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u_fine(iFine,nySub_fine-gridRatio,kFine)*ex(m)		! u . e
              ve	= v_fine(iFine,nySub_fine-gridRatio,kFine)*ey(m)		! v . e
              we	= w_fine(iFine,nySub_fine-gridRatio,kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(iFine,nySub_fine-gridRatio,kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_topXZ(m,i,k) = feq
              
           END DO
           
        END IF
        
     end do
  end do

  !Fill in the remaining points on the front and back planes
  do k = 1, nzSub
     kFine = 1 + (k-1)*gridRatio
     do j = 47, 55
        jFine = 1 + (j-45)*gridRatio
        IF(node_fine(1+gridRatio,jFine,kFine) .EQ. FLUID) THEN
           
           uu = u_fine(1+gridRatio,jFine,kFine)*u_fine(1+gridRatio,jFine,kFine) + v_fine(1+gridRatio,jFine,kFine)*v_fine(1+gridRatio,jFine,kFine) + w_fine(1+gridRatio,jFine,kFine)*w_fine(1+gridRatio,jFine,kFine)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u_fine(1+gridRatio,jFine,kFine)*ex(m)		! u . e
              ve	= v_fine(1+gridRatio,jFine,kFine)*ey(m)		! v . e
              we	= w_fine(1+gridRatio,jFine,kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(1+gridRatio,jFine,kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_frontYZ(m,j,k) = feq
              
           END DO
           
        END IF

        IF(node_fine(nxSub_fine-gridRatio, jFine, kFine) .EQ. FLUID) THEN
           
           uu = u_fine(nxSub_fine-gridRatio, jFine, kFine)*u_fine(nxSub_fine-gridRatio, jFine, kFine) + v_fine(nxSub_fine-gridRatio, jFine, kFine)*v_fine(nxSub_fine-gridRatio, jFine, kFine) + w_fine(nxSub_fine-gridRatio, jFine, kFine)*w_fine(nxSub_fine-gridRatio, jFine, kFine)						! u . u
           DO m=0,NumDistDirs
              
              ue	= u_fine(nxSub_fine-gridRatio, jFine, kFine)*ex(m)		! u . e
              ve	= v_fine(nxSub_fine-gridRatio, jFine, kFine)*ey(m)		! v . e
              we	= w_fine(nxSub_fine-gridRatio, jFine, kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(nxSub_fine-gridRatio, jFine, kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_backYZ(m,j,k) = feq
              
           END DO
           
        END IF
        
     end do
  end do

  
END SUBROUTINE ComputeEquilibriumForCoarseGrid


SUBROUTINE InterpolateToCoarseGrid      ! Interpolate required variables to coarse grid

  INTEGER :: i,j,k,m

  !Do the bottom and top x-z planes first
  do k=1,nzSub
     do i=46,56
        do m=0,NumDistDirs
           f(m,i,46,k) = feqFC_bottomXZ(m,i,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (f_fine(m,closestFineIindex(x(i)), closestFineJindex(y(46)), closestFineKindex(z(k))) -  feqFC_bottomXZ(m,i,k))
           f(m,i,56,k) = feqFC_topXZ(m,i,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (f_fine(m,closestFineIindex(x(i)), closestFineJindex(y(56)), closestFineKindex(z(k))) - feqFC_topXZ(m,i,k))
        end do
        rho(i,46,k) = rho_fine(closestFineIindex(x(i)), closestFineJindex(y(46)), closestFineKindex(z(k)))
        phi(i,46,k) = phi_fine(closestFineIindex(x(i)), closestFineJindex(y(46)), closestFineKindex(z(k)))
        rho(i,56,k) = rho_fine(closestFineIindex(x(i)), closestFineJindex(y(56)), closestFineKindex(z(k)))
        phi(i,56,k) = phi_fine(closestFineIindex(x(i)), closestFineJindex(y(56)), closestFineKindex(z(k)))        
     end do
  end do

  !Fill out the remaining points on the front and back y-z planes
  do k=1,nzSub
     do j=47,55
        do m=0,NumDistDirs
           f(m,46,j,k) = feqFC_frontYZ(m,j,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (f_fine(m,closestFineIindex(x(46)), closestFineJindex(y(j)), closestFineKindex(z(k))) - feqFC_frontYZ(m,j,k))
           f(m,56,j,k) = feqFC_backYZ(m,j,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (f_fine(m,closestFineIindex(x(56)), closestFineJindex(y(j)), closestFineKindex(z(k))) - feqFC_backYZ(m,j,k))
        end do
        rho(46,j,k) = rho_fine(closestFineIindex(x(46)), closestFineJindex(y(j)), closestFineKindex(z(k)))
        phi(46,j,k) = phi_fine(closestFineIindex(x(46)), closestFineJindex(y(j)), closestFineKindex(z(k)))
        rho(56,j,k) = rho_fine(closestFineIindex(x(56)), closestFineJindex(y(j)), closestFineKindex(z(k)))
        phi(56,j,k) = phi_fine(closestFineIindex(x(56)), closestFineJindex(y(j)), closestFineKindex(z(k)))
     end do
  end do
         
! temporalInterpolate(fFtoC_backYZ(1,j,k),fFtoC_backYZ(2,j,k), fFtoC_backYZ(3,j,k), desiredTime)

END SUBROUTINE InterpolateToCoarseGrid


FUNCTION temporalInterpolate(f1,f2,f3,n1,n2,n3,t)

  REAL(dbl) :: f1, f2, f3, t
  REAL(dbl)   :: n1,n2,n3
  REAL(dbl) :: temporalInterpolate

  if (n3 .eq. SOLID) then
     if (n2 .ne. SOLID) then
        temporalInterpolate = temporalExtrapolate_n1n2(f1,f2,t)
     else
        temporalInterpolate = 0.0
     end if
  else if (n1 .eq. SOLID) then
     if (n2 .ne. SOLID) then
        temporalInterpolate = temporalInterpolate_n2n3(f2,f3,t)
     end if
  else
     temporalInterpolate = temporalInterpolateAllThree(f1,f2,f3,t)
  end if
  
  RETURN

END FUNCTION temporalInterpolate

FUNCTION temporalExtrapolate_n1n2(f1,f2,t)

!!!Linear extrapolation 

  REAL(dbl) :: f1, f2, t
  REAL(dbl) :: temporalExtrapolate_n1n2

  temporalExtrapolate_n1n2 = t * (f2 - f1) + f2  

  RETURN
  
END FUNCTION temporalExtrapolate_n1n2

FUNCTION temporalInterpolate_n2n3(f2,f3,t)

  REAL(dbl) :: f2, f3, t
  REAL(dbl) :: temporalInterpolate_n2n3
  
  temporalInterpolate_n2n3 = f2 + (f3-f2) * t

  RETURN

END FUNCTION temporalInterpolate_n2n3

FUNCTION temporalInterpolateAllThree(f1,f2,f3,t)

  REAL(dbl) :: f1, f2, f3, t
  REAL(dbl) :: temporalInterpolateAllThree

  temporalInterpolateAllThree = f2 + ((f3-f1)*0.5 + ( (f1+f3)*0.5 - f2)*t) * t

  RETURN

END FUNCTION temporalInterpolateAllThree

FUNCTION spatialInterpolate(f1,f2,f3,f4,n1,n2,n3,n4,s)

!!!Symmetric Cubic spline temporal interpolation 

  REAL(dbl) :: f1, f2, f3, f4, s
  INTEGER   :: n1,n2,n3,n4
  REAL(dbl) :: spatialInterpolate
  REAL(dbl) :: aHat, bHat, cHat, dHat

  if (n3 .eq. SOLID) then
     if (n2 .ne. SOLID) then
        spatialInterpolate = spatialExtrapolate_n1n2(f1,f2,s)
     else
        spatialInterpolate = 0.0
     end if
     
  else if (n2 .eq. SOLID) then
     spatialInterpolate = spatialExtrapolate_n3n4(f3,f4,s)
     
  else if ( (n1 .ne. SOLID) .and. (n4 .eq. SOLID) ) then
     spatialInterpolate = spatialInterpolate_n1n2n3(f1,f2,f3,s)

  else if ( (n1 .eq. SOLID) .and. (n4 .ne. SOLID) ) then
     spatialInterpolate = spatialInterpolate_n2n3n4(f2,f3,f4,s)
     
  else if ( (n1 .eq. SOLID) .and. (n4 .eq. SOLID) ) then
     spatialInterpolate = spatialInterpolate_n2n3(f2,f3,s)

  else if ( (n1 .ne. SOLID) .and. (n2 .ne. SOLID) .and. (n3 .ne. SOLID) .and. (n4 .ne. SOLID) ) then
     spatialInterpolate = spatialInterpolateAllFour(f1,f2,f3,f4,s)
  end if
   
  RETURN
  
END FUNCTION spatialInterpolate

FUNCTION spatialExtrapolate_n1n2(f1,f2,s)

!!!Linear extrapolation 

  REAL(dbl) :: f1, f2, s
  REAL(dbl) :: spatialExtrapolate_n1n2

  spatialExtrapolate_n1n2 = s * (f2 - f1) + f2  

  RETURN
  
END FUNCTION spatialExtrapolate_n1n2

FUNCTION spatialInterpolate_n1n2n3(f1,f2,f3,s)

!!!Linear extrapolation 

  REAL(dbl) :: f1, f2, f3, s
  REAL(dbl) :: spatialInterpolate_n1n2n3

  spatialInterpolate_n1n2n3 = f2 + ((f3-f1)*0.5 + ( (f1+f3)*0.5 - f2)*s) * s

  RETURN
  
END FUNCTION spatialInterpolate_n1n2n3

FUNCTION spatialInterpolate_n2n3(f2,f3,s)

!!!Symmetric Linear interpolation 

  REAL(dbl) :: f2, f3, s
  REAL(dbl) :: spatialInterpolate_n2n3

  spatialInterpolate_n2n3 = s * (f3 - f2) + f2  ! s * f3 + (1.0_dbl-s) * f2

  RETURN
  
END FUNCTION spatialInterpolate_n2n3

FUNCTION spatialInterpolate_n2n3n4(f2,f3,f4,s)

!!!Symmetric Linear interpolation 

  REAL(dbl) :: f2, f3, f4, s
  REAL(dbl) :: spatialInterpolate_n2n3n4

  spatialInterpolate_n2n3n4 = f2 + ( (-1.5*f2 + 2.0_dbl*f3 - 0.5*f4) + (0.5*f2 - f3 + 0.5*f4)*s) * s

  RETURN
  
END FUNCTION spatialInterpolate_n2n3n4

FUNCTION spatialExtrapolate_n3n4(f3,f4,s)

!!!Symmetric Linear interpolation 

  REAL(dbl) :: f3, f4, s
  REAL(dbl) :: spatialExtrapolate_n3n4

  spatialExtrapolate_n3n4 = (2.0_dbl*f3 - f4) + (f4 - f3)*s 

  RETURN
  
END FUNCTION spatialExtrapolate_n3n4

FUNCTION spatialInterpolateAllFour(f1,f2,f3,f4,s)

!!!Symmetric Cubic spline temporal interpolation 

  REAL(dbl) :: f1, f2, f3, f4, s
  REAL(dbl) :: spatialInterpolateAllFour
  REAL(dbl) :: aHat, bHat, cHat, dHat

  aHat = (-f1 + 3*(f2 - f3) + f4)/6.0
  bHat = 0.5 * (f1 + f3) - f2
  dHat = f2
  cHat = f3 - aHat - bHat - dHat
  spatialInterpolateAllFour = dHat + s * (cHat + s * (bHat + s * aHat)) !Written in a weird way to save on multiplications and additions

  RETURN
  
END FUNCTION spatialInterpolateAllFour

SUBROUTINE testSpatialInterpolation

  write(31,*) 'Cubic spline interpolation test'
  write(31,*) 'Spatial co-ordinate runs from -1 to 2.0. The values at the 4 locations -1,0,1,2 are 1,2,3,4 respectively'
  write(31,*) 'The interpolated value at -1 is ', spatialInterpolate(1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl,1,1,1,1,-1.0_dbl)
  write(31,*) 'The interpolated value at 0 is ', spatialInterpolate(1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl,1,1,1,1,0.0_dbl)
  write(31,*) 'The interpolated value at 1 is ', spatialInterpolate(1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl,1,1,1,1,1.0_dbl)
  write(31,*) 'The interpolated value at 2 is ', spatialInterpolate(1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl,1,1,1,1,2.0_dbl)
  write(31,*) 'The interpolated value at 0.5 is ', spatialInterpolate(1.0_dbl,2.0_dbl,3.0_dbl,4.0_dbl,1,1,1,1,0.5_dbl)
  
END SUBROUTINE testSpatialInterpolation

SUBROUTINE InitializeAllTemporalInterpolation
!!! To begin with set all the 3 values for temporal interpolation to the latest available value

  INTEGER :: i,j,k,m !Counter variables
  
  do k=1,nzSub_fine
     do i=1,nxSub_fine
        do m=0,NumDistDirs
           fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the first time step
           fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
           fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the first time step
           fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
        end do
        dsCtoF_bottomXZ(:,1,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the first time step
        dsCtoF_bottomXZ(:,2,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
        dsCtoF_topXZ(:,1,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the first time step
        dsCtoF_topXZ(:,2,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
        node_fine_bottomXZ(:,i,k) = node_fine(i,1,k)
        node_fine_topXZ(:,i,k) = node_fine(i,ny_fine,k)

     end do
  end do

  do k=1,nzSub_fine
     do j=2,nySub_fine-1
        do m=0,NumDistDirs
           fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the first time step
           fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
           fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the first time step
           fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
        end do
        dsCtoF_frontYZ(:,1,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the first time step
        dsCtoF_frontYZ(:,2,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
        dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the first time step
        dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step        
        node_fine_frontYZ(:,j,k) = node_fine(1,j,k)
        node_fine_backYZ(:,j,k) = node_fine(nx_fine,j,k)
     end do
  end do
  
END SUBROUTINE InitializeAllTemporalInterpolation

!================================================
END MODULE LBM_FINE
!================================================
