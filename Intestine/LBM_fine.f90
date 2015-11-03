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

      IF(node_fine(i,j,k) .EQ. FLUID) THEN
      
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
          f_fine(m,i,j,k)		= f_fine(m,i,j,k) - oneOVERtau*(f_fine(m,i,j,k) - feq)							! collision
        
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

      IF(node(i,j,k) .EQ. FLUID) THEN
      
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
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 89: node_fine(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
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
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
            WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
            WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
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
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
            CALL BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f_fine(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node_fine(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
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
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
            CALL BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f_fine(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node_fine(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
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

        IF(rho(i,j,k) .NE. 0) THEN
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

      ELSE
      
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

SUBROUTINE InterpolateToFineGrid      ! Interpolate required variables to fine grid

  write(*,*) 'Nothing here for now'

END SUBROUTINE InterpolateToFineGrid



SUBROUTINE InterpolateToCoarseGrid    ! Interpolate required variable to coarse grid

  write(*,*) 'Nothing here for now'

END SUBROUTINE InterpolateToCoarseGrid

 


!================================================
END MODULE LBM_FINE
!================================================
