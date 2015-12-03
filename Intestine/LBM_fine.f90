!==================================================================================================
MODULE LBM_fine				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE Setup_fine
USE ICBC_fine
USE Output

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

FUNCTION lowerCoarseXindex(xf)
  !Returns the lower coarse X index
  REAL(dbl) :: xf
  INTEGER :: lowerCoarseXindex
!  write(*,*) 'xf = ', xf, ' x(1) = ' , x(1), ' xcf = ', xcf
  lowerCoarseXindex = FLOOR( ((xf-x(1))/xcf) ) + 1
  RETURN 
END FUNCTION lowerCoarseXindex

FUNCTION lowerCoarseYindex(yf)
  !Returns the lower coarse Y index
  REAL(dbl) :: yf
  INTEGER :: lowerCoarseYindex
!  write(*,*) 'yf = ', yf, ' y(1) = ' , y(1), ' ycf = ', ycf
  lowerCoarseYindex = FLOOR( ((yf-y(1))/ycf) ) + 1
  RETURN 
END FUNCTION lowerCoarseYindex

FUNCTION lowerCoarseZindex(zf)
  !Returns the lower coarse Z index
  REAL(dbl) :: zf
  INTEGER :: lowerCoarseZindex
!  write(*,*) 'zf = ', zf, ' z(1) = ' , z(1), ' zcf = ', zcf, ' kMin = ', kMin
  lowerCoarseZindex = FLOOR( ((zf-z(1))/zcf) ) + 1  
  RETURN 
END FUNCTION lowerCoarseZindex

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
        IF(node(i,46,k) .EQ. FLUID) THEN
           
           uu = u(i,46,k)*u(i,46,k) + v(i,46,k)*v(i,46,k) + w(i,46,k)*w(i,46,k)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u(i,46,k)*ex(m)		! u . e
              ve	= v(i,46,k)*ey(m)		! v . e
              we	= w(i,46,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(i,46,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function

              feqFF_bottomXZ(m,i,k) = feq
!              WRITE(31,*) 'feqFF_bottomXZ(', m, ',', i, ',', k, ') = ', feqFF_bottomXZ(m,i,k)
              
           END DO
           
        END IF

        IF(node(i,56,k) .EQ. FLUID) THEN
           
           uu = u(i,56,k)*u(i,56,k) + v(i,56,k)*v(i,56,k) + w(i,56,k)*w(i,56,k)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u(i,56,k)*ex(m)		! u . e
              ve	= v(i,56,k)*ey(m)		! v . e
              we	= w(i,56,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(i,56,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFF_topXZ(m,i,k) = feq
!              WRITE(31,*) 'feqFF_topXZ(', m, ',', i, ',', k, ') = ', feqFF_topXZ(m,i,k)
              
           END DO
           
        END IF
        
     end do
  end do

  CALL FLUSH(31)
  
  !Fill in the remaining points on the front and back planes
  do k = 1, nzSub
     do j = 44, 58
        IF(node(46,j,k) .EQ. FLUID) THEN
           
           uu = u(46,j,k)*u(46,j,k) + v(46,j,k)*v(46,j,k) + w(46,j,k)*w(46,j,k)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u(46,j,k)*ex(m)		! u . e
              ve	= v(46,j,k)*ey(m)		! v . e
              we	= w(46,j,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(46,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFF_frontYZ(m,j,k) = feq
              
           END DO
           
        END IF

        IF(node(56,j,k) .EQ. FLUID) THEN
           
           uu = u(56,j,k)*u(56,j,k) + v(56,j,k)*v(56,j,k) + w(56,j,k)*w(56,j,k)	! u . u
           DO m=1,NumDistDirs
              
              ue	= u(56,j,k)*ex(m)		! u . e
              ve	= v(56,j,k)*ey(m)		! v . e
              we	= w(56,j,k)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho(56,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
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

        lCxIndex = lowerCoarseXindex(x_fine(i))  ! Lower Coarse x Index
        lCzIndex = lowerCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        
        xInterp = dble( MODULO(i-1, gridRatio) ) / dble(gridRatio)

        do m=1,14	  

           fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_bottomXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex-1,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_bottomXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_bottomXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+1,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_bottomXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+2,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4,xInterp) !Interpolate the latest value to the last(third) time step

           write(31,*) 'tau_fine = ', tau_fine, ' tau = ', tau
           write(31,*) 'gridRatio = ', gridRatio
           WRITE(31,*) '(tau_fine - 1.0)/(gridRatio * (tau - 1.0)) = ', (tau_fine - 1.0)/(gridRatio * (tau - 1.0))
           WRITE(31,*) 'fPlus1 = ', fPlus(m,lCxIndex-1,46,lCzIndex), ' fPlus2 = ', fPlus(m,lCxIndex,46,lCzIndex), ' fPlus3 = ', fPlus(m,lCxIndex+1,46,lCzIndex), ' fPlus4 = ', fPlus(m,lCxIndex+2,46,lCzIndex)
           WRITE(31,*) 'feq1 = ', feqFF_bottomXZ(m,lCxIndex-1,lCzIndex), ' feq2 = ', feqFF_bottomXZ(m,lCxIndex,lCzIndex), ' feq3 = ', feqFF_bottomXZ(m,lCxIndex+1,lCzIndex), ' feq4 = ', feqFF_bottomXZ(m,lCxIndex+2,lCzIndex)
           WRITE(31,*) 'xInterp = ', xInterp, ' f1 = ', f1, ' f2 = ', f2, ' f3 = ', f3, ' f4 = ', f4
           WRITE(31,*) 'fCtoF_bottomXZ(', m, ',3', i, ',', k, ') = ', fCtoF_bottomXZ(m,3,i,k)
           
           fCtoF_bottomXZ(m,1,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_bottomXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_bottomXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_bottomXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_bottomXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex-1,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_bottomXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_bottomXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+1,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_bottomXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+2,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_bottomXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,xInterp) !Interpolate the latest value to the last(third) time step
!           WRITE(31,*) 'fCtoF_bottomXZ(', m, ',3', i, ',', nzSub_fine-gridRatio+1-k+1, ') = ', fCtoF_bottomXZ(m,3,i,nzSub_fine-gridRatio+1-k+1)

           fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_topXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex-1,56,lCzIndex) - feqFF_topXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_topXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex,56,lCzIndex) - feqFF_topXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_topXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+1,56,lCzIndex) - feqFF_topXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_topXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+2,56,lCzIndex) - feqFF_topXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_topXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4,xInterp) !Interpolate the latest value to the last(third) time step
!           WRITE(31,*) 'fCtoF_topXZ(', m, ',3', i, ',', k, ') = ', fCtoF_topXZ(m,3,i,k)

           fCtoF_topXZ(m,1,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_topXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_topXZ(m,2,i,nzSub_fine-gridRatio+1-k+1) = fCtoF_topXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_topXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex-1,56,lCzIndex) - feqFF_topXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_topXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex,56,lCzIndex) - feqFF_topXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_topXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+1,56,lCzIndex) - feqFF_topXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_topXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+2,56,lCzIndex) - feqFF_topXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_topXZ(m,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,xInterp) !Interpolate the latest value to the last(third) time step
!           WRITE(31,*) 'fCtoF_topXZ(', m, ',3', i, ',', nzSub_fine-gridRatio+1-k+1, ') = ', fCtoF_topXZ(m,3,i,nzSub_fine-gridRatio+1-k+1)

        end do
     end do
  end do
 
  !Fill out the remaining points on the front and back y-z planes
  !Only y-interpolation
  do k=1, gridRatio+1, gridRatio
     do j=2,nySub_fine-1

        lCyIndex = lowerCoarseYindex(y_fine(j))  ! Lower Coarse x Index
        lCzIndex = lowerCoarseZindex(z_fine(j))  ! Lower Coarse z Index - No interpolation in z
        
        yInterp = dble( MODULO(j-1, gridRatio) ) / dble(gridRatio)
        
        do m=1,14	  
           
           fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_frontYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex-1,lCzIndex) - feqFF_frontYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_frontYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex,lCzIndex) - feqFF_frontYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_frontYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex+1,lCzIndex) - feqFF_frontYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_frontYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex+2,lCzIndex) - feqFF_frontYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,yInterp) !Interpolate the latest value to the last(third) time step
           
           fCtoF_frontYZ(m,1,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_frontYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_frontYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_frontYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_frontYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex-1,lCzIndex) - feqFF_frontYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_frontYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex,lCzIndex) - feqFF_frontYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_frontYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex+1,lCzIndex) - feqFF_frontYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_frontYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex+2,lCzIndex) - feqFF_frontYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_frontYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,yInterp) !Interpolate the latest value to the last(third) time step

           fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_backYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,yInterp) !Interpolate the latest value to the last(third) time step

           fCtoF_backYZ(m,1,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_backYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
           fCtoF_backYZ(m,2,j,nzSub_fine-gridRatio+1-k+1) = fCtoF_backYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
           f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_backYZ(m,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,yInterp) !Interpolate the latest value to the last(third) time step
           
        end do
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
     do i=1,nxSub_fine

        lCxIndex = lowerCoarseXindex(x_fine(i))  ! Lower Coarse x Index
        lCzIndex = lowerCoarseZindex(z_fine(k))  ! Lower Coarse z Index
        
        xInterp = dble( MODULO(i-1, gridRatio) ) / dble(gridRatio)


        do m=1,14	  
           fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_bottomXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex-1,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_bottomXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_bottomXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+1,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_bottomXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+2,46,lCzIndex) - feqFF_bottomXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4,xInterp) !Interpolate the latest value to the last(third) time step
           
           fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
           fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
           f1 =  feqFF_topXZ(m,lCxIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex-1,56,lCzIndex) - feqFF_topXZ(m,lCxIndex-1,lCzIndex))
           f2 =  feqFF_topXZ(m,lCxIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex,56,lCzIndex) - feqFF_topXZ(m,lCxIndex,lCzIndex))
           f3 =  feqFF_topXZ(m,lCxIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+1,56,lCzIndex) - feqFF_topXZ(m,lCxIndex+1,lCzIndex))
           f4 =  feqFF_topXZ(m,lCxIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,lCxIndex+2,56,lCzIndex) - feqFF_topXZ(m,lCxIndex+2,lCzIndex))           
           fCtoF_topXZ(m,3,i,k) = spatialInterpolate(f1,f2,f3,f4,xInterp) !Interpolate the latest value to the last(third) time step

        end do
     end do
  end do
 
  !Fill out the remaining points on the front and back y-z planes
  !Only y-interpolation
  do k=2*gridRatio+1,nzSub_fine-3*gridRatio+1, gridRatio
     do j=2,nySub_fine-1

        lCyIndex = lowerCoarseYindex(y_fine(j))  ! Lower Coarse x Index
        lCzIndex = lowerCoarseZindex(z_fine(k))  ! Lower Coarse z Index - No interpolation in z
        
        yInterp = dble( MODULO(j-1, gridRatio) ) / dble(gridRatio)
        
        do m=1,14	  
           fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_frontYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex-1,lCzIndex) - feqFF_frontYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_frontYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex,lCzIndex) - feqFF_frontYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_frontYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex+1,lCzIndex) - feqFF_frontYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_frontYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,46,lCyIndex+2,lCzIndex) - feqFF_frontYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,yInterp) !Interpolate the latest value to the last(third) time step
           
           fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
           fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
           f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
           f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
           f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
           f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (fPlus(m,56,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
           fCtoF_backYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,yInterp) !Interpolate the latest value to the last(third) time step

        end do
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
           do m=1,4
              lFzIndex = k - MODULO(k-1, gridRatio)  ! Lower Fine z Index 
              
              zInterp = dble(MODULO(k-1, gridRatio)) / dble(gridRatio)

              fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,2,i,k) !Cycle the second time step to the first time step
              fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
              fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(fCtoF_bottomXZ(m,3,i,lFzIndex-gridRatio),fCtoF_bottomXZ(m,3,i,lFzIndex),fCtoF_bottomXZ(m,3,i,lFzIndex+gridRatio),fCtoF_bottomXZ(m,3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step

              fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
              fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
              fCtoF_topXZ(m,3,i,k) = spatialInterpolate(fCtoF_topXZ(m,3,i,lFzIndex-gridRatio),fCtoF_topXZ(m,3,i,lFzIndex),fCtoF_topXZ(m,3,i,lFzIndex+gridRatio),fCtoF_topXZ(m,3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
              
           end do
        end do
     END IF
  end do
  
  !Fill out the remaining points on the front and back y-z planes
  !z - interpolation
  do k=1,nzSub_fine
     IF ( MODULO(k-1, gridRatio) .gt. 0) THEN
        do j=2,nySub_fine-1
           do m=1,14
              lFzIndex = k - MODULO(k-1, gridRatio)  ! Lower Fine z Index 
              
              zInterp = dble( MODULO(k-1, gridRatio) ) / dble(gridRatio)
              
              fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,2,j,k) !Cycle the second time step to the first time step
              fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
              fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(fCtoF_frontYZ(m,3,j,lFzIndex-gridRatio),fCtoF_frontYZ(m,3,j,lFzIndex),fCtoF_frontYZ(m,3,j,lFzIndex+gridRatio),fCtoF_frontYZ(m,3,j,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
              
              fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
              fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
              fCtoF_backYZ(m,3,j,k) = spatialInterpolate(fCtoF_backYZ(m,3,j,lFzIndex-gridRatio),fCtoF_backYZ(m,3,j,lFzIndex),fCtoF_backYZ(m,3,j,lFzIndex+gridRatio),fCtoF_backYZ(m,3,j,lFzIndex+2*gridRatio),zInterp) 
              
           end do
        end do
     END IF
  end do
  
END SUBROUTINE ZSpatialInterpolateToFineGrid

SUBROUTINE TemporalInterpolateToFineGrid

  INTEGER :: i,j,k,m
  REAL(dbl) :: tInterp !The time to which the temporal interpolation has to be done - Non-dimensionalized by the coarse mesh time step.

  tInterp = dble(subIter/gridRatio)
  
  !Do the bottom and top x-z planes first
  do k=1,nzSub_fine
     do i=1,nxSub_fine
        do m=1,14
           fPlus_fine(m,i,1,k) = temporalInterpolate(fCtoF_bottomXZ(m,1,i,k),fCtoF_bottomXZ(m,1,i,k),fCtoF_bottomXZ(m,3,i,k),tInterp)
           fPlus_fine(m,i,ny_fine,k) = temporalInterpolate(fCtoF_topXZ(m,1,i,k),fCtoF_topXZ(m,1,i,k),fCtoF_topXZ(m,3,i,k),tInterp)
        end do
     end do
  end do

  !Fill out the remaining points on the front and back y-z planes
  do k=1,nzSub_fine
     do j=2,nySub_fine-1
        do m=1,14
           fPlus_fine(m,j,1,k) = temporalInterpolate(fCtoF_frontYZ(m,1,j,k),fCtoF_frontYZ(m,1,j,k),fCtoF_frontYZ(m,3,j,k),tInterp)
           fPlus_fine(m,j,ny_fine,k) = temporalInterpolate(fCtoF_backYZ(m,1,j,k),fCtoF_backYZ(m,1,j,k),fCtoF_backYZ(m,3,j,k),tInterp)
        end do
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
     kFine = 1 + k*gridRatio
     do i = 46, 56
        iFine = 1 + (i-46)*gridRatio
        IF(node_fine(iFine,1,kFine) .EQ. FLUID) THEN
           
           uu = u_fine(iFine,1,kFine)*u_fine(iFine,1,kFine) + v_fine(iFine,1,kFine)*v_fine(iFine,1,kFine) + w_fine(iFine,1,kFine)*w_fine(iFine,1,kFine)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u_fine(iFine,1,kFine)*ex(m)		! u . e
              ve	= v_fine(iFine,1,kFine)*ey(m)		! v . e
              we	= w_fine(iFine,1,kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(iFine,1,kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_bottomXZ(m,i,k) = feq
              
           END DO
           
        END IF

        IF(node_fine(iFine,nySub_fine,kFine) .EQ. FLUID) THEN
           
           uu = u_fine(iFine,nySub_fine,kFine)*u_fine(iFine,nySub_fine,kFine) + v_fine(iFine,nySub_fine,kFine)*v_fine(iFine,nySub_fine,kFine) + w_fine(iFine,nySub_fine,kFine)*w_fine(iFine,nySub_fine,kFine)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u_fine(iFine,nySub_fine,kFine)*ex(m)		! u . e
              ve	= v_fine(iFine,nySub_fine,kFine)*ey(m)		! v . e
              we	= w_fine(iFine,nySub_fine,kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(iFine,nySub_fine,kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_topXZ(m,i,k) = feq
              
           END DO
           
        END IF
        
     end do
  end do

  !Fill in the remaining points on the front and back planes
  do k = 1, nzSub
     kFine = 1 + k*gridRatio
     do j = 47, 55
        jFine = 1 + (j-46)*gridRatio
        IF(node_fine(1,jFine,kFine) .EQ. FLUID) THEN
           
           uu = u_fine(1,jFine,kFine)*u_fine(1,jFine,kFine) + v_fine(1,jFine,kFine)*v_fine(1,jFine,kFine) + w_fine(1,jFine,kFine)*w_fine(1,jFine,kFine)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u_fine(1,jFine,kFine)*ex(m)		! u . e
              ve	= v_fine(1,jFine,kFine)*ey(m)		! v . e
              we	= w_fine(1,jFine,kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(1,jFine,kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
              feqFC_frontYZ(m,j,k) = feq
              
           END DO
           
        END IF

        IF(node_fine(nxSub_fine, jFine, kFine) .EQ. FLUID) THEN
           
           uu = u_fine(nxSub_fine, jFine, kFine)*u_fine(nxSub_fine, jFine, kFine) + v_fine(nxSub_fine, jFine, kFine)*v_fine(nxSub_fine, jFine, kFine) + w_fine(nxSub_fine, jFine, kFine)*w_fine(nxSub_fine, jFine, kFine)						! u . u
           DO m=1,NumDistDirs
              
              ue	= u_fine(nxSub_fine, jFine, kFine)*ex(m)		! u . e
              ve	= v_fine(nxSub_fine, jFine, kFine)*ey(m)		! v . e
              we	= w_fine(nxSub_fine, jFine, kFine)*ez(m)		! w . e
              
              Usum	= ue + ve + we				! U . e
              
              feq = (wt(m)*rho_fine(nxSub_fine, jFine, kFine))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function
              
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
        do m=1,14
           fPlus(m,i,46,k) = feqFC_bottomXZ(m,i,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (fPlus_fine(m,closestFineIindex(x(i)), closestFineJindex(y(46)), closestFineKindex(z(k))) -  feqFC_bottomXZ(m,i,k))
           fPlus(m,i,56,k) = feqFC_topXZ(m,i,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (fPlus_fine(m,closestFineIindex(x(i)), closestFineJindex(y(56)), closestFineKindex(z(k))) - feqFC_topXZ(m,i,k))
        end do
     end do
  end do

  !Fill out the remaining points on the front and back y-z planes
  do k=1,nzSub
     do j=47,55
        do m=1,14
           fPlus(m,46,j,k) = feqFC_frontYZ(m,j,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (fPlus_fine(m,closestFineIindex(x(46)), closestFineJindex(y(j)), closestFineKindex(z(k))) - feqFC_frontYZ(m,j,k))
           fPlus(m,56,j,k) = feqFC_backYZ(m,j,k) + gridRatio * (tau - 1.0)/(tau_fine - 1.0) * (fPlus_fine(m,closestFineIindex(x(56)), closestFineJindex(y(j)), closestFineKindex(z(k))) - feqFC_backYZ(m,j,k))
        end do
     end do
  end do
         
! temporalInterpolate(fFtoC_backYZ(1,j,k),fFtoC_backYZ(2,j,k), fFtoC_backYZ(3,j,k), desiredTime)

END SUBROUTINE InterpolateToCoarseGrid


FUNCTION temporalInterpolate(f1,f2,f3,t)

  REAL(dbl) :: f1, f2, f3, t
  REAL(dbl) :: temporalInterpolate
!  write(*,*) 'Dummy temporal interpolation returning middle value f2 for now'
  
  temporalInterpolate = f2
  RETURN

END FUNCTION temporalInterpolate

FUNCTION spatialInterpolate(f1,f2,f3,f4,s)

!!!Symmetric Cubic spline temporal interpolation 

  REAL(dbl) :: f1, f2, f3, f4, s
  REAL(dbl) :: spatialInterpolate
  REAL(dbl) :: aHat, bHat, cHat, dHat

  aHat = (-f1 + 3*(f2 - f3) + f4)/6.0
  bHat = 0.5 * (f1 + f3) - f2
  dHat = f2
  cHat = f3 - aHat - bHat - dHat
  spatialInterpolate = dHat + s * (cHat + s * (bHat + s * aHat)) !Written in a weird way to save on multiplications and additions

  RETURN
  
END FUNCTION spatialInterpolate


!================================================
END MODULE LBM_FINE
!================================================
