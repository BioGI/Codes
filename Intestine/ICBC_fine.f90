!==================================================================================================
MODULE ICBC_fine		! Sets Initial and Boundary Conditions
					! Contains setup subroutines (ICs,BounceBack,IniParticles)
!==================================================================================================
USE SetPrecision
USE Setup
USE Setup_fine  
USE mpi

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE ICs_fine	! sets the initial conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE  

! Define local variables
INTEGER(lng) :: i,j,k,m									! index variables
REAL(dbl) :: feq
REAL(dbl) :: tmp

IF(restart) THEN											! restart from  file 
  
  OPEN(50,FILE='restart_fine.'//sub)						! open correct restart file
  
  DO k=0,nzSub_fine+1_lng
    DO j=0,nySub_fine+1_lng
      DO i=0,nxSub_fine+1_lng

        READ(50,*) node_fine(i,j,k)
        READ(50,*) u_fine(i,j,k)
        READ(50,*) v_fine(i,j,k)
        READ(50,*) w_fine(i,j,k)
        READ(50,*) rho_fine(i,j,k)
        READ(50,*) phi_fine(i,j,k)

        DO m=0,NumDistDirs
          READ(50,*) f_fine(m,i,j,k)
        END DO

      END DO
    END DO
  END DO

  READ(50,*) phiAbsorbed_fine
  READ(50,*) phiAbsorbedS_fine
  READ(50,*) phiAbsorbedV_fine
  READ(50,*) phiInOut_fine
 
  CLOSE(50)

ELSE															! clean start

  ! Initial conditions on velocity, density, and scalar
  DO k=0,nzSub_fine+1_lng
    DO j=0,nySub_fine+1_lng
      DO i=0,nxSub_fine+1_lng

        u_fine(i,j,k)   = 0.0_dbl							! x-velocity
        v_fine(i,j,k)        = -0.0_dbl
        tmp = sqrt( x_fine(i) * x_fine(i) + y_fine(j)*y_fine(j) )
!        w_fine(i,j,k)        = 1.0 - tmp/r_fine(k) !Linear velocity profile such that the velocity goes to zero at the wall! z-velocity
        !        w_fine(i,j,k) = 1.0
        w_fine(i,j,k) = 0.0
        rho_fine(i,j,k) = denL								! density
	! Balaji added
	! distribution functions (set to equilibrium)
	DO m=0,NumDistDirs
	  CALL Equilibrium_LOCAL_fine(m,rho_fine(i,j,k),u_fine(i,j,k),v_fine(i,j,k),w_fine(i,j,k),feq)	! distribution functions
	  f_fine(m,i,j,k) = feq
	END DO

      END DO
    END DO
  END DO

  ! Starting iteration
  iter0 = 1_lng
  !iter0 = 0_lng

  ! Initialize scalar values
  phiAbsorbed_fine	= 0.0_dbl								! total amount of scalar absorbed
  phiAbsorbedS_fine	= 0.0_dbl								! total amount of scalar absorbed through the macroscopic surface
  phiAbsorbedV_fine	= 0.0_dbl								! total amount of scalar absorbed through the villi
  phiInOut_fine	= 0.0_dbl								! total amount of scalar leaving the inlet/outlet
  delphi_particle_fine = 0.0_dbl								! Initialize the scalar contirbution from particles to 0.0. Once the particle
											! data is read, we can interpolate to get the delphi_particle. In any event, thi
										! variable is designed to store temporary data. 

END IF

!------------------------------------------------
END SUBROUTINE ICs_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarDistribution_fine		! Sets/Maintains initial distributions of scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,ii,jj		! lattice indices
	!write(*,*) iter
	!pause
! INTRODUCTION OF SCALAR
IF(iter .EQ. phiStart) THEN
	!write(*,*) sclrIC
	!pause
  SELECT CASE(sclrIC) 
      
    CASE(BLOB)							! blob of scalar at the center of the domain
      
      DO k=0,nzSub_fine+1
        DO j=0,nySub_fine+1
          DO i=0,nxSub_fine+1

             if ( (abs(x_fine(i)) .lt. 0.001) .and. (abs(y_fine(j)) .lt. 0.001) ) then
                phi_fine(i,j,k) = 1.0
             else
                phi_fine(i,j,k) = 0.0
             end if

          END DO
        END DO
      END DO

    CASE(LINE) 						! line of scalar along axis
  
      DO k=0,nzSub_fine+1
        DO j=0,nySub_fine+1
          DO i=0,nxSub_fine+1

            phi_fine(i,j,k) = phiIC*ee**(-((x_fine(i)**2 + y_fine(j)**2)/(2.0_dbl*sigma**2)))									! 2D Gaussian Distribution in x and y

          END DO
        END DO
      END DO

    CASE(INLET) 						! line of scalar at the inlet
  
      DO k=0,nzSub_fine+1
        DO j=0,nySub_fine+1
          DO i=0,nxSub_fine+1

            phi_fine(i,j,k) = phiIC*ee**(-((z_fine(k)**2)/(2.0_dbl*sigma**2)))													! 1D Gaussian Distribution in z

          END DO
        END DO
      END DO

    CASE(UNIFORM)						! uniform initial distribution
	    !write(*,*) "balaji"
	    !pause
      phi_fine(:,:,:) = phiIC			! set the full scalar field to phiIC

    CASE(LINEAR) 						! linear profile such that phi = 0 at the boundary
  
      DO k=0,nzSub_fine+1
        DO j=0,nySub_fine+1
          DO i=0,nxSub_fine+1

             phi_fine(i,j,k) = 1.0 - sqrt(x_fine(i)**2 + y_fine(j)**2)/r_fine(k)
             
          END DO
        END DO
      END DO

   CASE DEFAULT
   
      OPEN(1000,FILE="error.txt")
      WRITE(1000,*) "Error in ScalarIC in Setup.f90: sclrIC is not BLOB(1), LINE(2) or INLET(3)..."
      WRITE(1000,*) "sclrIC=", sclrIC
      CLOSE(1000)
      STOP

  END SELECT

  ! Calculate the intial amount of scalar
  phiTotal = 0.0_dbl
  DO k=1,nzSub_fine
    DO j=1,nySub_fine
      DO i=1,nxSub_fine

        IF(node_fine(i,j,k) .EQ. FLUID) THEN
          phiTotal = phiTotal + phi_fine(i,j,k)
        END IF

      END DO
    END DO
  END DO

ELSE

  ! MAINTAINENCE OF SCALAR
  SELECT CASE(sclrIC) 
      
    CASE(BLOB)							! blob of scalar at the center of the domain

      ! scalar is not maintained
  
    CASE(LINE) 						! line of scalar along axis

!       IF((SubID(2) .EQ. 0) .AND. (SubID(4) .EQ. 0)) THEN			! if no neighboring subdomains exist in the 2nd and 4th directions, then they lie at the centerline

! !        ! maintain scalar at centerline
! !        DO k=0,nzSub+1
! !          phi(1,1,k) = phiIC*ee**(-((x(1)**2 + y(1)**2)/(2.0_dbl*sigma**2)))
! !        END DO

!         DO k=0,nzSub+1
!           DO j=0,nySub+1
!             DO i=0,nxSub+1

!               IF((ABS(x(i)) .LE. 2.51_dbl*xcf) .AND. (ABS(y(j)) .LE. 2.51_dbl*ycf)) THEN				! (5.01 in case it is slightly higher than 5 due to round-off)
!                 phi(i,j,k) = phiIC																						! 2D Gaussian Distribution in x and y (maintain phi=1 at r=2.5*zcf)
!               END IF

!             END DO
!           END DO
!         END DO  

!       END IF

    CASE(INLET) 						! constant scalar at the inlet
 
      ! maintain scalar at inlet (1.0)
      IF(kMin .EQ. 1) THEN	
        DO j=1,ny
          DO i=1,nx
            DO k=0,1

	           phi(i,j,k) = phiIC*ee**(-((z(k)**2)/(2.0_dbl*sigma**2)))										! 1D Gaussian Distribution in z for the inlet nodes (maintain peak at z=0)

            END DO
          END DO
        END DO   
      END IF

    CASE(LINEAR)
      
    CASE(UNIFORM)						! uniform initial distribution 

      ! scalar is not maintained

    CASE DEFAULT
   
       OPEN(1000,FILE="error.txt")
       WRITE(1000,*) "Error in ScalarIC in Setup.f90: sclrIC is not BLOB(1), LINE(2) or INLET(3)..."
       WRITE(1000,*) "sclrIC=", sclrIC
       CLOSE(1000)
       STOP

  END SELECT

!  OPEN(6051,FILE='test-'//sub//'.dat')
!  WRITE(6051,'(A60)') 'VARIABLES = "x" "y" "z" "phi"'
!  WRITE(6051,'(A10,I4,A5,I4,A5,I4,A8)') 'ZONE I=',nxSub+2,' J=',nySub+2,' K=',nzSub+2,'F=POINT'
!  
!  DO k=0,nzSub+1
!    DO j=0,nySub+1
!      DO i=0,nxSub+1
!  
!        WRITE(6051,'(4E15.5)') x(i), y(j), z(k), phi(i,j,k)  
!  
!      END DO
!    END DO
!  END DO
!  
!  CLOSE(6051)

END IF	

!------------------------------------------------
END SUBROUTINE ScalarDistribution_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)			! implements the (moving) bounceback boundary conditions (1st order accurate - Ladd)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk													! radius of current node

rijk = SQRT(x_fine(im1)*x_fine(im1) + y_fine(jm1)*y_fine(jm1))				! radius at current location

cosTheta = x_fine(im1)/rijk											! COS(theta)
sinTheta = y_fine(jm1)/rijk											! SIN(theta)

ub = vel_fine(km1)*cosTheta											! x-component of the velocity at i,j,k
vb = vel_fine(km1)*sinTheta											! y-component of the velocity at i,j,k
wb = 0.0_dbl														! no z-component in this case			

fbb = fplus_fine(bb(m),i,j,k) + 6.0_dbl*wt(m)*1.0*(ub*ex(m) + vb*ey(m) + wb*ez(m))	! bounced back distribution function with added momentum

!------------------------------------------------
END SUBROUTINE BounceBackL_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBack2New_fine(m,i,j,k,im1,jm1,km1,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
! Implemented by Balaji 10/28/2014 using a method similar to Yanxing
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk													! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it			! loop index variables

ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction


IF((node_fine(ip1,jp1,kp1) .EQ. FLUID) .AND. (node_fine(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

!*****************************************************************************
   ! Initial fluid node guess
   x1=x_fine(i)
   y1=y_fine(j)
   z1=z_fine(k)
   
   ! Initial solid node guess
   x2=x_fine(im1)
   y2=y_fine(jm1)
   z2=z_fine(km1)
   
   IF (k.NE.km1) THEN
      DO it=1,10
         ! guess of boundary location 
         xt=(x1+x2)/2.0_dbl
         yt=(y1+y2)/2.0_dbl
         zt=(z1+z2)/2.0_dbl
         
         rt = SQRT(xt*xt + yt*yt)
         !Write(*,*) 'test'
         !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
         ht = ((zt-z_fine(k))*r_fine(km1)+(z_fine(km1)-zt)*r_fine(k))/(z_fine(km1)-z_fine(k))
         !ht = (r(km1)+r(k))/2.0_dbl
         
         IF(rt.GT.ht) then
            x2=xt
            y2=yt
            z2=zt
         ELSE
            x1=xt
            y1=yt
            z1=zt
         END IF
         
      END DO
      x1=x_fine(i)
      y1=y_fine(j)
      z1=z_fine(k)
      
      x2=x_fine(im1)
      y2=y_fine(jm1)
      z2=z_fine(km1)
      
      q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
   ELSE
      DO it=1,10
         ! guess of boundary location 
         xt=(x1+x2)/2.0_dbl
         yt=(y1+y2)/2.0_dbl
         zt=(z1+z2)/2.0_dbl
         
         rt = SQRT(xt*xt + yt*yt)
         !Write(*,*) 'test'
         !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
         !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
         ht = (r_fine(km1)+r_fine(k))/2.0_dbl
         
         IF(rt.GT.ht) then
            x2=xt
            y2=yt
            z2=zt
         ELSE
            x1=xt
            y1=yt
            z1=zt
         END IF
         
      END DO
      x1=x_fine(i)
      y1=y_fine(j)
      z1=z_fine(k)
      
      x2=x_fine(im1)
      y2=y_fine(jm1)
      z2=z_fine(km1)
      
      q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
   ENDIF
   cosTheta=xt/rt
   sinTheta=yt/rt
   IF (k.NE.km1) THEN
      !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
      vt = ((zt-z_fine(k))*vel_fine(km1)+(z_fine(km1)-zt)*vel_fine(k))/(z_fine(km1)-z_fine(k))
   ELSE
      vt = (vel_fine(k)+vel_fine(km1))*0.5_dbl
   ENDIF
   ub = vt*cosTheta										! x-component of the velocity at i,j,k
   vb = vt*sinTheta										! y-component of the velocity at i,j,k
   wb = 0.0_dbl											! no z-component in this case)
   
   ! bounced back distribution function with added momentum
   IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
      fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus_fine(bb(m),i,j,k) 															&
           + (1.0_dbl - 4.0_dbl*q*q)*fplus_fine(bb(m),ip1,jp1,kp1) 													& 
           - q*(1.0_dbl - 2.0_dbl*q)*fplus_fine(bb(m),ip2,jp2,kp2) 													&
           + 6.0_dbl*wt(m)*1.0*(ub*ex(m) + vb*ey(m) + wb*ez(m))
   ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
      fbb = fplus_fine(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
           + ((2.0_dbl*q - 1.0_dbl)*fplus_fine(m,i,j,k))/q																&
           - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus_fine(m,ip1,jp1,kp1)							&
           + (6.0_dbl*wt(m)*1.0*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
   ELSE
      OPEN(1000,FILE='error-'//sub//'.txt')
      WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
      WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
      CLOSE(1000)
      STOP
   END IF
   
ELSE
   
   CALL BounceBackL_fine(m,i,j,k,im1,jm1,km1,fbb)
   
END IF

!------------------------------------------------
END SUBROUTINE BounceBack2New_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalc_fine(m,i,j,k,im1,jm1,km1,q)			! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(OUT) :: q							! distance ratio
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: AB,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept			! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2						! terms used in calculation

! RAY
! point A (current node)
Ax = x_fine(i)
Ay = y_fine(j)
Az = z_fine(k)

! point B (solid node)
Bx = x_fine(im1)
By = y_fine(jm1)
Bz = z_fine(km1)

! distance from A to B
AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! unit vector (d) from point A to point B
dx = (x_fine(im1)-x_fine(i))/AB									! i direction
dy = (y_fine(jm1)-y_fine(j))/AB									! j direction
dz = (z_fine(km1)-z_fine(k))/AB									! k direction

! SURFACE
r1 = r_fine(k)													! radius at k (distance from CL)
r2 = r_fine(km1)													! radius at km1 (distance from CL)
z1 = z_fine(k)													! z-coordinate at k
z2 = z_fine(km1)													! z-coordinate at km1

IF(k .NE. km1) THEN
  slope = (r2-r1)/(z2-z1)								! approximate the surface as a conincal shell (linear between k values)
ELSE
  slope = 0.0_dbl
END IF

intercept = r1 - slope*z1								! z-intercept of the linearly approximated r-equation

! terms used in calculation
slope2 = slope*slope															! slope^2
term1 = Ax*dx + Ay*dy - Az*dz*slope2 - intercept*dz*slope		! reoccuring term
term2 = dx*dx + dy*dy - dz*dz*slope*slope								! reoccuring term

! calculate the distance from the current node (point A) to the wall (point P)
AP = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))

q = AP/AB													! distance ratio

! balaji added
!q=0.5

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "q=",q
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
  WRITE(1000,*) "r1=",r1,"r2=",r2
  WRITE(1000,*) "z1=",z1,"z2=",z2
  WRITE(1000,*) "slope=",slope
  WRITE(1000,*) "term1=",term1,"term2=",term2
  WRITE(1000,*) "intercept=",intercept
  WRITE(1000,*) "AB=",AB,"AP=",AP
  CLOSE(1000)
  STOP
END IF																																

!------------------------------------------------
END SUBROUTINE qCalc_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE PrintFieldsTEST_fine	! print velocity, density, and scalar to output files
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,ii,jj,kk,n		! index variables (local and global)
CHARACTER(7)	:: iter_char				! iteration stored as a character

  ! scale the iteration by 1/10 such that the numbers used in the output file aren't too large
  WRITE(iter_char(1:7),'(I7.7)') iter

  ! store the current iteration in "filenum"
  filenum(fileCount_fine) = iter
  fileCount_fine = fileCount_fine + 1_lng

  ! open the proper output file
  OPEN(60,FILE='out-TEST-'//iter_char//'-'//sub//'.dat')
  WRITE(60,*) 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
  WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nxSub,' J=',nySub,' K=',nzSub,'F=POINT'

  DO k=1,nzSub_fine
    DO j=1,nySub_fine
      DO i=1,nxSub_fine

         ! convert local i,j,k, to global ii,jj,kk
         ii = ((iMin - 1_lng) + i)
         jj = ((jMin - 1_lng) + j)
         kk = ((kMin - 1_lng) + k)

         WRITE(60,'(8E15.5,I6)') xx_fine(ii), yy_fine(jj), zz_fine(kk), u_fine(i,j,k)*vcf_fine, v_fine(i,j,k)*vcf_fine, w_fine(i,j,k)*vcf_fine, (rho_fine(i,j,k)-denL)*dcf_fine*pcf_fine,	&
                                 phi_fine(i,j,k), node_fine(i,j,k)

      END DO
    END DO
  END DO

  CLOSE(60)

!------------------------------------------------
END SUBROUTINE PrintFieldsTEST_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarBC_fine(m,i,j,k,im1,jm1,km1,phiBC)								! implements the scalar BCs 
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1								! index variables
REAL(dbl), INTENT(OUT) :: phiBC     											! scalar contribution from the boundary condition
INTEGER :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
REAL(dbl) :: q																			! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB																! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: feq_m																	! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m																! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: cosTheta, sinTheta													! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb																! wall velocity (x-, y-, z- components)
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing)
INTEGER   :: it

! Initial fluid node guess
x1=x_fine(i)
y1=y_fine(j)
z1=z_fine(k)

! Initial solid node guess
x2=x_fine(im1)
y2=y_fine(jm1)
z2=z_fine(km1)

IF (k.NE.km1) THEN
   DO it=1,10
      ! guess of boundary location 
      xt=(x1+x2)/2.0_dbl
      yt=(y1+y2)/2.0_dbl
      zt=(z1+z2)/2.0_dbl
      
      rt = SQRT(xt*xt + yt*yt)
      !Write(*,*) 'test'
      !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
      ht = ((zt-z_fine(k))*r_fine(km1)+(z_fine(km1)-zt)*r_fine(k))/(z_fine(km1)-z_fine(k))
      !ht = (r(km1)+r(k))/2.0_dbl
      
      IF(rt.GT.ht) then
         x2=xt
         y2=yt
         z2=zt
      ELSE
         x1=xt
         y1=yt
         z1=zt
      END IF
      
   END DO
   x1=x_fine(i)
   y1=y_fine(j)
   z1=z_fine(k)
   
   x2=x_fine(im1)
   y2=y_fine(jm1)
   z2=z_fine(km1)
   
   q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
   !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
ELSE
   DO it=1,10
      ! guess of boundary location 
      xt=(x1+x2)/2.0_dbl
      yt=(y1+y2)/2.0_dbl
      zt=(z1+z2)/2.0_dbl
      
      rt = SQRT(xt*xt + yt*yt)
      !Write(*,*) 'test'
      !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
      !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
      ht = (r_fine(km1)+r_fine(k))/2.0_dbl
      
      IF(rt.GT.ht) then
         x2=xt
         y2=yt
         z2=zt
      ELSE
         x1=xt
         y1=yt
         z1=zt
      END IF
      
   END DO
   x1=x_fine(i)
   y1=y_fine(j)
   z1=z_fine(k)
   
   x2=x_fine(im1)
   y2=y_fine(jm1)
   z2=z_fine(km1)
   
   q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
   !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
ENDIF
cosTheta=xt/rt
sinTheta=yt/rt
IF (k.NE.km1) THEN
   !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
   vt = ((zt-z_fine(k))*vel_fine(km1)+(z_fine(km1)-zt)*vel_fine(k))/(z_fine(km1)-z_fine(k))
ELSE
   vt = (vel_fine(k)+vel_fine(km1))*0.5_dbl
ENDIF
ub = vt*cosTheta										! x-component of the velocity at i,j,k
vb = vt*sinTheta										! y-component of the velocity at i,j,k
wb = 0.0_dbl											! no z-component in this case)

!---------------------------------------------------------------------------------------------------
!----- Computing phi at the wall in case of Dirichlet BC -------------------------------------------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .EQ. 0.0) then
   phiWall= coeffConst/coeffPhi
END IF

!---------------------------------------------------------------------------------------------------
!----- Estimating  phi at the wall in case of Neumann or Mixed BC ----------------------------------
!----- Two Fluid nodes are needed for extrapolation ------------------------------------------------
!----- 1st node; A, is adjacent to the boundary ----------------------------------------------------
!----- 2nd node; B, is a fluid neighbor of A, in the direction closest to geometry-normal ----------
!---------------------------------------------------------------------------------------------------
! IF (coeffGrad .NE. 0.0) then
!    Geom_norm_x= 1.0
!    Geom_norm_y= 0.0
!    Geom_norm_z= 0.0
!    n_prod_max= 0.0_dbl
!    !----- Finding the mth direction closest to normal vector
!    !----- Only one fluid neighboring node is needed for interpolation
!    DO mm=0,NumDistDirs
!       ip1 = i+ ex(mm)
!       jp1 = j+ ey(mm)
!       kp1 = k+ ez(mm)
!       IF (node_fine(ip1,jp1,kp1) .EQ. FLUID) THEN
!          n_prod= abs( Geom_norm_x * (ip1-i) + Geom_norm_y * (jp1-j) + Geom_norm_z * (kp1-k))
!          IF (n_prod_max .LT. n_prod) THEN
!             n_prod_max= n_prod
!             iB= ip1
!             jB= jp1
!             kB= kp1
!          END IF
!       END IF
!    END DO
!    phiWall= ((phiTemp_fine(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) -    &
!         (phiTemp_fine(iB,jB,kB)*q*q/(1.0+2.0*q)) -          &
!         (q*(1+q)/(1+2.0*q))* (coeffConst/coeffGrad)    )  &
!         /(1.0- (q*(1+q)/(1+2.0*q))*(coeffPhi/coeffGrad))
! END IF


! neighboring node (fluid side)	
ip1 = i + ex(m) 																		! i + 1
jp1 = j + ey(m)																		! j + 1
kp1 = k + ez(m)																		! k + 1

! if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node_fine(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

! assign values to boundary (density, scalar, f)
rhoB = (rho_fine(i,j,k) - rho_fine(ip1,jp1,kp1))*(1+q) + rho_fine(ip1,jp1,kp1)		! extrapolate the density
CALL Equilibrium_LOCAL_fine(m,rhoB,ub,vb,wb,feq_m)			        ! calculate the equibrium distribution function in the mth direction

!! Balaji added for sero flux BC. Otherwise set to constant value for Dirichlet BC
!phiWall = (phi(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - (phi(ip1,jp1,kp1)*q*q/(1.0+2.0*q)) 	! calculate phiWall for flux BC (eq. 28 in paper)

! find the contribution of scalar streamed from the wall to the current node (i,j,k), and from the current node to the next neighboring node (ip1,jp1,kp1)
phiB		= (feq_m/rhoB - wt(m)*Delta_fine)*phiWall								! contribution from the wall in the mth direction (zero if phiWall=0)
phiijk_m	= (fplus_fine(m,i,j,k)/rho_fine(i,j,k) - wt(m)*Delta_fine)*phiTemp_fine(i,j,k)	! contribution from the current node to the next node in the mth direction

! if q is too small, the extrapolation to phiBC can create a large error...
IF(q .LT. 0.25) THEN
  q = 0.25_dbl  																		! approximate the distance ratio as 0.25
END IF

! extrapolate using phiB and phijk_m to obtain contribution from the solid node to the current node
phiBC		= ((phiB - phiijk_m)/q) + phiB										! extrapolated scalar value at the solid node, using q

!------------------------------------------------
END SUBROUTINE ScalarBC_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Equilibrium_LOCAL_fine(m,rhoijk,uijk,vijk,wijk,feq_m)		! calculate and store the equilibrium distribution function
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN)	:: m											! distribution function direction		
REAL(dbl), INTENT(IN)		:: rhoijk,uijk,vijk,wijk				! density, and velocity components at the current node
REAL(dbl), INTENT(OUT)		:: feq_m										! density, and velocity components at the current node

REAL(dbl)						:: UU,ue,ve,we,Usum						! precalculated quantities for use in the feq equation

UU 	= uijk*uijk + vijk*vijk + wijk*wijk								! U . U
      
ue		= uijk*ex(m)															! u . e
ve		= vijk*ey(m)															! v . e
we		= wijk*ez(m)															! w . e

Usum	= ue + ve + we															! U . e
        
feq_m	= (wt(m)*rhoijk)*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function in the mth direction
        
!------------------------------------------------
END SUBROUTINE Equilibrium_LOCAL_fine
!------------------------------------------------

!================================================
END MODULE ICBC_FINE
!================================================
