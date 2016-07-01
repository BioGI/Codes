!==================================================================================================
MODULE PassiveScalar		! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs)
!==================================================================================================
USE SetPrecision
USE Setup
USE Setup_fine
USE ICBC
USE Geometry_fine

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_Setup	! sets up the passive scalar component
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! initialize arrays
phi    		= 0.0_dbl					! scalar
phiTemp		= 0.0_dbl					! temporary scalar
Over_Sat_Counter_l = 0
Over_Sat_Counter = 0
Largest_phi_l = 0.0_dbl
Largest_phi = 0.0_dbl
Negative_phi_Counter_l = 0
Negative_phi_Worst_l   = 0.0_dbl
Negative_phi_Counter   = 0
Negative_phi_Worst     = 0.0_dbl

! scalar parameters
Dm   			= nuL/Sc				! binary molecular diffusivity (scalar in fluid)
Dmcf			= (zcf*zcf)/tcf				! conversion factor for diffusivity
Delta			= 1.0_dbl - 6.0_dbl*Dm			! scalar diffusion parameter

! set scalar values at sources/sinks
phiWall		= 0.0_dbl								! value of scalar at the boundary

! set the scalar standard devation for gaussian distributions
sigma = 0.1_dbl*D										! 1/10 of the Diameter

! determine scalar starting iteration
phiStart = 0
!phiStart 	= NINT((phiPer*Tmix)/tcf)
!IF (phiPer.EQ.0.0) THEN
!	phiStart 	= NINT((phiPer*Tmix)/tcf) ! Balaji changed this to add 1 as for phiPer=0, phiSTart=0. But iter never has a value 0.
!ENDIF

CALL ScalarDistribution						! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]

!------------------------------------------------
END SUBROUTINE Scalar_Setup
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar								! calculates the evolution of scalar in the domain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: phiBC							! scalar contribution from boundary
REAL(dbl) :: zcf3
REAL(dbl) :: phiMaxInterface

Over_Sat_Counter_l = 0
Over_Sat_Counter = 0
Largest_phi_l = 0.0_dbl
Largest_phi = 0.0_dbl
Negative_phi_Counter_l = 0
Negative_phi_Worst_l   = 0.0_dbl
Negative_phi_Counter   = 0
Negative_phi_Worst     = 0.0_dbl

zcf3 = zcf * zcf * zcf

flagNodeIntersectCoarse = 0.0_dbl
CALL ScalarDistribution						! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]

! store the previous scalar values
phiTemp = phi

phiMaxInterface = 0.0
! Stream the scalar
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
      
      IF(node(i,j,k) .EQ. FLUID) THEN
      
        phi(i,j,k) = Delta*phiTemp(i,j,k)
	phi(i,j,k) = phi(i,j,k) + delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release

         DO m=0,NumDistDirs
      
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
             if( (i .eq. 46) .and. (j .eq. 50) .and. (k .eq. 1) ) then
!                write(31,*) m, phiTemp(i,j,k), fplus(m,im1,jm1,km1), rho(im1,jm1,km1), wt(m), Delta, phiTemp(im1,jm1,km1)
             end if
             phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1)
          ELSE IF(node(im1,jm1,km1) .EQ. FINEMESH) THEN
             if (phiTemp(im1,jm1,km1) .gt. phiMaxInterface) then
                phiMaxInterface = phiTemp(im1,jm1,km1)
             end if
             phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1)
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN															! macro- boundary
            CALL ScalarBC(m,i,j,k,im1,jm1,km1,phiBC)															! implement scalar boundary condition (using BB f's)	[MODULE: ICBC]
            phi(i,j,k) = phi(i,j,k) + phiBC     
            CALL FlagFineMeshNodesIntersectingWithCoarseMeshNodes(i,j,k)
            CALL AbsorbedScalarS(i,j,k,m,im1,jm1,km1,phiBC)	     ! measure the absorption rate
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter",iter
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

!	phi(i,j,k) = phi(i,j,k) + delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release

       	!  fix spurious oscillations in moment propagation method for high Sc #s
        IF(phi(i,j,k) .LT. 0.0_dbl) THEN
           Negative_phi_Counter_l = Negative_phi_Counter_l + 1.0
           Negative_phi_Total_l = Negative_phi_Total_l + phi(i,j,k) * (1.0-flagNodeIntersectFine(i,j,k)) * zcf3
           IF (phi(i,j,k) .LT. Negative_phi_Worst) THEN
              Negative_phi_Worst_l = phi(i,j,k)
           ENDIF
           phi(i,j,k) = 0.0_dbl
        ELSE IF (phi(i,j,k) .gt. Cs_mol) THEN
           Over_Sat_Counter_l = Over_Sat_Counter_l + 1
           IF (phi(i,j,k) .GT. Largest_Phi_l) THEN
              Largest_Phi_l = phi(i,j,k)
           ENDIF           
        END IF

      END IF

    END DO
  END DO
END DO


! Add the amount of scalar absorbed through the outer and villous surfaces
phiAbsorbed_coarse = 	phiAbsorbedS_coarse + phiAbsorbedV_coarse																		! total amount of scalar absorbed up to current time
      
!------------------------------------------------
END SUBROUTINE Scalar
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarS(i,j,k,m,im1,jm1,km1,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m,im1,jm1,km1				! index variables
REAL(dbl), INTENT(IN) :: phiBC     				! scalar contribution from the boundary condition

REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall
INTEGER(dbl) :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
REAL(dbl) :: q																			! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB																! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: feq_m																	! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m																! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: rijk													! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it			! loop index variables
REAL(dbl)    :: ubb,vbb,wbb
REAL(dbl)    :: feq_AO_u0
REAL(dbl)    :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar
REAL(dbl)    :: rhoA, PkA, feq_A
REAL(dbl)    :: fPlusBstar, rhoBstar, phiBstar, PkBstar


   ! Initial fluid node guess
   x1=x(i)
   y1=y(j)
   z1=z(k)
   
   ! Initial solid node guess
   x2=x(im1)
   y2=y(jm1)
   z2=z(km1)
   
   IF (k.NE.km1) THEN
      DO it=1,10
         ! guess of boundary location 
         xt=(x1+x2)/2.0_dbl
         yt=(y1+y2)/2.0_dbl
         zt=(z1+z2)/2.0_dbl
         
         rt = SQRT(xt*xt + yt*yt)
         !Write(*,*) 'test'
         !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
         ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
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
      x1=x(i)
      y1=y(j)
      z1=z(k)
      
      x2=x(im1)
      y2=y(jm1)
      z2=z(km1)
      
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
         ht = (r(km1)+r(k))/2.0_dbl
         
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
      x1=x(i)
      y1=y(j)
      z1=z(k)
      
      x2=x(im1)
      y2=y(jm1)
      z2=z(km1)
      
      q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
   ENDIF
   ubb= 0.0_dbl
   vbb= 0.0_dbl
   wbb= 0.0_dbl

!---------------------------------------------------------------------------------------------------
!----- Computing phiOUT ----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
CALL Equilibrium_LOCAL(bb(m),rho(i,j,k),ubb,vbb,wbb,feq_AO_u0)
phiOUT= (feq_AO_u0/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)

!---------------------------------------------------------------------------------------------------
!---- Conmputing phiIN------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!----- neighboring node (fluid side)
ip1 = i + ex(m)
jp1 = j + ey(m)
kp1 = k + ez(m)
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
   ip1 = i
   jp1 = j
   kp1 = k
END IF

!----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
rhoAstar= (rho(i,j,k)- rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoAstar,ubb,vbb,wbb,feq_Astar)! calculate the equibrium distribution function in the mth direction
phiAstar= phiWall! getting phi at the solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar! contribution from the wall in mth direction (0 if phiWall=0)

!---- Modification for moving boundary in case of using only A and A* for BC
rhoA= rho(i,j,k)
CALL Equilibrium_LOCAL(m,rhoA,ubb,vbb,wbb,feq_A)
PkA= (feq_A/rhoA - wt(m)*Delta)*phiTemp(i,j,k)
IF(q .LT. 0.25) THEN
   q = 0.25_dbl
END IF
phiIN   = ((PkAstar - PkA)/q) + PkAstar

phiAbsorbedS_coarse = phiAbsorbedS_coarse + (phiOUT-phiIN)! scalar absorbed at current location in mth direction


!------------------------------------------------
END SUBROUTINE AbsorbedScalarS
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarV(i,j,k,m,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m				! index variables
REAL(dbl), INTENT(IN) :: phiBC     				! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

phiIN 	= phiBC																						! contribution from the wall to the crrent node (in)
phiOUT	= (fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)		! contribution to the wall from the current node (out)

phiAbsorbedV_coarse = phiAbsorbedV_coarse + (phiOUT - phiIN)												! add the amount of scalar that has been absorbed at the current location in the current direction

!------------------------------------------------
END SUBROUTINE AbsorbedScalarV
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarInOut								! measure scalar that has left or entered the domain through the inlet or outlet
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1,iComm		! index variables
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

! XY Faces (z-direction)
DO iComm=5,6

  IF(SubID(iComm) .EQ. 0) THEN					

    k = XY_SendIndex(iComm)						! k index

    DO j=1,nySub
      DO i=1,nxSub

        IF(node(i,j,k) .NE. SOLID) THEN

          DO m=1,NumFs_face

            ! i,j,k location of neighboring node
            im1 = i - ex(bb(f_Comps(iComm,m)))
            jm1 = j - ey(bb(f_Comps(iComm,m)))
            km1 = k - ez(bb(f_Comps(iComm,m)))

            IF(node(im1,jm1,km1) .NE. SOLID) THEN
              phiIN = (fplus(bb(f_Comps(iComm,m)),im1,jm1,km1)/rho(im1,jm1,km1) - wt(bb(f_Comps(iComm,m)))*Delta)	&								! scalar contribution from inlet/outlet to current node
                      *phiTemp(im1,jm1,km1)		
              phiOUT	= (fplus(f_Comps(iComm,m),i,j,k)/rho(i,j,k) - wt(f_Comps(iComm,m))*Delta)*phiTemp(i,j,k)										! scalar contribution from current node to inlet/outlet
              phiInOut = phiInOut + (phiOUT - phiIN)
            END IF

          END DO
 
        END IF

      END DO
    END DO
  
  END IF

END DO

!------------------------------------------------
END SUBROUTINE ScalarInOut
!------------------------------------------------
!!------------------------------------------------
!SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
!! Called by Particle_Track (LBM.f90) to get delphi_particle
!!------------------------------------------------
!IMPLICIT NONE
!INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
!REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
!
!DO i=1,np
!	!ix0=FLOOR(xp(i))
!	!ix1=CEILING(xp(i))
!	!iy0=FLOOR(yp(i))
!	!iy1=CEILING(yp(i))
!	!iz0=FLOOR(zp(i))
!	!iz1=CEILING(zp(i))
!
!	ix0=FLOOR(xp(i))
!	ix1=FLOOR(xp(i))+1_lng
!	iy0=FLOOR(yp(i))
!	iy1=FLOOR(yp(i))+1_lng
!	iz0=FLOOR(zp(i))
!	iz1=FLOOR(zp(i))+1_lng
!
!	xd=(xp(i)-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
!	yd=(yp(i)-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
!	zd=(zp(i)-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
!
!	delphi_particle(ix0,iy0,iz0)=delNB(i)
!	
!!	!yd=0.0_dbl ! TEST: used to keep particle motion plainly 2-D ! Balaji added
!!
!!	! u-interpolation
!!	! Do first level linear interpolation in x-direction
!!	c00 = u(ix0,iy0,iz0)*(1.0_dbl-xd)+u(ix1,iy0,iz0)*xd	
!!	c01 = u(ix0,iy0,iz1)*(1.0_dbl-xd)+u(ix1,iy0,iz1)*xd	
!!	c10 = u(ix0,iy1,iz0)*(1.0_dbl-xd)+u(ix1,iy1,iz0)*xd	
!!	c11 = u(ix0,iy1,iz1)*(1.0_dbl-xd)+u(ix1,iy1,iz1)*xd	
!!	! Do second level linear interpolation in y-direction
!!	c0  = c00*(1.0_dbl-yd)+c10*yd
!!	c1  = c01*(1.0_dbl-yd)+c11*yd
!!	! Do third level linear interpolation in z-direction
!!	c   = c0*(1.0_dbl-zd)+c1*zd
!!        up(i)=c
!!
!!	! v-interpolation
!!	! Do first level linear interpolation in x-direction
!!	c00 = v(ix0,iy0,iz0)*(1.0_dbl-xd)+v(ix1,iy0,iz0)*xd	
!!	c01 = v(ix0,iy0,iz1)*(1.0_dbl-xd)+v(ix1,iy0,iz1)*xd	
!!	c10 = v(ix0,iy1,iz0)*(1.0_dbl-xd)+v(ix1,iy1,iz0)*xd	
!!	c11 = v(ix0,iy1,iz1)*(1.0_dbl-xd)+v(ix1,iy1,iz1)*xd	
!!	! Do second level linear interpolation in y-direction
!!	c0  = c00*(1.0_dbl-yd)+c10*yd
!!	c1  = c01*(1.0_dbl-yd)+c11*yd
!!	! Do third level linear interpolation in z-direction
!!	c   = c0*(1.0_dbl-zd)+c1*zd
!!        vp(i)=c
!!
!!	! w-interpolation
!!	! Do first level linear interpolation in x-direction
!!	c00 = w(ix0,iy0,iz0)*(1.0_dbl-xd)+w(ix1,iy0,iz0)*xd	
!!	c01 = w(ix0,iy0,iz1)*(1.0_dbl-xd)+w(ix1,iy0,iz1)*xd	
!!	c10 = w(ix0,iy1,iz0)*(1.0_dbl-xd)+w(ix1,iy1,iz0)*xd	
!!	c11 = w(ix0,iy1,iz1)*(1.0_dbl-xd)+w(ix1,iy1,iz1)*xd	
!!	! Do second level linear interpolation in y-direction
!!	c0  = c00*(1.0_dbl-yd)+c10*yd
!!	c1  = c01*(1.0_dbl-yd)+c11*yd
!!	! Do third level linear interpolation in z-direction
!!	c   = c0*(1.0_dbl-zd)+c1*zd
!!        wp(i)=c
!!
!!	!up(i)=0.5*(u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+u(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
!!	!vp(i)=0.5*(v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+v(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
!!	!wp(i)=0.5*(w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+w(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
!!	!up(i)=u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
!!	!vp(i)=v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
!!	!wp(i)=w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
!ENDDO
!!------------------------------------------------
!END SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
!!------------------------------------------------

!================================================
END MODULE PassiveScalar
!================================================
