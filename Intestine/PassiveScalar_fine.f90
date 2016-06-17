!==================================================================================================
MODULE PassiveScalar_fine		! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs)
!==================================================================================================
USE SetPrecision
USE Setup
USE ICBC
USE Setup_fine
USE ICBC_fine
USE PassiveScalar

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_Setup_fine	! sets up the passive scalar component
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER(lng) :: i,j,k  ! index variables
REAL(dbl) :: numFluids, numFluids_l! number of fluid nodes in the domain
REAL(dbl)    :: phiDomain, phiDomain_l, phiIC! current amount of scalar in the domain
INTEGER   :: mpierr


! initialize arrays
phi_fine  = 0.0_dbl		! scalar
phiTemp_fine  = 0.0_dbl		! temporary scalar

! scalar parameters
Dmcf_fine = (zcf_fine*zcf_fine)/tcf_fine		! conversion factor for diffusivity
Delta_fine = 1.0_dbl - gridRatio*(1.0_dbl - Delta)	! scalar diffusion parameter
write(31,*) 'nuL = ', nuL, ' Dm = ', Dm, ' Delta = ', Delta, ' Delta_fine = ', Delta_fine

CALL ScalarDistribution_fine			! sets/maintains initial distributions of scalar [MODULE: ICBC_fine.f90]

!----- Calculate the amount of scalar in the domain
numFluids = 0.0_dbl
phiDomain = 0.0_dbl
numFluids_l = 0.0_dbl
phiDomain_l = 0.0_dbl

DO k=1,nzSub_fine
   DO j=1,nySub_fine
      DO i=1,nxSub_fine
         
         IF(node_fine(i,j,k) .EQ. FLUID) THEN
            phiDomain_l = phiDomain_l + phi_fine(i,j,k)
            numFluids_l = numFluids_l + 1.0_dbl
         END IF

      END DO
   END DO
END DO

DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            phiDomain_l = phiDomain_l + (1.0-flagNodeIntersectFine(i,j,k)) * phi(i,j,k) * gridRatio * gridRatio * gridRatio
            numFluids_l = numFluids_l + (1.0-flagNodeIntersectFine(i,j,k)) * gridRatio * gridRatio * gridRatio
         END IF
      END DO
   END DO
END DO
CALL MPI_ALLREDUCE(phiDomain_l , phiDomain , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
Drug_Initial = phiDomain * zcf_fine * zcf_fine * zcf_fine

!------------------------------------------------
END SUBROUTINE Scalar_Setup_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_fine								! calculates the evolution of scalar in the domain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: phiBC				! scalar contribution from boundary
REAL(dbl) :: zcf3 
REAL(dbl) :: tmp
CALL ScalarDistribution_fine			! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]


zcf3 = zcf_fine * zcf_fine * zcf_fine
! store the previous scalar values
phiTemp_fine = phi_fine

tmp = sum(phi_fine(:,:,:)) * zcf3 

! Stream the scalar
DO k=1,nzSub_fine
  DO j=1,nySub_fine
    DO i=1,nxSub_fine
      
      IF(node_fine(i,j,k) .EQ. FLUID) THEN
      
        phi_fine(i,j,k) = Delta_fine*phiTemp_fine(i,j,k)
        phi_fine(i,j,k) = phi_fine(i,j,k) + delphi_particle_fine(i,j,k)
        
        DO m=0,NumDistDirs
      
!          i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

          IF(node_fine(im1,jm1,km1) .EQ. FLUID) THEN
            phi_fine(i,j,k) = phi_fine(i,j,k) + (fplus_fine(m,im1,jm1,km1)/rho_fine(im1,jm1,km1) - wt(m)*Delta_fine)*phiTemp_fine(im1,jm1,km1)
          ELSE IF (node_fine(im1,jm1,km1) .EQ. COARSEMESH) THEN
!             if ( abs(rho_fine(im1,jm1,km1) - 1.0) .gt. 0.5) then
!               rho_fine(im1,jm1,km1) = 1.0
!            end if
            phi_fine(i,j,k) = phi_fine(i,j,k) + (fplus_fine(m,im1,jm1,km1)/rho_fine(im1,jm1,km1) - wt(m)*Delta_fine)*phiTemp_fine(im1,jm1,km1)
          ELSE IF(node_fine(im1,jm1,km1) .EQ. SOLID) THEN ! macro- boundary
            CALL ScalarBC_fine(m,i,j,k,im1,jm1,km1,phiBC) ! implement scalar boundary condition (using BB f's)	[MODULE: ICBC]
            phi_fine(i,j,k) = phi_fine(i,j,k) + phiBC     
            CALL AbsorbedScalarS_fine(i,j,k,m,im1,jm1,km1,phiBC)	! measure the absorption rate
          ELSE
            OPEN(1000,FILE="error_fine.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar_fine.f90 at Line 66: node(im1,jm1,km1) is out of range"
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

       	!  fix spurious oscillations in moment propagation method for high Sc #s
        IF(phi_fine(i,j,k) .LT. 0.0_dbl) THEN
!           if(subIter .eq. gridRatio) then
              Negative_phi_Counter_l = Negative_phi_Counter_l + 1
              Negative_phi_Total_l   = Negative_phi_Total_l + phi_fine(i,j,k) * zcf3 
              IF (phi_fine(i,j,k) .LT. Negative_phi_Worst) THEN
                 Negative_phi_Worst_l = phi_fine(i,j,k)
              ENDIF
!           end if
           phi_fine(i,j,k) = 0.0_dbl
           
        ELSE IF (phi_fine(i,j,k) .gt. Cs_mol) THEN
           if(subIter .eq. gridRatio) then
              Over_Sat_Counter_l = Over_Sat_Counter_l + 1
              IF (phi_fine(i,j,k) .GT. Largest_Phi_l) THEN
                 Largest_Phi_l = phi_fine(i,j,k)
              ENDIF
           end if

        END IF
     END IF

     END DO
  END DO
END DO

! Add the amount of scalar absorbed through the outer and villous surfaces
phiAbsorbed_fine = phiAbsorbedS_fine + phiAbsorbedV_fine ! total amount of scalar absorbed up to current time
      
!------------------------------------------------
END SUBROUTINE Scalar_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarS_fine(i,j,k,m,im1,jm1,km1,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m,im1,jm1,km1	! index variables
REAL(dbl), INTENT(IN) :: phiBC ! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN  ! scalar values exchanged with the wall
INTEGER :: ip1,jp1,kp1  ! neighboring nodes (2 away from the wall)
REAL(dbl) :: q		! distance ratio from the current node to the solid node
REAL(dbl) :: feq_m	! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m	! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: cosTheta, sinTheta	 ! COS(theta), SIN(theta)
REAL(dbl) :: ubb, vbb, wbb	 ! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk	! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt  ! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it	! loop index variables
REAL(dbl)    :: feq_AO_u0
REAL(dbl)    :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar
REAL(dbl)    :: rhoA, PkA, feq_A
REAL(dbl)    :: fPlusBstar, rhoBstar, phiBstar, PkBstar


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
   ubb = 0.0
   vbb = 0.0
   wbb = 0.0
   
   CALL Equilibrium_LOCAL(bb(m),rho_fine(i,j,k),ubb,vbb,wbb,feq_AO_u0)
   phiOUT= (feq_AO_u0/rho_fine(i,j,k) - wt(bb(m))*Delta_fine)*phiTemp_fine(i,j,k)

   !---------------------------------------------------------------------------------------------------
   !---- Conmputing phiIN------------------------------------------------------------------------------
   !---------------------------------------------------------------------------------------------------
   !----- neighboring node (fluid side)
   ip1 = i + ex(m)
   jp1 = j + ey(m)
   kp1 = k + ez(m)
   IF(node_fine(ip1,jp1,kp1) .NE. FLUID) THEN
      ip1 = i
      jp1 = j
      kp1 = k
   END IF

   !----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
   rhoAstar= (rho_fine(i,j,k)- rho_fine(ip1,jp1,kp1))*(1+q)+ rho_fine(ip1,jp1,kp1)! extrapolate the density
   CALL Equilibrium_LOCAL_fine(m,rhoAstar,ubb,vbb,wbb,feq_Astar)! calculate the equibrium distribution function in the mth direction
   phiAstar= phiWall! getting phi at the solid surface
   PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta_fine)*phiAstar! contribution from the wall in mth direction (0 if phiWall=0)

   !---- Modification for moving boundary in case of using only A and A* for BC
   rhoA= rho_fine(i,j,k)
   CALL Equilibrium_LOCAL_fine(m,rhoA,ubb,vbb,wbb,feq_A)
   PkA= (feq_A/rhoA - wt(m)*Delta_fine)*phiTemp_fine(i,j,k)
   IF(q .LT. 0.25) THEN
      q = 0.25_dbl
   END IF
   phiIN   = ((PkAstar - PkA)/q) + PkAstar

   phiAbsorbedS_fine = phiAbsorbedS_fine + (phiOUT-phiIN)! scalar absorbed at current location in mth direction

!------------------------------------------------
END SUBROUTINE AbsorbedScalarS_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarInOut_fine		! measure scalar that has left or entered the domain through the inlet or outlet
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
              phiIN = (fplus(bb(f_Comps(iComm,m)),im1,jm1,km1)/rho(im1,jm1,km1) - wt(bb(f_Comps(iComm,m)))*Delta_fine)	&								! scalar contribution from inlet/outlet to current node
                      *phiTemp(im1,jm1,km1)		
              phiOUT	= (fplus(f_Comps(iComm,m),i,j,k)/rho(i,j,k) - wt(f_Comps(iComm,m))*Delta_fine)*phiTemp(i,j,k)										! scalar contribution from current node to inlet/outlet
              phiInOut = phiInOut + (phiOUT - phiIN)
            END IF

          END DO
 
        END IF

      END DO
    END DO
  
  END IF

END DO

!------------------------------------------------
END SUBROUTINE ScalarInOut_fine
!------------------------------------------------

!================================================
END MODULE PassiveScalar_fine
!================================================
