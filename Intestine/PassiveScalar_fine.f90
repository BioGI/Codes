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

! initialize arrays
phi_fine  = 0.0_dbl		! scalar
phiTemp_fine  = 0.0_dbl		! temporary scalar

! scalar parameters
Dmcf_fine = (zcf_fine*zcf_fine)/tcf_fine		! conversion factor for diffusivity
Delta_fine = 1.0_dbl - gridRatio*(1.0_dbl - Delta)	! scalar diffusion parameter
write(31,*) 'nuL = ', nuL, ' Dm = ', Dm, ' Delta = ', Delta, ' Delta_fine = ', Delta_fine

CALL ScalarDistribution_fine			! sets/maintains initial distributions of scalar [MODULE: ICBC_fine.f90]

!------------------------------------------------
END SUBROUTINE Scalar_Setup_fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_fine								! calculates the evolution of scalar in the domain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: phiBC				! scalar contribution from boundary
REAL(dbl) :: tmp
REAL(dbl) :: zcf3
CALL ScalarDistribution_fine			! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]

zcf3 = zcf_fine * zcf_fine * zcf_fine

! store the previous scalar values
phiTemp_fine = phi_fine


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
            phi_fine(i,j,k) = phi_fine(i,j,k) + (fplus_fine(m,im1,jm1,km1)/(rho_fine(im1,jm1,km1)+1e-10) - wt(m)*Delta_fine)*phiTemp_fine(im1,jm1,km1)
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

!	phi_fine(i,j,k) = phi_fine(i,j,k) + delphi_particle_fine(i,j,k) ! Balaji added to introduce drug concentration release

       	!  fix spurious oscillations in moment propagation method for high Sc #s
        IF(phi_fine(i,j,k) .LT. 0.0_dbl) THEN

           Negative_phi_Counter = Negative_phi_Counter + 1.0
           Negative_phi_Total   = Negative_phi_Total + phi_fine(i,j,k) * zcf3 
           IF (phi_fine(i,j,k) .LT. Negative_phi_Worst) THEN
              Negative_phi_Worst = phi_fine(i,j,k)
           ENDIF
           
           phi_fine(i,j,k) = 0.0_dbl
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

INTEGER(lng), INTENT(IN) :: i,j,k,im1,jm1,km1,m				! index variables
REAL(dbl), INTENT(IN) :: phiBC     				! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

REAL(dbl)    :: feq_AO_u0
REAL(dbl)    :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar
REAL(dbl)    :: rhoA, PkA, feq_A
REAL(dbl)    :: fPlusBstar, rhoBstar, phiBstar, PkBstar
REAL(dbl)    :: ub,vb,wb, ubb,vbb,wbb
REAL(dbl)    :: q
REAL(dbl)    :: xt,yt,zt,rt			             	! boundary coordinate

CALL  qCalc_iter_fine(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)

ubb= 0.0_dbl
vbb= 0.0_dbl
wbb= 0.0_dbl

!---------------------------------------------------------------------------------------------------
!----- Computing phiOUT ----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
CALL Equilibrium_LOCAL(bb(m),rho_fine(i,j,k),ubb,vbb,wbb,feq_AO_u0)
phiOUT= (feq_AO_u0/rho_fine(i,j,k) - wt(bb(m))*Delta)*phiTemp_fine(i,j,k)

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
rhoAstar= (rho_fine(i,j,k)- rho_fine(ip1,jp1,kp1))*(1+q)+ rho_fine(ip1,jp1,kp1)	! extrapolate the density
CALL Equilibrium_LOCAL_fine(m,rhoAstar,ubb,vbb,wbb,feq_Astar)		! calculate the equibrium distribution function in the mth direction
phiAstar= phiWall							! getting phi at the solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! contribution from the wall in mth direction (0 if phiWall=0)

!------ Computing values at B* & scalar streamed from B* (Chpter 3 paper)
!rhoBstar=   (1-q)*rho_fine(ip1,jp1,kp1)     + q*rho_fine(i,j,k)
!CALL Equilibrium_LOCAL(m,rhoBstar,ubb,vbb,wbb,feq_Bstar)
!phiBstar=   (1-q)*phiTemp_fine(ip1,jp1,kp1) + q*phiTemp_fine(i,j,k)
!PkBstar=    (feq_Bstar/rhoBstar - wt(m)*Delta)*phiBstar

!phiIN= PkAstar+ (PkAstar- PkBstar)*(1-q)


!---- Modification for moving boundary in case of using only A and A* for BC
rhoA= rho_fine(i,j,k)
CALL Equilibrium_LOCAL_fine(m,rhoA,ubb,vbb,wbb,feq_A) 
PkA= (feq_A/rhoA - wt(m)*Delta)*phiTemp_fine(i,j,k) 
IF(q .LT. 0.25) THEN
  q = 0.25_dbl
END IF 
phiIN   = ((PkAstar - PkA)/q) + PkAstar

!--- No Modifications in book-keeping for moving boundaries
!phiIN= phiBC                                                    	 ! contribution from wall to crrent node _fine(in)
!phiOUT= (fplus(bb(m),i,j,k)/rho_fine(i,j,k)-wt(bb(m))*Delta)*phiTemp_fine(i,j,k)

phiAbsorbedS_fine = phiAbsorbedS_fine + (phiOUT-phiIN)				! scalar absorbed at current location in mth direction
!===================================================================================================


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
