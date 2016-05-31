!==================================================================================================
MODULE LBM_fine				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs, ParticleTracking)
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
  
  !===================================================================================================
  SUBROUTINE Particle_Setup_fine
    !===================================================================================================
    
    IMPLICIT NONE
    
    IF (restart) THEN
    ELSE
       CALL Interp_Parvel_fine
    ENDIF
    
    !===================================================================================================
  END SUBROUTINE Particle_Setup_Fine
  !===================================================================================================
  
  !------------------------------------------------
  SUBROUTINE Particle_Track_fine
    !------------------------------------------------
    IMPLICIT NONE
    INTEGER(lng)   		 :: i,ipartition,ii,jj,kk, CaseNo
    REAL(dbl)      		 :: xpold(1:np),ypold(1:np),zpold(1:np) 	! old particle coordinates (working coordinates are stored in xp,yp,zp)
    REAL(dbl)      		 :: upold(1:np),vpold(1:np),wpold(1:np) 	! old particle velocity components (new vales are stored in up, vp, wp)
    REAL(dbl)                :: Cb_Domain, Cb_Local, Cb_Hybrid, V_eff_Ratio
    TYPE(ParRecord), POINTER :: current
    TYPE(ParRecord), POINTER :: next
    
    ParticleTransfer_fine = .FALSE. 						! AT this time we do not know if any particles need to be transferred.
    delphi_particle = 0.0_dbl 						! set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
    delphi_particle_fine = 0.0_dbl 						! set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
    
    tausgs_particle_x_fine = 0.0_dbl
    tausgs_particle_y_fine = 0.0_dbl
    tausgs_particle_z_fine = 0.0_dbl
    
    !--Second order interpolation in time
    !--Backup particle data from previous time step using a linked list of particle records
    
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next 						! copy pointer of next node
       
       IF ( flagParticleCF(current%pardata%parid) )  THEN  !Check if particle is in fine mesh
          IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++
             current%pardata%xpold = current%pardata%xp
             current%pardata%ypold = current%pardata%yp
             current%pardata%zpold = current%pardata%zp
             
             current%pardata%upold = current%pardata%up
             current%pardata%vpold = current%pardata%vp
             current%pardata%wpold = current%pardata%wp
             
             current%pardata%xp=current%pardata%xpold+current%pardata%up * tcf_fine
             current%pardata%yp=current%pardata%ypold+current%pardata%vp * tcf_fine
             current%pardata%zp=current%pardata%zpold+current%pardata%wp * tcf_fine
          END IF
          
       END IF
       
       current => next
    ENDDO

!    write(31,*) 'Calling Interp_Parvel_fine'
    CALL Interp_Parvel_fine
    
    !--Using a linked list of particle records
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next 						! copy pointer of next node
       
       IF ( flagParticleCF(current%pardata%parid) ) THEN  !Check if particle is in fine mesh
          
          IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++         
             current%pardata%xp=current%pardata%xpold+0.5*(current%pardata%up+current%pardata%upold) * tcf_fine
             current%pardata%yp=current%pardata%ypold+0.5*(current%pardata%vp+current%pardata%vpold) * tcf_fine
             current%pardata%zp=current%pardata%zpold+0.5*(current%pardata%wp+current%pardata%wpold) * tcf_fine
             
          END IF
       END IF
       current => next
    ENDDO

!    write(31,*) 'Calling Interp_Parvel_fine for the second time'    
    CALL Interp_Parvel_fine						! interpolate final particle velocities after the final position is ascertained. 
    CALL Particle_Transfer_fine
    
    !   CALL Interp_bulkconc(Cb_Local)  					! interpolate final bulk_concentration after the final position is ascertained.
    !   CALL Calc_Global_Bulk_Scalar_Conc(Cb_Domain)
!    write(31,*) 'Calling Compute_Cb_fine'    
    CALL Compute_Cb_fine(V_eff_Ratio,CaseNo,Cb_Hybrid)  
    
    open(172,file='Cb-history.dat',position='append')
    write(172,*) iter, V_eff_Ratio, CaseNo, Cb_Local, Cb_Domain, Cb_Hybrid
    
    !   CALL Update_Sh 							! Update the Sherwood number for each particle depending on the shear rate at the particle location. 
    CALL Calc_Scalar_Release_fine					! Updates particle radius, calculates new drug conc release rate delNB. 
    CALL Interp_ParToNodes_Conc_Fine 					! distributes released drug concentration to neightbouring nodes 
    !drug molecules released by the particle at this new position
    
    
    
    !---- Now update tausgs only for those cells that have non-zero values of tausgs
    DO kk=0,nzSub_fine+1
       DO jj=0,nySub_fine+1
          DO ii=0,nxSub_fine+1
             if (tausgs_particle_x_fine(ii,jj,kk).ne.0.0_dbl) then
                tausgs_particle_x_fine(ii,jj,kk) = u_fine(ii,jj,kk)*phi_fine(ii,jj,kk)
             endif
             if (tausgs_particle_y_fine(ii,jj,kk).ne.0.0_dbl) then
                tausgs_particle_y_fine(ii,jj,kk) = v_fine(ii,jj,kk)*phi_fine(ii,jj,kk)
             endif
             if (tausgs_particle_z_fine(ii,jj,kk).ne.0.0_dbl) then
                tausgs_particle_z_fine(ii,jj,kk) = w_fine(ii,jj,kk)*phi_fine(ii,jj,kk)
             endif
          ENDDO
       ENDDO
    ENDDO
    
    
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next ! copy pointer of next node
       
       IF ( flagParticleCF(current%pardata%parid) ) THEN  !Check if particle is in coarse mesh
          
          if (current%pardata%cur_part .eq. mySub) then
          SELECT CASE(current%pardata%parid)             
          CASE(1_lng)
             open(72,file='particle1-history.dat',position='append')
             write(72,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up*vcf,current%pardata%vp*vcf,current%pardata%wp*vcf,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(72)
          CASE(3_lng)
             open(73,file='particle3-history.dat',position='append')
             write(73,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(73) 
          CASE(5_lng)
             open(74,file='particle5-history.dat',position='append')
             write(74,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(74)
          CASE(7_lng)
             open(75,file='particle7-history.dat',position='append')
             write(75,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(75)
          CASE(9_lng)
             open(76,file='particle9-history.dat',position='append')
             write(76,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(76) 
          CASE(10_lng)
             open(77,file='particle10-history.dat',position='append')
             write(77,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(77)
          CASE(8_lng)
             open(78,file='particle8-history.dat',position='append')
             write(78,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(78)
          CASE(6_lng)
             open(79,file='particle6-history.dat',position='append')
             write(79,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(79)
          CASE(4_lng)
             open(80,file='particle4-history.dat',position='append')
             write(80,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(80)
          CASE(2_lng)
             open(81,file='particle2-history.dat',position='append')
             write(81,*) iter,(iter-1)*tcf+subIter*tcf_fine,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNB,current%pardata%cur_part,current%pardata%new_part
             close(81)
             
          END SELECT
       end if
          
       END IF
       
       current => next   				! point to next node in the list
    ENDDO
    
    !------------------------------------------------
  END SUBROUTINE Particle_Track_fine
  !------------------------------------------------
  
  
  !-----------------------------------------------!===================================================================================================
  SUBROUTINE Particle_Transfer_fine
    !===================================================================================================
    IMPLICIT NONE
    
    INTEGER(lng)   		 :: i,ipartition,ii,jj,kk
    INTEGER             :: RANK
    INTEGER(lng)             :: mpierr
    TYPE(ParRecord), POINTER :: current
    TYPE(ParRecord), POINTER :: next
    
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next ! copy pointer of next node
       IF (mySub .EQ.current%pardata%cur_part) THEN
          
          IF ( flagParticleCF(current%pardata%parid) ) THEN  !Check if particle is in coarse mesh
             
             !------- Wrappign around in z-direction for periodic BC in z
             IF (current%pardata%zp.GE.REAL(L-zcf_fine,dbl)) THEN
                current%pardata%zp = current%pardata%zp - L !MOD(current%pardata%zp,REAL(L,dbl))
             ELSE IF (current%pardata%zp .LE. (-1.0 * zcf_fine) ) THEN
                current%pardata%zp = current%pardata%zp+REAL(L,dbl)
             ENDIF
             
             !------- Estimate to which partition the updated position belongs to.
             DO ipartition = 1_lng,NumSubsTotal
                
                IF (( ((current%pardata%xp - xx_fine(1))/xcf_fine + 1) .GE. REAL(iMinDomain_fine(ipartition),dbl)-1.0_dbl).AND.&
                     ( ((current%pardata%xp - xx_fine(1))/xcf_fine + 1) .LT. (REAL(iMaxDomain_fine(ipartition),dbl)+0.0_dbl)).AND. &
                     ( ((current%pardata%yp - yy_fine(1))/ycf_fine + 1) .GE. REAL(jMinDomain_fine(ipartition),dbl)-1.0_dbl).AND. &
                     ( ((current%pardata%yp - yy_fine(1))/ycf_fine + 1) .LT. (REAL(jMaxDomain_fine(ipartition),dbl)+0.0_dbl)).AND. &
                     ( ((current%pardata%zp - zz_fine(1))/zcf_fine + 1) .GE. REAL(kMinDomain_fine(ipartition),dbl)-1.0_dbl).AND. &
                     ( ((current%pardata%zp - zz_fine(1))/zcf_fine + 1) .LT. (REAL(kMaxDomain_fine(ipartition),dbl)+0.0_dbl))) THEN
                   
                   current%pardata%new_part = ipartition
                END IF
                
             END DO
             
          END IF
       END IF
       current => next
    ENDDO
    
    !---- Parallel communication between all processors
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next 
       RANK= current%pardata%cur_part - 1
       current%pardata%cur_part = current%pardata%new_part 
       CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%xp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%yp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%zp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%up,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%vp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%wp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%rp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%sh,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%xpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%ypold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%zpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%upold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%vpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%wpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%rpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%delNB,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%par_conc,  1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%bulk_conc, 1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%gamma_cont,1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%Nbj,       1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%S,         1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%Sst,       1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%cur_part,  1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       CALL MPI_BCast(current%pardata%new_part,  1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
       current => next  
    ENDDO
    !===================================================================================================
  END SUBROUTINE Particle_Transfer_fine
  !===================================================================================================
  
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
  
  !===================================================================================================
  SUBROUTINE Interp_Parvel_fine ! Using Trilinear interpolation
    !===================================================================================================
    
    IMPLICIT NONE
    INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
    REAL(dbl)     :: xp,yp,zp,c00,c01,c10,c11,c0,c1,c,xd,yd,zd
    TYPE(ParRecord), POINTER :: current
    TYPE(ParRecord), POINTER :: next
    
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next ! copy pointer of next node
       
       IF ( flagParticleCF(current%pardata%parid) ) THEN  !Check if particle is in fine mesh
          IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++
             
             xp = (current%pardata%xp-xx_fine(1))/xcf_fine + 1 - REAL(iMin_fine-1_lng,dbl)
             yp = (current%pardata%yp-yy_fine(1))/ycf_fine + 1 - REAL(jMin_fine-1_lng,dbl)
             zp = (current%pardata%zp-zz_fine(1))/zcf_fine + 1 - REAL(kMin_fine-1_lng,dbl)
             
             ix0=FLOOR(xp)
             ix1=CEILING(xp)
             iy0=FLOOR(yp)
             iy1=CEILING(yp)
             iz0=FLOOR(zp)
             iz1=CEILING(zp)
!!!!!! MAKE SURE THE ABOVE NODES ARE FLUID NODES
             
             IF (ix1 /= ix0) THEN 
                xd=(xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
             ELSE
                xd = 0.0_dbl
             END IF
             
             IF (iy1 /= iy0) THEN 
                yd=(yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
             ELSE
                yd = 0.0_dbl
             END IF
             
             
             IF (iz1 /= iz0) THEN 
                zd=(zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
             ELSE
                zd = 0.0_dbl
             END IF
             
             ! u-interpolation
             ! Do first level linear interpolation in x-direction
             c00 = u_fine(ix0,iy0,iz0)*(1.0_dbl-xd)+u_fine(ix1,iy0,iz0)*xd	
             c01 = u_fine(ix0,iy0,iz1)*(1.0_dbl-xd)+u_fine(ix1,iy0,iz1)*xd	
             c10 = u_fine(ix0,iy1,iz0)*(1.0_dbl-xd)+u_fine(ix1,iy1,iz0)*xd	
             c11 = u_fine(ix0,iy1,iz1)*(1.0_dbl-xd)+u_fine(ix1,iy1,iz1)*xd	
             
             ! Do second level linear interpolation in y-direction
             c0  = c00*(1.0_dbl-yd)+c10*yd
             c1  = c01*(1.0_dbl-yd)+c11*yd
             
             ! Do third level linear interpolation in z-direction
             c   = c0*(1.0_dbl-zd)+c1*zd
             current%pardata%up=c * vcf_fine
             
             
             ! v-interpolation
             ! Do first level linear interpolation in x-direction
             c00 = v_fine(ix0,iy0,iz0)*(1.0_dbl-xd)+v_fine(ix1,iy0,iz0)*xd
             c01 = v_fine(ix0,iy0,iz1)*(1.0_dbl-xd)+v_fine(ix1,iy0,iz1)*xd
             c10 = v_fine(ix0,iy1,iz0)*(1.0_dbl-xd)+v_fine(ix1,iy1,iz0)*xd
             c11 = v_fine(ix0,iy1,iz1)*(1.0_dbl-xd)+v_fine(ix1,iy1,iz1)*xd	
             
             ! Do second level linear interpolation in y-direction
             c0  = c00*(1.0_dbl-yd)+c10*yd
             c1  = c01*(1.0_dbl-yd)+c11*yd
             
             ! Do third level linear interpolation in z-direction
             c   = c0*(1.0_dbl-zd)+c1*zd
             current%pardata%vp=c * vcf_fine
             
             ! w-interpolation
             ! Do first level linear interpolation in x-direction

             if( (w_fine(ix0,iy0,iz0) .lt. 1e-18) .or. (w_fine(ix1,iy0,iz0) .lt. 1e-18) .or. (w_fine(ix0,iy1,iz0) .lt. 1e-18) .or. (w_fine(ix0,iy1,iz0) .lt. 1e-18) .or. (w_fine(ix0,iy0,iz1) .lt. 1e-18) .or. (w_fine(ix1,iy0,iz1) .lt. 1e-18) .or. (w_fine(ix0,iy1,iz1) .lt. 1e-18) .or. (w_fine(ix0,iy1,iz1) .lt. 1e-18) ) then
                
                write(31,*) 'parid = ', current%pardata%parid
                write(31,*) 'ix0,iy0,iz0 = ', ix0, iy0, iz0
                write(31,*) 'w_fine(ix0,iy0,iz0) = ',w_fine(ix0,iy0,iz0)
                write(31,*) 'w_fine(ix1,iy0,iz0) = ',w_fine(ix1,iy0,iz0)
                write(31,*) 'w_fine(ix0,iy1,iz0) = ',w_fine(ix0,iy1,iz0)
                write(31,*) 'w_fine(ix1,iy1,iz0) = ',w_fine(ix1,iy1,iz0)            
                write(31,*) 'w_fine(ix0,iy0,iz1) = ',w_fine(ix0,iy0,iz1)
                write(31,*) 'w_fine(ix1,iy0,iz1) = ',w_fine(ix1,iy0,iz1)
                write(31,*) 'w_fine(ix0,iy1,iz1) = ',w_fine(ix0,iy1,iz1)
                write(31,*) 'w_fine(ix1,iy1,iz1) = ',w_fine(ix1,iy1,iz1)
             end if
             
             c00 = w_fine(ix0,iy0,iz0)*(1.0_dbl-xd)+w_fine(ix1,iy0,iz0)*xd	
             c01 = w_fine(ix0,iy0,iz1)*(1.0_dbl-xd)+w_fine(ix1,iy0,iz1)*xd	
             c10 = w_fine(ix0,iy1,iz0)*(1.0_dbl-xd)+w_fine(ix1,iy1,iz0)*xd	
             c11 = w_fine(ix0,iy1,iz1)*(1.0_dbl-xd)+w_fine(ix1,iy1,iz1)*xd	
             
             ! Do second level linear interpolation in y-direction
             c0  = c00*(1.0_dbl-yd)+c10*yd
             c1  = c01*(1.0_dbl-yd)+c11*yd
             
             ! Do third level linear interpolation in z-direction
             c   = c0*(1.0_dbl-zd)+c1*zd
             current%pardata%wp=c * vcf_fine
             
          END IF
       END IF
       ! point to next node in the list
       current => next
       
    ENDDO
    
    !===================================================================================================
  END SUBROUTINE Interp_Parvel_fine ! Using Trilinear interpolation
  !===================================================================================================
  
  !===================================================================================================
  SUBROUTINE Compute_Cb_fine(V_eff_Ratio,CaseNo,Cb_Hybrid) ! Computes the mesh-independent bulk concentration
    !===================================================================================================
    
    IMPLICIT NONE
    
    INTEGER(lng)  		:: i,j,k, kk, CaseNo
    INTEGER(lng)  		:: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11			! Trilinear interpolation parameters
    REAL(dbl)			:: fluids_Veff, fluids_Veff_l
    INTEGER,DIMENSION(2)   	:: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
    INTEGER,DIMENSION(2)   	:: NEP_x, NEP_y, NEP_z                  	! Lattice Nodes Surronding the particle
    
    REAL(dbl)     		:: c00,c01,c10,c11,c0,c1,c,xd,yd,zd		! Trilinear interpolation parameters
    REAL(dbl)  	   		:: xp,yp,zp! Local particle or node indices
    REAL(dbl)                 :: zp_fine, zp_coarse ! Local particle location after adjusting for periodicity
    
    REAL(dbl)			:: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
    REAL(dbl)       	        :: N_b         					! Modeling parameter to extend the volume of influence  
    REAL(dbl)    	        :: R_P, Sh_P, delta_P
    REAl(dbl)               	:: R_influence_p, L_influence_p			! Parameters related to particle's volume of influence
    REAl(dbl)               	:: V_influence_P, V_eff_Ratio			! Parameters related to particle's volume of influence
    REAL(dbl)			:: Cb_Total_Veff, Cb_Total_Veff_l
    REAL(dbl),DIMENSION(2)   	:: VIB_x, VIB_y, VIB_z	 			! Volume of Influence's Borders
    REAL(dbl),DIMENSION(2)      :: VIB_zc 			! Volume of Influence's Borders on coarse mesh
    REAL(dbl),DIMENSION(2)      :: VIB_zf 			! Volume of Influence's Borders on fine mesh
    REAL(dbl),DIMENSION(2)    	:: NVB_x, NVB_y, NVB_z				! Node Volume's Borders
    REAL(dbl)                   :: Delta_X, Delta_Y, Delta_Z
    REAL(dbl)                   :: x_DP, y_DP, z_DP				! Coordinates of "Discretized Point" (DP)
    REAL(dbl)                   :: Cb_Hybrid
    REAL(dbl)                   :: checkEffVolumeOverlapFineMesh                ! Check for how much of the effective volume of the particle overlaps with the fine mesh
    REAL(dbl)                   :: overlapCoarseProc, overlapFineProc
    LOGICAL                     :: hardCheckCoarseMesh                          ! Flag to ffigure out whether a given point is in the coarse or fine mesh
    
    TYPE(ParRecord), POINTER  	:: current
    TYPE(ParRecord), POINTER  	:: next
    INTEGER(lng)  		:: mpierr
    
    delta_mesh = 1.0_dbl
    
    zcf3 = xcf_fine*ycf_fine*zcf_fine
    
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       
       !------ Copy pointer of next node
       next => current%next
       
       IF ( flagParticleCF(current%pardata%parid) ) THEN  !Check if particle is in coarse mesh
          !------ Calculate length scale for jth particle:  delta = R / Sh
          !------ Calculate effective radius: R_influence_P = R + (N_b *delta)
          !------ Note: need to convert this into Lattice units and not use the physical length units
          !------ Then compute equivalent cubic mesh length scale
          N_b = 1.0
          R_P  = current%pardata%rp
          Sh_P = current%pardata%sh
          delta_P = R_P / Sh_P
          R_influence_P = (R_P + N_b * delta_P) 
          
          !------ Computing equivalent cubic mesh length scale
          L_influence_P = ( (4.0_dbl*PI/3.0_dbl)**(1.0_dbl/3.0_dbl) ) * R_influence_P
          V_influence_P= (4.0_dbl/3.0_dbl) * PI * R_influence_P**3.0_dbl
          V_eff_Ratio = (L_influence_P / zcf_fine)**3 					! Ratio of the effective volume to cell size
!          write(31,*) 'L_influence_P = ', L_influence_P, ' in fine mesh units = ', L_influence_P/xcf_fine
          
!          write(31,*) 'V_eff_Ratio = ', V_eff_Ratio
!          flush(31)
          
          !------ Finding particle location in this processor
          xp = (current%pardata%xp-xx_fine(1))/xcf_fine + 1 - REAL(iMin_fine-1_lng,dbl)
          yp = (current%pardata%yp-yy_fine(1))/ycf_fine + 1 - REAL(jMin_fine-1_lng,dbl)
          zp = (current%pardata%zp-zz_fine(1))/zcf_fine + 1 - REAL(kMin_fine-1_lng,dbl)
          
          !----------------------------------------------------------------------------------------------------------------------
          !------ Veff is smaller than the mesh volume --> Cb = Trilinear interpolation of the concentration at particle location
          !----------------------------------------------------------------------------------------------------------------------
          IF (V_eff_Ratio .LE. 1.0) THEN 					
             CaseNo = 1
!             write(31,*) 'CaseNo = ', CaseNo
             IF (mySub .EQ.current%pardata%cur_part) THEN
                ix0 =FLOOR(xp)
                ix1 =CEILING(xp)
                iy0 =FLOOR(yp)
                iy1 =CEILING(yp)
                iz0 =FLOOR(zp)
                iz1 =CEILING(zp)
                
                !----------TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
                
                IF (ix1 /= ix0) THEN
                   xd=(xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
                ELSE
                   xd = 0.0_dbl
                END IF
                
                IF (iy1 /= iy0) THEN
                   yd=(yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
                ELSE
                   yd = 0.0_dbl
                END IF
                
                IF (iz1 /= iz0) THEN
                   zd=(zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                ELSE
                   zd = 0.0_dbl
                END IF
                
                !--------- Concentration Trilinear Iinterpolation
                !--------- Interpolation in x-direction
!                write(31,*) 'ix0,iy0,iz0 = ', ix0, iy0, iz0
!                write(31,*) 'ix1,ix1,iz1 = ', ix1, iy1, iz1
!                write(31,*) phi_fine(ix0,iy0,iz0), phi_fine(ix1,iy0,iz0)
!                write(31,*) phi_fine(ix0,iy0,iz1), phi_fine(ix1,iy0,iz1)
!                write(31,*) phi_fine(ix0,iy1,iz0), phi_fine(ix1,iy1,iz0)
!                write(31,*) phi_fine(ix0,iy1,iz1), phi_fine(ix1,iy1,iz1)
                flush(31)
                c00 = phi_fine(ix0,iy0,iz0) * (1.0_dbl-xd) + phi_fine(ix1,iy0,iz0) * xd
                c01 = phi_fine(ix0,iy0,iz1) * (1.0_dbl-xd) + phi_fine(ix1,iy0,iz1) * xd
                c10 = phi_fine(ix0,iy1,iz0) * (1.0_dbl-xd) + phi_fine(ix1,iy1,iz0) * xd
                c11 = phi_fine(ix0,iy1,iz1) * (1.0_dbl-xd) + phi_fine(ix1,iy1,iz1) * xd
                !--------- Interpolation in y-direction
                c0  = c00 * (1.0_dbl-yd) + c10 * yd
                c1  = c01 * (1.0_dbl-yd) + c11 * yd
                !--------- Interpolation in z-direction
                c   = c0 * (1.0_dbl-zd) + c1 * zd
                
                Cb_Hybrid= c 
!                write(31,*) 'Cb_Hybrid = ', Cb_Hybrid
                flush(31)
                
             END IF
             !----------------------------------------------------------------------------------------------------------------------
             !------ Veff is slightly larger than mesh volume --> Volume of influence is discretized 
             !------ Cb= Average of concentration interpolated on each of the descritized nodes inside volume of influence 
             !----------------------------------------------------------------------------------------------------------------------
          ELSE IF ( (V_eff_Ratio .GT. 1.0) .AND. (V_eff_Ratio .LT. 27.0 ) ) THEN		
             CaseNo = 2
!             write(31,*) 'CaseNo = ', CaseNo
             
             !Volume of Influence Border (VIB) for this particle
             VIB_x(1)= current%pardata%xp - 0.5_dbl * L_influence_P
             VIB_x(2)= current%pardata%xp + 0.5_dbl * L_influence_P
             VIB_y(1)= current%pardata%yp - 0.5_dbl * L_influence_P
             VIB_y(2)= current%pardata%yp + 0.5_dbl * L_influence_P
             VIB_z(1)= current%pardata%zp - 0.5_dbl * L_influence_P
             VIB_z(2)= current%pardata%zp + 0.5_dbl * L_influence_P
             
             !Check if Volume of Influence Border overlaps with the current processor domain
             
             if (MAX ( MIN(VIB_z(2)+L,kMaxDomain_fine(mySub)*zcf_fine) - MAX(VIB_z(1)+L,(kMinDomain_fine(mySub)-1)*zcf_fine), 0.0_dbl) .gt. 0) then
                
                VIB_z(1) = VIB_z(1)+L
                VIB_z(2) = VIB_z(2)+L
                
             else  if (MAX ( MIN(VIB_z(2)-L,kMaxDomain_fine(mySub)*zcf_fine) - MAX(VIB_z(1)-L,(kMinDomain_fine(mySub)-1)*zcf_fine), 0.0_dbl) .gt. 0) then
                
                VIB_z(1) = VIB_z(1)-L
                VIB_z(2) = VIB_z(2)-L
                
             end if
             
             overlapFineProc = MAX ( MIN(VIB_z(2),kMaxDomain_fine(mySub)*zcf_fine) - MAX(VIB_z(1),(kMinDomain_fine(mySub)-1)*zcf_fine) , 0.0_dbl)
             
             Cb_Total_Veff = 0
             fluids_Veff = 0
             Cb_Total_Veff_l = 0
             fluids_Veff_l = 0
             
             if (overlapFineProc .gt. 0) then
                
                checkEffVolumeOverlapFineMesh = ( MAX ( MIN(VIB_x(2), fractionDfine * D * 0.5 + xcf) - MAX(VIB_x(1), -fractionDfine * D * 0.5 - xcf), 0.0_dbl) / L_influence_P ) * &
                     ( MAX ( MIN(VIB_y(2),fractionDfine * D * 0.5 + xcf) - MAX(VIB_y(1), -fractionDfine * D * 0.5 - ycf), 0.0_dbl) / L_influence_P ) - 0.99
                
                !--------- Discretizing the volume of influence to  make sure at least 64 points are available
                Delta_X = (VIB_x(2)-VIB_x(1)) / 3.0 
                Delta_Y = (VIB_y(2)-VIB_y(1)) / 3.0
                Delta_Z = (VIB_z(2)-VIB_z(1)) / 3.0 
                
                
                
                !--------- Loop over discretized points and averaging the concentration
                DO i= 0, 3
                   DO j= 0, 3
                      DO k= 0, 3
                         xp = VIB_x(1) + (i * Delta_X) 
                         yp = VIB_y(1) + (j * Delta_Y)
                         zp = VIB_z(1) + (k * Delta_Z)
                         
                         
                         hardCheckCoarseMesh = ( (xp - fractionDfine * D * 0.5 - xcf) * (xp + fractionDfine * D * 0.5 + xcf) > 0 ) .or. ( (yp - fractionDfine * D * 0.5 - ycf) * (yp + fractionDfine * D * 0.5 + ycf) > 0 )
                         
                         if (hardCheckCoarseMesh) then
!                            write(31,*) 'Point in coarse mesh', xp, yp, zp
                            
                            if( (zp - kMaxDomain(mySub)*zcf) * (zp - (kMinDomain(mySub)-1)*zcf) .lt. 0)  then !If this point is in this processor
                               
                               
                               !------------------ Finding Lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                               x_DP = (xp - xx(1))/xcf + 1 - REAL(iMin-1_lng,dbl)
                               y_DP = (yp - yy(1))/ycf + 1 - REAL(jMin-1_lng,dbl)
                               z_DP = (zp - zz(1))/zcf + 1 - REAL(kMin-1_lng,dbl)
                               
                               ix0 = FLOOR(x_DP)
                               ix1 = CEILING(x_DP)
                               iy0 = FLOOR(y_DP)
                               iy1 = CEILING(y_DP)
                               iz0 = FLOOR(z_DP)
                               iz1 = CEILING(z_DP)
                               
                               !------------------ TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
                               IF (ix1 /= ix0) THEN
                                  xd=(x_DP-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
                               ELSE
                                  xd = 0.0_dbl
                               END IF
                               
                               IF (iy1 /= iy0) THEN
                                  yd=(y_DP-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
                               ELSE
                                  yd = 0.0_dbl
                               END IF
                               
                               IF (iz1 /= iz0) THEN
                                  zd=(z_DP-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                               ELSE
                                  zd = 0.0_dbl
                               END IF
                               
                               !------------------ Taking care of the periodic BC in Z-dir
                               iz00 = iz0
                               IF (iz0 .gt. nz) THEN
                                  iz00 = iz0 - (nz - 1)
                               ELSE IF (iz0 .lt. 1) THEN
                                  iz00 = iz0 + (nz-1)
                               END IF
                               
                               iz11 = iz1         
                               IF (iz1 .gt. nz) THEN
                                  iz11 = iz1 - (nz-1)
                               ELSE IF (iz1 .lt. 1) THEN
                                  iz11 = iz1 + (nz-1)
                               END IF
                               
                               !------------------ Concentration Trilinear Iinterpolation
                               !------------------ Interpolation in x-direction
                               c00 = phi(ix0,iy0,iz00) * (1.0_dbl-xd) + phi(ix1,iy0,iz00) * xd
                               c01 = phi(ix0,iy0,iz11) * (1.0_dbl-xd) + phi(ix1,iy0,iz11) * xd
                               c10 = phi(ix0,iy1,iz00) * (1.0_dbl-xd) + phi(ix1,iy1,iz00) * xd
                               c11 = phi(ix0,iy1,iz11) * (1.0_dbl-xd) + phi(ix1,iy1,iz11) * xd
                               !------------------ Interpolation in y-direction
                               c0  = c00 * (1.0_dbl-yd) + c10 * yd
                               c1  = c01 * (1.0_dbl-yd) + c11 * yd
                               !------------------ Interpolation in z-direction
                               c   = c0 * (1.0_dbl-zd) + c1 * zd
                               
                               Cb_Total_Veff_l  = Cb_Total_Veff_l  + c
                               fluids_Veff_l = fluids_Veff_l + 1.0_dbl
                               
                            end if
                         else
                            
!                            write(31,*) 'Point in fine mesh', xp, yp, zp
                            if( (zp - kMaxDomain_fine(mySub)*zcf_fine) * (zp - (kMinDomain_fine(mySub)-1)*zcf_fine) .lt. 0)  then !If this point is in this processor                    
                               
                               
                               !------------------ Finding Lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                               x_DP = (xp - xx_fine(1))/xcf_fine + 1 - REAL(iMin_fine-1_lng,dbl)
                               y_DP = (yp - yy_fine(1))/ycf_fine + 1 - REAL(jMin_fine-1_lng,dbl)
                               z_DP = (zp - zz_fine(1))/zcf_fine + 1 - REAL(kMin_fine-1_lng,dbl)
                               
                               !------------------ Finding Lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                               ix0 = MAX(2,FLOOR(x_DP))
                               ix1 = MIN(nx_fine-1,CEILING(x_DP))
                               iy0 = MAX(2,FLOOR(y_DP))
                               iy1 = MIN(ny_fine-1,CEILING(y_DP))
                               iz0 = MAX(1,FLOOR(z_DP))
                               iz1 = MIN(nzSub_fine,CEILING(z_DP))
                               
                               !------------------ TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
                               IF (ix1 /= ix0) THEN
                                  xd=(x_DP-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
                               ELSE
                                  xd = 0.0_dbl
                               END IF
                               
                               IF (iy1 /= iy0) THEN
                                  yd=(y_DP-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
                               ELSE
                                  yd = 0.0_dbl
                               END IF
                               
                               IF (iz1 /= iz0) THEN
                                  zd=(z_DP-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                               ELSE
                                  zd = 0.0_dbl
                               END IF
                               
                               !------------------ Taking care of the periodic BC in Z-dir
                               iz00 = iz0
                               IF (iz0 .gt. nz_fine) THEN
                                  iz00 = iz0 - (nz_fine - 1)
                               ELSE IF (iz0 .lt. 1) THEN
                                  iz00 = iz0 + (nz_fine-1)
                               END IF
                               
                               iz11 = iz1         
                               IF (iz1 .gt. nz_fine) THEN
                                  iz11 = iz1 - (nz_fine-1)
                               ELSE IF (iz1 .lt. 1) THEN
                                  iz11 = iz1 + (nz_fine-1)
                               END IF
                               
                               !------------------ Concentration Trilinear Iinterpolation
                               !------------------ Interpolation in x-direction
                               c00 = phi_fine(ix0,iy0,iz00) * (1.0_dbl-xd) + phi_fine(ix1,iy0,iz00) * xd
                               c01 = phi_fine(ix0,iy0,iz11) * (1.0_dbl-xd) + phi_fine(ix1,iy0,iz11) * xd
                               c10 = phi_fine(ix0,iy1,iz00) * (1.0_dbl-xd) + phi_fine(ix1,iy1,iz00) * xd
                               c11 = phi_fine(ix0,iy1,iz11) * (1.0_dbl-xd) + phi_fine(ix1,iy1,iz11) * xd
                               !------------------ Interpolation in y-direction
                               c0  = c00 * (1.0_dbl-yd) + c10 * yd
                               c1  = c01 * (1.0_dbl-yd) + c11 * yd
                               !------------------ Interpolation in z-direction
                               c   = c0 * (1.0_dbl-zd) + c1 * zd
                               
                               Cb_Total_Veff_l  = Cb_Total_Veff_l  + c
                               fluids_Veff_l = fluids_Veff_l + 1.0_dbl
                               
                            end if
                         end if
                         
                      END DO
                   END DO
                END DO
                
             end if
             
             !----Communication with other processors for V_eff greater than 1
             CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
             CALL MPI_ALLREDUCE(Cb_Total_Veff_l , Cb_Total_Veff , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
             CALL MPI_ALLREDUCE(fluids_Veff_l, fluids_Veff, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
             Cb_Hybrid= Cb_Total_Veff / fluids_Veff
             
             !----------------------------------------------------------------------------------------------------------------------
             !------ Veff is much larger than mesh volume --> Cb= total number of moles in volume of influence / volume of influence 
             !----------------------------------------------------------------------------------------------------------------------
          ELSE IF (V_eff_Ratio .GE. 27.0) THEN                             
             CaseNo = 3
!             write(31,*) 'CaseNo = ', CaseNo
             flush(31)
             !Volume of Influence Border (VIB) for this particle
             VIB_x(1)= current%pardata%xp - 0.5_dbl * L_influence_P 
             VIB_x(2)= current%pardata%xp + 0.5_dbl * L_influence_P
             VIB_y(1)= current%pardata%yp - 0.5_dbl * L_influence_P
             VIB_y(2)= current%pardata%yp + 0.5_dbl * L_influence_P
             VIB_z(1)= current%pardata%zp - 0.5_dbl * L_influence_P
             VIB_z(2)= current%pardata%zp + 0.5_dbl * L_influence_P
             
             checkEffVolumeOverlapFineMesh = ( MAX ( MIN(VIB_x(2), fractionDfine * D * 0.5 + xcf) - MAX(VIB_x(1), -fractionDfine * D * 0.5 - xcf), 0.0_dbl) / L_influence_P ) * &
                  ( MAX ( MIN(VIB_y(2),fractionDfine * D * 0.5 + xcf) - MAX(VIB_y(1), -fractionDfine * D * 0.5 - ycf), 0.0_dbl) / L_influence_P ) - 0.99
             
             !Check if Volume of Influence Border overlaps with the current processor domain
             
!             write(31,*) 'mySub = ', mySub
!             write(31,*) 'VIB_z(1) = ', VIB_z(1)
!             write(31,*) 'VIB_z(2) = ', VIB_z(2)
!             write(31,*) 'VIB_x(1) = ', VIB_x(1)
!             write(31,*) 'VIB_x(2) = ', VIB_x(2)
!             write(31,*) 'VIB_y(1) = ', VIB_y(1)
!             write(31,*) 'VIB_y(2) = ', VIB_y(2)
!             write(31,*) 'L = ', L
!             write(31,*) 'z end = ', kMaxDomain(mySub)*zcf
!             write(31,*) 'z begin = ', (kMinDomain(mySub)-1)*zcf
!             write(31,*) 'MIN(VIB_z(2)+L,kMaxDomain(mySub)*zcf) = ', MIN(VIB_z(2)+L,kMaxDomain(mySub)*zcf)
!             write(31,*) 'MAX(VIB_z(1)+L,(kMinDomain(mySub)-1)*zcf) = ', MAX(VIB_z(1)+L,(kMinDomain(mySub)-1)*zcf)
!             write(31,*) 'MIN(VIB_z(2)-L,kMaxDomain(mySub)*zcf) = ', MIN(VIB_z(2)-L,kMaxDomain(mySub)*zcf)
!             write(31,*) 'MAX(VIB_z(1)-L,(kMinDomain(mySub)-1)*zcf) = ', MAX(VIB_z(1)-L,(kMinDomain(mySub)-1)*zcf)
!             write(31,*) 'MIN(VIB_z(2),kMaxDomain(mySub)*zcf) = ', MIN(VIB_z(2),kMaxDomain(mySub)*zcf)
!             write(31,*) 'MAX(VIB_z(1),(kMinDomain(mySub)-1)*zcf) = ', MAX(VIB_z(1),(kMinDomain(mySub)-1)*zcf)


             if (MAX ( MIN(VIB_z(2)+L,(kMaxDomain(mySub)-0.5)*zcf) - MAX(VIB_z(1)+L,(kMinDomain(mySub)-1.5)*zcf) , 0.0_dbl) .gt. 0) then
                
                zp_coarse = current%pardata%zp + L           
                VIB_zc(1) = VIB_z(1)+L
                VIB_zc(2) = VIB_z(2)+L
                
             else  if (MAX ( MIN(VIB_z(2)-L,(kMaxDomain(mySub)-0.5)*zcf) - MAX(VIB_z(1)-L,(kMinDomain(mySub)-1.5)*zcf) , 0.0_dbl) .gt. 0) then
                
                zp_coarse = current%pardata%zp - L
                VIB_zc(1) = VIB_z(1)-L
                VIB_zc(2) = VIB_z(2)-L
                
             else

                zp_coarse = current%pardata%zp                
                VIB_zc(1) = VIB_z(1)
                VIB_zc(2) = VIB_z(2)
                
             end if
             
!             write(31,*) 'MIN(VIB_x(2), iMaxDomain(mySub)*xcf ) - MAX(VIB_x(1), (iMinDomain(mySub)-1)*xcf) = ', MIN(VIB_x(2), iMaxDomain(mySub)*xcf ) - MAX(VIB_x(1), (iMinDomain(mySub)-1)*xcf)
!             write(31,*) 'MIN(VIB_y(2),jMaxDomain(mySub)*ycf ) - MAX(VIB_y(1),(jMinDomain(mySub)-1)*ycf) = ', MIN(VIB_y(2),jMaxDomain(mySub)*ycf ) - MAX(VIB_y(1),(jMinDomain(mySub)-1)*ycf)
!             write(31,*) 'MIN(VIB_zc(2),kMaxDomain(mySub)*zcf) - MAX(VIB_zc(1),(kMinDomain(mySub)-1)*zcf) = ', MIN(VIB_zc(2),kMaxDomain(mySub)*zcf) - MAX(VIB_zc(1),(kMinDomain(mySub)-1)*zcf)

!             write(31,*) 'VIB_zc = ', VIB_zc(1), VIB_zc(2), ' zMinMax  = ', (kMaxDomain(mySub)-0.5)*zcf, (kMinDomain(mySub)-1.5)*zcf
             overlapCoarseProc = MAX ( MIN(VIB_zc(2),(kMaxDomain(mySub)-0.5)*zcf) - MAX(VIB_zc(1),(kMinDomain(mySub)-1.5)*zcf) , 0.0_dbl)
!             write(31,*) 'overlapCoarseProc = ', overlapCoarseProc, ' Ratio = ', overlapCoarseProc/(VIB_zc(2)-VIB_zc(1))
!             flush(31)
             
             if (MAX ( MIN(VIB_z(2)+L,(kMaxDomain_fine(mySub)-0.5)*zcf_fine) - MAX(VIB_z(1)+L,(kMinDomain_fine(mySub)-1.5)*zcf_fine), 0.0_dbl) .gt. 0) then
                
                zp_fine = current%pardata%zp + L                           
                VIB_zf(1) = VIB_z(1)+L
                VIB_zf(2) = VIB_z(2)+L
                
             else  if (MAX ( MIN(VIB_z(2)-L,(kMaxDomain_fine(mySub)-0.5)*zcf_fine) - MAX(VIB_z(1)-L,(kMinDomain_fine(mySub)-1.5)*zcf_fine), 0.0_dbl) .gt. 0) then
                
                zp_fine = current%pardata%zp - L
                VIB_zf(1) = VIB_z(1)-L
                VIB_zf(2) = VIB_z(2)-L
                
             else
                
                zp_fine = current%pardata%zp
                VIB_zf(1) = VIB_z(1)
                VIB_zf(2) = VIB_z(2)
                
             end if

!             write(31,*) 'VIB_zf = ', VIB_zf(1), VIB_zf(2)
!             write(31,*) 'zMinMax = ', (kMinDomain_fine(mySub)-1.5)*zcf_fine, (kMaxDomain_fine(mySub)-0.5)*zcf_fine
             
             overlapFineProc = MAX ( MIN(VIB_zf(2),(kMaxDomain_fine(mySub)-0.5)*zcf_fine) - MAX(VIB_zf(1),(kMinDomain_fine(mySub)-1.5)*zcf_fine) , 0.0_dbl)
!             write(31,*) 'overlapFineProc = ', overlapFineProc, 'Ratio = ', overlapFineProc/(VIB_zf(2)-VIB_zf(1))
!             flush(31)
             
             Cb_Total_Veff = 0.0_dbl
             fluids_Veff = 0.0_dbl
             Cb_Total_Veff_l = 0.0_dbl
             fluids_Veff_l = 0.0_dbl
             
             if (overlapFineProc .gt. 0) then
                
                !------ Finding particle location in this processor
                xp = (current%pardata%xp-xx_fine(1))/xcf_fine + 1 - REAL(iMin_fine-1_lng,dbl)
                yp = (current%pardata%yp-yy_fine(1))/ycf_fine + 1 - REAL(jMin_fine-1_lng,dbl)
                zp = (zp_fine-zz_fine(1))/zcf_fine + 1 - REAL(kMin_fine-1_lng,dbl)
                
                !Finding the lattice "Nodes Effected by Particle"
                NEP_x(1)= MAX(2,CEILING(xp - 0.5*L_influence_P/xcf_fine))
                NEP_x(2)= MIN(nxSub_fine-1, FLOOR(xp + 0.5*L_influence_P/xcf_fine))
                NEP_y(1)= MAX(2,CEILING(yp - 0.5*L_influence_P/ycf_fine))
                NEP_y(2)= MIN(nySub_fine-1, FLOOR(yp + 0.5*L_influence_P/ycf_fine))
                NEP_z(1)= MAX(1,CEILING(zp - 0.5*L_influence_P/zcf_fine))
                NEP_z(2)= MIN(nzSub_fine, FLOOR(zp + 0.5*L_influence_P/zcf_fine))
                
                ! write(31,*) 'Fine mesh'
                ! write(31,*) 'NEP_x = ', NEP_x(1), NEP_x(2)
                ! write(31,*) 'NEP_y = ', NEP_y(1), NEP_y(2)
                ! write(31,*) 'zp, NEP_z = ', zp, NEP_z(1), NEP_z(2)
                ! write(31,*) 'estimate of fluids_Veff_l = ', (NEP_x(2)-NEP_x(1)+1)*(NEP_y(2)-NEP_y(1)+1)*(NEP_z(2)-NEP_z(1)+1)
                ! flush(31)

                ! write(31,*) 'fluids_Veff_l before = ', fluids_Veff_l

                DO i= NEP_x(1),NEP_x(2) 
                   DO j= NEP_y(1),NEP_y(2)
                      DO k= NEP_z(1),NEP_z(2)
                         
                         !---- Taking care of the periodic BC in Z-dir
                         kk = k
                         
                         IF (node_fine(i,j,kk) .EQ. FLUID) THEN
                            Cb_Total_Veff_l  = Cb_Total_Veff_l  + phi_fine(i,j,kk)
                            fluids_Veff_l = fluids_Veff_l + 1.0_dbl
                         END IF
                      END DO
                   END DO
                END DO
                
             end if
!             write(31,*) 'fluids_Veff_l = ', fluids_Veff_l
!             write(31,*) 'Overlap fine mesh', fluids_Veff_l*xcf_fine*xcf_fine*xcf_fine/(L_influence_P*L_influence_P*L_influence_P)
!             flush(31)
             
             if (overlapCoarseProc .gt. 0) then
                if (checkEffVolumeOverlapFineMesh .lt. 0) then
                   
                   !------ Finding particle location in this processor
                   xp = (current%pardata%xp-xx(1))/xcf + 1 - REAL(iMin-1_lng,dbl)
                   yp = (current%pardata%yp-yy(1))/ycf + 1 - REAL(jMin-1_lng,dbl)
                   zp = (zp_fine-zz(1))/zcf + 1 - REAL(kMin-1_lng,dbl)
                   
                   NEP_x(1)= MAX(1,CEILING(xp - 0.5*L_influence_P/xcf))
                   NEP_x(2)= MIN(nxSub, FLOOR(xp + 0.5*L_influence_P/xcf))
                   NEP_y(1)= MAX(1,CEILING(yp - 0.5*L_influence_P/ycf))
                   NEP_y(2)= MIN(nySub, FLOOR(yp + 0.5*L_influence_P/ycf))
                   NEP_z(1)= MAX(1,CEILING(zp - 0.5*L_influence_P/zcf))
                   NEP_z(2)= MIN(nzSub, FLOOR(zp + 0.5*L_influence_P/zcf))
                   
                  ! write(31,*) 'Coarse mesh'
                  ! write(31,*) 'NEP_x = ', NEP_x(1), NEP_x(2)
                  ! write(31,*) 'NEP_y = ', NEP_y(1), NEP_y(2)
                  ! write(31,*) 'zp = ', zp, ' zp - 0.5*L_influence_P/zcf = ', zp - 0.5*L_influence_P/zcf, ' zp + 0.5*L_influence_P/zcf = ', zp + 0.5*L_influence_P/zcf
                  ! write(31,*) 'NEP_z = ', NEP_z(1), NEP_z(2)

                  !  flush(31)                   
                   
                   DO i= NEP_x(1),NEP_x(2) 
                      DO j= NEP_y(1),NEP_y(2)
                         DO k= NEP_z(1),NEP_z(2)
                            
                            !---- Taking care of the periodic BC in Z-dir
                            kk = k
                            
                            IF (node(i,j,kk) .EQ. FLUID) THEN
                               Cb_Total_Veff_l  = Cb_Total_Veff_l  + phi(i,j,kk) * (1.0 - flagNodeIntersectFine(i,j,kk)) * (gridRatio * gridRatio * gridRatio)
                               fluids_Veff_l = fluids_Veff_l + (1.0 - flagNodeIntersectFine(i,j,kk)) * (gridRatio * gridRatio * gridRatio)
                            END IF
                         END DO
                      END DO
                   END DO
                end if
             end if
!             write(31,*) 'Overlap fine+coarse mesh', fluids_Veff_l*xcf_fine*xcf_fine*xcf_fine/(L_influence_P*L_influence_P*L_influence_P)             
             
!             write(31,*) 'Bulk conc calculation coarse mesh ', Cb_Total_Veff_l, fluids_Veff_l
             !----Communication with other processors for V_eff greater than 1
             CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
             CALL MPI_ALLREDUCE(Cb_Total_Veff_l , Cb_Total_Veff , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
             CALL MPI_ALLREDUCE(fluids_Veff_l, fluids_Veff, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
             Cb_Hybrid= Cb_Total_Veff / fluids_Veff

!             write(31,*) 'Overlap total', fluids_Veff*xcf_fine*xcf_fine*xcf_fine/(L_influence_P*L_influence_P*L_influence_P)             
             
             current%pardata%bulk_conc = Cb_Hybrid
!             write(31,*) 'current%pardata%bulk_conc = ', current%pardata%bulk_conc
!             flush(31)
             
          END IF
       END IF
       current => next
       
    END DO
    
    !===================================================================================================
  END SUBROUTINE Compute_Cb_fine
  !===================================================================================================
  
  !===================================================================================================
  SUBROUTINE Calc_Scalar_Release_fine! Calculate rate of scalar release at every time step  
    !===================================================================================================
    
    ! Called by Particle_Track (LBM.f90) to get delNB, update particle radius,
    ! Sh(t)- sherwood number
    
    IMPLICIT NONE
    INTEGER(lng)  :: numFluids,i,j,k
    REAL(dbl)     :: deltaR,bulkVolume,temp,cbt,zcf3,bulkconc
    TYPE(ParRecord), POINTER :: current
    TYPE(ParRecord), POINTER :: next
    INTEGER(lng)  :: RANK,mpierr  
    
    
    !bulkVolume=xcf*ycf_fine*zcf_fine*Cb_numFluids/num_particles
    zcf3=xcf_fine*ycf_fine*zcf_fine
    
    !calculate delNB for each particle in the domain
    current => ParListHead%next
    DO WHILE (ASSOCIATED(current))
       next => current%next ! copy pointer of next node
       
       IF ( flagParticleCF(current%pardata%parid) )  THEN  !Check if particle is in fine mesh
          
          IF (mySub .EQ.current%pardata%cur_part) THEN 
             
              current%pardata%rpold=current%pardata%rp
             
              bulkconc = current%pardata%bulk_conc
             
              temp = current%pardata%rpold**2.0_dbl-4.0_dbl*tcf_fine*molarvol*diffm*current%pardata%sh*max((current%pardata%par_conc-bulkconc),0.0_dbl)
              ! write(31,*) '4.0_dbl*tcf*molarvol*diffm = ', 4.0_dbl*tcf*molarvol*diffm
              ! write(31,*) 'current%pardata%sh = ', current%pardata%sh
              ! write(31,*) 'max((current%pardata%par_conc-bulkconc),0.0_dbl) = ', max((current%pardata%par_conc-bulkconc),0.0_dbl)
              ! write(31,*) 'current%pardata%par_conc, bulkconc', current%pardata%par_conc, bulkconc
             
              IF (temp.GE.0.0_dbl) THEN
                 current%pardata%rp=0.5_dbl*(current%pardata%rpold+sqrt(temp))
              ELSE
                 temp = 0.0_dbl
                 current%pardata%rp=0.5_dbl*(current%pardata%rpold+sqrt(temp))
              END IF
             
             deltaR=current%pardata%rpold-current%pardata%rp
             
             current%pardata%delNB=(4.0_dbl/3.0_dbl)*PI*(current%pardata%rpold*current%pardata%rpold*current%pardata%rpold &
                  -current%pardata%rp*current%pardata%rp*current%pardata%rp) &
                  /molarvol
             
!             write(*,*) 'Not allowing particle radius to change in the fine mesh'
!             current%pardata%delNB= (4.0* PI* current%pardata%rp) * current%pardata%sh* diffm * max((current%pardata%par_conc-current%pardata%bulk_conc) , 0.0_dbl) * tcf_fine
             
!             write(31,*) 'current%pardata%rpold = ', current%pardata%rpold
!             write(31,*) 'current%pardata%rp = ', current%pardata%rp         
!             write(31,*) 'current%pardata%delNB = ', current%pardata%delNB
!             write(31,*) 'molarvol, bulkVolume ', molarvol, bulkVolume
!             flush(31)
             
             IF (associated(current,ParListHead%next)) THEN
                write(9,*) iter*tcf_fine,current%pardata%parid,current%pardata%rp,current%pardata%Sh,Cb_global*zcf3*Cb_numFluids,current%pardata%delNB,Cb_global,Cb_numFluids
                CALL FLUSH(9)
             ENDIF
             
             
          END IF
          
       END IF
       ! point to next node in the list
       
       CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
       RANK= current%pardata%cur_part - 1
       CALL MPI_BCast(current%pardata%delNB, 1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD, mpierr)
       CALL MPI_BCast(current%pardata%rp, 1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD, mpierr)     
       
       current => next
    ENDDO
    
    
    IF (myid .EQ. 0) THEN
       open(799,file='testoutput.dat',position='append')
       write(799,*) iter*tcf_fine,Cb_global*zcf3*Cb_numFluids,Cb_global,Cb_numFluids
       close(799)
    ENDIF
    
    !===================================================================================================
  END SUBROUTINE Calc_Scalar_Release_Fine
  !===================================================================================================
  
  !===================================================================================================
  SUBROUTINE Interp_ParToNodes_Conc_Fine
    !===================================================================================================
    
    !--- Interpolate Particle concentration release to node locations 
    !--- Called by Particle_Track (LBM.f90) to get delphi_particle
    
    IMPLICIT NONE
    INTEGER(lng) 	      :: i,j,k,kk
    REAL(dbl)     	      :: xp,yp,zp ! Local particle or node indices
    REAL(dbl)                 :: zp_fine, zp_coarse ! Local particle location after adjusting for periodicity
    REAL(dbl)		      :: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
    REAL(dbl)                 :: N_d         				! Modeling parameter to extend the volume of influence around 
    REAL(dbl)                 :: R_P, Sh_P, delta_P
    REAL(dbl)                 :: R_influence_P, L_influence_P
    REAL(dbl),DIMENSION(2)    :: VIB_x, VIB_y, VIB_z 			! Volume of Influence's Borders
    REAL(dbl),DIMENSION(2)    :: VIB_zc 			! Volume of Influence's Borders on coarse mesh
    REAL(dbl),DIMENSION(2)    :: VIB_zf 			! Volume of Influence's Borders on fine mesh
    REAL(dbl),DIMENSION(2)    :: NVB_x, NVB_y, NVB_z			! Node Volume's Borders
    INTEGER  ,DIMENSION(2)    :: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
    INTEGER  ,DIMENSION(2)    :: NEP_x, NEP_y, NEP_z                        ! Lattice Nodes Surronding the particle
    REAL(dbl)		  :: tmp, Overlap_sum, Overlap_sum_l, Overlap_sum_coarse, Overlap_sum_fine, checkEffVolumeOverlapFineMesh
    REAL(dbl)                 :: overlapCoarseProc, overlapFineProc
    TYPE(ParRecord), POINTER  :: current
    TYPE(ParRecord), POINTER  :: next
    INTEGER(lng)  		  :: mpierr
    
    
    delta_mesh = 1.0_dbl
    zcf3 = xcf_fine*ycf_fine*zcf_fine
    current => ParListHead%next
    
    DO WHILE (ASSOCIATED(current))
       
       !------ Copy pointer of next node
       next => current%next 
       
       IF ( flagParticleCF(current%pardata%parid) )  THEN  !Check if particle is in fine mesh
          
          
          !------ Calculate length scale for jth particle:  delta = R / Sh
          !------ Calculate effective radius: R_influence_P = R + (N_d *delta)
          !------ Note: need to convert this into Lattice units and not use the physical length units
          !------ Then compute equivalent cubic mesh length scale
          N_d = 1.0
          R_P  = current%pardata%rp
          Sh_P = current%pardata%sh
          delta_P = R_P / Sh_P
          R_influence_P = (R_P + N_d * delta_P)

          !------ Computing equivalent cubic mesh length scale
          L_influence_P = ( (4.0_dbl*PI/3.0_dbl)**(1.0_dbl/3.0_dbl) ) * R_influence_P
          
          !------ NEW: Volume of Influence Border (VIB) for this particle
          VIB_x(1)= current%pardata%xp - 0.5_dbl * L_influence_P 
          VIB_x(2)= current%pardata%xp + 0.5_dbl * L_influence_P
          VIB_y(1)= current%pardata%yp - 0.5_dbl * L_influence_P
          VIB_y(2)= current%pardata%yp + 0.5_dbl * L_influence_P
          VIB_z(1)= current%pardata%zp - 0.5_dbl * L_influence_P
          VIB_z(2)= current%pardata%zp + 0.5_dbl * L_influence_P
          
          !Check if Volume of Influence Border overlaps with the current processor domain
          
          if (MAX ( MIN(VIB_z(2)+L,(kMaxDomain(mySub)-0.5)*zcf) - MAX(VIB_z(1)+L,(kMinDomain(mySub)-1.5)*zcf), 0.0_dbl) .gt. 0) then
             zp_coarse = current%pardata%zp + L
             VIB_zc(1) = VIB_z(1)+L
             VIB_zc(2) = VIB_z(2)+L
             
          else  if (MAX ( MIN(VIB_z(2)-L,(kMaxDomain(mySub)-0.5)*zcf) - MAX(VIB_z(1)-L,(kMinDomain(mySub)-1.5)*zcf), 0.0_dbl) .gt. 0) then
             zp_coarse = current%pardata%zp - L             
             VIB_zc(1) = VIB_z(1)-L
             VIB_zc(2) = VIB_z(2)-L

          else
             zp_coarse = current%pardata%zp
             VIB_zc(1) = VIB_z(1)
             VIB_zc(2) = VIB_z(2)             
             
          end if
          
!          write(31,*) 'VIB_zc = ', VIB_zc(1), VIB_zc(2), ' zMinMax  = ', (kMaxDomain(mySub)-0.5)*zcf, (kMinDomain(mySub)-1.5)*zcf
          overlapCoarseProc = MAX ( MIN(VIB_zc(2),(kMaxDomain(mySub)-0.5)*zcf) - MAX(VIB_zc(1),(kMinDomain(mySub)-1.5)*zcf) , 0.0_dbl)
!          write(31,*) 'overlapCoarseProc = ', overlapCoarseProc, ' Ratio = ', overlapCoarseProc/(VIB_zc(2)-VIB_zc(1))
!          flush(31)
          
          
          if (MAX ( MIN(VIB_z(2)+L,(kMaxDomain_fine(mySub)-0.5)*zcf_fine) - MAX(VIB_z(1)+L,(kMinDomain_fine(mySub)-1.5)*zcf_fine), 0.0_dbl) .gt. 0) then
             zp_fine = current%pardata%zp + L             
             VIB_zf(1) = VIB_z(1)+L
             VIB_zf(2) = VIB_z(2)+L
             
          else  if (MAX ( MIN(VIB_z(2)-L,(kMaxDomain_fine(mySub)-0.5)*zcf_fine) - MAX(VIB_z(1)-L,(kMinDomain_fine(mySub)-1.5)*zcf_fine), 0.0_dbl) .gt. 0) then
             zp_fine = current%pardata%zp - L             
             VIB_zf(1) = VIB_z(1)-L
             VIB_zf(2) = VIB_z(2)-L

          else
             zp_fine = current%pardata%zp
             VIB_zf(1) = VIB_z(1)
             VIB_zf(2) = VIB_z(2)
             
          end if
          
          overlapFineProc = MAX ( MIN(VIB_zf(2),(kMaxDomain_fine(mySub)-0.5)*zcf_fine) - MAX(VIB_zf(1),(kMinDomain_fine(mySub)-1)*zcf_fine) , 0.0_dbl)
!          write(31,*) 'overlapFineProc = ', overlapFineProc, ' Ratio = ', overlapFineProc/(VIB_zf(2)-VIB_zf(1))
          
          Overlap_sum_coarse = 0.0_dbl
          Overlap_sum_fine = 0.0_dbl
          Overlap_sum_l = 0.0_dbl
          Overlap_sum = 0.0_dbl
          
          !If volume of influence extends into coarse mesh, also calculate scalar release for coarse mesh
          checkEffVolumeOverlapFineMesh = ( MAX ( MIN(VIB_x(2), fractionDfine * D * 0.5 + xcf) - MAX(VIB_x(1), -fractionDfine * D * 0.5 - xcf), 0.0_dbl) / L_influence_P ) * &
               ( MAX ( MIN(VIB_y(2),fractionDfine * D * 0.5 + xcf) - MAX(VIB_y(1), -fractionDfine * D * 0.5 - ycf), 0.0_dbl) / L_influence_P ) - 0.99
!          write(31,*) 'checkEffVolumeOverlapFineMesh = ', checkEffVolumeOverlapFineMesh
          
          if (overlapCoarseProc .gt. 0) then
             
             if (checkEffVolumeOverlapFineMesh .lt. 0) then
                
 !               write(31,*) 'Particle also overlaps with coarse mesh'
                !------ Finding particle location in this processor
                xp = (current%pardata%xp-xx(1))/xcf + 1 - REAL(iMin-1_lng,dbl)
                yp = (current%pardata%yp-yy(1))/ycf + 1 - REAL(jMin-1_lng,dbl)
                zp = (zp_coarse-zz(1))/zcf + 1 - REAL(kMin-1_lng,dbl)
                
                !------ NEW: Finding the lattice "Nodes Effected by Particle" 
                !------ NEW: Finding the lattice "Nodes Effected by Particle" 
                NEP_x(1)= MAX(1,NINT(xp - 0.5*L_influence_P/xcf))
                NEP_x(2)= MIN(nxSub,NINT(xp + 0.5*L_influence_P/xcf))
                NEP_y(1)= MAX(1,NINT(yp - 0.5*L_influence_P/xcf))
                NEP_y(2)= MIN(nySub,NINT(yp + 0.5*L_influence_P/ycf))
                NEP_z(1)= MAX(1,NINT(zp - 0.5*L_influence_P/xcf))
                NEP_z(2)= MIN(nzSub,NINT(zp + 0.5*L_influence_P/zcf))
                
                !------ NEW: Finding the volume overlapping between particle-effetive-volume and the volume around each lattice node
                
                Overlap= 0.0
                
                DO i= NEP_x(1),NEP_x(2) 
                   DO j= NEP_y(1),NEP_y(2)
                      DO k= NEP_z(1),NEP_z(2)
                         NVB_x(1) = x(i) - 0.5_dbl*xcf
                         NVB_x(2) = x(i) + 0.5_dbl*xcf
                         NVB_y(1) = y(j) - 0.5_dbl*ycf
                         NVB_y(2) = y(j) + 0.5_dbl*ycf
                         NVB_z(1) = z(k) - 0.5_dbl*zcf
                         NVB_z(2) = z(k) + 0.5_dbl*zcf
                         
                         tmp = MAX ( MIN(VIB_x(2),NVB_x(2)) - MAX(VIB_x(1),NVB_x(1)), 0.0_dbl) * &
                              MAX ( MIN(VIB_y(2),NVB_y(2)) - MAX(VIB_y(1),NVB_y(1)), 0.0_dbl) * &
                              MAX ( MIN(VIB_zc(2),NVB_z(2)) - MAX(VIB_zc(1),NVB_z(1)), 0.0_dbl)
                         
                         !---- Taking care of the periodic BC in Z-dir
                         kk = k 
                         ! IF (k .gt. nz) THEN
                         !    kk = k - (nz - 1)
                         ! ELSE IF (k .lt. 1) THEN
                         !    kk = k + (nz-1)
                         ! END IF
                         
!                         write(31,*) 'i,j,k = ', i,j,kk, ' Overlap_coarse = ', tmp
                         flush(31)
                         IF (node(i,j,kk) .EQ. FLUID) THEN
                            Overlap(i,j,kk)= tmp * (1.0-flagNodeIntersectFine(i,j,kk)) 
                            Overlap_sum_coarse = Overlap_sum_coarse + Overlap(i,j,kk)
                         END IF
                      END DO
                   END DO
                END DO
                
             else
!                write(31,*) 'Particle overlaps only with fine mesh and not with coarse mesh'
                
             end if
          end if
          
          if (overlapFineProc .gt. 0) then
             
             !------ Finding particle location in this processor
             xp = (current%pardata%xp-xx_fine(1))/xcf_fine + 1 - REAL(iMin_fine-1_lng,dbl)
             yp = (current%pardata%yp-yy_fine(1))/ycf_fine + 1 - REAL(jMin_fine-1_lng,dbl)
             zp = (zp_fine-zz_fine(1))/zcf_fine + 1 - REAL(kMin_fine-1_lng,dbl)
!             write(31,*) 'xp, yp, zp = ', xp, yp, zp
             
             !------ NEW: Finding the lattice "Nodes Effected by Particle" 
             NEP_x(1)= MAX(2,NINT(xp - 0.5*L_influence_P/xcf_fine))
             NEP_x(2)= MIN(nxSub_fine-1,NINT(xp + 0.5*L_influence_P/xcf_fine))
             NEP_y(1)= MAX(2,NINT(yp - 0.5*L_influence_P/xcf_fine))
             NEP_y(2)= MIN(nySub_fine-1,NINT(yp + 0.5*L_influence_P/ycf_fine))
             NEP_z(1)= MAX(1,NINT(zp - 0.5*L_influence_P/xcf_fine))
             NEP_z(2)= MIN(nzSub_fine,NINT(zp + 0.5*L_influence_P/zcf_fine))
!             write(31,*) 'NEP_x(1) = ', NEP_x(1), 'NEP_x(2) = ', NEP_x(2)
!             write(31,*) 'NEP_y(1) = ', NEP_y(1), 'NEP_y(2) = ', NEP_y(2)
!             write(31,*) 'NEP_z(1) = ', NEP_z(1), 'NEP_z(2) = ', NEP_z(2)
             
             !------ NEW: Finding the volume overlapping between particle-effetive-volume and the volume around each lattice node
             Overlap_sum_fine = 0.0_dbl
             Overlap_fine= 0.0
             
             DO i= NEP_x(1),NEP_x(2) 
                DO j= NEP_y(1),NEP_y(2)
                   DO k= NEP_z(1),NEP_z(2)
                      NVB_x(1) = x_fine(i) - 0.5_dbl*xcf_fine
                      NVB_x(2) = x_fine(i) + 0.5_dbl*xcf_fine
                      NVB_y(1) = y_fine(j) - 0.5_dbl*ycf_fine
                      NVB_y(2) = y_fine(j) + 0.5_dbl*ycf_fine
                      NVB_z(1) = z_fine(k) - 0.5_dbl*zcf_fine
                      NVB_z(2) = z_fine(k) + 0.5_dbl*zcf_fine
                      
                      tmp = MAX ( MIN(VIB_x(2),NVB_x(2)) - MAX(VIB_x(1),NVB_x(1)), 0.0_dbl) * &
                           MAX ( MIN(VIB_y(2),NVB_y(2)) - MAX(VIB_y(1),NVB_y(1)), 0.0_dbl) * &
                           MAX ( MIN(VIB_zf(2),NVB_z(2)) - MAX(VIB_zf(1),NVB_z(1)), 0.0_dbl)
                      
                      !---- Taking care of the periodic BC in Z-dir
                      kk = k 
                      ! IF (k .gt. nz_fine) THEN
                      !    kk = k - (nz_fine - 1)
                      ! ELSE IF (k .lt. 1) THEN
                      !    kk = k + (nz_fine-1)
                      ! END IF
                      
!                      write(31,*) 'i,j,k = ', i,j,kk, ' Overlap_fine = ', tmp
                      IF (node_fine(i,j,kk) .EQ. FLUID) THEN
                         Overlap_fine(i,j,kk)= tmp
                         Overlap_sum_fine= Overlap_sum_fine + Overlap_fine(i,j,kk)
                      END IF
                   END DO
                END DO
             END DO
             
          end if
          
          Overlap_sum_l = Overlap_sum_coarse + Overlap_sum_fine !Add the overlap in the coarse and the fine meshes.
          
!          write(31,*) 'Overlap_sum_coarse = ', Overlap_sum_coarse, ' Ratio = ', Overlap_sum_coarse/(L_influence_P * L_influence_P * L_influence_P)
!          write(31,*) 'Overlap_sum_fine = ', Overlap_sum_fine, ' Ratio = ', Overlap_sum_fine/(L_influence_P * L_influence_P * L_influence_P)
          CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
          CALL MPI_ALLREDUCE(Overlap_sum_l, Overlap_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

!          write(31,*) 'Overlap = ', Overlap_sum, ' Ratio = ', Overlap_sum/(L_influence_P * L_influence_P * L_influence_P)
          
          if (overlapFineProc .gt. 0) then
             
             !------ Computing particle release contribution to scalar field at each lattice node
             DO i= NEP_x(1),NEP_x(2)
                DO j= NEP_y(1),NEP_y(2)
                   DO k= NEP_z(1),NEP_z(2)
                      
                      !----- Taking care of the periodic BC in Z-dir
                      kk = k 
                      ! IF (k .gt. nz_fine) THEN
                      !    kk = k - (nz_fine-1)
                      ! ELSE IF (k .lt. 1) THEN
                      !    kk = k + (nz_fine-1)
                      ! END IF
                      
                      IF (node_fine(i,j,kk) .EQ. FLUID) THEN                 
                         
                         !----- Overlap_sum going to zero when the particle is disapearing
                         IF (Overlap_sum .gt. 1e-40) THEN
                            Overlap_fine(i,j,kk) = Overlap_fine(i,j,kk) / Overlap_sum
                         ELSE
                            Overlap_fine(i,j,kk) = 0.0
                         END IF
                         
                         delphi_particle_fine(i,j,kk)  = delphi_particle_fine(i,j,kk)  + current%pardata%delNB * Overlap_fine(i,j,kk)/ zcf3
!                         write(31,*) 'i,j,kk = ', i,j,kk, ' delphi_fine = ', current%pardata%delNB, Overlap_fine(i,j,kk), current%pardata%delNB * Overlap_fine(i,j,kk)
!                         flush(31)
                         
                         
                      END IF
                   END DO
                END DO
             END DO
             
          end if
          
          if (overlapCoarseProc .gt. 0) then
             
             !Now the coarse mesh
             if (checkEffVolumeOverlapFineMesh .lt. 0) then
                
                !------ Finding particle location in this processor
                xp = (current%pardata%xp-xx(1))/xcf + 1 - REAL(iMin-1_lng,dbl)
                yp = (current%pardata%yp-yy(1))/ycf + 1 - REAL(jMin-1_lng,dbl)
                zp = (zp_coarse-zz(1))/zcf + 1 - REAL(kMin-1_lng,dbl)
                
                !------ NEW: Finding the lattice "Nodes Effected by Particle" 
                !------ NEW: Finding the lattice "Nodes Effected by Particle" 
                NEP_x(1)= MAX(1,NINT(xp - 0.5*L_influence_P/xcf))
                NEP_x(2)= MIN(nxSub,NINT(xp + 0.5*L_influence_P/xcf))
                NEP_y(1)= MAX(1,NINT(yp - 0.5*L_influence_P/xcf))
                NEP_y(2)= MIN(nySub,NINT(yp + 0.5*L_influence_P/ycf))
                NEP_z(1)= MAX(1,NINT(zp - 0.5*L_influence_P/xcf))
                NEP_z(2)= MIN(nzSub,NINT(zp + 0.5*L_influence_P/zcf))
                
                !------ Computing particle release contribution to scalar field at each lattice node
                DO i= NEP_x(1),NEP_x(2)
                   DO j= NEP_y(1),NEP_y(2)
                      DO k= NEP_z(1),NEP_z(2)
                         
                         !----- Taking care of the periodic BC in Z-dir
                         kk = k 
                         ! IF (k .gt. nz) THEN
                         !    kk = k - (nz-1)
                         ! ELSE IF (k .lt. 1) THEN
                         !    kk = k + (nz-1)
                         ! END IF
                         
                         IF (node(i,j,kk) .EQ. FLUID) THEN                 
                            
                            !----- Overlap_sum going to zero when the particle is disapearing
                            IF (Overlap_sum .gt. 1e-40) THEN
                               Overlap(i,j,kk) = Overlap(i,j,kk) / Overlap_sum
                            ELSE
                               Overlap(i,j,kk) = 0.0
                            END IF
                            
                            delphi_particle(i,j,kk)  = delphi_particle(i,j,kk)  + current%pardata%delNB * Overlap(i,j,kk) / (xcf * ycf * zcf * (1.0-flagNodeIntersectFine(i,j,k)) )
!                            write(31,*) 'i,j,kk = ', i,j,kk, ' delphi = ', current%pardata%delNB, Overlap(i,j,kk), current%pardata%delNB * Overlap(i,j,kk)
                         END IF
                      END DO
                   END DO
                END DO
                
             end if
             
          end if
       END IF
       
       
       !------ point to next node in the list
       current => next
    ENDDO
    
    !===================================================================================================
  END SUBROUTINE Interp_ParToNodes_Conc_Fine
  !===================================================================================================
  
  
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
          velCtoF_bottomXZ(:,1,i,k) = velCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
          velCtoF_bottomXZ(:,2,i,k) = velCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
          f1 =  u(lCxIndex-1,45,lCzIndex)
          f2 =  u(lCxIndex,45,lCzIndex) 
          f3 =  u(lCxIndex+1,45,lCzIndex)
          f4 =  u(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(lCxIndex-1,45,lCzIndex)
          f2 =  v(lCxIndex,45,lCzIndex) 
          f3 =  v(lCxIndex+1,45,lCzIndex)
          f4 =  v(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(lCxIndex-1,45,lCzIndex)
          f2 =  w(lCxIndex,45,lCzIndex) 
          f3 =  w(lCxIndex+1,45,lCzIndex)
          f4 =  w(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(3,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step


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
          velCtoF_bottomXZ(:,1,i,nzSub_fine-gridRatio+1-k+1) = velCtoF_bottomXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
          velCtoF_bottomXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) = velCtoF_bottomXZ(:,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
          f1 =  u(lCxIndex-1,45,lCzIndex)
          f2 =  u(lCxIndex,45,lCzIndex)
          f3 =  u(lCxIndex+1,45,lCzIndex)
          f4 =  u(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(1,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(lCxIndex-1,45,lCzIndex)
          f2 =  v(lCxIndex,45,lCzIndex)
          f3 =  v(lCxIndex+1,45,lCzIndex)
          f4 =  v(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(2,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time ste        
          f1 =  w(lCxIndex-1,45,lCzIndex)
          f2 =  w(lCxIndex,45,lCzIndex)
          f3 =  w(lCxIndex+1,45,lCzIndex)
          f4 =  w(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(3,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time ste        


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
          velCtoF_topXZ(:,1,i,k) = velCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
          velCtoF_topXZ(:,2,i,k) = velCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
          f1 =  u(lCxIndex-1,57,lCzIndex) 
          f2 =  u(lCxIndex,57,lCzIndex) 
          f3 =  u(lCxIndex+1,57,lCzIndex) 
          f4 =  u(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(lCxIndex-1,57,lCzIndex) 
          f2 =  v(lCxIndex,57,lCzIndex) 
          f3 =  v(lCxIndex+1,57,lCzIndex) 
          f4 =  v(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(lCxIndex-1,57,lCzIndex) 
          f2 =  w(lCxIndex,57,lCzIndex) 
          f3 =  w(lCxIndex+1,57,lCzIndex) 
          f4 =  w(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(3,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          
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
          velCtoF_topXZ(:,1,i,nzSub_fine-gridRatio+1-k+1) = velCtoF_topXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
          velCtoF_topXZ(:,2,i,nzSub_fine-gridRatio+1-k+1) = velCtoF_topXZ(:,3,i,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
          f1 =  u(lCxIndex-1,57,lCzIndex) 
          f2 =  u(lCxIndex,57,lCzIndex) 
          f3 =  u(lCxIndex+1,57,lCzIndex) 
          f4 =  u(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(1,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(lCxIndex-1,57,lCzIndex) 
          f2 =  v(lCxIndex,57,lCzIndex) 
          f3 =  v(lCxIndex+1,57,lCzIndex) 
          f4 =  v(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(2,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(lCxIndex-1,57,lCzIndex) 
          f2 =  w(lCxIndex,57,lCzIndex) 
          f3 =  w(lCxIndex+1,57,lCzIndex) 
          f4 =  w(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(3,3,i,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),xInterp) !Interpolate the latest value to the last(third) time step
          
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
          velCtoF_frontYZ(:,1,j,k) = velCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
          velCtoF_frontYZ(:,2,j,k) = velCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
          f1 =  u(45,lCyIndex-1,lCzIndex) 
          f2 =  u(45,lCyIndex,lCzIndex) 
          f3 =  u(45,lCyIndex+1,lCzIndex) 
          f4 =  u(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(45,lCyIndex-1,lCzIndex) 
          f2 =  v(45,lCyIndex,lCzIndex) 
          f3 =  v(45,lCyIndex+1,lCzIndex) 
          f4 =  v(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(45,lCyIndex-1,lCzIndex) 
          f2 =  w(45,lCyIndex,lCzIndex) 
          f3 =  w(45,lCyIndex+1,lCzIndex) 
          f4 =  w(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(3,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
          
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
          velCtoF_frontYZ(:,1,j,nzSub_fine-gridRatio+1-k+1) = velCtoF_frontYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
          velCtoF_frontYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) = velCtoF_frontYZ(:,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
          f1 =  u(45,lCyIndex-1,lCzIndex) 
          f2 =  u(45,lCyIndex,lCzIndex) 
          f3 =  u(45,lCyIndex+1,lCzIndex) 
          f4 =  u(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(1,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(45,lCyIndex-1,lCzIndex) 
          f2 =  v(45,lCyIndex,lCzIndex) 
          f3 =  v(45,lCyIndex+1,lCzIndex) 
          f4 =  v(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(2,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(45,lCyIndex-1,lCzIndex) 
          f2 =  w(45,lCyIndex,lCzIndex) 
          f3 =  w(45,lCyIndex+1,lCzIndex) 
          f4 =  w(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(3,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
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
          velCtoF_backYZ(:,1,j,k) = velCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
          velCtoF_backYZ(:,2,j,k) = velCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
          f1 =  u(57,lCyIndex-1,lCzIndex) 
          f2 =  u(57,lCyIndex,lCzIndex) 
          f3 =  u(57,lCyIndex+1,lCzIndex) 
          f4 =  u(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(57,lCyIndex-1,lCzIndex) 
          f2 =  v(57,lCyIndex,lCzIndex) 
          f3 =  v(57,lCyIndex+1,lCzIndex) 
          f4 =  v(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(57,lCyIndex-1,lCzIndex) 
          f2 =  w(57,lCyIndex,lCzIndex) 
          f3 =  w(57,lCyIndex+1,lCzIndex) 
          f4 =  w(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(3,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
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
          velCtoF_backYZ(:,1,j,nzSub_fine-gridRatio+1-k+1) = velCtoF_backYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) !Cycle the second time step to the first time step
          velCtoF_backYZ(:,2,j,nzSub_fine-gridRatio+1-k+1) = velCtoF_backYZ(:,3,j,nzSub_fine-gridRatio+1-k+1) !Cycle the last time step to the second time step
          f1 =  u(57,lCyIndex-1,lCzIndex) 
          f2 =  u(57,lCyIndex,lCzIndex) 
          f3 =  u(57,lCyIndex+1,lCzIndex) 
          f4 =  u(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(1,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(57,lCyIndex-1,lCzIndex) 
          f2 =  v(57,lCyIndex,lCzIndex) 
          f3 =  v(57,lCyIndex+1,lCzIndex) 
          f4 =  v(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(2,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(57,lCyIndex-1,lCzIndex) 
          f2 =  w(57,lCyIndex,lCzIndex) 
          f3 =  w(57,lCyIndex+1,lCzIndex) 
          f4 =  w(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(3,3,j,nzSub_fine-gridRatio+1-k+1) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
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
          velCtoF_bottomXZ(:,1,i,k) = velCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
          velCtoF_bottomXZ(:,2,i,k) = velCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
          f1 =  u(lCxIndex-1,45,lCzIndex)
          f2 =  u(lCxIndex,45,lCzIndex) 
          f3 =  u(lCxIndex+1,45,lCzIndex)
          f4 =  u(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(lCxIndex-1,45,lCzIndex)
          f2 =  v(lCxIndex,45,lCzIndex) 
          f3 =  v(lCxIndex+1,45,lCzIndex)
          f4 =  v(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(lCxIndex-1,45,lCzIndex)
          f2 =  w(lCxIndex,45,lCzIndex) 
          f3 =  w(lCxIndex+1,45,lCzIndex)
          f4 =  w(lCxIndex+2,45,lCzIndex)
          velCtoF_bottomXZ(3,3,i,k) = spatialInterpolate(f1,f2,f3,f4,node(lCxIndex-1,45,lCzIndex),node(lCxIndex,45,lCzIndex),node(lCxIndex+1,45,lCzIndex),node(lCxIndex+2,45,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
          
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
          velCtoF_topXZ(:,1,i,k) = velCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
          velCtoF_topXZ(:,2,i,k) = velCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
          f1 =  u(lCxIndex-1,57,lCzIndex) 
          f2 =  u(lCxIndex,57,lCzIndex) 
          f3 =  u(lCxIndex+1,57,lCzIndex) 
          f4 =  u(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(1,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(lCxIndex-1,57,lCzIndex) 
          f2 =  v(lCxIndex,57,lCzIndex) 
          f3 =  v(lCxIndex+1,57,lCzIndex) 
          f4 =  v(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(2,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(lCxIndex-1,57,lCzIndex) 
          f2 =  w(lCxIndex,57,lCzIndex) 
          f3 =  w(lCxIndex+1,57,lCzIndex) 
          f4 =  w(lCxIndex+2,57,lCzIndex) 
          velCtoF_topXZ(3,3,i,k) = spatialInterpolate(f1,f2,f3,f4, node(lCxIndex-1,57,lCzIndex),node(lCxIndex,57,lCzIndex),node(lCxIndex+1,57,lCzIndex),node(lCxIndex+2,57,lCzIndex), xInterp) !Interpolate the latest value to the last(third) time step
          
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
             fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCyIndex-1,45,lCzIndex),node(lCyIndex,45,lCzIndex),node(lCyIndex+1,45,lCzIndex),node(lCyIndex+2,45,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
             
             fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
             fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
             f1 =  feqFF_backYZ(m,lCyIndex-1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex-1,lCzIndex) - feqFF_backYZ(m,lCyIndex-1,lCzIndex))
             f2 =  feqFF_backYZ(m,lCyIndex,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex,lCzIndex) - feqFF_backYZ(m,lCyIndex,lCzIndex))
             f3 =  feqFF_backYZ(m,lCyIndex+1,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+1,lCzIndex) - feqFF_backYZ(m,lCyIndex+1,lCzIndex))
             f4 =  feqFF_backYZ(m,lCyIndex+2,lCzIndex) + (tau_fine - 1.0)/(gridRatio * (tau - 1.0)) * (f(m,57,lCyIndex+2,lCzIndex) - feqFF_backYZ(m,lCyIndex+2,lCzIndex))           
             fCtoF_backYZ(m,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(57,lCyIndex-1,lCzIndex),node(57,lCyIndex,lCzIndex),node(57,lCyIndex+1,lCzIndex),node(57,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
             
          end do
          
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
          velCtoF_frontYZ(:,1,j,k) = velCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
          velCtoF_frontYZ(:,2,j,k) = velCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
          f1 =  u(45,lCyIndex-1,lCzIndex) 
          f2 =  u(45,lCyIndex,lCzIndex) 
          f3 =  u(45,lCyIndex+1,lCzIndex) 
          f4 =  u(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(45,lCyIndex-1,lCzIndex) 
          f2 =  v(45,lCyIndex,lCzIndex) 
          f3 =  v(45,lCyIndex+1,lCzIndex) 
          f4 =  v(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(45,lCyIndex-1,lCzIndex) 
          f2 =  w(45,lCyIndex,lCzIndex) 
          f3 =  w(45,lCyIndex+1,lCzIndex) 
          f4 =  w(45,lCyIndex+2,lCzIndex) 
          velCtoF_frontYZ(3,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(45,lCyIndex-1,lCzIndex),node(45,lCyIndex,lCzIndex),node(45,lCyIndex+1,lCzIndex),node(45,lCyIndex+2,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
          dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
          dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
          f1 =  rho(57,lCyIndex-1,lCzIndex) 
          f2 =  rho(57,lCyIndex,lCzIndex) 
          f3 =  rho(57,lCyIndex+1,lCzIndex) 
          f4 =  rho(57,lCyIndex+2,lCzIndex) 
          dsCtoF_backYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCyIndex-1,57,lCzIndex),node(lCyIndex,57,lCzIndex),node(lCyIndex+1,57,lCzIndex),node(lCyIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  phi(57,lCyIndex-1,lCzIndex) 
          f2 =  phi(57,lCyIndex,lCzIndex) 
          f3 =  phi(57,lCyIndex+1,lCzIndex) 
          f4 =  phi(57,lCyIndex+2,lCzIndex) 
          dsCtoF_backYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCyIndex-1,57,lCzIndex),node(lCyIndex,57,lCzIndex),node(lCyIndex+1,57,lCzIndex),node(lCyIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          velCtoF_backYZ(:,1,j,k) = velCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
          velCtoF_backYZ(:,2,j,k) = velCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
          f1 =  u(57,lCyIndex-1,lCzIndex) 
          f2 =  u(57,lCyIndex,lCzIndex) 
          f3 =  u(57,lCyIndex+1,lCzIndex) 
          f4 =  u(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(1,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCyIndex-1,57,lCzIndex),node(lCyIndex,57,lCzIndex),node(lCyIndex+1,57,lCzIndex),node(lCyIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  v(57,lCyIndex-1,lCzIndex) 
          f2 =  v(57,lCyIndex,lCzIndex) 
          f3 =  v(57,lCyIndex+1,lCzIndex) 
          f4 =  v(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(2,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCyIndex-1,57,lCzIndex),node(lCyIndex,57,lCzIndex),node(lCyIndex+1,57,lCzIndex),node(lCyIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          f1 =  w(57,lCyIndex-1,lCzIndex) 
          f2 =  w(57,lCyIndex,lCzIndex) 
          f3 =  w(57,lCyIndex+1,lCzIndex) 
          f4 =  w(57,lCyIndex+2,lCzIndex) 
          velCtoF_backYZ(3,3,j,k) = spatialInterpolate(f1,f2,f3,f4,node(lCyIndex-1,57,lCzIndex),node(lCyIndex,57,lCzIndex),node(lCyIndex+1,57,lCzIndex),node(lCyIndex+2,57,lCzIndex),yInterp) !Interpolate the latest value to the last(third) time step
          
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
                fCtoF_bottomXZ(m,3,i,k) = spatialInterpolate(fCtoF_bottomXZ(m,3,i,lFzIndex-gridRatio),fCtoF_bottomXZ(m,3,i,lFzIndex),fCtoF_bottomXZ(m,3,i,lFzIndex+gridRatio),fCtoF_bottomXZ(m,3,i,lFzIndex+2*gridRatio), node_fine_bottomXZ(3,i,lFzIndex-gridRatio), node_fine_bottomXZ(3,i,lFzIndex), node_fine_bottomXZ(3,i,lFzIndex+gridRatio), node_fine_bottomXZ(3,i,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
                
                fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,2,i,k) !Cycle the second time step to the first time step
                fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
                fCtoF_topXZ(m,3,i,k) = spatialInterpolate(fCtoF_topXZ(m,3,i,lFzIndex-gridRatio),fCtoF_topXZ(m,3,i,lFzIndex),fCtoF_topXZ(m,3,i,lFzIndex+gridRatio),fCtoF_topXZ(m,3,i,lFzIndex+2*gridRatio), node_fine_topXZ(3,i,lFzIndex-gridRatio), node_fine_topXZ(3,i,lFzIndex), node_fine_topXZ(3,i,lFzIndex+gridRatio), node_fine_topXZ(3,i,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
                
             end do
             
             dsCtoF_bottomXZ(:,1,i,k) = dsCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
             dsCtoF_bottomXZ(:,2,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
             dsCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(dsCtoF_bottomXZ(1,3,i,lFzIndex-gridRatio),dsCtoF_bottomXZ(1,3,i,lFzIndex),dsCtoF_bottomXZ(1,3,i,lFzIndex+gridRatio),dsCtoF_bottomXZ(1,3,i,lFzIndex+2*gridRatio),node_fine_bottomXZ(3,i,lFzIndex-gridRatio), node_fine_bottomXZ(3,i,lFzIndex), node_fine_bottomXZ(3,i,lFzIndex+gridRatio), node_fine_bottomXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             dsCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(dsCtoF_bottomXZ(2,3,i,lFzIndex-gridRatio),dsCtoF_bottomXZ(2,3,i,lFzIndex),dsCtoF_bottomXZ(2,3,i,lFzIndex+gridRatio),dsCtoF_bottomXZ(2,3,i,lFzIndex+2*gridRatio),node_fine_bottomXZ(3,i,lFzIndex-gridRatio), node_fine_bottomXZ(3,i,lFzIndex), node_fine_bottomXZ(3,i,lFzIndex+gridRatio), node_fine_bottomXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             
             velCtoF_bottomXZ(:,1,i,k) = velCtoF_bottomXZ(:,2,i,k) !Cycle the second time step to the first time step
             velCtoF_bottomXZ(:,2,i,k) = velCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
             velCtoF_bottomXZ(1,3,i,k) = spatialInterpolate(velCtoF_bottomXZ(1,3,i,lFzIndex-gridRatio),velCtoF_bottomXZ(1,3,i,lFzIndex),velCtoF_bottomXZ(1,3,i,lFzIndex+gridRatio),velCtoF_bottomXZ(1,3,i,lFzIndex+2*gridRatio),node_fine_bottomXZ(3,i,lFzIndex-gridRatio), node_fine_bottomXZ(3,i,lFzIndex), node_fine_bottomXZ(3,i,lFzIndex+gridRatio), node_fine_bottomXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             velCtoF_bottomXZ(2,3,i,k) = spatialInterpolate(velCtoF_bottomXZ(2,3,i,lFzIndex-gridRatio),velCtoF_bottomXZ(2,3,i,lFzIndex),velCtoF_bottomXZ(2,3,i,lFzIndex+gridRatio),velCtoF_bottomXZ(2,3,i,lFzIndex+2*gridRatio),node_fine_bottomXZ(3,i,lFzIndex-gridRatio), node_fine_bottomXZ(3,i,lFzIndex), node_fine_bottomXZ(3,i,lFzIndex+gridRatio), node_fine_bottomXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             velCtoF_bottomXZ(3,3,i,k) = spatialInterpolate(velCtoF_bottomXZ(3,3,i,lFzIndex-gridRatio),velCtoF_bottomXZ(3,3,i,lFzIndex),velCtoF_bottomXZ(3,3,i,lFzIndex+gridRatio),velCtoF_bottomXZ(3,3,i,lFzIndex+2*gridRatio),node_fine_bottomXZ(3,i,lFzIndex-gridRatio), node_fine_bottomXZ(3,i,lFzIndex), node_fine_bottomXZ(3,i,lFzIndex+gridRatio), node_fine_bottomXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step             

             dsCtoF_topXZ(:,1,i,k) = dsCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
             dsCtoF_topXZ(:,2,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
             dsCtoF_topXZ(1,3,i,k) = spatialInterpolate(dsCtoF_topXZ(1,3,i,lFzIndex-gridRatio),dsCtoF_topXZ(1,3,i,lFzIndex),dsCtoF_topXZ(1,3,i,lFzIndex+gridRatio),dsCtoF_topXZ(1,3,i,lFzIndex+2*gridRatio), node_fine_topXZ(3,i,lFzIndex-gridRatio), node_fine_topXZ(3,i,lFzIndex), node_fine_topXZ(3,i,lFzIndex+gridRatio), node_fine_topXZ(3,i,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
             dsCtoF_topXZ(2,3,i,k) = spatialInterpolate(dsCtoF_topXZ(2,3,i,lFzIndex-gridRatio),dsCtoF_topXZ(2,3,i,lFzIndex),dsCtoF_topXZ(2,3,i,lFzIndex+gridRatio),dsCtoF_topXZ(2,3,i,lFzIndex+2*gridRatio), node_fine_topXZ(3,i,lFzIndex-gridRatio), node_fine_topXZ(3,i,lFzIndex), node_fine_topXZ(3,i,lFzIndex+gridRatio), node_fine_topXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             
             velCtoF_topXZ(:,1,i,k) = velCtoF_topXZ(:,2,i,k) !Cycle the second time step to the first time step
             velCtoF_topXZ(:,2,i,k) = velCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
             velCtoF_topXZ(1,3,i,k) = spatialInterpolate(velCtoF_topXZ(1,3,i,lFzIndex-gridRatio),velCtoF_topXZ(1,3,i,lFzIndex),velCtoF_topXZ(1,3,i,lFzIndex+gridRatio),velCtoF_topXZ(1,3,i,lFzIndex+2*gridRatio), node_fine_topXZ(3,i,lFzIndex-gridRatio), node_fine_topXZ(3,i,lFzIndex), node_fine_topXZ(3,i,lFzIndex+gridRatio), node_fine_topXZ(3,i,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
             velCtoF_topXZ(2,3,i,k) = spatialInterpolate(velCtoF_topXZ(2,3,i,lFzIndex-gridRatio),velCtoF_topXZ(2,3,i,lFzIndex),velCtoF_topXZ(2,3,i,lFzIndex+gridRatio),velCtoF_topXZ(2,3,i,lFzIndex+2*gridRatio), node_fine_topXZ(3,i,lFzIndex-gridRatio), node_fine_topXZ(3,i,lFzIndex), node_fine_topXZ(3,i,lFzIndex+gridRatio), node_fine_topXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             velCtoF_topXZ(3,3,i,k) = spatialInterpolate(velCtoF_topXZ(3,3,i,lFzIndex-gridRatio),velCtoF_topXZ(3,3,i,lFzIndex),velCtoF_topXZ(3,3,i,lFzIndex+gridRatio),velCtoF_topXZ(3,3,i,lFzIndex+2*gridRatio), node_fine_topXZ(3,i,lFzIndex-gridRatio), node_fine_topXZ(3,i,lFzIndex), node_fine_topXZ(3,i,lFzIndex+gridRatio), node_fine_topXZ(3,i,lFzIndex+2*gridRatio),zInterp) !Interpolate the latest value to the last(third) time step
             
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
                fCtoF_frontYZ(m,3,j,k) = spatialInterpolate(fCtoF_frontYZ(m,3,j,lFzIndex-gridRatio),fCtoF_frontYZ(m,3,j,lFzIndex),fCtoF_frontYZ(m,3,j,lFzIndex+gridRatio),fCtoF_frontYZ(m,3,j,lFzIndex+2*gridRatio), node_fine_frontYZ(3,j,lFzIndex-gridRatio), node_fine_frontYZ(3,j,lFzIndex), node_fine_frontYZ(3,j,lFzIndex+gridRatio), node_fine_frontYZ(3,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
                
                fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,2,j,k) !Cycle the second time step to the first time step
                fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
                fCtoF_backYZ(m,3,j,k) = spatialInterpolate(fCtoF_backYZ(m,3,j,lFzIndex-gridRatio),fCtoF_backYZ(m,3,j,lFzIndex),fCtoF_backYZ(m,3,j,lFzIndex+gridRatio),fCtoF_backYZ(m,3,j,lFzIndex+2*gridRatio), node_fine_backYZ(3,j,lFzIndex-gridRatio), node_fine_backYZ(3,j,lFzIndex), node_fine_backYZ(3,j,lFzIndex+gridRatio), node_fine_backYZ(3,j,lFzIndex+2*gridRatio), zInterp)              
                
             end do
             
             dsCtoF_frontYZ(:,1,j,k) = dsCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
             dsCtoF_frontYZ(:,2,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
             dsCtoF_frontYZ(1,3,j,k) = spatialInterpolate(dsCtoF_frontYZ(1,3,j,lFzIndex-gridRatio),dsCtoF_frontYZ(1,3,j,lFzIndex),dsCtoF_frontYZ(1,3,j,lFzIndex+gridRatio),dsCtoF_frontYZ(1,3,j,lFzIndex+2*gridRatio), node_fine_frontYZ(3,j,lFzIndex-gridRatio), node_fine_frontYZ(3,j,lFzIndex), node_fine_frontYZ(3,j,lFzIndex+gridRatio), node_fine_frontYZ(3,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
             dsCtoF_frontYZ(2,3,j,k) = spatialInterpolate(dsCtoF_frontYZ(2,3,j,lFzIndex-gridRatio),dsCtoF_frontYZ(2,3,j,lFzIndex),dsCtoF_frontYZ(2,3,j,lFzIndex+gridRatio),dsCtoF_frontYZ(2,3,j,lFzIndex+2*gridRatio), node_fine_frontYZ(3,j,lFzIndex-gridRatio), node_fine_frontYZ(3,j,lFzIndex), node_fine_frontYZ(3,j,lFzIndex+gridRatio), node_fine_frontYZ(3,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step

             velCtoF_frontYZ(:,1,j,k) = velCtoF_frontYZ(:,2,j,k) !Cycle the second time step to the first time step
             velCtoF_frontYZ(:,2,j,k) = velCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
             velCtoF_frontYZ(1,3,j,k) = spatialInterpolate(velCtoF_frontYZ(1,3,j,lFzIndex-gridRatio),velCtoF_frontYZ(1,3,j,lFzIndex),velCtoF_frontYZ(1,3,j,lFzIndex+gridRatio),velCtoF_frontYZ(1,3,j,lFzIndex+2*gridRatio), node_fine_frontYZ(3,j,lFzIndex-gridRatio), node_fine_frontYZ(3,j,lFzIndex), node_fine_frontYZ(3,j,lFzIndex+gridRatio), node_fine_frontYZ(3,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
             velCtoF_frontYZ(2,3,j,k) = spatialInterpolate(velCtoF_frontYZ(2,3,j,lFzIndex-gridRatio),velCtoF_frontYZ(2,3,j,lFzIndex),velCtoF_frontYZ(2,3,j,lFzIndex+gridRatio),velCtoF_frontYZ(2,3,j,lFzIndex+2*gridRatio), node_fine_frontYZ(3,j,lFzIndex-gridRatio), node_fine_frontYZ(3,j,lFzIndex), node_fine_frontYZ(3,j,lFzIndex+gridRatio), node_fine_frontYZ(3,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
             velCtoF_frontYZ(3,3,j,k) = spatialInterpolate(velCtoF_frontYZ(3,3,j,lFzIndex-gridRatio),velCtoF_frontYZ(3,3,j,lFzIndex),velCtoF_frontYZ(3,3,j,lFzIndex+gridRatio),velCtoF_frontYZ(3,3,j,lFzIndex+2*gridRatio), node_fine_frontYZ(3,j,lFzIndex-gridRatio), node_fine_frontYZ(3,j,lFzIndex), node_fine_frontYZ(3,j,lFzIndex+gridRatio), node_fine_frontYZ(3,j,lFzIndex+2*gridRatio), zInterp) !Interpolate the latest value to the last(third) time step
             
             dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
             dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
             dsCtoF_backYZ(1,3,j,k) = spatialInterpolate(dsCtoF_backYZ(1,3,j,lFzIndex-gridRatio),dsCtoF_backYZ(1,3,j,lFzIndex),dsCtoF_backYZ(1,3,j,lFzIndex+gridRatio),dsCtoF_backYZ(1,3,j,lFzIndex+2*gridRatio),  node_fine_backYZ(3,j,lFzIndex-gridRatio), node_fine_backYZ(3,j,lFzIndex), node_fine_backYZ(3,j,lFzIndex+gridRatio), node_fine_backYZ(3,j,lFzIndex+2*gridRatio), zInterp) 
             dsCtoF_backYZ(2,3,j,k) = spatialInterpolate(dsCtoF_backYZ(2,3,j,lFzIndex-gridRatio),dsCtoF_backYZ(2,3,j,lFzIndex),dsCtoF_backYZ(2,3,j,lFzIndex+gridRatio),dsCtoF_backYZ(2,3,j,lFzIndex+2*gridRatio),  node_fine_backYZ(3,j,lFzIndex-gridRatio), node_fine_backYZ(3,j,lFzIndex), node_fine_backYZ(3,j,lFzIndex+gridRatio), node_fine_backYZ(3,j,lFzIndex+2*gridRatio), zInterp)

             velCtoF_backYZ(:,1,j,k) = velCtoF_backYZ(:,2,j,k) !Cycle the second time step to the first time step
             velCtoF_backYZ(:,2,j,k) = velCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step
             velCtoF_backYZ(1,3,j,k) = spatialInterpolate(velCtoF_backYZ(1,3,j,lFzIndex-gridRatio),velCtoF_backYZ(1,3,j,lFzIndex),velCtoF_backYZ(1,3,j,lFzIndex+gridRatio),velCtoF_backYZ(1,3,j,lFzIndex+2*gridRatio),  node_fine_backYZ(3,j,lFzIndex-gridRatio), node_fine_backYZ(3,j,lFzIndex), node_fine_backYZ(3,j,lFzIndex+gridRatio), node_fine_backYZ(3,j,lFzIndex+2*gridRatio), zInterp) 
             velCtoF_backYZ(2,3,j,k) = spatialInterpolate(velCtoF_backYZ(2,3,j,lFzIndex-gridRatio),velCtoF_backYZ(2,3,j,lFzIndex),velCtoF_backYZ(2,3,j,lFzIndex+gridRatio),velCtoF_backYZ(2,3,j,lFzIndex+2*gridRatio),  node_fine_backYZ(3,j,lFzIndex-gridRatio), node_fine_backYZ(3,j,lFzIndex), node_fine_backYZ(3,j,lFzIndex+gridRatio), node_fine_backYZ(3,j,lFzIndex+2*gridRatio), zInterp)
             velCtoF_backYZ(3,3,j,k) = spatialInterpolate(velCtoF_backYZ(3,3,j,lFzIndex-gridRatio),velCtoF_backYZ(3,3,j,lFzIndex),velCtoF_backYZ(3,3,j,lFzIndex+gridRatio),velCtoF_backYZ(3,3,j,lFzIndex+2*gridRatio),  node_fine_backYZ(3,j,lFzIndex-gridRatio), node_fine_backYZ(3,j,lFzIndex), node_fine_backYZ(3,j,lFzIndex+gridRatio), node_fine_backYZ(3,j,lFzIndex+2*gridRatio), zInterp)
             
          end do
       END IF
    end do
    
  END SUBROUTINE ZSpatialInterpolateToFineGrid
  
  SUBROUTINE TemporalInterpolateToFineGrid
    
    INTEGER :: i,j,k,m
    REAL(dbl) :: tInterp !The time to which the temporal interpolation has to be done - Non-dimensionalized by the coarse mesh time step.
    REAL(dbl) :: tmp
    
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
          u_fine(i,1,k) = temporalInterpolate(velCtoF_bottomXZ(1,1,i,k),velCtoF_bottomXZ(1,2,i,k),velCtoF_bottomXZ(1,3,i,k), node_fine_bottomXZ(1,i,k), node_fine_bottomXZ(2,i,k), node_fine_bottomXZ(3,i,k), tInterp)
          v_fine(i,1,k) = temporalInterpolate(velCtoF_bottomXZ(2,1,i,k),velCtoF_bottomXZ(2,2,i,k),velCtoF_bottomXZ(2,3,i,k), node_fine_bottomXZ(1,i,k), node_fine_bottomXZ(2,i,k), node_fine_bottomXZ(3,i,k), tInterp)
          w_fine(i,1,k) = temporalInterpolate(velCtoF_bottomXZ(3,1,i,k),velCtoF_bottomXZ(3,2,i,k),velCtoF_bottomXZ(3,3,i,k), node_fine_bottomXZ(1,i,k), node_fine_bottomXZ(2,i,k), node_fine_bottomXZ(3,i,k), tInterp)          
          
          rho_fine(i,ny_fine,k) = temporalInterpolate(dsCtoF_topXZ(1,1,i,k),dsCtoF_topXZ(1,2,i,k),dsCtoF_topXZ(1,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
          phi_fine(i,ny_fine,k) = temporalInterpolate(dsCtoF_topXZ(2,1,i,k),dsCtoF_topXZ(2,2,i,k),dsCtoF_topXZ(2,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
          u_fine(i,ny_fine,k) = temporalInterpolate(velCtoF_topXZ(1,1,i,k),velCtoF_topXZ(1,2,i,k),velCtoF_topXZ(1,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
          v_fine(i,ny_fine,k) = temporalInterpolate(velCtoF_topXZ(2,1,i,k),velCtoF_topXZ(2,2,i,k),velCtoF_topXZ(2,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
          w_fine(i,ny_fine,k) = temporalInterpolate(velCtoF_topXZ(3,1,i,k),velCtoF_topXZ(3,2,i,k),velCtoF_topXZ(3,3,i,k), node_fine_topXZ(1,i,k), node_fine_topXZ(2,i,k), node_fine_topXZ(3,i,k), tInterp)
          
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
          u_fine(1,j,k) = temporalInterpolate(velCtoF_frontYZ(1,1,j,k),velCtoF_frontYZ(1,2,j,k),velCtoF_frontYZ(1,3,j,k), node_fine_frontYZ(1,j,k), node_fine_frontYZ(2,j,k), node_fine_frontYZ(3,j,k), tInterp)
          v_fine(1,j,k) = temporalInterpolate(velCtoF_frontYZ(2,1,j,k),velCtoF_frontYZ(2,2,j,k),velCtoF_frontYZ(2,3,j,k), node_fine_frontYZ(1,j,k), node_fine_frontYZ(2,j,k), node_fine_frontYZ(3,j,k), tInterp)
          w_fine(1,j,k) = temporalInterpolate(velCtoF_frontYZ(3,1,j,k),velCtoF_frontYZ(3,2,j,k),velCtoF_frontYZ(3,3,j,k), node_fine_frontYZ(1,j,k), node_fine_frontYZ(2,j,k), node_fine_frontYZ(3,j,k), tInterp)

          rho_fine(nx_fine,j,k) = temporalInterpolate(dsCtoF_backYZ(1,1,j,k),dsCtoF_backYZ(1,2,j,k),dsCtoF_backYZ(1,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)
          phi_fine(nx_fine,j,k) = temporalInterpolate(dsCtoF_backYZ(2,1,j,k),dsCtoF_backYZ(2,2,j,k),dsCtoF_backYZ(2,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)        
          u_fine(nx_fine,j,k) = temporalInterpolate(velCtoF_backYZ(1,1,j,k),velCtoF_backYZ(1,2,j,k),velCtoF_backYZ(1,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)
          v_fine(nx_fine,j,k) = temporalInterpolate(velCtoF_backYZ(2,1,j,k),velCtoF_backYZ(2,2,j,k),velCtoF_backYZ(2,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)
          w_fine(nx_fine,j,k) = temporalInterpolate(velCtoF_backYZ(3,1,j,k),velCtoF_backYZ(3,2,j,k),velCtoF_backYZ(3,3,j,k), node_fine_backYZ(1,j,k), node_fine_backYZ(2,j,k), node_fine_backYZ(3,j,k), tInterp)        
          
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
    INTEGER   :: n1,n2,n3
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
       else
          temporalInterpolate = f3
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
    
    do k=-gridRatio+1,nzSub_fine+gridRatio+1
       do i=1,nxSub_fine
          do m=0,NumDistDirs
             fCtoF_bottomXZ(m,1,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the first time step
             fCtoF_bottomXZ(m,2,i,k) = fCtoF_bottomXZ(m,3,i,k) !Cycle the last time step to the second time step
             fCtoF_topXZ(m,1,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the first time step
             fCtoF_topXZ(m,2,i,k) = fCtoF_topXZ(m,3,i,k) !Cycle the last time step to the second time step
          end do
          dsCtoF_bottomXZ(:,1,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the first time step
          dsCtoF_bottomXZ(:,2,i,k) = dsCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
          velCtoF_bottomXZ(:,1,i,k) = velCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the first time step
          velCtoF_bottomXZ(:,2,i,k) = velCtoF_bottomXZ(:,3,i,k) !Cycle the last time step to the second time step
          
          dsCtoF_topXZ(:,1,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the first time step
          dsCtoF_topXZ(:,2,i,k) = dsCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step
          velCtoF_topXZ(:,1,i,k) = velCtoF_topXZ(:,3,i,k) !Cycle the last time step to the first time step
          velCtoF_topXZ(:,2,i,k) = velCtoF_topXZ(:,3,i,k) !Cycle the last time step to the second time step

          node_fine_bottomXZ(1,i,k) = node_fine_bottomXZ(3,i,k)
          node_fine_bottomXZ(2,i,k) = node_fine_bottomXZ(3,i,k)
          node_fine_topXZ(1,i,k) = node_fine_topXZ(3,i,k)
          node_fine_topXZ(2,i,k) = node_fine_topXZ(3,i,k)
       end do
    end do
    
    do k=-gridRatio+1,nzSub_fine+gridRatio+1
       do j=2,nySub_fine-1
          do m=0,NumDistDirs
             fCtoF_frontYZ(m,1,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the first time step
             fCtoF_frontYZ(m,2,j,k) = fCtoF_frontYZ(m,3,j,k) !Cycle the last time step to the second time step
             fCtoF_backYZ(m,1,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the first time step
             fCtoF_backYZ(m,2,j,k) = fCtoF_backYZ(m,3,j,k) !Cycle the last time step to the second time step
          end do
          dsCtoF_frontYZ(:,1,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the first time step
          dsCtoF_frontYZ(:,2,j,k) = dsCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step
          velCtoF_frontYZ(:,1,j,k) = velCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the first time step
          velCtoF_frontYZ(:,2,j,k) = velCtoF_frontYZ(:,3,j,k) !Cycle the last time step to the second time step

          dsCtoF_backYZ(:,1,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the first time step
          dsCtoF_backYZ(:,2,j,k) = dsCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step        
          velCtoF_backYZ(:,1,j,k) = velCtoF_backYZ(:,3,j,k) !Cycle the last time step to the first time step
          velCtoF_backYZ(:,2,j,k) = velCtoF_backYZ(:,3,j,k) !Cycle the last time step to the second time step        

          node_fine_frontYZ(1,j,k) = node_fine_frontYZ(3,j,k)
          node_fine_frontYZ(2,j,k) = node_fine_frontYZ(3,j,k)
          node_fine_backYZ(1,j,k) = node_fine_backYZ(3,j,k)
          node_fine_backYZ(2,j,k) = node_fine_backYZ(3,j,k)
       end do
    end do
    
  END SUBROUTINE InitializeAllTemporalInterpolation
  
  !================================================
END MODULE LBM_FINE
!================================================
