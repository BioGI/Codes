!==================================================================================================
PROGRAM LBM3D	! 3D Parallelized LBM Simulation
               ! Gino Banco (2008-20010)
!==================================================================================================
  USE SetPrecision 
  USE Setup
  USE Setup_fine
  USE Parallel    
  USE Parallel_fine   
  USE LBM      
  USE LBM_fine      
  USE Geometry
  USE Geometry_fine
  USE PassiveScalar
  USE PassiveScalar_fine
  USE ICBC	
  USE ICBC_fine	
  USE Output
  USE Output_fine
  
  IMPLICIT NONE
  
  INTEGER(lng) :: mpierr					! MPI standard error variable
  INTEGER :: i,j,k,m
  REAL(dbl) :: feq
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MPI Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  CALL MPI_INIT(mpierr)					! initialize parallelization [Intrinsic]
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,mpierr)	! get the size of the parallel "world" (number of processing units) [Intrinsic]
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,mpierr)		! assign each prossessing unit a number (myid) for identification [Intrinsic]
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  CALL Global_Setup_fine					! set up the simulation {MODULE: Setup]
  CALL Global_Setup					! set up the simulation {MODULE: Setup]
  
  CALL MPI_Setup			! set up MPI component of the simulation [MODULE: Parallel]
  CALL MPI_Setup_Fine		! set up MPI component of the simulation [MODULE: Parallel_fine]
  CALL LBM_Setup			! set up LBM simulation [MODULE: LBM]
  CALL LBM_Setup_Fine		! set up LBM simulation [MODULE: LBM_fine]
  CALL Geometry_Setup		! set up the geometry of the physical simulation [MODULE: Geometry]
  CALL Geometry_Setup_Fine	! set up the geometry of the fine mesh in the physical simulation [MODULE: Geometry_fine]
  CALL Output_Setup		! set up the output [MODULE: Output]
  CALL Output_Setup_fine		! set up the output [MODULE: Output_fine]
  CALL OpenOutputFiles		! opens output files for writing [MODULE: Output.f90]
  CALL OpenOutputFiles_fine	! opens output files for writing [MODULE: Output_fine.f90]
  CALL Scalar_Setup		! set up the passive scalar component of the simluation [MODULE: Scalar]
  CALL Scalar_Setup_fine		! set up the passive scalar component of the simluation [MODULE: Scalar_fine]
  CALL FlagCoarseMeshNodesIntersectingWithFineMeshNodes
  CALL ICs			! set initial conditions [MODULE: ICBC]
  CALL ICs_fine			! set initial conditions [MODULE: ICBC_fine]       
  
  ! Setup interpolation 
  iter = 0  
  write(*,*) 'Using iter=0 for initializing SetNodesInterface_nPlus1_fine. Has to be modified for restart'
  CALL SetNodesInterface_nPlus1_fine               ! Calculate the node values on the fine mesh boundary at the next time step for temporal interpolation
  CALL ComputeEquilibriumForFineGrid               ! Compute the equilibrium distribution function at the coarse grid interface for the fine grid 
  CALL XYSpatialInterpolateBufferToFineGrid        ! Do the XY spatial interpolation on the buffer nodes for required variables to fine grid
  CALL PackAndSendDataBufferInterpolation          ! Send the data on the buffer nodes
  CALL XYSpatialInterpolateInternalNodesToFineGrid ! Do the XY spatial interpolation on the internal nodes for required variables to fine grid
  CALL ReceiveAndUnpackDataBufferInterpolation     ! Receive the buffer data
  CALL ZSpatialInterpolateToFineGrid               ! Do the Z spatial interpolation for required variables to fine grid
  CALL InitializeAllTemporalInterpolation          ! To begin with set all the 3 values for temporal interpolation to the latest available value
  
  iter=0
  subIter=0
  CALL temporalInterpolateToFineGrid !Using the spatial interpolation at the three time points, n-1, n and n+1, perform temporal interpolation to the current sub Iteration                   
  
  CALL PrintParams		! print simulation info [MODULE: Output]
  CALL PrintFields		! output the velocity, density, and scalar fields [MODULE: Output]
  CALL PrintParams_fine		! print simulation info [MODULE: Output]
  CALL PrintFields_fine		! output the velocity, density, and scalar fields [MODULE: Output]
  CALL PrintStatus		! Start simulation timer, print status [MODULE: Output]
  
  IF(ParticleTrack.EQ.ParticleOn) THEN 				! If particle tracking is 'on' then do the following
     CALL IniParticles
     CALL Particle_Setup
     CALL Particle_Setup_fine           
  ENDIF
  
  IF(restart) THEN		! calculate the villous locations/angles at iter0-1 [MODULE: Geometry]
     write(*,*) "Cannot restart yet. Implementation not complete"
     !CALL AdvanceGeometry
  END IF
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processes before starting simulation [Intrinsic]
  
  
  !	CALL PrintTime 				! print time (scalability) information [MODULE: Output]
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  DO iter=iter0,nt
     
     CALL AdvanceGeometry		! advance the geometry to the next time step [MODULE: Geometry]
     !         fPlus = f
     CALL Stream			! perform the streaming operation (with Lallemand 2nd order BB) [MODULE: Algorithm]
     CALL Macro			! calcuate the macroscopic quantities [MODULE: Algorithm]
     
     IF(ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN 	! If particle tracking is 'on' then do the following
        !	   CALL Calc_Global_Bulk_Scalar_Conc				! Estimate bluk	scalar concentration in each partition
        !	   CALL Collect_Distribute_Global_Bulk_Scalar_Conc		! Collect Cb_Global from different processors, average it and distribute it to all the processors.
        
        CALL Particle_Track
        !	   CALL Particle_MPI_Transfer
     ENDIF
     
     IF(iter .GE. phiStart) THEN
        CALL Scalar		! calcuate the evolution of scalar in the domain [MODULE: Algorithm]
     END IF
     phi_fine = phi_fine + delphi_particle_fine !Add the drug release corresponding to any particle in the coarse mesh whose effective volume interfaces with the fine mesh
     
     CALL Collision			! collision step [MODULE: Algorithm]
     CALL MPI_Transfer		! transfer the data (distribution functions, density, scalar) [MODULE: Parallel]
     
     CALL SetNodesInterface_nPlus1_fine               ! Calculate the node values on the fine mesh boundary at the next time step for temporal interpolation
     CALL ComputeEquilibriumForFineGrid               ! Compute the equilibrium distribution function at the coarse grid interface for the fine grid 
     CALL XYSpatialInterpolateBufferToFineGrid        ! Do the XY spatial interpolation on the buffer nodes for required variables to fine grid
     CALL PackAndSendDataBufferInterpolation          ! Send the data on the buffer nodes
     CALL XYSpatialInterpolateInternalNodesToFineGrid ! Do the XY spatial interpolation on the internal nodes for required variables to fine grid
     CALL ReceiveAndUnpackDataBufferInterpolation     ! Receive the buffer data
     CALL ZSpatialInterpolateToFineGrid               ! Do the Z spatial interpolation for required variables to fine grid
     
     DO subiter=1,gridRatio
        CALL AdvanceGeometry_Fine   ! Advance the geometry on the fine grid
        CALL temporalInterpolateToFineGrid !Using the spatial interpolation at the three time points, n-1, n and n+1, perform temporal interpolation to the current sub Iteration
        !           fPlus_fine = f_fine
        CALL Stream_Fine            ! Stream fine grid
        CALL Macro_Fine             ! Calculate Macro properties on fine grid
        
        IF(ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN 	! If particle tracking is 'on' then do the following
           CALL Particle_Track_fine
           !              CALL Particle_MPI_Transfer !The first time this is called, it should technically do the work for any particles in the coarse mesh that have to be transferred across processors as well.
        ENDIF
        CALL Scalar_Fine       ! Calculate Scalar stuff on fine grid
        phi = phi + delphi_particle
        
        CALL Collision_Fine     ! Collision step on the fine grid
        CALL MPI_Transfer_Fine  ! Transfer the data across processor boundaries on the fine grid
        
     END DO
     
     CALL PrintFields			! output the velocity, density, and scalar fields [MODULE: Output]
     CALL PrintFields_fine		! output the velocity, density, and scalar fields [MODULE: Output]                       
     CALL ComputeEquilibriumForCoarseGrid ! Compute the equilibrium distribution function at the fine grid interface for the coarse grid 
     CALL InterpolateToCoarseGrid    ! Interpolate required variable to coarse grid
     
     IF(ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN 	! If particle tracking is 'on' then do the following
        CALL PrintParticles						! output the particle velocity, radius, position and con. [MODULE: Output]
     ENDIF
     
     CALL PrintDrugConservation		! print the total absorbed/entering/leaving scalar as a function of time [MODULE: Output]
     CALL PrintMass			! print the total mass in the system (TEST)
     
     !	  CALL PrintPeriodicRestart	! print periodic restart files (SAFE GUARD) [MODULE: Output]
     
     CALL PrintStatus 		! print current status [MODULE: Output]
     
     CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processing units before next loop [Intrinsic]
     
  END DO
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Simulation Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !	CALL PrintTime 			! print time (scalability) information [MODULE: Output]
  
  CALL PrintFinalRestart		! print a final set of restart files to continue if desired [MODULE: Output]
  
  CALL DEAllocateArrays		! clean up the memory [MODULE: Setup]
  CALL DEAllocateArrays_fine		! clean up the memory [MODULE: Setup_fine]
  
  CALL CloseOutputFiles		! closes output files [MODULE: Output.f90]
  CALL CloseOutputFiles_fine		! closes output files [MODULE: Output_fine.f90]
  
  !   CALL MergeOutput		! combine the subdomain output into an output file for the entire computational domain [MODULE: Output]
  
  CALL MPI_FINALIZE(mpierr)	! end the MPI simulation [Intrinsic]
  
  !================================================ 
END PROGRAM LBM3D
!================================================
