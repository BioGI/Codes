!==================================================================================================
MODULE Parallel_fine		! Defines Parallel (MPI) Variables
							! Contains Parallel (MPI) Subroutines (MPI_Sub_Info, MPI_Initialize, FillSendArrays, MPI_Transfer, SendData, RecvData)
							! Written by Yanxing Wang (2008)
			 				! Modified by Gino Banco (2008-2009)
!================================================================================================== 
USE SetPrecision
USE Setup
USE Setup_fine
USE Setup
USE LBM
USE LBM_fine
USE PassiveScalar
USE MPI					! Intrinsic MPI definitions module

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE MPI_Setup_fine	! setup the MPI (parallel) component
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE


! Initialize variables and arrays
Corner_SendIndex_Fine	= 0_lng		! i, j, and k indices for each corner
Corner_RecvIndex_Fine	= 0_lng		! i, j, and k indices for each corner (phantom node for recieving data)
Z_SendIndex_Fine			= 0_lng		! i and j indices for each Z side 
Z_RecvIndex_Fine			= 0_lng		! i and j indices for each Z side (phantom node for recieving data)
X_SendIndex_Fine			= 0_lng		! j and k indices for each X side 
X_RecvIndex_Fine			= 0_lng		! j and k indices for each X side (phantom node for recieving data)
Y_SendIndex_Fine			= 0_lng		! i and k indices for each Y side 
Y_RecvIndex_Fine			= 0_lng		! i and k indices for each Y side (phantom node for recieving data)
YZ_SendIndex_Fine		= 0_lng		! i index for each YZ face 
YZ_RecvIndex_Fine		= 0_lng		! i index for each YZ face (phantom node for recieving data)
ZX_SendIndex_Fine		= 0_lng		! j index for each ZX face 
ZX_RecvIndex_Fine		= 0_lng		! j index for each ZX face (phantom node for recieving data)
XY_SendIndex_Fine		= 0_lng		! k index for each XY face 
XY_RecvIndex_Fine		= 0_lng		! k index for each XY face (phantom node for recieving data)
! OppCommDir - Same as coarse mesh
CommDataStart_f_Fine	= 0_lng		! array of starting indices in the send arrays for the distribution functions from each communication direction 
CommDataStart_rho_Fine	= 0_lng		! array of starting indices in the send arrays for the density from each communication direction
CommDataStart_phi_Fine	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_u_Fine	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_v_Fine	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_w_Fine	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
fSize_Fine				= 0_lng		! array of the number of elements sent for each communication direction (distribution functions)
dsSize_Fine				= 0_lng		! array of the number of elements sent for each communication direction (density and scalar)
uvwSize_Fine				= 0_lng		! array of the number of elements sent for each communication direction (density and scalar)

! Fill out the MPI arrays
CALL MPI_Initialize_Fine

!------------------------------------------------
END SUBROUTINE MPI_Setup_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE MPI_Initialize_Fine	! initialize the MPI arrays and variables
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: iComm												! index variable
INTEGER(lng) :: YZ_FaceSize_fine, ZX_FaceSize_fine, XY_FaceSize_fine		! number of nodes on a subdomain face oriented on the respective faces
INTEGER(lng) :: f_SendSize_fine, ds_SendSize_fine, uvw_SendSize_fine, total_SendSize_fine	! sizes of the distribution function, density, velocity, and scalar data transfer arrays respectively
!INTEGER(lng) :: msgSze												! size of the message send/recv arrays (maximum size)

! OppCommDir - Same as coarse mesh

! Distribution function components transferred
! f_Comps - Same as coarse mesh


! Fill out the size arrays
YZ_FaceSize_fine		= nySub_fine*nzSub_fine	 ! number of nodes on a subdomain face oriented in the ZY plane
ZX_FaceSize_fine		= nzSub_fine*nxSub_fine	 ! number of nodes on a subdomain face oriented in the ZX plane
XY_FaceSize_fine		= nxSub_fine*nySub_fine	 ! number of nodes on a subdomain face oriented in the XY plane

fSize_fine		= 0_lng			! initialize array
dsSize_fine		= 0_lng			! initialize array
uvwSize_fine		= 0_lng			! initialize array

fSize_fine(1:2)		= YZ_FaceSize_fine*NumFs_face	! YZ faces
fSize_fine(3:4)		= ZX_FaceSize_fine*NumFs_face			! ZX faces
fSize_fine(5:6)		= XY_FaceSize_fine*NumFs_face			! XY faces
fSize_fine(7:10)		= nzSub_fine*NumFs_side					! Z sides
fSize_fine(11:14)	= nxSub_fine*NumFs_side					! X sides
fSize_fine(15:18)	= nySub_fine*NumFs_side					! Y sides
fSize_fine(19:26)	= 1_lng*NumFs_corner					! corners

dsSize_fine(1:2)		= YZ_FaceSize_fine							! YZ faces
dsSize_fine(3:4)		= ZX_FaceSize_fine							! ZX faces
dsSize_fine(5:6)		= XY_FaceSize_fine							! XY faces
dsSize_fine(7:10)	= nzSub_fine									! Z sides
dsSize_fine(11:14)	= nxSub_fine									! X sides
dsSize_fine(15:18)	= nySub_fine									! Y sides
dsSize_fine(19:26)	= 1_lng									! corners


uvwSize_fine(1:2)		= YZ_FaceSize_fine							! YZ faces
uvwSize_fine(3:4)		= ZX_FaceSize_fine							! ZX faces
uvwSize_fine(5:6)		= XY_FaceSize_fine							! XY faces
uvwSize_fine(7:10)	= nzSub_fine									! Z sides
uvwSize_fine(11:14)	= nxSub_fine									! X sides
uvwSize_fine(15:18)	= nySub_fine									! Y sides
uvwSize_fine(19:26)	= 1_lng									! corners

msgSize_fine(:)		= fSize_fine(:) + 2_lng*(dsSize_fine(:))+3_lng*(uvwSize_fine(:))	! total message sizes

f_SendSize_fine	= SUM(fSize_fine)
ds_SendSize_fine	= SUM(dsSize_fine)
uvw_SendSize_fine	= SUM(dsSize_fine)
total_SendSize_fine  = SUM(msgSize_fine)

ALLOCATE(msgSend_fine(total_SendSize_fine))						
ALLOCATE(msgRecv_fine(total_SendSize_fine))

! Fill out the '3D_Index' arrays (for converting back from 1D array)
! Faces
YZ_SendIndex_fine(1)	= nxSub_fine
YZ_RecvIndex_fine(1)	= nxSub_fine+1_lng		

YZ_SendIndex_fine(2)	= 1_lng	
YZ_RecvIndex_fine(2)	= 0_lng	

ZX_SendIndex_fine(3)	= nySub_fine		
ZX_RecvIndex_fine(3)	= nySub_fine+1_lng		

ZX_SendIndex_fine(4)	= 1_lng	
ZX_RecvIndex_fine(4)	= 0_lng	

XY_SendIndex_fine(5)	= nzSub_fine
XY_RecvIndex_fine(5)	= nzSub_fine+1_lng		
		
XY_SendIndex_fine(6)	= 1_lng
XY_RecvIndex_fine(6)	= 0_lng	


! Sides
! Z Sides
Z_SendIndex_fine(7,1)	= nxSub_fine
Z_RecvIndex_fine(7,1)	= nxSub_fine+1_lng

Z_SendIndex_fine(7,2)	= nySub_fine
Z_RecvIndex_fine(7,2)	= nySub_fine+1_lng


Z_SendIndex_fine(8,1)	= 1_lng
Z_RecvIndex_fine(8,1)	= 0_lng

Z_SendIndex_fine(8,2)	= 1_lng
Z_RecvIndex_fine(8,2)	= 0_lng


Z_SendIndex_fine(9,1)	= nxSub_fine
Z_RecvIndex_fine(9,1)	= nxSub_fine+1_lng

Z_SendIndex_fine(9,2)	= 1_lng
Z_RecvIndex_fine(9,2)	= 0_lng


Z_SendIndex_fine(10,1)	= 1_lng
Z_RecvIndex_fine(10,1)	= 0_lng

Z_SendIndex_fine(10,2)	= nySub_fine
Z_RecvIndex_fine(10,2)	= nySub_fine+1_lng


! X Sides
X_SendIndex_fine(11,1)	= nySub_fine
X_RecvIndex_fine(11,1)	= nySub_fine+1_lng

X_SendIndex_fine(11,2)	= nzSub_fine
X_RecvIndex_fine(11,2)	= nzSub_fine+1_lng


X_SendIndex_fine(12,1)	= 1_lng
X_RecvIndex_fine(12,1)	= 0_lng

X_SendIndex_fine(12,2)	= 1_lng
X_RecvIndex_fine(12,2)	= 0_lng


X_SendIndex_fine(13,1)	= nySub_fine
X_RecvIndex_fine(13,1)	= nySub_fine+1_lng

X_SendIndex_fine(13,2)	= 1_lng
X_RecvIndex_fine(13,2)	= 0_lng


X_SendIndex_fine(14,1)	= 1_lng
X_RecvIndex_fine(14,1)	= 0_lng

X_SendIndex_fine(14,2)	= nzSub_fine
X_RecvIndex_fine(14,2)	= nzSub_fine+1_lng


! Y Sides
Y_SendIndex_fine(15,1)	= nxSub_fine
Y_RecvIndex_fine(15,1)	= nxSub_fine+1_lng

Y_SendIndex_fine(15,2)	= nzSub_fine
Y_RecvIndex_fine(15,2)	= nzSub_fine+1_lng


Y_SendIndex_fine(16,1)	= 1_lng
Y_RecvIndex_fine(16,1)	= 0_lng

Y_SendIndex_fine(16,2)	= 1_lng
Y_RecvIndex_fine(16,2)	= 0_lng


Y_SendIndex_fine(17,1)	= 1
Y_RecvIndex_fine(17,1)	= 0_lng

Y_SendIndex_fine(17,2)	= nzSub_fine
Y_RecvIndex_fine(17,2)	= nzSub_fine+1_lng


Y_SendIndex_fine(18,1)	= nxSub_fine
Y_RecvIndex_fine(18,1)	= nxSub_fine+1_lng

Y_SendIndex_fine(18,2)	= 1_lng
Y_RecvIndex_fine(18,2)	= 0_lng


! Corners
Corner_SendIndex_fine(19,1) = nxSub_fine
Corner_SendIndex_fine(19,2) = nySub_fine
Corner_SendIndex_fine(19,3) = nzSub_fine

Corner_RecvIndex_fine(19,1) = nxSub_fine+1_lng
Corner_RecvIndex_fine(19,2) = nySub_fine+1_lng
Corner_RecvIndex_fine(19,3) = nzSub_fine+1_lng


Corner_SendIndex_fine(20,1) = 1_lng
Corner_SendIndex_fine(20,2) = 1_lng
Corner_SendIndex_fine(20,3) = 1_lng

Corner_RecvIndex_fine(20,1) = 0_lng
Corner_RecvIndex_fine(20,2) = 0_lng
Corner_RecvIndex_fine(20,3) = 0_lng


Corner_SendIndex_fine(21,1) = nxSub_fine
Corner_SendIndex_fine(21,2) = nySub_fine
Corner_SendIndex_fine(21,3) = 1_lng

Corner_RecvIndex_fine(21,1) = nxSub_fine+1_lng
Corner_RecvIndex_fine(21,2) = nySub_fine+1_lng
Corner_RecvIndex_fine(21,3) = 0_lng


Corner_SendIndex_fine(22,1) = 1_lng
Corner_SendIndex_fine(22,2) = 1_lng
Corner_SendIndex_fine(22,3) = nzSub_fine

Corner_RecvIndex_fine(22,1) = 0_lng
Corner_RecvIndex_fine(22,2) = 0_lng
Corner_RecvIndex_fine(22,3) = nzSub_fine+1_lng


Corner_SendIndex_fine(23,1) = 1_lng
Corner_SendIndex_fine(23,2) = nySub_fine
Corner_SendIndex_fine(23,3) = nzSub_fine

Corner_RecvIndex_fine(23,1) = 0_lng
Corner_RecvIndex_fine(23,2) = nySub_fine+1_lng
Corner_RecvIndex_fine(23,3) = nzSub_fine+1_lng


Corner_SendIndex_fine(24,1) = nxSub_fine
Corner_SendIndex_fine(24,2) = 1_lng
Corner_SendIndex_fine(24,3) = 1_lng

Corner_RecvIndex_fine(24,1) = nxSub_fine+1_lng
Corner_RecvIndex_fine(24,2) = 0_lng
Corner_RecvIndex_fine(24,3) = 0_lng


Corner_SendIndex_fine(25,1) = nxSub_fine
Corner_SendIndex_fine(25,2) = 1_lng
Corner_SendIndex_fine(25,3) = nzSub_fine

Corner_RecvIndex_fine(25,1) = nxSub_fine+1_lng
Corner_RecvIndex_fine(25,2) = 0_lng
Corner_RecvIndex_fine(25,3) = nzSub_fine+1_lng


Corner_SendIndex_fine(26,1) = 1_lng
Corner_SendIndex_fine(26,2) = nySub_fine
Corner_SendIndex_fine(26,3) = 1_lng

Corner_RecvIndex_fine(26,1) = 0_lng
Corner_RecvIndex_fine(26,2) = nySub_fine+1_lng
Corner_RecvIndex_fine(26,3) = 0_lng


! Fill out the 'CommDataStart' arrays
! Initialize arrays
CommDataStart_f_fine  		= 0_lng					! distribution functions
CommDataStart_rho_fine  	= 0_lng					! density
CommDataStart_phi_fine  	= 0_lng					! scalar
CommDataStart_u_fine  	= 0_lng					! velocity
CommDataStart_v_fine  	= 0_lng					! velocity
CommDataStart_w_fine  	= 0_lng					! velocity

CommDataStart_f_fine(1)	= 1_lng	
CommDataStart_rho_fine(1) = CommDataStart_f_fine(1)	+ fSize_fine(1)
CommDataStart_phi_fine(1) = CommDataStart_rho_fine(1)	+ dsSize_fine(1)
CommDataStart_u_fine(1) = CommDataStart_phi_fine(1)	+ dsSize_fine(1)
CommDataStart_v_fine(1) = CommDataStart_u_fine(1)		+ uvwSize_fine(1)
CommDataStart_w_fine(1) = CommDataStart_v_fine(1)		+ uvwSize_fine(1)

DO iComm=2,NumCommDirs							! fill out for communication directions 2-NumCommDirs
  CommDataStart_f_fine(iComm)	= CommDataStart_w_fine(iComm-1) 	+ uvwSize_fine(iComm-1)
  CommDataStart_rho_fine(iComm)	= CommDataStart_f_fine(iComm)	+ fSize_fine(iComm)
  CommDataStart_phi_fine(iComm)	= CommDataStart_rho_fine(iComm) 	+ dsSize_fine(iComm)
  CommDataStart_u_fine(iComm)	= CommDataStart_phi_fine(iComm) 	+ dsSize_fine(iComm)
  CommDataStart_v_fine(iComm)	= CommDataStart_u_fine(iComm) 	+ uvwSize_fine(iComm)
  CommDataStart_w_fine(iComm)	= CommDataStart_v_fine(iComm) 	+ uvwSize_fine(iComm)
END DO

!WRITE(6678,*) 'f_SendSize_fine', total_SendSize_fine
!WRITE(6678,*) 'ds_SendSize_fine', total_SendSize_fine
!WRITE(6678,*) 'total_SendSize_fine', total_SendSize_fine
!DO iComm=1,NumCommDirs
!  WRITE(6678,*) 'CommDataStart_f_fine(iComm)', CommDataStart_f_fine(iComm)
!  WRITE(6678,*) 'CommDataStart_rho_fine(iComm)', CommDataStart_rho_fine(iComm) 
!  WRITE(6678,*) 'CommDataStart_phi_fine(iComm)', CommDataStart_phi_fine(iComm) 
!END DO
!STOP

! Allocate the MPI_WAITALL status array
ALLOCATE(waitStat_fine(MPI_STATUS_SIZE,2*NumCommDirs))

ALLOCATE(waitStatInterpolationBuffer(MPI_STATUS_SIZE,16))

!------------------------------------------------
END SUBROUTINE MPI_Initialize_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE PackData_Fine	! transfer the distribution functions between neighboring initialize the MPI arrays and variables
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: i,j,k,m,iComm,ii	! index variables

! Initialize array
msgSend_fine		= 0.0_dbl

! Fill Arrays
! FACES
! YZ Faces
DO iComm=1,2

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)

    i = YZ_SendIndex_fine(iComm)								! i index for assigning the proper component of the distribution function to the 1D send array
  
    DO k=1,nzSub_fine
      DO j=1,nySub_fine

        ! density and scalar
        ii = CommDataStart_rho_fine(iComm)			&		! start location
           + (j - 1_lng)							&		! convert 3D-coordinate into proper 1D-array coordinate
           + (k - 1_lng)*nySub_fine							! convert 3D-coordinate into proper 1D-array coordinate

        msgSend_fine(ii) 						= rho_fine(i,j,k)	! store the proper density in the send array
        msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)	! store the proper scalar quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array

        ! distribution functions
        DO m=1,NumFs_face
     
          ii = CommDataStart_f_fine(iComm)			&		! start location
             + (m - 1_lng)							&		! current distribution function location
             + (j - 1_lng)*NumFs_face			&		! convert 3D-coordinate into proper 1D-array coordinate
             + (k - 1_lng)*NumFs_face*nySub_fine			! convert 3D-coordinate into proper 1D-array coordinate

          msgSend_fine(ii) = f_fine(f_Comps(iComm,m),i,j,k)	! store the proper component of the distribution function in the send array

        END DO   

      END DO
    END DO

  END IF

END DO

! ZX Faces
DO iComm=3,4

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)

    j = ZX_SendIndex_fine(iComm)								! j index for assigning the proper component of the distribution function to the 1D send array
  
    DO i=1,nxSub_fine
      DO k=1,nzSub_fine

        ! density and scalar
        ii = CommDataStart_rho_fine(iComm)			&		! start location
           + (k - 1_lng)							&		! convert 3D-coordinate into proper 1D-array coordinate
           + (i - 1_lng)*nzSub_fine							! convert 3D-coordinate into proper 1D-array coordinate
    
        msgSend_fine(ii) 						= rho_fine(i,j,k)	! store the proper density in the send array
        msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)	! store the proper scalar quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array

        ! distribution functions
        DO m=1,NumFs_face
    
          ii = CommDataStart_f_fine(iComm)			&		! start location
             + (m - 1_lng)							&		! current distribution function location
             + (k - 1_lng)*NumFs_face			&		! convert 3D-coordinate into proper 1D-array coordinate
             + (i - 1_lng)*NumFs_face*nzSub_fine			! convert 3D-coordinate into proper 1D-array coordinate

          msgSend_fine(ii) = f_fine(f_Comps(iComm,m),i,j,k)	! store the proper component of the distribution function in the send array

        END DO      

      END DO
    END DO

  END IF  

END DO

! XY Faces
DO iComm=5,6

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)  

    k = XY_SendIndex_fine(iComm)								! j index for assigning the proper component of the distribution function to the 1D send array
  
    DO j=1,nySub_fine
      DO i=1,nxSub_fine

        ! density and scalar
        ii = CommDataStart_rho_fine(iComm)			&		! start location
           + (i - 1_lng)							&		! convert 3D-coordinate into proper 1D-array coordinate
           + (j - 1_lng)*nxSub_fine							! convert 3D-coordinate into proper 1D-array coordinate     
    
        msgSend_fine(ii) 						= rho_fine(i,j,k)	! store the proper density in the send array
        msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)	! store the proper scalar quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array
      
        DO m=1,NumFs_face
      
          ii = CommDataStart_f_fine(iComm)		&			! start location
             + (m - 1_lng)						&			! current distribution function location
             + (i - 1_lng)*NumFs_face		&			! convert 3D-coordinate into proper 1D-array coordinate
             + (j - 1_lng)*NumFs_face*nxSub_fine			! convert 3D-coordinate into proper 1D-array coordinate

          msgSend_fine(ii) = f_fine(f_Comps(iComm,m),i,j,k)	! store the proper component of the distribution function in the send array
  
        END DO
      
      END DO
    END DO
  
  END IF

END DO


! SIDES
! Z Sides
DO iComm=7,10

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)    

    i = Z_SendIndex_fine(iComm,1)								! i index for assigning the proper component of the distribution function to the 1D send array        
    j = Z_SendIndex_fine(iComm,2)								! j index for assigning the proper component of the distribution function to the 1D send array
  
    DO k=1,nzSub_fine

      	! density and scalar
	ii = CommDataStart_rho_fine(iComm)			&			! start location
           + (k - 1_lng)										! convert 3D-coordinate into proper 1D-array coordinate
    
	msgSend_fine(ii) 					= rho_fine(i,j,k)		! store the proper density in the send array
	msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)		! store the proper scalar quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array
    
      DO m=1,NumFs_side
  
        ! distribution functions   
        ii = CommDataStart_f_fine(iComm)			&			! start location
           + (m - 1_lng)						&			! current distribution function location
           + (k - 1_lng)*NumFs_side						! convert 3D-coordinate into proper 1D-array coordinate

        msgSend_fine(ii) = f_fine(f_Comps(iComm,m),i,j,k)		! store the proper component of the distribution function in the send array
 
      END DO
  
    END DO

  END IF  
  
END DO

! X Sides
DO iComm=11,14

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)     

    j = X_SendIndex_fine(iComm,1)								! j index for assigning the proper component of the distribution function to the 1D send array        
    k = X_SendIndex_fine(iComm,2)								! k index for assigning the proper component of the distribution function to the 1D send array
  
    DO i=1,nxSub_fine

      ! density and scalar
      ii = CommDataStart_rho_fine(iComm)			&			! start location
         + (i - 1_lng)										! convert 3D-coordinate into proper 1D-array coordinate
    
      msgSend_fine(ii) 					= rho_fine(i,j,k)		! store the proper density in the send array
      msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)		! store the proper scalar quantity in the send array 
      msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
      msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
      msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array
    
      DO m=1,NumFs_side

        ! distribution functions      
        ii = CommDataStart_f_fine(iComm)			&			! start location
           + (m - 1_lng)						&			! current distribution function location
           + (i - 1_lng)*NumFs_side						! convert 3D-coordinate into proper 1D-array coordinate

        msgSend_fine(ii) = f_fine(f_Comps(iComm,m),i,j,k)		! store the proper component of the distribution function in the send array

      END DO
    
    END DO

  END IF
  
END DO

! Y Sides
DO iComm=15,18

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)   

    i = Y_SendIndex_fine(iComm,1)								! i index for assigning the proper component of the distribution function to the 1D send array        
    k = Y_SendIndex_fine(iComm,2)								! k index for assigning the proper component of the distribution function to the 1D send array
  
    DO j=1,nySub_fine

      ! density and scalar
      ii = CommDataStart_rho_fine(iComm)			&			! start location
         + (j - 1_lng)										! convert 3D-coordinate into proper 1D-array coordinate
    
        msgSend_fine(ii) 						= rho_fine(i,j,k)	! store the proper density in the send array
        msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)	! store the proper scalar quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
        msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array
    
      DO m=1,NumFs_side

        ! distribution functions
        ii = CommDataStart_f_fine(iComm)			&			! start location
           + (m - 1_lng)						&			! current distribution function locationcurrent distribution function location
           + (j - 1_lng)*NumFs_side						! convert 3D-coordinate into proper 1D-array coordinate

        msgSend_fine(ii) = f_fine(f_Comps(iComm,m),i,j,k)		! store the proper component of the distribution function in the send array

      END DO
     
    END DO
  
  END IF

END DO


! CORNERS
DO iComm=19,26
  
  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)   
  
    i = Corner_SendIndex_fine(iComm,1)						! j index for assigning the proper component of the distribution function to the 1D send array 	
    j = Corner_SendIndex_fine(iComm,2)						! j index for assigning the proper component of the distribution function to the 1D send array        
    k = Corner_SendIndex_fine(iComm,3)						! k index for assigning the proper component of the distribution function to the 1D send array

    ii = CommDataStart_rho_fine(iComm)						! start location
    msgSend_fine(ii) 					= rho_fine(i,j,k)			! store the proper density in the send array
    msgSend_fine(ii+dsSize_fine(iComm))	= phi_fine(i,j,k)			! store the proper scalar quantity in the send array
    msgSend_fine(ii+2_lng*dsSize_fine(iComm))	= u_fine(i,j,k)	! store the proper velocity quantity in the send array
    msgSend_fine(ii+2_lng*dsSize_fine(iComm)+uvwSize_fine(icomm))= v_fine(i,j,k)	! store the proper velocity quantity in the send array
    msgSend_fine(ii+2_lng*dsSize_fine(iComm)+2_lng*uvwSize_fine(icomm))= w_fine(i,j,k)	! store the proper velocity quantity in the send array
  
    ii = CommDataStart_f_fine(iComm)							! start location
    msgSend_fine(ii) = f_fine(f_Comps(iComm,1),i,j,k)			! store the proper component of the distribution function in the send array

  END IF
  	
END DO

!------------------------------------------------
END SUBROUTINE PackData_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE MPI_Transfer_Fine	! transfer the distribution functions between neighboring initialize the MPI arrays and variables
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: iComm											! index variable
INTEGER(lng) :: numReqs											! number of send/recv requests
INTEGER(lng) :: mpierr											! MPI standard error object

CALL PackData_Fine														! fill out the arrays to be transferred 

numReqs = 0_lng

! Post the receives
DO iComm = 1,NumCommDirs

  IF(SubID(OppCommDir(iComm)) .NE. 0) THEN 
    numReqs = numReqs + 1_lng   
    CALL PostRecv_Fine(iComm,numReqs)										! receive data
  END IF

END DO

! Send the data
DO iComm = 1,NumCommDirs

  IF(SubID(iComm) .NE. 0) THEN
    numReqs = numReqs + 1_lng 
    CALL SendData_Fine(iComm,numReqs)						  				! send data
  END IF

END DO

CALL MPI_WAITALL(numReqs,req(1:numReqs),waitStat_fine,mpierr)

! Store the Data
DO iComm = 1,NumCommDirs

  IF(SubID(OppCommDir(iComm)) .NE. 0) THEN    
    CALL UnPackData_Fine(iComm)										! store the sent data in the proper location at the destination
  END IF

END DO

!------------------------------------------------
END SUBROUTINE MPI_Transfer_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE PostRecv_Fine(iComm,numReqs)	! receives information from a neighboring subdomain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: iComm								! communication direction from which to receive
INTEGER(lng), INTENT(IN) :: numReqs								! number of send/recv requests
INTEGER(lng) :: msgStart,msgEnd									! start and end indices of the message
INTEGER(lng) :: src													! source processing unit
INTEGER(lng) :: tag													! message tag
INTEGER(lng) :: mpierr												! MPI standard error variable

! starting/ending indices of the message in the mgsRecv array
msgStart = CommDataStart_f_fine(OppCommDir(iComm))
msgEnd	= msgStart + msgSize_fine(OppCommDir(iComm))

src	= SubID(OppCommDir(iComm)) - 1_lng						! rank of processing unit sending message TO this processing unit
tag	= iComm + 100_lng												! message tag 

CALL MPI_IRECV(msgRecv_fine(msgStart:msgEnd),msgSize_fine(OppCommDir(iComm)),MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,req(numReqs),mpierr)		! receive data

!------------------------------------------------
END SUBROUTINE PostRecv_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SendData_Fine(iComm,numReqs)											! sends information to a neighboring subdomain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: iComm								! communication direction
INTEGER(lng), INTENT(IN) :: numReqs								! number of send/recv requests
INTEGER(lng) :: msgStart,msgEnd									! start and end indices of the distribution function data
INTEGER(lng) :: dest													! rank of destination processing unit
INTEGER(lng) :: tag													! message tag
INTEGER(lng) :: mpierr												! MPI standard error variable 

! starting/ending indices of the message in the mgsSend array
msgStart		= CommDataStart_f_fine(iComm)	
msgEnd		= msgStart + msgSize_fine(iComm)

dest 			= SubID(iComm) - 1_lng								! rank of processing unit receiving message from the current processing unit (-1 to correspond to rank (myid))
tag			= iComm + 100_lng										! message tag 

CALL MPI_ISEND(msgSend_fine(msgStart:msgEnd),msgSize_fine(iComm),MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,req(numReqs),mpierr)					! send data

!------------------------------------------------
END SUBROUTINE SendData_Fine
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE UnPackData_Fine(iComm)	! store the recieved data in the proper locations
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: iComm																				! communication direction
INTEGER(lng) :: i,j,k,m,ii																							! index variables
INTEGER(lng) :: stat(MPI_STATUS_SIZE)																			! MPI status object
INTEGER(lng) :: mpierr																								! MPI standard error variable

SELECT CASE(OppCommDir(iComm))
 
  CASE(1,2)			! YZ Faces

    i = YZ_RecvIndex_fine(OppCommDir(iComm))																		! i index for obtaining the proper information from the 1D transfer array 	

    DO k=1,nzSub_fine
      DO j=1,nySub_fine

          ii = j + (k-1)*nySub_fine																					! location of density and scalar function to recieve

          rho_fine(i,j,k) = msgRecv_fine((CommDataStart_rho_fine(OppCommDir(iComm))-1) + ii)						! store recieved density in proper place
          phi_fine(i,j,k) = msgRecv_fine((CommDataStart_phi_fine(OppCommDir(iComm))-1) + ii)						! store recieved scalar quantity in proper place          
          u_fine(i,j,k) = msgRecv_fine((CommDataStart_u_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          v_fine(i,j,k) = msgRecv_fine((CommDataStart_v_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          w_fine(i,j,k) = msgRecv_fine((CommDataStart_w_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
        
        DO m=1,NumFs_Face
          
          ii = m 							&																			! location of distribution function to recieve
             + (j-1)*NumFs_Face 		&
             + (k-1)*nySub_fine*NumFs_Face
             
          f_fine(f_Comps(iComm,m),i,j,k) = msgRecv_fine((CommDataStart_f_fine(OppCommDir(iComm))-1) + ii)	! store the recieved distribution function component in the proper place in the f array

        END DO
        
      END DO
    END DO

  CASE(3,4)			! ZX Faces

    j = ZX_RecvIndex_fine(OppCommDir(iComm))																		! j index for obtaining the proper information from the 1D transfer array 

    DO i=1,nxSub_fine
      DO k=1,nzSub_fine

          ii = k + (i-1)*nzSub_fine																					! location of density and scalar function to recieve

          rho_fine(i,j,k) = msgRecv_fine((CommDataStart_rho_fine(OppCommDir(iComm))-1) + ii)						! store recieved density in proper place
          phi_fine(i,j,k) = msgRecv_fine((CommDataStart_phi_fine(OppCommDir(iComm))-1) + ii)						! store recieved scalar quantity in proper place              
          u_fine(i,j,k) = msgRecv_fine((CommDataStart_u_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          v_fine(i,j,k) = msgRecv_fine((CommDataStart_v_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          w_fine(i,j,k) = msgRecv_fine((CommDataStart_w_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
        
        DO m=1,NumFs_Face
          
          ii = m 							&																			! location of distribution function to recieve
             + (k-1)*NumFs_Face 		&																			! location of density and scalar function to recieve
             + (i-1)*nzSub_fine*NumFs_Face
             
          f_fine(f_Comps(iComm,m),i,j,k) = msgRecv_fine((CommDataStart_f_fine(OppCommDir(iComm))-1) + ii)	! store the recieved distribution function component in the proper place in the f array

        END DO
        
      END DO
    END DO

  CASE(5,6)			! XY Faces

    k = XY_RecvIndex_fine(OppCommDir(iComm))																		! k index for obtaining the proper information from the 1D transfer array 

    DO j=1,nySub_fine
      DO i=1,nxSub_fine

        ii = i + (j-1)*nxSub_fine																						! location of density and scalar function to recieve
             	
          rho_fine(i,j,k) = msgRecv_fine((CommDataStart_rho_fine(OppCommDir(iComm))-1) + ii)						! store recieved density in proper place
          phi_fine(i,j,k) = msgRecv_fine((CommDataStart_phi_fine(OppCommDir(iComm))-1) + ii)						! store recieved scalar quantity in proper place         
          u_fine(i,j,k) = msgRecv_fine((CommDataStart_u_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          v_fine(i,j,k) = msgRecv_fine((CommDataStart_v_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          w_fine(i,j,k) = msgRecv_fine((CommDataStart_w_fine(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
        
        DO m=1,NumFs_Face
          
          ii = m 							&																			! location of distribution function to recieve
             + (i-1)*NumFs_Face 		&																			! location of density and scalar function to recieve
             + (j-1)*nxSub_fine*NumFs_Face
             
          f_fine(f_Comps(iComm,m),i,j,k) = msgRecv_fine((CommDataStart_f_fine(OppCommDir(iComm))-1) + ii)	! store the recieved distribution function component in the proper place in the f array
  
        END DO
        
      END DO
    END DO

  CASE(7:10)		! Z Sides

    i = Z_RecvIndex_fine(OppCommDir(iComm),1)																		! i index for obtaining the proper information from the 1D transfer array 	
    j = Z_RecvIndex_fine(OppCommDir(iComm),2)																		! j index for obtaining the proper information from the 1D transfer array 

    DO k=1,nzSub_fine

      rho_fine(i,j,k) = msgRecv_fine((CommDataStart_rho_fine(OppCommDir(iComm))-1) + k)							! store recieved density in proper place
      phi_fine(i,j,k) = msgRecv_fine((CommDataStart_phi_fine(OppCommDir(iComm))-1) + k)							! store recieved scalar quantity in proper place            
      u_fine(i,j,k) = msgRecv_fine((CommDataStart_u_fine(OppCommDir(iComm))-1) + k)						! store recieved velocity quantity in proper place          
      v_fine(i,j,k) = msgRecv_fine((CommDataStart_v_fine(OppCommDir(iComm))-1) + k)						! store recieved velocity quantity in proper place          
      w_fine(i,j,k) = msgRecv_fine((CommDataStart_w_fine(OppCommDir(iComm))-1) + k)						! store recieved velocity quantity in proper place          
        
      DO m=1,NumFs_Side
          
        ii = m + (k-1)*NumFs_Side																				! location of distribution function to recieve
             
        f_fine(f_Comps(iComm,m),i,j,k) = msgRecv_fine((CommDataStart_f_fine(OppCommDir(iComm))-1) + ii)		! store the recieved distribution function component in the proper place in the f array

      END DO

    END DO

  CASE(11:14)		! X Sides

    j = X_RecvIndex_fine(OppCommDir(iComm),1)																		! j index for obtaining the proper information from the 1D transfer array 	
    k = X_RecvIndex_fine(OppCommDir(iComm),2)																		! k index for obtaining the proper information from the 1D transfer array 

    DO i=1,nxSub_fine

      rho_fine(i,j,k) = msgRecv_fine((CommDataStart_rho_fine(OppCommDir(iComm))-1) + i)							! store recieved density in proper place
      phi_fine(i,j,k) = msgRecv_fine((CommDataStart_phi_fine(OppCommDir(iComm))-1) + i)							! store recieved scalar quantity in proper place  
      u_fine(i,j,k) = msgRecv_fine((CommDataStart_u_fine(OppCommDir(iComm))-1) + i)						! store recieved velocity quantity in proper place          
      v_fine(i,j,k) = msgRecv_fine((CommDataStart_v_fine(OppCommDir(iComm))-1) + i)						! store recieved velocity quantity in proper place          
      w_fine(i,j,k) = msgRecv_fine((CommDataStart_w_fine(OppCommDir(iComm))-1) + i)						! store recieved velocity quantity in proper place          
        
      DO m=1,NumFs_Side
          
        ii = m + (i-1)*NumFs_Side																				! location of distribution function to recieve
             
        f_fine(f_Comps(iComm,m),i,j,k) = msgRecv_fine((CommDataStart_f_fine(OppCommDir(iComm))-1) + ii)		! store the recieved distribution function component in the proper place in the f array

      END DO

    END DO

  CASE(15:18)		! Y Sides

    i = Y_RecvIndex_fine(OppCommDir(iComm),1)																		! i index for obtaining the proper information from the 1D transfer array 	
    k = Y_RecvIndex_fine(OppCommDir(iComm),2)																		! k index for obtaining the proper information from the 1D transfer array 

    DO j=1,nySub_fine

      rho_fine(i,j,k) = msgRecv_fine((CommDataStart_rho_fine(OppCommDir(iComm))-1) + j)							! store recieved density in proper place
      phi_fine(i,j,k) = msgRecv_fine((CommDataStart_phi_fine(OppCommDir(iComm))-1) + j)							! store recieved scalar quantity in proper place  
      u_fine(i,j,k) = msgRecv_fine((CommDataStart_u_fine(OppCommDir(iComm))-1) + j)						! store recieved velocity quantity in proper place          
      v_fine(i,j,k) = msgRecv_fine((CommDataStart_v_fine(OppCommDir(iComm))-1) + j)						! store recieved velocity quantity in proper place          
      w_fine(i,j,k) = msgRecv_fine((CommDataStart_w_fine(OppCommDir(iComm))-1) + j)						! store recieved velocity quantity in proper place          
        
      DO m=1,NumFs_Side
          
        ii = m + (j-1)*NumFs_Side																				! location of distribution function to recieve
             
        f_fine(f_Comps(iComm,m),i,j,k) = msgRecv_fine((CommDataStart_f_fine(OppCommDir(iComm))-1) + ii)		! store the recieved distribution function component in the proper place in the f array

      END DO

    END DO

  CASE(19:26)		! Corners

    i = Corner_RecvIndex_fine(OppCommDir(iComm),1) 																! i index for obtaining the proper information from the 1D transfer array 
    j = Corner_RecvIndex_fine(OppCommDir(iComm),2)															 	! j index for obtaining the proper information from the 1D transfer array 	
    k = Corner_RecvIndex_fine(OppCommDir(iComm),3)																! k index for obtaining the proper information from the 1D transfer array    

    rho_fine(i,j,k) = msgRecv_fine(CommDataStart_rho_fine(OppCommDir(iComm)))											! store recieved density in proper place
    phi_fine(i,j,k) = msgRecv_fine(CommDataStart_phi_fine(OppCommDir(iComm)))											! store recieved scalar quantity in proper place  
    u_fine(i,j,k) = msgRecv_fine(CommDataStart_u_fine(OppCommDir(iComm)))						! store recieved velocity quantity in proper place          
    v_fine(i,j,k) = msgRecv_fine(CommDataStart_v_fine(OppCommDir(iComm)))						! store recieved velocity quantity in proper place          
    w_fine(i,j,k) = msgRecv_fine(CommDataStart_w_fine(OppCommDir(iComm)))						! store recieved velocity quantity in proper place          
    f_fine(f_Comps(iComm,1),i,j,k) = msgRecv_fine(CommDataStart_f_fine(OppCommDir(iComm)))						! store the recieved distribution function component in the proper place in the f array

  CASE DEFAULT

    OPEN(1000,FILE="error.txt")
    WRITE(1000,*) "Error in UnPackData in Parallel.f90: iComm is not 1-26..."
    WRITE(1000,*) "iComm",iComm
    CLOSE(1000)
    STOP

END SELECT

!------------------------------------------------
END SUBROUTINE UnPackData_Fine
!------------------------------------------------

SUBROUTINE PackAndSendDataBufferInterpolation
  
!Packs the data for the buffers for the cubic spline interpolation for the multiGrid algorithm

  integer :: m,i,j,k       !Iteration variables
  integer :: iComm         !Communication direction
  integer :: dest	   ! rank of destination processing unit
  integer :: tag	   ! message tag
  integer :: mpierr	   ! MPI standard error variable 
  integer :: bufferSize    ! Size of message being transferred

  fC_bufferSendLeft_bottomXZ(1:NumDistDirs,1,:) = fCtoF_bottomXZ(:,3,:,1)
  fC_bufferSendLeft_bottomXZ(1:NumDistDirs,2,:) = fCtoF_bottomXZ(:,3,:,1+gridRatio)
  fC_bufferSendLeft_topXZ(1:NumDistDirs,1,:) = fCtoF_topXZ(:,3,:,1)
  fC_bufferSendLeft_topXZ(1:NumDistDirs,2,:) = fCtoF_topXZ(:,3,:,1+gridRatio)

  fC_bufferSendLeft_bottomXZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_bottomXZ(:,3,:,1)
  fC_bufferSendLeft_bottomXZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_bottomXZ(:,3,:,1+gridRatio)  
  fC_bufferSendLeft_topXZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_topXZ(:,3,:,1)
  fC_bufferSendLeft_topXZ(NumDistDirs+1:NumDistDirs+2,2,:) = dsCtoF_topXZ(:,3,:,1+gridRatio)
  
  fC_bufferSendRight_bottomXZ(1:NumDistDirs,1,:) = fCtoF_bottomXZ(:,3,:,nzSub_fine-gridRatio+1)
  fC_bufferSendRight_topXZ(1:NumDistDirs,1,:) = fCtoF_topXZ(:,3,:,nzSub_fine-gridRatio+1)

  fC_bufferSendRight_bottomXZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_bottomXZ(:,3,:,nzSub_fine-gridRatio+1)
  fC_bufferSendRight_topXZ(NumDistDirs+1:NumDistDirs+2,1,:) = fCtoF_topXZ(:,3,:,nzSub_fine-gridRatio+1)

  fC_bufferSendLeft_frontYZ(1:NumDistDirs,1,:) = fCtoF_frontYZ(:,3,:,1)
  fC_bufferSendLeft_frontYZ(1:NumDistDirs,2,:) = fCtoF_frontYZ(:,3,:,1+gridRatio)
  fC_bufferSendLeft_backYZ(1:NumDistDirs,1,:) = fCtoF_backYZ(:,3,:,1)
  fC_bufferSendLeft_backYZ(1:NumDistDirs,2,:) = fCtoF_backYZ(:,3,:,1+gridRatio)

  fC_bufferSendLeft_frontYZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_frontYZ(:,3,:,1)
  fC_bufferSendLeft_frontYZ(NumDistDirs+1:NumDistDirs+2,2,:) = dsCtoF_frontYZ(:,3,:,1+gridRatio)
  fC_bufferSendLeft_backYZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_backYZ(:,3,:,1)
  fC_bufferSendLeft_backYZ(NumDistDirs+1:NumDistDirs+2,2,:) = dsCtoF_backYZ(:,3,:,1+gridRatio)

  fC_bufferSendRight_frontYZ(1:NumDistDirs,1,:) = fCtoF_frontYZ(:,3,:,nzSub_fine-gridRatio+1)
  fC_bufferSendRight_backYZ(1:NumDistDirs,1,:) = fCtoF_backYZ(:,3,:,nzSub_fine-gridRatio+1)

  fC_bufferSendRight_frontYZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_frontYZ(:,3,:,nzSub_fine-gridRatio+1)
  fC_bufferSendRight_backYZ(NumDistDirs+1:NumDistDirs+2,1,:) = dsCtoF_backYZ(:,3,:,nzSub_fine-gridRatio+1)

  iComm = 6                      !Send to processor on the left in the z direction
  dest = SubID(iComm) - 1_lng    ! rank of processing unit receiving message from the current processing unit (-1 to correspond to rank (myid))
  bufferSize = (NumDistDirs+2) * 2 * nxSub_fine
  tag = iComm*100_lng		  
  CALL MPI_ISEND(fC_bufferSendLeft_bottomXZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(1),mpierr)
  tag = iComm*100_lng + 1	  
  CALL MPI_ISEND(fC_bufferSendLeft_topXZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(2),mpierr)
  bufferSize = (NumDistDirs+2) * 2 * (nySub_fine-2)
  tag = iComm*100_lng + 2
  CALL MPI_ISEND(fC_bufferSendLeft_frontYZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(3),mpierr)
  tag = iComm*100_lng + 3
  CALL MPI_ISEND(fC_bufferSendLeft_backYZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(4),mpierr)
  
  iComm = 5                      !Send to processor on the right in the z direction
  dest = SubID(iComm) - 1_lng    ! rank of processing unit receiving message from the current processing unit (-1 to correspond to rank (myid))
  bufferSize = (NumDistDirs+2) * 1 * nxSub_fine
  tag = iComm*100_lng		  
  CALL MPI_ISEND(fC_bufferSendRight_bottomXZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(5),mpierr)
  tag = iComm*100_lng + 1	  
  CALL MPI_ISEND(fC_bufferSendRight_topXZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(6),mpierr)
  bufferSize = (NumDistDirs+2) * 1 * (nySub_fine-2)
  tag = iComm*100_lng + 2
  CALL MPI_ISEND(fC_bufferSendRight_frontYZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(7),mpierr)
  tag = iComm*100_lng + 3
  CALL MPI_ISEND(fC_bufferSendRight_backYZ,bufferSize,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,reqInterpolationBuffer(8),mpierr)

END SUBROUTINE PackAndSendDataBufferInterpolation


SUBROUTINE ReceiveAndUnpackDataBufferInterpolation
  
!Packs the data for the buffers for the cubic spline interpolation for the multiGrid algorithm

  integer :: m,i,j,k       ! Iteration variables
  integer :: iComm         ! Communication direction
  integer :: src	   ! rank of source processing unit
  integer :: tag	   ! message tag
  integer :: mpierr	   ! MPI standard error variable 
  integer :: bufferSize    ! Size of message being transferred
  
  iComm = 6                ! Send to processor on the left in the z direction
  src = SubID(OppCommDir(iComm)) - 1_lng  
  bufferSize = (NumDistDirs+2) * 2 * nxSub_fine
  tag = iComm*100_lng		  
  CALL MPI_IRECV(fC_bufferRecvRight_bottomXZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(1),mpierr)
  tag = iComm*100_lng + 1	  
  CALL MPI_IRECV(fC_bufferRecvRight_topXZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(2),mpierr)
  bufferSize = (NumDistDirs+2) * 2 * (nySub_fine - 2)
  tag = iComm*100_lng + 2
  CALL MPI_IRECV(fC_bufferRecvRight_frontYZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(3),mpierr)
  tag = iComm*100_lng + 3
  CALL MPI_IRECV(fC_bufferRecvRight_backYZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(4),mpierr)

  iComm = 5                      !Send to processor on the right in the z direction
  src = SubID(OppCommDir(iComm)) - 1_lng    ! rank of processing unit receiving message from the current processing unit (-1 to correspond to rank (myid))
  bufferSize = (NumDistDirs+2) * 1 * nxSub_fine
  tag = iComm*100_lng		  
  CALL MPI_IRECV(fC_bufferRecvLeft_bottomXZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(5),mpierr)
  tag = iComm*100_lng + 1	  
  CALL MPI_IRECV(fC_bufferRecvLeft_topXZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(6),mpierr)
  bufferSize = (NumDistDirs+2) * 1 * (nySub_fine - 2)
  tag = iComm*100_lng + 2
  CALL MPI_IRECV(fC_bufferRecvLeft_frontYZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(7),mpierr)
  tag = iComm*100_lng + 3
  CALL MPI_IRECV(fC_bufferRecvLeft_backYZ,bufferSize,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,reqInterpolationBuffer(8),mpierr)

  CALL MPI_WAITALL(8,reqInterpolationBuffer,waitStatInterpolationBuffer,mpierr)  

  fCtoF_bottomXZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_bottomXZ(1:NumDistDirs,1,:)
  dsCtoF_bottomXZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_bottomXZ(NumDistDirs+1:NumDistDirs+2,1,:)

  fCtoF_topXZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_topXZ(1:NumDistDirs,1,:)
  dsCtoF_topXZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_topXZ(NumDistDirs+1:NumDistDirs+2,1,:)
  
  fCtoF_bottomXZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_bottomXZ(1:NumDistDirs,1,:)
  dsCtoF_bottomXZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_bottomXZ(NumDistDirs+1:NumDistDirs+2,1,:)
  fCtoF_bottomXZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_bottomXZ(1:NumDistDirs,2,:)
  dsCtoF_bottomXZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_bottomXZ(NumDistDirs+1:NumDistDirs+2,2,:)
  
  fCtoF_topXZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_topXZ(1:NumDistDirs,1,:)
  dsCtoF_topXZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_topXZ(NumDistDirs+1:NumDistDirs+2,1,:)
  fCtoF_topXZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_topXZ(1:NumDistDirs,2,:)
  dsCtoF_topXZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_topXZ(NumDistDirs+1:NumDistDirs+2,2,:)

  fCtoF_frontYZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_frontYZ(1:NumDistDirs,1,:)
  dsCtoF_frontYZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_frontYZ(NumDistDirs+1:NumDistDirs+2,1,:)
  
  fCtoF_backYZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_backYZ(1:NumDistDirs,1,:)
  dsCtoF_backYZ(:,3,:,-gridRatio+1) = fC_bufferRecvLeft_backYZ(NumDistDirs+1:NumDistDirs+2,1,:)

  fCtoF_frontYZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_frontYZ(1:NumDistDirs,1,:)
  dsCtoF_frontYZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_frontYZ(NumDistDirs+1:NumDistDirs+2,1,:)
  fCtoF_frontYZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_frontYZ(1:NumDistDirs,2,:)
  dsCtoF_frontYZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_frontYZ(NumDistDirs+1:NumDistDirs+2,2,:)

  fCtoF_backYZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_backYZ(1:NumDistDirs,1,:)
  dsCtoF_backYZ(:,3,:,nzSub_fine+1) = fC_bufferRecvRight_backYZ(NumDistDirs+1:NumDistDirs+2,1,:)
  fCtoF_backYZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_backYZ(1:NumDistDirs,2,:)
  dsCtoF_backYZ(:,3,:,nzSub_fine+1+gridRatio) = fC_bufferRecvRight_backYZ(NumDistDirs+1:NumDistDirs+2,2,:)

END SUBROUTINE ReceiveAndUnpackDataBufferInterpolation

!================================================
END MODULE Parallel_Fine
!================================================
