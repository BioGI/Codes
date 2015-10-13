IMPLICIT NONE
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(13,307)	! floating point precision (double)
INTEGER, PARAMETER :: lng = KIND(10000000)					! maximum integer value ("long")
INTEGER, allocatable :: seed1(:)
INTEGER :: seed_size1,seed_date1(8),nbins,nbinmiddle,i,np,j,npindex,k
REAL(dbl), ALLOCATABLE	:: randnomono(:),randno(:),v0R(:),Q0R(:),Q0RdR(:),v0RdR(:),Rbins(:),Q0RdRint(:),Radlist(:)
REAl(dbl) :: R0,Rstar,sigR,sigmax,vptotal,fourbythreepi,xmax,xmin,ymax,ymin,zmax,zmin,deltaR,sumvolume	

CALL DATE_AND_TIME(VALUES=seed_date1)
CALL RANDOM_SEED(size=seed_size1)
ALLOCATE(seed1(seed_size1))
!write(*,*) 'seed_size_________',seed_size
CALL RANDOM_SEED(GET=seed1)
seed1=972
!seed=seed*(seed_date(8)-500_lng)
CALL RANDOM_SEED(put=seed1)
!write(*,*) 'seed_________',seed
DEALLOCATE(seed1)

!z0=10.0_dbl
!x0=30.0_dbl
!x1=50.0_dbl
!y0=30.0_dbl
!y1=50.0_dbl
!npx=10_lng
!npy=10_lng
!
!open(50,file='particle-a.txt')
!write(50,*) npx*npy
!  do i=1,npx
!	do j=1,npy
!  		write(50,*) (x1-x0)/(npx-1)*(i-1)+x0,(y1-y0)/(npy-1)*(j-1)+y0,z0
!	enddo
!  end do

!z0=10.0_dbl
!x0=16.0_dbl
!x1=26.0_dbl
!y0=21.0_dbl
!np=10_lng
!open(50,file='particle-a.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,y0,z0
!end do

!z0=10.0_dbl
!x0=22.0_dbl
!x1=32.0_dbl
!y0=42.0_dbl
!np=100_lng
!open(50,file='particle-a-100.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,y0,z0
!end do


!z0=10.0_dbl
!x0=22.0_dbl
!x1=32.0_dbl
!y0=36.0_dbl
!y1=46.0_dbl
!np=1000_lng
!open(50,file='particle-a-1000.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,(y1-y0)/(np-1)*(i-1)+y0,z0
!end do

!z0=10.0_dbl
!x0=10.0_dbl
!x1=30.0_dbl
!y0=21.0_dbl
!y1=21.0_dbl
!np=14_lng
!open(50,file='particle-a-14.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,(y1-y0)/(np-1)*(i-1)+y0,z0
!end do

!z0=10.0_dbl
!x0=30.0_dbl
!x1=35.0_dbl
!y0=36.0_dbl
!y1=36.0_dbl
!np=14_lng
!open(50,file='particle-a-14.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,(y1-y0)/(np-1)*(i-1)+y0,z0
!end do

!z0=10.0_dbl
!x0=40.0_dbl
!x1=80.0_dbl
!y0=61.0_dbl
!y1=61.0_dbl
!np=14_lng
!open(50,file='particle-a-14.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,(y1-y0)/(np-1)*(i-1)+y0,z0
!end do

!z0=10.0_dbl
!x0=40.0_dbl
!x1=45.0_dbl
!y0=61.0_dbl
!y1=61.0_dbl
!np=14_lng
!R0 = 0.0026 ! cm
!open(50,file='particle-a-14.txt')
!write(50,*) np
!do i=1,np
!  write(50,*) i,(x1-x0)/(np-1)*(i-1)+x0,(y1-y0)/(np-1)*(i-1)+y0,z0,R0
!end do

! Monodisperse Collection
zmin=1.0_dbl+2.0_dbl
xmin=0.2_dbl*240_dbl+2.0_dbl
ymin=0.2_dbl*240_dbl+2.0_dbl
zmax=240.0_dbl - zmin
xmax=240.0_dbl - xmin
ymax=240.0_dbl - ymin
np=14_lng!250_lng!200_lng
ALLOCATE(randnomono(3_lng*np))
CALL RANDOM_NUMBER(randnomono)
!R0 = 0.0026_dbl ! cm
!R0 = 0.00100623293815_dbl ! cm
R0 = 0.00263008138299_dbl ! cm
open(52,file='particle-a-14.txt')
write(52,*) np
do i=1,np
  write(52,*) i,xmin+(xmax-xmin)*randnomono(3*(i-1)+1),ymin+(ymax-ymin)*randnomono(3*(i-1)+2),zmin+(zmax-zmin)*randnomono(3*(i-1)+3),R0
end do
close(52)
write(*,*) np*(88.0_dbl/21.0_dbl)*(R0**3.0_dbl)

!!*****************************************************************************************************
!! Polydisperse Collection
!zmin=1.0_dbl+2.0_dbl
!xmin=0.2_dbl*240_dbl+2.0_dbl
!ymin=0.2_dbl*240_dbl+2.0_dbl
!zmax=240.0_dbl-2_dbl
!xmax=240.0_dbl - xmin
!ymax=240.0_dbl - ymin
!np=200_lng
!R0 = 0.0026_dbl ! cm
!fourbythreepi = 88.0_dbl/21.0_dbl
!vptotal = fourbythreepi*np*(R0**3_dbl)
!sigR = 0.2
!sigmax = 0.4
!nbins = 51_lng
!Rstar = R0
!nbinmiddle = ANINT(0.5_dbl*(nbins-1_lng))+1_lng
!deltaR = (1.0_dbl/nbinmiddle)*(sigmax*Rstar)
!ALLOCATE(Rbins(nbins))
!ALLOCATE(v0R(nbins))
!ALLOCATE(Q0R(nbins))
!ALLOCATE(Q0RdR(nbins))
!ALLOCATE(Q0RdRint(nbins))
!do i=1_lng,nbinmiddle
!	Rbins(i) = Rstar*(1.0_dbl+(sigmax)*(i-nbinmiddle)/(nbinmiddle))
!	Rbins(nbins-(i-1_lng)) = Rstar*(1.0_dbl+(sigmax)*(nbins-(i-1_lng)-nbinmiddle)/(nbinmiddle))
!	!write(*,*) Rbins(i),Rstar
!end do
!
!do i=1_lng,nbins
!	v0R(i) = sqrt(1.0_dbl/(44.0_dbl*(sigR**2_dbl)*(Rstar**2_dbl)/7_dbl))*exp(-0.5_dbl*((Rbins(i)-Rstar)/(sigR*Rstar))**2_dbl)
!	Q0R(i) = v0R(i)*vptotal/(fourbythreepi*(Rbins(i)**3_dbl))
!	Q0RdR(i) = (Q0R(i)*deltaR)
!	Q0RdRint(i) = NINT(Q0R(i)*deltaR)
!end do
!
!open(51,file='polydisperse_dist.txt')
!do i=1,nbins 
!  write(51,*) i,Rbins(i),v0R(i),Q0R(i),Q0RdR(i),Q0RdRint(i)
!end do
!close(51)
!np =sum(Q0RdRint)
!write(*,*) np
!ALLOCATE(randno(3_lng*np))
!CALL RANDOM_NUMBER(randno)
!open(50,file='particle-a-polydisperse.txt')
!write(50,*) np
!npindex = 0
!do i=1,nbins
!	do j = 1,INT(Q0RdRint(i),lng)
!		npindex=npindex+1
! 		write(50,*) npindex,xmin+(xmax-xmin)*randno(3*(npindex-1)+1),ymin+(ymax-ymin)*randno(3*(npindex-1)+2),zmin+(zmax-zmin)*randno(3*(npindex-1)+3),Rbins(i)
!	enddo
!enddo
!!*****************************************************************************************************

!*****************************************************************************************************
! Polydisperse Collection From Yanxing
zmin=1.0_dbl+2.0_dbl
xmin=0.2_dbl*240_dbl+2.0_dbl
ymin=0.2_dbl*240_dbl+2.0_dbl
zmax=240.0_dbl-2_dbl
xmax=240.0_dbl - xmin
ymax=240.0_dbl - ymin
np=250_lng
fourbythreepi = 88.0_dbl/21.0_dbl
nbins = 20_lng
ALLOCATE(Rbins(nbins))
ALLOCATE(v0R(nbins))
ALLOCATE(v0RdR(nbins))
ALLOCATE(Q0RdR(nbins))
ALLOCATE(Q0RdRint(nbins))
ALLOCATE(Radlist(np))
open(51,file='np250-nb20.txt')
do i=1,nbins 
  read(51,*) Rbins(i),v0R(i),v0RdR(i),Q0RdR(i),Q0RdRint(i)
end do
close(51)


!do i=1_lng,nbinmiddle
!	Rbins(i) = Rstar*(1.0_dbl+(sigmax)*(i-nbinmiddle)/(nbinmiddle))
!	Rbins(nbins-(i-1_lng)) = Rstar*(1.0_dbl+(sigmax)*(nbins-(i-1_lng)-nbinmiddle)/(nbinmiddle))
!	!write(*,*) Rbins(i),Rstar
!end do
!
!do i=1_lng,nbins
!	v0R(i) = sqrt(1.0_dbl/(44.0_dbl*(sigR**2_dbl)*(Rstar**2_dbl)/7_dbl))*exp(-0.5_dbl*((Rbins(i)-Rstar)/(sigR*Rstar))**2_dbl)
!	Q0R(i) = v0R(i)*vptotal/(fourbythreepi*(Rbins(i)**3_dbl))
!	Q0RdR(i) = (Q0R(i)*deltaR)
!	Q0RdRint(i) = NINT(Q0R(i)*deltaR)
!end do
!open(51,file='polydisperse_dist.txt')
!do i=1,nbins 
!  write(51,*) i,Rbins(i),v0R(i),Q0R(i),Q0RdR(i),Q0RdRint(i)
!end do
!close(51)

np =sum(Q0RdRint)
k = 0
!sumvolume = 0.0_dbl
do i = 1,nbins
	!sumvolume = sumvolume + INT(Q0RdRint(i))*(fourbythreepi/8.0_dbl)*((Rbins(i)*0.0001_dbl)**3.0)
	do j = 1,INT(Q0RdRint(i))
		k = k + 1
		Radlist(k) = 0.5_dbl*Rbins(i)*0.0001_dbl
		!sumvolume = sumvolume + 1.0_dbl*(fourbythreepi)*((0.5_dbl*Rbins(i)*0.0001_dbl)**3.0)
		write(*,*) k,Radlist(k)
	enddo
enddo
!write(*,*) sumvolume 
sumvolume = 0.0_dbl
do i = 1,np
	sumvolume = sumvolume +  fourbythreepi*(Radlist(i)**3.0)
enddo 

write(*,*) k,np,sumvolume
	


ALLOCATE(randno(3_lng*np))
CALL RANDOM_NUMBER(randno)
open(50,file='particle-a-polydisperse.txt')
write(50,*) np
do i=1,np
 write(50,*) i,xmin+(xmax-xmin)*randno(3*(i-1)+1),ymin+(ymax-ymin)*randno(3*(i-1)+2),zmin+(zmax-zmin)*randno(3*(i-1)+3),Radlist(i)
enddo
!*****************************************************************************************************


!z0=10.0_dbl
!x0=22.0_dbl
!x1=32.0_dbl
!y0=36.0_dbl
!y1=46.0_dbl
!np=100_lng
!open(50,file='particle-a-10000.txt')
!write(50,*) np*np
!do j=1,np
!	do i = 1,np
!  		write(50,*) i+(j-1)*np,(x1-x0)/(np-1)*(i-1)+x0,(y1-y0)/(np-1)*(j-1)+y0,z0
!	enddo
!end do

!x1=40.5
!x2=80.5
!y0=43.0
!y1=58.0
!y2=78.0
!np=10
!open(50,file='particle-a.dat')
!write(50,*) np
!write(50,*) 2*np
!do i=1,np
!  write(50,*) x1,(y1-y0)/(np-1)*(i-1)+y0,0.0
!end do
!do i=1,np
!  write(50,*) x2,(y2-y0)/(np-1)*(i-1)+y0,0.0
!end do

close(50)

stop
end
