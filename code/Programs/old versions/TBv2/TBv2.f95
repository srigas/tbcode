program TB
	implicit none
	real*8 :: a_1(3), a_2(3), a_3(3), RMAX, R0, TTPRIME(3), E0, KPOINT(3), RPOINT(3), k_1(3), k_2(3)
	real*8, allocatable, dimension(:) :: W, RWORK
	real*8, allocatable, dimension(:,:) :: KPTS, TPTS, RLATT, K, RESULTS
	complex*16, allocatable, dimension(:) :: WORK
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN
	complex*16 :: CI
	character :: answer*1
	integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, N_k, LWORK, INFO, NMAX, NCELLS

	CI = (0.0,1.0) !setting the imaginary unit

	print *, 'If everything has been setup correctly in config.dat, press y to continue, otherwise press n to exit.'
	read *, answer

	if (answer == 'n') then
		call exit(123)
	endif

	open(10, file = 'config.dat', action = 'read')
	read(10,*) a_1
	read(10,*) a_2
	read(10,*) a_3
	read(10,*) NUMKX,NUMKY,NUMKZ
	read(10,*) RMAX
	read(10,*) R0
	read(10,*) E0
	read(10,*) NCELLS
	close(10)

	if (DOT_PRODUCT(a_1,a_2) /= 0 .or. DOT_PRODUCT(a_1,a_3) /= 0 .or. DOT_PRODUCT(a_2,a_3) /= 0 .or. RMAX < 0 .or. R0 <= 0) then
		print *, 'A value inserted in config.dat is incorrect. Please try again after everything has been corrected.'
		call exit(123)
	endif

	call BZ(a_1,a_2,a_3,NUMKX,NUMKY,NUMKZ,KPTS,NUMK)
	
	!The following reads the number of basis vectors from a file named basisvectors.dat
	NUMT = 0
	open (1, file = 'basisvectors.dat', action = 'read')
	do
		read(1,*,iostat=io)
		if (io/=0) exit
		NUMT = NUMT + 1
	end do
	close(1)

	allocate(TPTS(3,NUMT))

	open (1, file = 'basisvectors.dat', action = 'read')
	do i = 1,NUMT
		read(1,*) TPTS(:,i)
	end do
	close(1)
	
	allocate(HAMILTONIAN(NUMT,NUMT))

	IRLATTMAX =  (2*NCELLS+1)**3

	call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) !Here we construct the RLATT Matrix

	print *, 'Please insert the vector k_1 as X,X,X.'
	read *, k_1(1),k_1(2),k_1(3)
	
	print *, 'Please insert the vector k_2 as X,X,X.'
	read *, k_2(1),k_2(2),k_2(3)

	print *, 'Please insert the number of partitions N_k.'
	read *, N_k
	
	!The following are for the configuration of zheev.
	allocate(W(NUMT))
	allocate(RWORK(3*NUMT - 2))
	LWORK = NUMT*(NUMT+1)
	allocate(WORK(LWORK))
	allocate(RESULTS(3+NUMT,N_k+1))

	allocate(K(3,N_k+1))

	do i = 0,N_k
		K(1:3,i+1) = (k_2-k_1)*i/N_k
	end do

        do kcounter = 1,N_k+1
		KPOINT = K(1:3,kcounter)

		do i = 1, NUMT
			do j = 1, NUMT
				HAMILTONIAN(i,j) = (0.0,0.0)
				TTPRIME = TPTS(1:3,j) - TPTS(1:3,i)
				do IRLATT = 1, IRLATTMAX
					RPOINT = RLATT(1:3,IRLATT)
					HAMILTONIAN(i,j) = HAMILTONIAN(i,j) + exp(CI*DOT_PRODUCT(KPOINT,RPOINT))*HOP(RPOINT + TTPRIME,E0,R0,RMAX)
				end do
			end do
		end do

		call zheev ('N', 'U', NUMT, HAMILTONIAN, NUMT, W, WORK, LWORK, RWORK, INFO) !calculates the eigenvalues
		RESULTS(1:3,kcounter) = KPOINT
        RESULTS(4:3+NUMT,kcounter) = W

	end do

	NMAX = min(NUMT+3,13) !The maximum number of NUMT is 10 in this case.

    open(11, file = 'results.txt', action = 'write')
	do kcounter = 1,N_k+1
        write(11,100) (RESULTS(i,kcounter), i = 1, NMAX)
    end do

	100 format(3F10.4,10E16.8)
		
	!------------------------------------------------------------------------------------------------------------------

	contains

	subroutine RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX)
		implicit none

		integer :: i,j,k,IRLATTMAX,counter,NCELLS
		real*8 :: a_1(3), a_2(3), a_3(3)
		real*8, allocatable, dimension(:,:) :: RLATT
		counter = 1

		allocate(RLATT(3,IRLATTMAX))

		do i = -NCELLS,NCELLS
			do j = -NCELLS,NCELLS
				do k = -NCELLS,NCELLS
					RLATT(1:3,counter) = i*a_1 + j*a_2 + k*a_3
					counter = counter + 1
				end do
			end do
		end do

	end subroutine RPTS

	subroutine BZ(a,b,c,N_x,N_y,N_z,KPTS,Ntot)
		implicit none

		integer :: N_x, N_y, N_z, Ntot, counter, c_1, c_2, c_3
		real*8 :: a(3), b(3), c(3), b_1(3), b_2(3), b_3(3)
		real*8, allocatable, dimension(:,:) :: KPTS
		real*8 :: pi, volume
		pi = 4.D0*atan(1.D0)
		counter = 1
		
		volume = DOT_PRODUCT(a, CROSS_PRODUCT(b,c)) !The volume of the unit cell

		!At this point we calculate the reciprocal space vectors.
		b_1 = 2*pi*CROSS_PRODUCT(b,c)/volume
		b_2 = 2*pi*CROSS_PRODUCT(c,a)/volume
		b_3 = 2*pi*CROSS_PRODUCT(a,b)/volume

		Ntot = N_x*N_y*N_z !The total number of different wavevectors in reciprocal space.

		allocate(KPTS(3,Ntot))
		open(12, file = 'BZ.txt', action = 'write')

		do c_1 = 1, N_x
			do c_2 = 1, N_y
				do c_3 = 1, N_z
					KPTS(1, counter) = (b_1(1)*c_1)/dfloat(N_x) + (b_2(1)*c_2)/dfloat(N_y) + (b_3(1)*c_3)/dfloat(N_z)
					KPTS(2, counter) = (b_1(2)*c_1)/dfloat(N_x) + (b_2(2)*c_2)/dfloat(N_y) + (b_3(2)*c_3)/dfloat(N_z)
					KPTS(3, counter) = (b_1(3)*c_1)/dfloat(N_x) + (b_2(3)*c_2)/dfloat(N_y) + (b_3(3)*c_3)/dfloat(N_z)
					counter = counter + 1
				end do
			end do
		end do

		write (12,'(3F12.6)') KPTS !'(3F10.4)' formats the way it is printed on the .txt file
		close(12)

	end subroutine BZ

	function CROSS_PRODUCT(x,y) result(cross)
		implicit none
		real*8, dimension(3), intent(in) :: x, y
		real*8, dimension(3) :: cross

		cross(1) = x(2)*y(3) - x(3)*y(2)
		cross(2) = x(3)*y(1) - x(1)*y(3)
		cross(3) = x(1)*y(2) - x(2)*y(1)
		
	end function CROSS_PRODUCT

	function HOP(TTPRIME,E0,R0,RMAX) result(hopping) !TTPRIME = t - t'
		implicit none
		real*8, intent(in) :: TTPRIME(3)
		real*8 :: hopping, E0, R0, RMAX

		if (norm2(TTPRIME) == 0.0) then !This case corresponds to t = t', R = 0
			hopping = E0
		else if (norm2(TTPRIME) < RMAX) then
			hopping = exp((-1)*norm2(TTPRIME)/R0)
		else
			hopping = 0.0
		endif
	end function HOP

end program TB