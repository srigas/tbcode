program TB
	implicit none
	real*8 :: a_1(3), a_2(3), a_3(3), RMAX, R0, TTPRIME(3), KPOINT(3), RPOINT(3), TEMPVAL
	real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, rho, rhozero, EIGENVALUES, EIGENVALUESEXTRA
	real*8, allocatable, dimension(:,:) :: KPTS, TPTS, RLATT
	complex*16, allocatable, dimension(:) :: WORK, TEMPVEC
	complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, EIGENVECTORS
	complex*16 :: CI
	character :: answer*1
	integer, allocatable, dimension(:) :: INUMEL, ALPHAS
	integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, LWORK, INFO, NCELLS, NUMEL, ini, fin, minl(1)

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
	allocate(E0(NUMT))
	allocate(INUMEL(NUMT)) ! NUMT x 1 column with number of electrons of each basis atom
	allocate(ULCN(NUMT)) ! NUMT x 1 column with the U_LCN constant for each basis atom
	allocate(rhozero(NUMT)) ! NUMT x 1 column with the charges rho_0 of each basis atom
	allocate(rho(NUMT)) ! NUMT x 1 column with the charges rho of each basis atom

	open (1, file = 'basisvectors.dat', action = 'read')
	do i = 1,NUMT
		read(1,*) TPTS(1:3,i), E0(i), INUMEL(i), ULCN(i), rhozero(i)
	end do
	close(1)
	
	!calculates the number of electrons in the unit cell
	NUMEL = 0
	do i = 1, NUMT
		NUMEL = NUMEL + INUMEL(i)
	end do

	allocate(EIGENVECTORS(NUMT,NUMT*NUMK))
	allocate(EIGENVALUES(NUMT*NUMK))
	allocate(HAMILTONIAN(NUMT,NUMT))
	allocate(EIGENVALUESEXTRA(NUMT*NUMK))
	allocate(ALPHAS(NUMT*NUMK))

	IRLATTMAX =  (2*NCELLS+1)**3

	call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) !Here we construct the RLATT Matrix

	!Sets a set of initial values for rho which will be corrected in the self consistent run of the algorithm later
	do i = 1, NUMT
		print *, 'Please insert the initial value for the charge of atom number', i
		read *, rho(i)
	end do
	
	!The following are for the configuration of zheev.
	!-------------------------------------------------
	allocate(W(NUMT))
	allocate(RWORK(3*NUMT - 2))
	LWORK = NUMT*(NUMT+1)
	allocate(WORK(LWORK))
	!-------------------------------------------------

	do kcounter = 1, NUMK
		KPOINT = KPTS(1:3,kcounter)
		do i = 1, NUMT
			do j = 1, NUMT
				HAMILTONIAN(i,j) = (0.0,0.0)
				TTPRIME = TPTS(1:3,j) - TPTS(1:3,i)
				do IRLATT = 1, IRLATTMAX
					RPOINT = RLATT(1:3,IRLATT)
					HAMILTONIAN(i,j) = HAMILTONIAN(i,j) + exp(CI*DOT_PRODUCT(KPOINT,RPOINT))*HOP(RPOINT + TTPRIME,NUMT,E0,j,R0,RMAX)
				end do
			end do
		end do

		call zheev ('V', 'U', NUMT, HAMILTONIAN, NUMT, W, WORK, LWORK, RWORK, INFO) !calculates the eigenvalues and eigenvectors

		ini = (kcounter-1)*NUMT + 1
		fin = kcounter*NUMT
		EIGENVALUES(ini:fin) = W
		EIGENVECTORS(:,ini:fin) = HAMILTONIAN
	end do

	open(11, file = 'FINALS1.txt', action = 'write')
	do kcounter = 1,NUMK*NUMT
        write(11,100) (EIGENVECTORS(i,kcounter), i = 1, NUMT)
    end do

	100 format(10E16.8)

	close(11)

	open(12, file = 'FINALS2.txt', action = 'write')
	write(12,*) EIGENVALUES
	close(12)

	allocate(TEMPVEC(NUMT))
	
	EIGENVALUESEXTRA = EIGENVALUES

	do i = 1, NUMEL
		minl = MINLOC(EIGENVALUESEXTRA)
		alphas(minl) = 1
		EIGENVALUESEXTRA(minl) = MAXVAL(EIGENVALUESEXTRA)
	end do

	do i = 1, NUMK*NUMT
		if (alphas(i) /= 1) then
			alphas(i) = 0
		endif
	end do

	open(15, file = 'alphas.txt', action = 'write')
	write(15,*) alphas
	close(15)

	!At that point all the eigenvalues are in the form (.,.,.,.,...) and all the eigenvectors are NUMT*NUMK columns of NUMT rows
	!Now we want to find the N_e LOWEST eigenvalues and their corresponding eigenvectors.
	do i = 1, NUMT*NUMK
		do j = i, NUMT*NUMK
			if (EIGENVALUES(i) > EIGENVALUES(j)) then
				TEMPVAL = EIGENVALUES(j)
				EIGENVALUES(j) = EIGENVALUES(i)
				EIGENVALUES(i) = TEMPVAL

				TEMPVEC = EIGENVECTORS(:,j)
				EIGENVECTORS(:,j) = EIGENVECTORS(:,i)
				EIGENVECTORS(:,i) = TEMPVEC
			endif
		end do
	end do


	open(13, file = 'FINALS3.txt', action = 'write')
	do kcounter = 1,NUMK*NUMT
        write(13,102) (EIGENVECTORS(i,kcounter), i = 1, NUMT)
    end do

	102 format(10E16.8)

	close(13)

	open(14, file = 'FINALS4.txt', action = 'write')
	write(14,*) EIGENVALUES
	close(14)


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

	function HOP(TTPRIME,NUMT,E0,i,R0,RMAX) result(hopping) !TTPRIME = t - t'
		implicit none
		real*8, intent(in) :: TTPRIME(3)
		real*8 :: hopping, R0, RMAX
		integer :: NUMT, i
		real*8, dimension(NUMT) :: E0

		if (norm2(TTPRIME) == 0.0) then !This case corresponds to t = t', R = 0
			hopping = E0(i) !+ ULCN(i)*(rho(i)-rhozero(i))
		else if (norm2(TTPRIME) < RMAX) then
			hopping = exp((-1)*norm2(TTPRIME)/R0)
		else
			hopping = 0.0
		endif
	end function HOP

end program TB