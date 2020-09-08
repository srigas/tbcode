program TB
	implicit none
	real*8 :: a_1(3), a_2(3), a_3(3), RMAX, R0, TTPRIME(3), KPOINT(3), RPOINT(3), epsilon, NUMEL, min_val, max_val, mixfactor
	real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, rho, newrho, rhozero, EIGENVALUES, diff, TEMPEIGENVALUES, SORTEDEIGVALS, &
	& UNIQUEEIGVALS, INUMEL
	real*8, allocatable, dimension(:,:) :: KPTS, TPTS, RLATT, RHOCHAIN
	complex*16, allocatable, dimension(:) :: WORK
	complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, EIGENVECTORS
	complex*16 :: CI
	character :: answer*1
	integer, allocatable, dimension(:) :: alphas, multiplicity
	real*8, allocatable, dimension(:,:) :: intnumdensity, numdensity
	integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, LWORK, INFO, NCELLS, ini, fin, reps, &
	& maxreps, minl(1), uniquecounter

	CI = (0.0,1.0) ! setting the imaginary unit

	print *, 'If everything has been setup correctly in config.dat, press y to continue, otherwise press n to exit.'
	read *, answer

	if (answer /= 'y') then
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

	allocate(TPTS(3,NUMT))
	allocate(E0(NUMT))
	allocate(INUMEL(NUMT)) ! NUMT x 1 column with number of electrons of each basis atom
	allocate(ULCN(NUMT)) ! NUMT x 1 column with the U_LCN constant for each basis atom
	allocate(rhozero(NUMT)) ! NUMT x 1 column with the charges rho_0 of each basis atom
	allocate(rho(NUMT)) ! NUMT x 1 column with the charges rho of each basis atom
	allocate(newrho(NUMT)) ! NUMT x 1 column with the charges rho of each basis atom
	allocate(diff(NUMT))

	NUMT = 1000
    
    do i = 1,NUMT
        TPTS(1,i) = i-1 !These create a NUMT-site chain
        TPTS(2,i) = 0.0
        TPTS(3,i) = 0.0
        E0(i) = 0.0
        ULCN(i) = 0.0
	end do

    print *, 'Please enter the center atom`s on site energy.'
    read *, E0(INT((NUMT+1)/2))
    print *, 'Please enter the center atom`s U.'
    read *, ULCN(INT((NUMT+1)/2))

    print *, 'Please enter the charge per atom'
    read *, INUMEL(1)

    do i = 1, NUMT
        INUMEL(i) = INUMEL(1)
        rhozero(i) = INUMEL(1)
        rho(i) = 0.3 !The initial value for the charges
    end do

	! Calculates the total number of electrons
	NUMEL = 0
	do i = 1, NUMT
		NUMEL = NUMEL + INUMEL(i)*NUMK
	end do

	allocate(EIGENVECTORS(NUMT,NUMT*NUMK))
	allocate(EIGENVALUES(NUMT*NUMK))
	allocate(HAMILTONIAN(NUMT,NUMT))
	allocate(TEMPEIGENVALUES(NUMT*NUMK))
	allocate(alphas(NUMT*NUMK))

	IRLATTMAX =  (2*NCELLS+1)**3

	call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) !Here we construct the RLATT Matrix

	
	!The following are for the configuration of zheev.
	!-------------------------------------------------
	allocate(W(NUMT))
	allocate(RWORK(3*NUMT - 2))
	LWORK = NUMT*(NUMT+1)
	allocate(WORK(LWORK))
	!-------------------------------------------------

	! These configurations ensure that the following while loop is initiated
	print *, 'Please insert the maximum difference epsilon for convergence.'
	read *, epsilon
	reps = 0
	do i = 1, NUMT
		diff(i) = 1.0
	end do
	
	print *, 'Please insert the maximum number of runs for the procedure, even if convergence is not achieved.'
	read *, maxreps
	print *, 'Please insert the mixing factor for the calculation of the new charges after every run of the self-consistent cycle.'
	read *, mixfactor

	print *, 'Initiating charge densities calculation...'
    do while (MAXVAL(diff) > epsilon .and. reps < maxreps) ! Check for convergence
        
        do i = 1, NUMT*NUMK
            alphas(i) = 0
        end do

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

			call zheev ('V', 'U', NUMT, HAMILTONIAN, NUMT, W, WORK, LWORK, RWORK, INFO)

			ini = (kcounter-1)*NUMT + 1
			fin = kcounter*NUMT
			EIGENVALUES(ini:fin) = W
			EIGENVECTORS(:,ini:fin) = HAMILTONIAN
		end do

		! At that point all the eigenvalues are in the form (.,.,.,.,...) and all the eigenvectors are NUMT*NUMK columns of NUMT rows
		! Now we want to find the N_e LOWEST eigenvalues and their corresponding eigenvectors. We thus assign a row with 0's and 1's depending
		! on whether the corresponding value is amongst the N_e LOWEST or not.
		TEMPEIGENVALUES = EIGENVALUES

		do i = 1, NINT(NUMEL)
			minl = MINLOC(TEMPEIGENVALUES)
			alphas(minl) = 1
			TEMPEIGENVALUES(minl) = MAXVAL(TEMPEIGENVALUES)
		end do

		! Now we can ensure that all calculations are performed for E < E_Fermi
		do i = 1, NUMT
			newrho(i) = 0.0
			do j = 1, NUMT*NUMK
				if (alphas(j) == 1) then
					newrho(i) = newrho(i) + abs(EIGENVECTORS(i,j))**2
				endif
			end do
			newrho(i) = newrho(i)/NUMK
			diff(i) = abs(newrho(i) - rho(i))
		end do

        rho = (1.0 - mixfactor)*rho + mixfactor*newrho
		reps = reps + 1

	end do
	print *, 'Charge densities evaluated.'

	deallocate(TEMPEIGENVALUES)

	! This prints the charges for the chain model test
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	allocate(RHOCHAIN(2,NUMT))
	do i = 1, NUMT
		RHOCHAIN(1,i) = i
		RHOCHAIN(2,i) = newrho(i)
	end do
	open(16, file = 'rhochain.txt', action = 'write')
		do j = 1, NUMT
			write (16,100) (RHOCHAIN(i,j), i = 1,2)
		end do
		100 format(3F17.8)
	close(16)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!At that point we have a good approximation for the charges. We move on to calculate once again the Hamiltonian and then the states density.
	!do kcounter = 1, NUMK
		!KPOINT = KPTS(1:3,kcounter)
		!do i = 1, NUMT
			!do j = 1, NUMT
			!	HAMILTONIAN(i,j) = (0.0,0.0)
			!	TTPRIME = TPTS(1:3,j) - TPTS(1:3,i)
			!	do IRLATT = 1, IRLATTMAX
			!		RPOINT = RLATT(1:3,IRLATT)
			!		HAMILTONIAN(i,j) = HAMILTONIAN(i,j) + exp(CI*DOT_PRODUCT(KPOINT,RPOINT))*HOP(RPOINT + TTPRIME,NUMT,E0,j,R0,RMAX)
			!	end do
			!end do
		!end do

		!call zheev ('V', 'U', NUMT, HAMILTONIAN, NUMT, W, WORK, LWORK, RWORK, INFO)

!		ini = (kcounter-1)*NUMT + 1
!		fin = kcounter*NUMT
!		EIGENVALUES(ini:fin) = W
!		EIGENVECTORS(:,ini:fin) = HAMILTONIAN
!	end do

	!This takes the different eigenvalues and organizes them in ascending order
!	allocate(UNIQUEEIGVALS(NUMT*NUMK))
!	uniquecounter = 0
!	min_val = MINVAL(EIGENVALUES) - 1.0 ! -1.0 Is inserted for a case of complete degeneracy
!	max_val = MAXVAL(EIGENVALUES)
!	do while (min_val < max_val)
!		uniquecounter = uniquecounter + 1
!		min_val = MINVAL(EIGENVALUES, mask = EIGENVALUES > min_val)
!		UNIQUEEIGVALS(uniquecounter) = min_val
!	enddo
!	allocate(SORTEDEIGVALS(uniquecounter))
!	SORTEDEIGVALS = UNIQUEEIGVALS(1:uniquecounter)

!	call INT_NUM_DEN(uniquecounter,EIGENVALUES,NUMT,NUMK,SORTEDEIGVALS,multiplicity,intnumdensity)
	
!	call NUM_DEN(uniquecounter,SORTEDEIGVALS,EIGENVALUES,numdensity,NUMT,NUMK)

	!------------------------------------------------------------------------------------------------------------------

	contains

	subroutine INT_NUM_DEN(uniquecounter,EIGENVALUES,NUMT,NUMK,SORTEDEIGVALS,multiplicity,intnumdensity)
		implicit none

		integer, allocatable, dimension (:) :: multiplicity
		real*8, allocatable, dimension(:,:) :: intnumdensity
		integer :: uniquecounter, NUMT, NUMK, i, j
		real*8 :: EIGENVALUES(NUMT*NUMK), SORTEDEIGVALS(uniquecounter)

		allocate(multiplicity(uniquecounter))
		allocate(intnumdensity(2,uniquecounter))
		intnumdensity(1,:) = SORTEDEIGVALS

		do i = 1, uniquecounter
			multiplicity(i) = 0
			do j = 1, NUMT*NUMK
				if (SORTEDEIGVALS(i) == EIGENVALUES(j)) then
					multiplicity(i) = multiplicity(i) + 1
				endif
			end do
			if (i == 1) then
				intnumdensity(2,i) = multiplicity(i)
			else
				intnumdensity(2,i) = intnumdensity(2,i-1) + multiplicity(i)
			endif
		end do

		open(13, file = 'intnumdensity.txt', action = 'write')
		do j = 1, uniquecounter
			write (13,100) (intnumdensity(i,j), i = 1,2)
		end do
		100 format(3F17.8)
		close(13)

		! This calculates and prints the Fermi Energy
		do i = 1, uniquecounter
			if (intnumdensity(2,i) > NINT(NUMEL)) then
				print *, 'The Fermi Energy is', intnumdensity(1,i-1)
				exit
			endif
		end do
		
	end subroutine INT_NUM_DEN

	subroutine NUM_DEN(uniquecounter,SORTEDEIGVALS,EIGENVALUES,numdensity,NUMT,NUMK)
		implicit none

		integer :: uniquecounter, i, NUMT, NUMK, numenergyintervals
		real*8, allocatable, dimension(:,:) :: numdensity
		real*8 :: pi, delta, SORTEDEIGVALS(uniquecounter), EIGENVALUES(NUMT*NUMK), energyintervals

		pi = 4.D0*atan(1.D0)
		print *, 'Please insert the lorentzian delta factor.'
		read *, delta

		print *, 'Please enter the number of intervals for the plot of n(E).'
		read *, numenergyintervals

		allocate(numdensity(2,numenergyintervals+1))

		energyintervals = (MAXVAL(SORTEDEIGVALS) - MINVAL(SORTEDEIGVALS))/numenergyintervals

		do i = 0, numenergyintervals
			numdensity(1,i+1) = MINVAL(SORTEDEIGVALS) + energyintervals*i
		end do

		do i = 0, numenergyintervals
			numdensity(2,i+1) = 0.0
			do j = 1, NUMT*NUMK
				numdensity(2,i+1) = numdensity(2,i+1) + (1.0/pi)*delta/((numdensity(1,i+1) - EIGENVALUES(j))**2 + delta**2)
			end do
		end do

		open(15, file = 'numdensity.txt', action = 'write')
		do j = 1, numenergyintervals
			write (15,100) (numdensity(i,j), i = 1,2)
		end do
		100 format(3F17.8)
		close(15)

	end subroutine NUM_DEN

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
		!It is not necesary to have these on a file
		!open(12, file = 'BZ.txt', action = 'write')

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

		!write (12,'(3F12.6)') KPTS !'(3F10.4)' formats the way it is printed on the .txt file
		!close(12)

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

		if (norm2(TTPRIME) == 0.0) then ! This case corresponds to t = t', R = 0
			hopping = E0(i) + ULCN(i)*(rho(i) - rhozero(i))
		else if (norm2(TTPRIME) < RMAX) then
			hopping = exp((-1)*norm2(TTPRIME)/R0)
		else
			hopping = 0.0
		endif
	end function HOP

end program TB