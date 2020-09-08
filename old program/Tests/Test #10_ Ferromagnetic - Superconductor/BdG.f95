program TB
    implicit none
    real*8 :: a_1(3), a_2(3), a_3(3), RMAX, R0, TTPRIME(3), KPOINT(3), RPOINT(3), epsilon, min_val, max_val, mixfactorN, &
	&chempot, mixfactorD, readcharge
    real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, nu, newnu, nuzero, EIGENVALUES, SORTEDEIGVALS, &
	& UNIQUEEIGVALS, BETA, magnet, USUPCOND, nuup, nudown, diffN, diffD
    real*8, allocatable, dimension(:,:) :: KPTS, TPTS, RLATT
    complex*16, allocatable, dimension(:) :: WORK, DELTA, newDELTA
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, EIGENVECTORS
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), positionhamiltonian(2,2), expon, deltaterm, &
    & readdelta
    integer, allocatable, dimension(:) :: multiplicity
    real*8, allocatable, dimension(:,:) :: intnumdensity, numdensity, numdensityperatom
    integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, LWORK, INFO, NCELLS, ini, fin, reps, &
	& maxreps, uniquecounter

    CI = (0.0,1.0) ! setting the imaginary unit
	
    chempot = 0.0 ! Set to zero by default because we want to study superconductivity

    ! Setting the Pauli matrices !
    ! -------------------------------------------------
    IdentityPauli(1,1) = (1.0,0.0)
    IdentityPauli(1,2) = (0.0,0.0)
    IdentityPauli(2,1) = (0.0,0.0)
    IdentityPauli(2,2) = (1.0,0.0)
    xPauli(1,1) = (0.0,0.0)
    xPauli(1,2) = (1.0,0.0)
    xPauli(2,1) = (1.0,0.0)
    xPauli(2,2) = (0.0,0.0)
    yPauli(1,1) = (0.0,0.0)
    yPauli(1,2) = (0.0,-1.0)
    yPauli(2,1) = (0.0,1.0)
    yPauli(2,2) = (0.0,0.0)
    zPauli(1,1) = (1.0,0.0)
    zPauli(1,2) = (0.0,0.0)
    zPauli(2,1) = (0.0,0.0)
    zPauli(2,2) = (-1.0,0.0)
    ! -------------------------------------------------

    open(10, file = 'config.dat', action = 'read')
    read(10,*) a_1
    read(10,*) a_2
    read(10,*) a_3
    read(10,*) NUMKX,NUMKY,NUMKZ
    read(10,*) RMAX
    read(10,*) R0
    read(10,*) NCELLS
    close(10)

    if (DOT_PRODUCT(a_1,a_2) /= 0 .or. DOT_PRODUCT(a_1,a_3) /= 0 .or. DOT_PRODUCT(a_2,a_3) /= 0 .or. RMAX < 0 .or. R0 <= 0 &
	& .or. NCELLS < 0) then
        print *, 'A value inserted in config.dat is incorrect. Please try again after everything has been corrected.'
        call exit(123)
    endif

    call BZ(a_1,a_2,a_3,NUMKX,NUMKY,NUMKZ,KPTS,NUMK) ! We create the Brillouin Zone's k-points to be used later on
	
	! The following reads the number of basis vectors from a file named basisvectors.dat
    NUMT = 0
    open (1, file = 'basisvectors.dat', action = 'read')
    do
        read(1,*,iostat=io)
        if (io/=0) exit
        NUMT = NUMT + 1
    end do
    close(1)

    NUMT = NUMT - 2 ! To ignore the final 2 lines in basisvectors.dat which are the configuration settings.

    allocate(TPTS(3,NUMT))
    allocate(E0(NUMT))
    allocate(ULCN(NUMT)) ! NUMT x 1 column with the U_LCN constant for each basis atom
    allocate(nuzero(NUMT)) ! NUMT x 1 column with the charges n_0 of each basis atom
    allocate(nu(NUMT)) ! NUMT x 1 column with the charges n of each basis atom
    allocate(newnu(NUMT)) ! NUMT x 1 column with the charges n of each basis atom
    allocate(BETA(NUMT)) ! NUMT x 1 column with the B for each basis atom
    allocate(diffN(NUMT))
    allocate(diffD(NUMT))
    allocate(magnet(NUMT))
    allocate(DELTA(NUMT))
    allocate(newDELTA(NUMT))
    allocate(USUPCOND(NUMT))
    allocate(nuup(NUMT))
    allocate(nudown(NUMT))

    open (1, file = 'basisvectors.dat', action = 'read')
    do i = 1,NUMT
        read(1,*) TPTS(1:3,i), E0(i), ULCN(i), nuzero(i), BETA(i), USUPCOND(i)
    end do
    close(1)

    ! The *4 factors are now due to spin and particle-hole
    allocate(EIGENVECTORS(4*NUMT,4*NUMT*NUMK))
    allocate(EIGENVALUES(4*NUMT*NUMK))
    allocate(HAMILTONIAN(4*NUMT,4*NUMT))

    IRLATTMAX =  (2*NCELLS+1)**3 ! Configures how many neighbouring cells are taken into account

    call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) ! Here we construct the RLATT Matrix consisting of the lattice sites

	! Sets a set of initial values for n and D which will be corrected in the self consistent run of the algorithm later
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    print *, 'Please insert an initial value for all charges.'
    read *, readcharge
    print *, 'Please insert an initial value for all D.'
    read *, readdelta
    do i = 1, NUMT
        nu(i) = readcharge
        DELTA(i) = readdelta
    end do
	!do i = 1, NUMT
		!print *, 'Please insert the initial value for the charge of atom number', i
		!read *, nu(i)
	!end do
	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	!The following are for the configuration of zheev. The *4 factors are due to spin and particle-hole
	!-------------------------------------------------
    allocate(W(4*NUMT))
    allocate(RWORK(3*(4*NUMT) - 2))
    LWORK = 4*NUMT*(4*NUMT+1)
    allocate(WORK(LWORK))
	!-------------------------------------------------

	! These configurations ensure that the following while loop is initiated
    print *, 'Please insert the maximum difference epsilon for convergence.'
    read *, epsilon
    reps = 0
    do i = 1, NUMT
        diffN(i) = 1.0
        diffD(i) = 1.0
    end do

    print *, 'Please insert the maximum number of runs for the procedure, even if convergence is not achieved.'
    read *, maxreps
    print *, 'Please insert the mixing factor for the calculation of the new charges after every run of the self-consistent cycle.'
    read *, mixfactorN
    print *, 'Please insert the mixing factor for the calculation of D after every run of the self-consistent cycle.'
    read *, mixfactorD

    print *, 'Initiating charge densities calculation...'
    do while ((MAXVAL(diffN) > epsilon .or. MAXVAL(diffD) > epsilon) .and. reps < maxreps) ! Check for convergence or maxreps

        do kcounter = 1, NUMK
            KPOINT = KPTS(1:3,kcounter)

			! This loop sets all initial values of H to zero, so that the sum afterwords can work
            do i = 1, 4*NUMT
                do j = 1, 4*NUMT
                    HAMILTONIAN(i,j) = (0.0,0.0)
                end do
            end do

			! This begins the construction of H(k) for the specific k-point
            do i = 1, NUMT
                do j = 1, NUMT
                    
                    TTPRIME = TPTS(1:3,j) - TPTS(1:3,i) ! This calculates t - t' for the basis atoms
                    do IRLATT = 1, IRLATTMAX ! This is the summation over all lattice points
                        RPOINT = RLATT(1:3,IRLATT)
                        call HAM(zPauli,IdentityPauli,chempot,RPOINT + TTPRIME,NUMT,E0,j,&
                        &R0,RMAX,ULCN,nu,nuzero,BETA,DELTA,positionhamiltonian,deltaterm) ! Calls h(r-r')
                        expon = exp(-CI*DOT_PRODUCT(KPOINT,RPOINT)) ! The e^(-ikr) factor

						! This is practically the Fourier transform of h(r-r') to H(k)

                        HAMILTONIAN(i,j) = HAMILTONIAN(i,j) + expon*positionhamiltonian(1,1) ! 1-1 block
                        HAMILTONIAN(i,j+NUMT) = HAMILTONIAN(i,j+NUMT) + expon*positionhamiltonian(1,2) ! 1-2 block, for spinup-spindown interactions
						! HAMILTONIAN(i,j+2*NUMT), i.e. the 1-3 block, remains zero in any case
                        HAMILTONIAN(i,j+3*NUMT) = HAMILTONIAN(i,j+3*NUMT) + expon*deltaterm ! 1-4 block, the superconductivity pairing

                        HAMILTONIAN(i+NUMT,j) = HAMILTONIAN(i+NUMT,j) + expon*positionhamiltonian(2,1) ! 2-1 block, for spindown-spinup interactions
                        HAMILTONIAN(i+NUMT,j+NUMT) = HAMILTONIAN(i+NUMT,j+NUMT) + expon*positionhamiltonian(2,2) ! 2-2 block
                        HAMILTONIAN(i+NUMT,j+2*NUMT) = HAMILTONIAN(i+NUMT,j+2*NUMT) + expon*deltaterm ! 2-3 block, the superconductivity pairing
						! HAMILTONIAN(i+NUMT,j+3*NUMT), i.e. the 2-4 block, remains zero in any case

						! HAMILTONIAN(i+2*NUMT,j), i.e. the 3-1 block, remains zero in any case
                        HAMILTONIAN(i+2*NUMT,j+NUMT) = HAMILTONIAN(i+2*NUMT,j+NUMT) + expon*CONJG(deltaterm) ! 3-2 block, the complex conjugate of the superconductivity pairing
                        HAMILTONIAN(i+2*NUMT,j+2*NUMT) = HAMILTONIAN(i+2*NUMT,j+2*NUMT) - expon*CONJG(positionhamiltonian(1,1)) ! 3-3 block
                        HAMILTONIAN(i+2*NUMT,j+3*NUMT) = HAMILTONIAN(i+2*NUMT,j+3*NUMT) + expon*CONJG(positionhamiltonian(1,2)) ! 3-4 block, for spinup-spindown interactions

                        HAMILTONIAN(i+3*NUMT,j) = HAMILTONIAN(i+3*NUMT,j) + expon*CONJG(deltaterm) ! 4-1 block, the complex conjugate of the superconductivity pairing
						! HAMILTONIAN(i+3*NUMT,j+NUMT), i.e. the 4-2 block, remains zero in any case
                        HAMILTONIAN(i+3*NUMT,j+2*NUMT) = HAMILTONIAN(i+3*NUMT,j+2*NUMT) + expon*CONJG(positionhamiltonian(2,1)) ! 4-3 block, for spindown-spinup interactions
                        HAMILTONIAN(i+3*NUMT,j+3*NUMT) = HAMILTONIAN(i+3*NUMT,j+3*NUMT) - expon*CONJG(positionhamiltonian(2,2)) ! 4-4 block
                    end do
                end do
            end do

            call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO) ! Don't forget to reconfigure those whenever the dimensions change!

			! The eigenvectors are in the form of 4-spinors: (u-up, u-down, v-up, v-down)

            ini = (kcounter-1)*4*NUMT + 1 ! The *4 factors are due spin and particle-hole
            fin = kcounter*4*NUMT ! The *4 factors are due spin and particle-hole
            EIGENVALUES(ini:fin) = W
            EIGENVECTORS(:,ini:fin) = HAMILTONIAN
        end do

		! At that point all the eigenvalues are in the form (.,.,.,.,...) and all the eigenvectors are 4*NUMT*NUMK columns of 4*NUMT rows

		! Calculates the charges n-up, n-down and n as well as Delta
        do i = 1, NUMT
            nuup(i) = 0.0
            nudown(i) = 0.0
            newDELTA(i) = (0.0,0.0)

            do j = 1, 4*NUMT*NUMK
                if (EIGENVALUES(j) > chempot) then ! This is now the "Fermi condition"
                    nuup(i) = nuup(i) + abs(EIGENVECTORS(2*NUMT+i,j))**2

                    nudown(i) = nudown(i) + abs(EIGENVECTORS(3*NUMT+i,j))**2

                    newDELTA(i) = newDELTA(i) + 0.5*USUPCOND(i)*( EIGENVECTORS(i,j)*CONJG(EIGENVECTORS(3*NUMT+i,j)) +&
					& EIGENVECTORS(NUMT+i,j)*CONJG(EIGENVECTORS(2*NUMT+i,j)) )
                endif
            end do

            newnu(i) = nuup(i) + nudown(i) ! This is the density of the i-th atom
            newnu(i) = newnu(i)/(NUMK) ! Since the eigenvectors come normalized, this ensures the normalization of n
            diffN(i) = abs(newnu(i) - nu(i))
			
            newDELTA(i) = newDELTA(i)/NUMK ! Similarly, the normalization of D
            diffD(i) = abs(abs(newDELTA(i)) - abs(DELTA(i)))

        end do
		
        DELTA = (1.0 - mixfactorD)*DELTA + mixfactorD*newDELTA
        nu = (1.0 - mixfactorN)*nu + mixfactorN*newnu
        reps = reps + 1
    end do
    print *, 'Charge densities evaluated.'

	!At that point we have a good approximation for the charges. We move on to calculate once again the Hamiltonian and then the states density.    
    do kcounter = 1, NUMK
        KPOINT = KPTS(1:3,kcounter)

		! This loop sets all initial values of H to zero, so that the sum afterwords can work
        do i = 1, 4*NUMT
            do j = 1, 4*NUMT
                HAMILTONIAN(i,j) = (0.0,0.0)
            end do
        end do

        do i = 1, NUMT
            do j = 1, NUMT
				
                TTPRIME = TPTS(1:3,j) - TPTS(1:3,i) ! This calculates t - t' for the basis atoms
                do IRLATT = 1, IRLATTMAX ! This is the summation over all lattice points
                    RPOINT = RLATT(1:3,IRLATT)
                    call HAM(zPauli,IdentityPauli,chempot,RPOINT + TTPRIME,NUMT,E0,j,&
                    &R0,RMAX,ULCN,nu,nuzero,BETA,DELTA,positionhamiltonian,deltaterm) ! Calls h(r-r')
                    expon = exp(-CI*DOT_PRODUCT(KPOINT,RPOINT)) ! The e^(-ikr) factor

					! This is practically the Fourier transform of h(r-r') to H(k)

                    HAMILTONIAN(i,j) = HAMILTONIAN(i,j) + expon*positionhamiltonian(1,1) ! 1-1 block
                    HAMILTONIAN(i,j+NUMT) = HAMILTONIAN(i,j+NUMT) + expon*positionhamiltonian(1,2) ! 1-2 block, for spinup-spindown interactions
                    ! HAMILTONIAN(i,j+2*NUMT), i.e. the 1-3 block, remains zero in any case
                    HAMILTONIAN(i,j+3*NUMT) = HAMILTONIAN(i,j+3*NUMT) + expon*deltaterm ! 1-4 block, the superconductivity pairing

                    HAMILTONIAN(i+NUMT,j) = HAMILTONIAN(i+NUMT,j) + expon*positionhamiltonian(2,1) ! 2-1 block, for spindown-spinup interactions
                    HAMILTONIAN(i+NUMT,j+NUMT) = HAMILTONIAN(i+NUMT,j+NUMT) + expon*positionhamiltonian(2,2) ! 2-2 block
                    HAMILTONIAN(i+NUMT,j+2*NUMT) = HAMILTONIAN(i+NUMT,j+2*NUMT) + expon*deltaterm ! 2-3 block, the superconductivity pairing
                    ! HAMILTONIAN(i+NUMT,j+3*NUMT), i.e. the 2-4 block, remains zero in any case

                    ! HAMILTONIAN(i+2*NUMT,j), i.e. the 3-1 block, remains zero in any case
                    HAMILTONIAN(i+2*NUMT,j+NUMT) = HAMILTONIAN(i+2*NUMT,j+NUMT) + expon*CONJG(deltaterm) ! 3-2 block, the complex conjugate of the superconductivity pairing
                    HAMILTONIAN(i+2*NUMT,j+2*NUMT) = HAMILTONIAN(i+2*NUMT,j+2*NUMT) - expon*CONJG(positionhamiltonian(1,1)) ! 3-3 block
                    HAMILTONIAN(i+2*NUMT,j+3*NUMT) = HAMILTONIAN(i+2*NUMT,j+3*NUMT) + expon*CONJG(positionhamiltonian(1,2)) ! 3-4 block, for spinup-spindown interactions

                    HAMILTONIAN(i+3*NUMT,j) = HAMILTONIAN(i+3*NUMT,j) + expon*CONJG(deltaterm) ! 4-1 block, the complex conjugate of the superconductivity pairing
                    ! HAMILTONIAN(i+3*NUMT,j+NUMT), i.e. the 4-2 block, remains zero in any case
                    HAMILTONIAN(i+3*NUMT,j+2*NUMT) = HAMILTONIAN(i+3*NUMT,j+2*NUMT) + expon*CONJG(positionhamiltonian(2,1)) ! 4-3 block, for spindown-spinup interactions
                    HAMILTONIAN(i+3*NUMT,j+3*NUMT) = HAMILTONIAN(i+3*NUMT,j+3*NUMT) - expon*CONJG(positionhamiltonian(2,2)) ! 4-4 block
                end do
            end do
        end do

        call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO) ! Don't forget to reconfigure those whenever the dimensions change!

        ini = (kcounter-1)*4*NUMT + 1 ! The *4 factors are due spin and particle-hole
        fin = kcounter*4*NUMT ! The *4 factors are due spin and particle-hole
        EIGENVALUES(ini:fin) = W
        EIGENVECTORS(:,ini:fin) = HAMILTONIAN
        
    end do

	! At that point, we calculate the density, magnetization and D for the final Hamiltonian using the converged values.
    do i = 1, NUMT
        nuup(i) = 0.0
        nudown(i) = 0.0
        newDELTA(i) = (0.0,0.0)
        magnet(i) = 0.0

        do j = 1, 4*NUMT*NUMK
            if (EIGENVALUES(j) > chempot) then
                nuup(i) = nuup(i) + abs(EIGENVECTORS(2*NUMT+i,j))**2
                    
                nudown(i) = nudown(i) + abs(EIGENVECTORS(3*NUMT+i,j))**2

                newDELTA(i) = newDELTA(i) + 0.5*USUPCOND(i)*( EIGENVECTORS(i,j)*CONJG(EIGENVECTORS(3*NUMT+i,j)) +&
				& EIGENVECTORS(NUMT+i,j)*CONJG(EIGENVECTORS(2*NUMT+i,j)) )
            endif
        end do

        newnu(i) = (nuup(i) + nudown(i))/NUMK ! Final Density per atom
        newDELTA(i) = newDELTA(i)/NUMK ! Final D per atom
        magnet(i) = (nuup(i) - nudown(i))/NUMK ! Final magnetization

        print *, 'n = ', newnu(i)
        print *, 'D = ', newDELTA(i)
        print *, 'M = ', magnet(i)
    end do

    !open(16, file = 'name.txt', action = 'write')
		!do j = 1, NUMT
			!write (16,100) (ARRAY(i,j), i = 1,2)
		!end do
		!100 format(3F17.8)
	!close(16)

	!This takes the DIFFERENT eigenvalues and organizes them in ascending order
	!---------------------------------------------------------------------------------------------
    allocate(UNIQUEEIGVALS(4*NUMT*NUMK))
    uniquecounter = 0
    min_val = MINVAL(EIGENVALUES) - 1.0 ! -1.0 Is inserted for a case of complete degeneracy
    max_val = MAXVAL(EIGENVALUES)
    do while (min_val < max_val)
        uniquecounter = uniquecounter + 1
        min_val = MINVAL(EIGENVALUES, mask = EIGENVALUES > min_val)
        UNIQUEEIGVALS(uniquecounter) = min_val
    enddo
    allocate(SORTEDEIGVALS(uniquecounter))
    SORTEDEIGVALS = UNIQUEEIGVALS(1:uniquecounter)
	!---------------------------------------------------------------------------------------------

    call INT_NUM_DEN(uniquecounter,EIGENVALUES,NUMT,NUMK,SORTEDEIGVALS,multiplicity,intnumdensity)
	
    call NUM_DEN(uniquecounter,SORTEDEIGVALS,EIGENVALUES,EIGENVECTORS,numdensityperatom,numdensity,NUMT,NUMK)

	!------------------------------------------------------------------------------------------------------------------

    contains

	! Now, after the addition of spin, the result becomes a 2x2 matrix and the previous function is now a subroutine.
	! After the addition of particle-hole duality, we don't need to increase its dimensions once again, we will only
	! manually handle how the 4*NUMT*4*NUMT Hamiltonian comes out. If we want a non-hand-fixed result, we should
	! introduce 4x4 matrices, which are direct products of Pauli matrices in spin and particle-hole space.

    subroutine HAM(zPauli,IdentityPauli,chempot,TTPRIME,NUMT,E0,j,R0,RMAX,ULCN,nu,nuzero,BETA,DELTA,positionhamiltonian,deltaterm)
        implicit none
        real*8, intent(in) :: TTPRIME(3)
        real*8 :: R0, RMAX, chempot
        integer :: NUMT, j
        real*8 :: E0(NUMT), ULCN(NUMT), nu(NUMT), nuzero(NUMT), BETA(NUMT)
        complex*16 :: zPauli(2,2),IdentityPauli(2,2), positionhamiltonian(2,2), deltaterm, DELTA(NUMT)

        if (norm2(TTPRIME) == 0.0) then ! This case corresponds to t = t', R = 0
            positionhamiltonian = (E0(j) - chempot + ULCN(j)*(nu(j) - nuzero(j)))*IdentityPauli - BETA(j)*zPauli
            deltaterm = DELTA(j) ! This ensures on-site superconducting pairing (s-wave superconductivity)
        else if (norm2(TTPRIME) < RMAX) then
            positionhamiltonian = (exp((-1)*norm2(TTPRIME)/R0))*IdentityPauli
            deltaterm = (0.0,0.0) ! This ensures on-site superconducting pairing (s-wave superconductivity)
        else
            positionhamiltonian = (0.0,0.0)
            deltaterm = (0.0,0.0) ! This ensures on-site superconducting pairing (s-wave superconductivity)
        endif
    end subroutine HAM

    subroutine INT_NUM_DEN(uniquecounter,EIGENVALUES,NUMT,NUMK,SORTEDEIGVALS,multiplicity,intnumdensity)
        implicit none

        integer, allocatable, dimension (:) :: multiplicity
        real*8, allocatable, dimension(:,:) :: intnumdensity
        integer :: uniquecounter, NUMT, NUMK, i, j
        real*8 :: EIGENVALUES(4*NUMT*NUMK), SORTEDEIGVALS(uniquecounter)

        allocate(multiplicity(uniquecounter))
        allocate(intnumdensity(2,uniquecounter))
        intnumdensity(1,:) = SORTEDEIGVALS

        do i = 1, uniquecounter
            multiplicity(i) = 0
            do j = 1, 4*NUMT*NUMK
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

    end subroutine INT_NUM_DEN

    ! --------------------------------------------------------------------------------------------------------------------------------
    subroutine NUM_DEN(uniquecounter,SORTEDEIGVALS,EIGENVALUES,EIGENVECTORS,numdensityperatom,numdensity,NUMT,NUMK)
        implicit none

        integer :: uniquecounter, i, NUMT, NUMK, numenergyintervals,j,k
        real*8, allocatable, dimension(:,:) :: numdensityperatom, numdensity
        real*8 :: pi, delta, SORTEDEIGVALS(uniquecounter), EIGENVALUES(4*NUMT*NUMK), energyintervals
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT*NUMK)

        pi = 4.D0*atan(1.D0)
        print *, 'Please insert the lorentzian delta factor for the number density.'
        read *, delta

        print *, 'Please enter the number of intervals for the plot of n(E).'
        read *, numenergyintervals

        allocate(numdensityperatom(1+NUMT,numenergyintervals+1))
        allocate(numdensity(2,numenergyintervals+1))

        energyintervals = (MAXVAL(SORTEDEIGVALS) - MINVAL(SORTEDEIGVALS))/numenergyintervals

        ! These are the energies E
        do i = 0, numenergyintervals
            numdensityperatom(1,i+1) = MINVAL(SORTEDEIGVALS) + energyintervals*i
            numdensity(1,i+1) = MINVAL(SORTEDEIGVALS) + energyintervals*i
            numdensity(2,i+1) = 0.0
        end do

        ! Calculation of the DoS PER ATOM
        do k = 1, NUMT
            do i = 0, numenergyintervals
                numdensityperatom(1+k,i+1) = 0.0
                do j = 1, 4*NUMT*NUMK
                    if (EIGENVALUES(j) > chempot) then
                        numdensityperatom(1+k,i+1) = numdensityperatom(1+k,i+1) + (delta/pi)*((abs(EIGENVECTORS(k,j))**2 +&
                        & abs(EIGENVECTORS(k+NUMT,j))**2)/((numdensityperatom(1,i+1) - EIGENVALUES(j))**2 + delta**2) +&
                        & (abs(EIGENVECTORS(k+2*NUMT,j))**2 + abs(EIGENVECTORS(k+3*NUMT,j))**2)/((numdensityperatom(1,i+1) +&
                        & EIGENVALUES(j))**2 + delta**2))
                    else if (EIGENVALUES(j) == chempot) then
                        ! the 0.5 is inserted to avoid double counting in the absence of superconducting effects for zero eigenvalues
                        numdensityperatom(1+k,i+1) = numdensityperatom(1+k,i+1) + 0.5*(delta/pi)*((abs(EIGENVECTORS(k,j))**2 +&
                        & abs(EIGENVECTORS(k+NUMT,j))**2)/((numdensityperatom(1,i+1) - EIGENVALUES(j))**2 + delta**2) +&
                        & (abs(EIGENVECTORS(k+2*NUMT,j))**2 + abs(EIGENVECTORS(k+3*NUMT,j))**2)/((numdensityperatom(1,i+1) +&
                        & EIGENVALUES(j))**2 + delta**2))
                    endif
                end do
            end do
        end do

        open(16, file = 'numdensityperatom.txt', action = 'write')
        do j = 1, numenergyintervals
            write (16,100) (numdensityperatom(i,j), i = 1,1+NUMT)
        end do
        100 format(41F17.8)
        close(16)

        ! Calculation of the full density of states
        do k = 1, NUMT
            numdensity(2,:) = numdensity(2,:) + numdensityperatom(1+k,:)
        end do

        open(15, file = 'numdensity.txt', action = 'write')
        do j = 1, numenergyintervals
            write (15,102) (numdensity(i,j), i = 1,2)
        end do
        102 format(3F17.8)
        close(15)

    end subroutine NUM_DEN
    ! --------------------------------------------------------------------------------------------------------------------------------
    
    ! ---------------------------- These have not changed after the addition of spin or particle-hole duality

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

    function CROSS_PRODUCT(x,y) result(cross)
        implicit none
        real*8, dimension(3), intent(in) :: x, y
        real*8, dimension(3) :: cross

        cross(1) = x(2)*y(3) - x(3)*y(2)
        cross(2) = x(3)*y(1) - x(1)*y(3)
        cross(3) = x(1)*y(2) - x(2)*y(1)
		
    end function CROSS_PRODUCT

end program TB