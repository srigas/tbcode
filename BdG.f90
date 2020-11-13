program TB
    implicit none
    integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, LWORK, INFO, NCELLS, ini, fin, reps, &
    & maxreps, uniquecounter, NUMIMP, imppointer, NUMCHEMTYPES, bandpointer, intnumdenpointer, NUME, metalorno, JATOM, rcheck,&
    & MAXNEIGHB, FTIMO
    integer, allocatable, dimension(:) :: multiplicity, CHEMTYPE, NEIGHBNUM
    integer, allocatable, dimension(:,:) :: IMPPTSVAR, JTYPE
    real*8 :: ALAT, a_1(3), a_2(3), a_3(3), RMAX, R0, KPOINT(3), epsilon, min_val, max_val, mixfactorN, RPOINT(3), TTPRIME(3),&
	&chempot, mixfactorD, inicharge, T, PI, KB, b_1(3), b_2(3), b_3(3), DETCHECK, lorentzbroad, diffchem, newchempot, lambda, TOL
    real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, nu, newnu, nuzero, EIGENVALUES, SORTEDEIGVALS, &
	& UNIQUEEIGVALS, magnet, VSUPCOND, nuup, nudown, diffN, diffD, DOSATMU
    real*8, allocatable, dimension(:,:) :: KPTS, TPTS, RLATT, BETA, LHOPS, PREFACTORS, NNDIS, HOPPVALS
    real*8, allocatable, dimension(:,:,:) :: RCONNECT
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), inidelta
    complex*16, allocatable, dimension(:) :: WORK, DELTA, newDELTA, METALDELTA
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, EIGENVECTORS, HAMILTONIANPREP
    complex*16, allocatable, dimension(:,:,:) :: EXPONS

    ! Important notice about write format. Example: format(5F17.8)
    ! 5F -> 5 = number of entries before changing row, F is because we want numbers
    ! 17 = Total digits of entry. E.g. if entry = PI, then 3. are the first two digits, therefore 17-2=15 are remaining.
    ! Of these 15, the 8 (that's what the .8 stands for) are reserved as decimals. The 15-8=7 remaining are spaces before the next entry.

    TOL = 0.0001 ! The fault tolerance for lattice vectors' norms
    
    call CONSTANTS(IdentityPauli,xPauli,yPauli,zPauli,CI,PI,KB) ! Sets some universal constants.

    open(1, file = 'config.dat', action = 'read')
    read(1,*) ALAT
    read(1,*) a_1
    read(1,*) a_2
    read(1,*) a_3
    read(1,*) NUMKX,NUMKY,NUMKZ
    read(1,*) RMAX ! To determine nearest neighbours
    read(1,*) R0 ! Constant at the exponential of the hopping element
    read(1,*) NCELLS ! NCELLS creates a (2NCELLS+1)^3 mini-cube from there the lattice points are used to determine neighbours
    read(1,*) T ! the system's temperature. Default: 0.0
    read(1,*) chempot ! initial value for chemical potential
    read(1,*) inicharge ! initial value for charges
    read(1,*) inidelta ! initial value for Delta
    read(1,*) mixfactorN ! mixing factor for N
    read(1,*) mixfactorD ! mixing factor for Delta
    read(1,*) NUME ! Number of points for DoS diagrams
    read(1,*) lorentzbroad ! Lorentz Gamma broadening
    close(1)

    DETCHECK = a_1(1)*(a_2(2)*a_3(3)-a_2(3)*a_3(2)) - a_2(1)*(a_1(2)*a_3(3)-a_3(2)*a_1(3)) + &
    &a_3(1)*(a_1(2)*a_2(3)-a_1(3)*a_2(2))

    if (DETCHECK <= 0.0001 .or. RMAX < 0 .or. R0 <= 0 .or. NCELLS < 0 .or. T < 0.0) then
        print *, 'A value inserted in config.dat is incorrect. Please try again after everything has been corrected.'
        call exit(123)
    endif

    ! This is because the lattice vectors are inserted through Bravais coordinates
    a_1 = ALAT*a_1
    a_2 = ALAT*a_2
    a_3 = ALAT*a_3

    IRLATTMAX = (2*NCELLS+1)**3 ! Configures how many neighbouring cells are taken into account

    call BZ(a_1,a_2,a_3,NUMKX,NUMKY,NUMKZ,PI,KPTS,NUMK,b_1,b_2,b_3) ! We create the Brillouin Zone's k-points to be used later on
	
	! The following reads the number of basis vectors from a file named basisvectors.dat
    NUMT = 0
    open (2, file = 'basisvectors.dat', action = 'read')
    do
        read(2,*,iostat=io)
        if (io/=0) exit
        NUMT = NUMT + 1
    end do
    close(2)

    NUMT = NUMT - 2 ! To ignore the final 2 lines in basisvectors.dat which are the configuration settings.

    allocate(TPTS(3,NUMT))
    allocate(E0(NUMT))
    allocate(ULCN(NUMT)) ! NUMT x 1 column with the U_LCN constant for each basis atom
    allocate(nuzero(NUMT)) ! NUMT x 1 column with the charges n_0 of each basis atom
    allocate(nu(NUMT)) ! NUMT x 1 column with the charges n of each basis atom
    allocate(newnu(NUMT)) ! NUMT x 1 column with the charges n of each basis atom
    allocate(BETA(3,NUMT)) ! NUMT x 3 matrix with the B_x, B_y, B_z for each basis atom
    allocate(diffN(NUMT))
    allocate(diffD(NUMT))
    allocate(magnet(NUMT))
    allocate(DELTA(NUMT))
    allocate(METALDELTA(NUMT))
    allocate(DOSATMU(NUMT))
    allocate(newDELTA(NUMT))
    allocate(VSUPCOND(NUMT))
    allocate(nuup(NUMT))
    allocate(nudown(NUMT))
    allocate(CHEMTYPE(NUMT)) ! Disciminates between different chemical elements

    open (1, file = 'basisvectors.dat', action = 'read')
    do i = 1,NUMT
        read(1,*) TPTS(1:3,i), CHEMTYPE(i), E0(i), ULCN(i), nuzero(i), BETA(1,i), BETA(2,i), BETA(3,i), VSUPCOND(i)
    end do
    close(1)

    do i = 1, NUMT
        nu(i) = inicharge
        DELTA(i) = inidelta
        METALDELTA(i) = (0.0,0.0)
    end do

    NUMCHEMTYPES = MAXVAL(CHEMTYPE)
    allocate(LHOPS(NUMCHEMTYPES,NUMCHEMTYPES))
    allocate(PREFACTORS(NUMCHEMTYPES,NUMCHEMTYPES))
    allocate(NNDIS(NUMCHEMTYPES,NUMCHEMTYPES))

    ! The *4 factors are now due to spin and particle-hole
    allocate(EIGENVECTORS(4*NUMT,4*NUMT*NUMK))
    allocate(EIGENVALUES(4*NUMT*NUMK))
    allocate(HAMILTONIAN(4*NUMT,4*NUMT))
    allocate(HAMILTONIANPREP(4*NUMT,4*NUMT))

    call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) ! Here we construct the RLATT Matrix consisting of the lattice sites

    call HOPPS(RLATT,IRLATTMAX,R0,NUMCHEMTYPES,LHOPS,NUMT,CHEMTYPE,TPTS,PREFACTORS)
	
	!The following are for the configuration of zheev. The *4 factors are due to spin and particle-hole
	!-------------------------------------------------
    allocate(W(4*NUMT))
    allocate(RWORK(3*(4*NUMT) - 2))
    LWORK = 4*NUMT*(4*NUMT+1)
    allocate(WORK(LWORK))
	!-------------------------------------------------

    ! That part calculates the number of neighbours that belong in the ith atom's cluster, excluding itself
    allocate(NEIGHBNUM(NUMT))
    NEIGHBNUM(1:NUMT) = 0

    do i = 1, NUMT
        do j = 1, NUMT

            TTPRIME = TPTS(1:3,j) - TPTS(1:3,i)

            do IRLATT = 1, IRLATTMAX
                RPOINT = RLATT(1:3,IRLATT)

                if (norm2(RPOINT+TTPRIME) < (RMAX+TOL) .and. norm2(RPOINT+TTPRIME) > TOL) then
                    NEIGHBNUM(i) = NEIGHBNUM(i) + 1
                end if
            end do
                
        end do
    end do

    ! This part catalogues all the nearest neighbour distances for each atom i, while also listing
    ! what type of neighbours they are, excluding the on-site interaction
    MAXNEIGHB = MAXVAL(NEIGHBNUM)
    allocate(JTYPE(MAXNEIGHB,NUMT))
    allocate(RCONNECT(3,MAXNEIGHB,NUMT))

    RCONNECT(1:3,:,:) = 0.0
    JTYPE(:,:) = 0
    
    do i = 1, NUMT
        rcheck = 1
        do j = 1, NUMT

            TTPRIME = TPTS(1:3,j) - TPTS(1:3,i)

            do IRLATT = 1, IRLATTMAX
                RPOINT = RLATT(1:3,IRLATT)

                if (norm2(RPOINT+TTPRIME) < (RMAX+TOL) .and. norm2(RPOINT+TTPRIME) > TOL) then
                    RCONNECT(1:3,rcheck,i) = RPOINT + TTPRIME
                    JTYPE(rcheck,i) = j
                    rcheck = rcheck + 1
                end if
            end do
                
        end do
    end do

    ! We finally calculate the hopping elements, as well as the fourier exponentials, so that we don't have to calculate
    ! them after every self-consistency cycle.
    allocate(EXPONS(NUMK,MAXNEIGHB,NUMT))
    allocate(HOPPVALS(MAXNEIGHB,NUMT))
    do kcounter = 1, NUMK
        KPOINT = KPTS(1:3,kcounter)

        do i = 1, NUMT
            do j = 1, NEIGHBNUM(i)

                RPOINT = RCONNECT(1:3,j,i)
                EXPONS(kcounter,j,i) = exp(-CI*DOT_PRODUCT(KPOINT,RPOINT))

                JATOM = JTYPE(j,i)
                lambda = PREFACTORS(CHEMTYPE(i),CHEMTYPE(JATOM))
                
                HOPPVALS(j,i) = -lambda*exp(-norm2(RPOINT)/R0)

            end do
        end do
    end do

	! These configurations ensure that the following while loop is initiated
    epsilon = 0.001
    metalorno = 1
    diffchem = 1.0
    reps = 0
    do i = 1, NUMT
        diffN(i) = 1.0
    end do

    print *, 'Please insert the self-consistency cycles for the metal calculations.'
    read *, maxreps

    print *, 'Initiating METAL self-consistency procedure...'
    do while ((MAXVAL(diffN) > epsilon .or. diffchem > epsilon) .and. reps < maxreps)

        call HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,METALDELTA,HAMILTONIANPREP)

        do kcounter = 1, NUMK

            call FOURIERHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)

            call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO)

            ini = (kcounter-1)*4*NUMT + 1
            fin = kcounter*4*NUMT
            EIGENVALUES(ini:fin) = W
            EIGENVECTORS(:,ini:fin) = HAMILTONIAN
        end do

        newchempot = chempot
        do i = 1, NUMT
            nuup(i) = 0.0
            nudown(i) = 0.0
            DOSATMU(i) = 0.0

            FTIMO = 4*(i-1)

            do j = 1, 4*NUMT*NUMK
                nuup(i) = nuup(i) + FERMI(EIGENVALUES(j),T,KB)*abs(EIGENVECTORS(1+FTIMO,j))**2
                nudown(i) = nudown(i) + FERMI(-EIGENVALUES(j),T,KB)*abs(EIGENVECTORS(4+FTIMO,j))**2

                DOSATMU(i) = DOSATMU(i) + (lorentzbroad/pi)*(1.0/NUMK)*((abs(EIGENVECTORS(1+FTIMO,j))**2 +&
                &abs(EIGENVECTORS(2+FTIMO,j))**2)/((chempot - EIGENVALUES(j))**2 + lorentzbroad**2))
            end do

            newnu(i) = (nuup(i) + nudown(i))/NUMK ! Normalization
            diffN(i) = abs(newnu(i) - nu(i))

            newchempot = newchempot - (newnu(i)-nuzero(i))*DOSATMU(i)
        end do
        diffchem = abs(chempot-newchempot)

        nu = (1.0 - mixfactorN)*nu + mixfactorN*newnu
        chempot = newchempot ! new chempot for next run
        reps = reps + 1

    end do
    print *, 'METAL self-consistency finished.'

    ! We now perform a DoS calculation for the metal, in order to then be able to compare it to the SC.
    call NUM_DEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,NUME,lorentzbroad,metalorno)

    ! Using the chemical potential that we obtained, as well as the charges (nu), as initial values, we proceed to use them for the
    ! self-consistency cycle of the superconductor, i.e. Delta =/= 0.
    metalorno = 0
    reps = 0
    diffchem = 1.0
    do i = 1, NUMT
        diffN(i) = 1.0
        diffD(i) = 1.0
    end do

    print *, 'Please insert the self-consistency cycles for the superconductor calculations.'
    read *, maxreps

    print *, 'Initiating SC self-consistency procedure...'
    do while ((MAXVAL(diffN) > epsilon .or. MAXVAL(diffD) > epsilon .or. diffchem > epsilon) .and. reps < maxreps) ! Check for convergence or maxreps

        call HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,DELTA,HAMILTONIANPREP)

        do kcounter = 1, NUMK

            call FOURIERHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)

            call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO)

            ini = (kcounter-1)*4*NUMT + 1
            fin = kcounter*4*NUMT
            EIGENVALUES(ini:fin) = W
            EIGENVECTORS(:,ini:fin) = HAMILTONIAN
        end do

		! At that point all the eigenvalues are in the form (.,.,.,.,...) and all the eigenvectors are 4*NUMT*NUMK columns of 4*NUMT rows

        newchempot = chempot
		! Calculates the charges n-up, n-down and n as well as Delta and chempot
        do i = 1, NUMT
            nuup(i) = 0.0
            nudown(i) = 0.0
            newDELTA(i) = (0.0,0.0)
            DOSATMU(i) = 0.0

            FTIMO = 4*(i-1)

            do j = 1, 4*NUMT*NUMK
                nuup(i) = nuup(i) + FERMI(EIGENVALUES(j),T,KB)*abs(EIGENVECTORS(1+FTIMO,j))**2

                nudown(i) = nudown(i) + FERMI(-EIGENVALUES(j),T,KB)*abs(EIGENVECTORS(4+FTIMO,j))**2

                newDELTA(i) = newDELTA(i) -&
                & FERMI(EIGENVALUES(j),T,KB)*VSUPCOND(i)*EIGENVECTORS(1+FTIMO,j)*CONJG(EIGENVECTORS(4+FTIMO,j))

                DOSATMU(i) = DOSATMU(i) + (lorentzbroad/pi)*(1.0/NUMK)*((abs(EIGENVECTORS(1+FTIMO,j))**2 +&
                &abs(EIGENVECTORS(2+FTIMO,j))**2)/((chempot - EIGENVALUES(j))**2 + lorentzbroad**2))
            end do

            newnu(i) = (nuup(i) + nudown(i))/NUMK ! Normalization
            diffN(i) = abs(newnu(i) - nu(i))
			
            newDELTA(i) = newDELTA(i)/NUMK ! Normalization
            diffD(i) = abs(abs(newDELTA(i)) - abs(DELTA(i)))

            newchempot = newchempot - (newnu(i)-nuzero(i))*DOSATMU(i)
        end do
        diffchem = abs(chempot-newchempot)

        DELTA = (1.0 - mixfactorD)*DELTA + mixfactorD*newDELTA
        nu = (1.0 - mixfactorN)*nu + mixfactorN*newnu
        chempot = newchempot ! new chempot for next run
        reps = reps + 1
    end do
    print *, 'SC self-consistency finished.'

	!At that point we have a good approximation for the charges. We move on to calculate once again the Hamiltonian and then the states density.    
    call HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,DELTA,HAMILTONIANPREP)

    do kcounter = 1, NUMK

        call FOURIERHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)

        call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO)

        ini = (kcounter-1)*4*NUMT + 1
        fin = kcounter*4*NUMT
        EIGENVALUES(ini:fin) = W
        EIGENVECTORS(:,ini:fin) = HAMILTONIAN
    end do

	! At that point, we calculate the density, magnetization and D for the final Hamiltonian using the converged values.
    newchempot = chempot
    do i = 1, NUMT
        nuup(i) = 0.0
        nudown(i) = 0.0
        newDELTA(i) = (0.0,0.0)
        magnet(i) = 0.0
        DOSATMU(i) = 0.0

        FTIMO = 4*(i-1)

        do j = 1, 4*NUMT*NUMK
            nuup(i) = nuup(i) + FERMI(EIGENVALUES(j),T,KB)*abs(EIGENVECTORS(1+FTIMO,j))**2

            nudown(i) = nudown(i) + FERMI(-EIGENVALUES(j),T,KB)*abs(EIGENVECTORS(4+FTIMO,j))**2

            newDELTA(i) = newDELTA(i) -&
            & FERMI(EIGENVALUES(j),T,KB)*VSUPCOND(i)*EIGENVECTORS(1+FTIMO,j)*CONJG(EIGENVECTORS(4+FTIMO,j))

            DOSATMU(i) = DOSATMU(i) + (lorentzbroad/pi)*(1.0/NUMK)*((abs(EIGENVECTORS(1+FTIMO,j))**2 +&
            &abs(EIGENVECTORS(2+FTIMO,j))**2)/((chempot - EIGENVALUES(j))**2 + lorentzbroad**2))
        end do

        newnu(i) = (nuup(i) + nudown(i))/NUMK ! Final Density per atom
        newDELTA(i) = newDELTA(i)/NUMK ! Final D per atom
        magnet(i) = (nuup(i) - nudown(i))/NUMK ! Final magnetization
        newchempot = newchempot - (newnu(i)-nuzero(i))*DOSATMU(i)

        print *, 'n = ', newnu(i), 'for atom No. ', i
        print *, 'D = ', newDELTA(i), 'for atom No. ', i
        print *, 'M = ', magnet(i), 'for atom No. ', i
    end do
    
    print *, 'chempot = ', newchempot

    ! Routine for making band diagrams
    print *, 'To export band diagram data press 0, otherwise press any other number.'
    read *, bandpointer
    if (bandpointer == 0) then
        call BANDS(NUMT,W,WORK,LWORK,RWORK,INFO,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0, &
        &ULCN,nu,nuzero,BETA,DELTA,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,RCONNECT,ALAT)
    endif

	!This takes the DIFFERENT eigenvalues and organizes them in ascending order
	!---------------------------------------------------------------------------------------------
    print *, 'To calculate integrated number density press 0, otherwise press any other number.'
    read *, intnumdenpointer
    if (intnumdenpointer == 0) then
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

        call INT_NUM_DEN(uniquecounter,EIGENVALUES,NUMT,NUMK,SORTEDEIGVALS,multiplicity)
    endif
	!---------------------------------------------------------------------------------------------

    call NUM_DEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,NUME,lorentzbroad,metalorno)

	!------------------------------------------------------------------------------------------------------------------

    print *, 'Choose the input way for the impurity sites.'
    print *, '0 to input via vectors, 1 to input via integers, other number to skip Greens functions calculations.'
    read *, imppointer

    if (imppointer == 0 .or. imppointer == 1) then
        call IDENTIFIER(imppointer,b_1,b_2,b_3,PI,a_1,a_2,a_3,TPTS,NUMT,NUMIMP,IMPPTSVAR)

        print *, 'Initiating Green functions calculations.'
        call GREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NUME,lorentzbroad)
    endif

    ! Prepares two config files to be used by the impurity program
    open (1, file = 'impconfig.dat', action = 'write')
        write (1,*) NUMIMP, '! Number of impurities.'
        write (1,*) NUME, '! Number of energy values.'
    close(1)

    open (1, file = 'impatoms.dat', action = 'write')
    do i = 1, NUMIMP
        j = IMPPTSVAR(4,i)
        write (1,'(6F15.7, A, I)') E0(j), BETA(1,j), BETA(2,j), BETA(3,j), REAL(newDELTA(j)), &
        &AIMAG(newDELTA(j)), '    ! Impurity No. ', i
        write (1,*)
    end do
    write (1,'(A)') '------------------------------------------------------------'
    write (1,'(A)') 'Format: E_0,           B_x,           B_y,           B_z,           Re(D),         Im(D)'
    write (1,'(A)') 'Please insert the corresponding impurity value next to each element.'
    close(1)

    contains

    subroutine HOPPS(RLATT,IRLATTMAX,R0,NUMCHEMTYPES,LHOPS,NUMT,CHEMTYPE,TPTS,PREFACTORS)
        implicit none

        integer :: NUMT, i, j, ITYPE, JTYPE, IRLATT, IRLATTMAX, CHEMTYPE(NUMT), NUMCHEMTYPES
        real*8 :: TTPRIME(3), RPOINT(3), TPTS(3,NUMT), RLATT(3,IRLATTMAX), LHOPS(NUMCHEMTYPES,NUMCHEMTYPES), &
        &PREFACTORS(NUMCHEMTYPES,NUMCHEMTYPES), NNDIS(NUMCHEMTYPES,NUMCHEMTYPES), expon, R0, TOL

        TOL = 0.0001

        open (1, file = 'hoppings.dat', action = 'read')
        do i = 1, NUMCHEMTYPES
            read(1,*) (LHOPS(i,j), j = 1, NUMCHEMTYPES)
        end do
        close(1)
        
        do i = 1, NUMCHEMTYPES
            do j = 1, NUMCHEMTYPES
                if (LHOPS(i,j) /= LHOPS(j,i)) then
                    print *, 'Wrong value inserted as hopping element.'
                    print *, 'The', i,j, 'value is different from the', j,i, 'value.'
                    call exit(123)
                endif
            end do
        end do

        ! This calculates the nearest-neighbour distance for each interaction
        do ITYPE = 1, NUMCHEMTYPES
            do JTYPE = 1, NUMCHEMTYPES
                NNDIS(JTYPE,ITYPE) = 100000.0

                do i = 1, NUMT
                    do j = 1, NUMT

                        if (CHEMTYPE(i) == ITYPE .and. CHEMTYPE(j) == JTYPE) then
                            TTPRIME = TPTS(1:3,j) - TPTS(1:3,i)
                            do IRLATT = 1, IRLATTMAX
                                RPOINT = RLATT(1:3,IRLATT)

                                if (NNDIS(JTYPE,ITYPE) > norm2(RPOINT+TTPRIME) .and. norm2(RPOINT+TTPRIME) > TOL) then
                                    NNDIS(JTYPE,ITYPE) = norm2(RPOINT+TTPRIME)
                                endif

                            end do
                        endif

                    end do
                end do
    
                expon = exp(-NNDIS(JTYPE,ITYPE)/R0)
                PREFACTORS(JTYPE,ITYPE) = LHOPS(JTYPE,ITYPE)/expon ! This sets the prefactors to be entered in the Hamiltonian
            end do
        end do

    end subroutine HOPPS

    subroutine HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,DELTA,HAMILTONIANPREP)
        implicit none
        integer :: NUMT, i, FTIMO
        real*8 :: chempot, E0(NUMT), ULCN(NUMT), nu(NUMT), nuzero(NUMT), BETA(3,NUMT)
        complex*16 :: xPauli(2,2), yPauli(2,2), zPauli(2,2), IdentityPauli(2,2), helperham(2,2), DELTA(NUMT),&
        & HAMILTONIANPREP(4*NUMT,4*NUMT)

        HAMILTONIANPREP(:,:) = (0.0,0.0)

        do i = 1, NUMT

            FTIMO = 4*(i-1)

            helperham = (E0(i) - chempot + ULCN(i)*(nu(i) - nuzero(i)))*IdentityPauli -&
            &BETA(1,i)*xPauli - BETA(2,i)*yPauli - BETA(3,i)*zPauli

            HAMILTONIANPREP(1 + FTIMO, 1 + FTIMO) = helperham(1,1)
            HAMILTONIANPREP(1 + FTIMO, 2 + FTIMO) = helperham(1,2)
            !HAMILTONIANPREP(1 + FTIMO, 3 + FTIMO) = (0.0,0.0)
            HAMILTONIANPREP(1 + FTIMO, 4 + FTIMO) = DELTA(i)

            HAMILTONIANPREP(2 + FTIMO, 1 + FTIMO) = helperham(2,1)
            HAMILTONIANPREP(2 + FTIMO, 2 + FTIMO) = helperham(2,2)
            HAMILTONIANPREP(2 + FTIMO, 3 + FTIMO) = DELTA(i)
            !HAMILTONIANPREP(2 + FTIMO, 4 + FTIMO) = (0.0,0.0)

            !HAMILTONIANPREP(3 + FTIMO, 1 + FTIMO) = (0.0,0.0)
            HAMILTONIANPREP(3 + FTIMO, 2 + FTIMO) = CONJG(DELTA(i))
            HAMILTONIANPREP(3 + FTIMO, 3 + FTIMO) = -CONJG(helperham(1,1))
            HAMILTONIANPREP(3 + FTIMO, 4 + FTIMO) = CONJG(helperham(1,2))

            HAMILTONIANPREP(4 + FTIMO, 1 + FTIMO) = CONJG(DELTA(i))
            !HAMILTONIANPREP(4 + FTIMO, 2 + FTIMO) = (0.0,0.0)
            HAMILTONIANPREP(4 + FTIMO, 3 + FTIMO) = CONJG(helperham(2,1))
            HAMILTONIANPREP(4 + FTIMO, 4 + FTIMO) = -CONJG(helperham(2,2))

        end do

    end subroutine HAMPREP

    subroutine FOURIERHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)
        implicit none
        integer, intent(in) :: kcounter
        integer :: NUMT, NUMK, i, jneighb, NEIGHBNUM(NUMT), MAXNEIGHB, JATOM, JTYPE(MAXNEIGHB,NUMT), FTIMO, FTJMO
        real*8 :: HOPPVALS(MAXNEIGHB,NUMT)
        complex*16 :: HAMILTONIAN(4*NUMT,4*NUMT), HAMILTONIANPREP(4*NUMT,4*NUMT), EXPONS(NUMK,MAXNEIGHB,NUMT), term

        HAMILTONIAN(:,:) = HAMILTONIANPREP(:,:)

        do i = 1, NUMT
            FTIMO = 4*(i-1)
            do jneighb = 1, NEIGHBNUM(i)

                JATOM = JTYPE(jneighb,i)
                FTJMO = 4*(JATOM-1)

                term = EXPONS(kcounter,jneighb,i)*HOPPVALS(jneighb,i)

                HAMILTONIAN(1 + FTIMO,1 + FTJMO) = HAMILTONIAN(1 + FTIMO,1 + FTJMO) + term
                !HAMILTONIAN(1 + FTIMO,2 + FTJMO) = HAMILTONIAN(1 + FTIMO,2 + FTJMO) + term
                !HAMILTONIAN(2 + FTIMO,1 + FTJMO) = HAMILTONIAN(2 + FTIMO,1 + FTJMO) + term
                HAMILTONIAN(2 + FTIMO,2 + FTJMO) = HAMILTONIAN(2 + FTIMO,2 + FTJMO) + term
                HAMILTONIAN(3 + FTIMO,3 + FTJMO) = HAMILTONIAN(3 + FTIMO,3 + FTJMO) - term
                !HAMILTONIAN(3 + FTIMO,4 + FTJMO) = HAMILTONIAN(3 + FTIMO,4 + FTJMO) + term
                !HAMILTONIAN(4 + FTIMO,3 + FTJMO) = HAMILTONIAN(4 + FTIMO,3 + FTJMO) + term
                HAMILTONIAN(4 + FTIMO,4 + FTJMO) = HAMILTONIAN(4 + FTIMO,4 + FTJMO) - term

            end do
        end do

    end subroutine FOURIERHAM

    subroutine RANDOMKHAM(KPOINT,HAMILTONIANPREP,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,RCONNECT,NUMT,HAMILTONIAN)
        implicit none
        real*8, intent(in) :: KPOINT(3)
        integer :: NUMT, i, jneighb, MAXNEIGHB, NEIGHBNUM(NUMT), JATOM, JTYPE(MAXNEIGHB,NUMT), FTIMO, FTJMO
        real*8 :: hopping, HOPPVALS(MAXNEIGHB,NUMT), RPOINT(3), lambda, RCONNECT(3,MAXNEIGHB,NUMT)
        complex*16 :: expon, HAMILTONIAN(4*NUMT,4*NUMT), HAMILTONIANPREP(4*NUMT,4*NUMT)

        HAMILTONIAN(:,:) = HAMILTONIANPREP(:,:)

        do i = 1, NUMT
            FTIMO = 4*(i-1)
            do jneighb = 1, NEIGHBNUM(i)

                RPOINT = RCONNECT(1:3,jneighb,i)
                JATOM = JTYPE(jneighb,i)
                FTJMO = 4*(JATOM-1)

                hopping = HOPPVALS(jneighb,i)

                expon = exp(-CI*DOT_PRODUCT(KPOINT,RPOINT))

                HAMILTONIAN(1 + FTIMO,1 + FTJMO) = HAMILTONIAN(1 + FTIMO,1 + FTJMO) + expon*hopping
                !HAMILTONIAN(1 + FTIMO,2 + FTJMO) = HAMILTONIAN(1 + FTIMO,2 + FTJMO) + expon*hopping
                !HAMILTONIAN(2 + FTIMO,1 + FTJMO) = HAMILTONIAN(2 + FTIMO,1 + FTJMO) + expon*hopping
                HAMILTONIAN(2 + FTIMO,2 + FTJMO) = HAMILTONIAN(2 + FTIMO,2 + FTJMO) + expon*hopping
                HAMILTONIAN(3 + FTIMO,3 + FTJMO) = HAMILTONIAN(3 + FTIMO,3 + FTJMO) - CONJG(expon*hopping)
                !HAMILTONIAN(3 + FTIMO,4 + FTJMO) = HAMILTONIAN(3 + FTIMO,4 + FTJMO) + CONJG(expon*hopping)
                !HAMILTONIAN(4 + FTIMO,3 + FTJMO) = HAMILTONIAN(4 + FTIMO,3 + FTJMO) + CONJG(expon*hopping)
                HAMILTONIAN(4 + FTIMO,4 + FTJMO) = HAMILTONIAN(4 + FTIMO,4 + FTJMO) - CONJG(expon*hopping)

            end do
        end do        

    end subroutine RANDOMKHAM

    subroutine GREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NUME,lorentzbroad)
        implicit none

        integer :: NUMT, NUMK, NUME, IE, i, j, k, n, ini, fin, NUMIMP, IMPPTSVAR(4,NUMIMP), &
        &l, m, a, aprime, dosorno, FTIMO, FTJMO, NUMTKONE
        real*8, allocatable, dimension(:,:) :: greendensityperatom, greendensity
        complex*16, allocatable, dimension(:) :: energies
        real*8 :: PI, EIGENVALUES(4*NUMT*NUMK), energyintervals, lorentzbroad, &
        &a_1(3), a_2(3), a_3(3), RPOINT(3), RPRIMEPOINT(3), FOURIERVEC(3), TPTS(3,NUMT), KPTS(3,NUMK), KPOINT(3)
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT*NUMK), GMATRIX(4*NUMT,4*NUMT), EZ, ENFRAC, GK(4*NUMT,4*NUMT*NUMK),&
        &GREENR(4*NUMIMP,4*NUMIMP), expon

        ! This part sets up the E points to be used in plots concerning the Green's function
        ! At the moment it simply gives NUME real energies.

        print *, 'If you want to perform a DoS calculation using the Green functions press 0, otherwise press any other number.'
        read *, dosorno

        allocate(energies(NUME+1))
        if (dosorno == 0) then
            allocate(greendensityperatom(1+NUMT,NUME+1)) ! The first column is for the energies, the rest for the values
            allocate(greendensity(2,NUME+1)) ! The first column is for the energies, the other for the values
        endif

        energyintervals = (MAXVAL(EIGENVALUES) - MINVAL(EIGENVALUES))/NUME

        ! These are the "complex" energies E
        do i = 1, NUME+1
            energies(i) = cmplx(MINVAL(EIGENVALUES) + energyintervals*(i-1), lorentzbroad)
        end do

        ! This part writes all the Fouriered Green function elements per energy at this text file
        open(1, file = 'greenimp.txt', action = 'write')

        do IE = 1, NUME+1 ! Begins a loop over the energies, in order to find G(E) for each E

            EZ = energies(IE)

            if (dosorno == 0) then
                ! Startup values for the densities
                greendensityperatom(1,IE) = REAL(EZ)
                greendensity(1,IE) = REAL(EZ)
                greendensity(2,IE) = 0.0
                do i = 1, NUMT
                    greendensityperatom(1+i,IE) = 0.0
                end do
            endif

            ! This initiates the calculation of the Green's function matrix G(α,α';E) per k-point
            do k = 1, NUMK ! k
                NUMTKONE = 4*NUMT*(k-1)

                ! Set all G-matrix values equal to zero, so that the following summation can work
                GMATRIX(:,:) = (0.0,0.0)

                do i = 1, NUMT ! α
                    FTIMO = 4*(i-1)
                    do j = 1, NUMT ! α'
                        FTJMO = 4*(j-1)

                        do n = 1, 4*NUMT ! This is the sum over all eigenenergies per k

                            ENFRAC = (1.0/(EZ-EIGENVALUES(n + NUMTKONE)))
                            
                            GMATRIX(1 + FTIMO, 1 + FTJMO) = GMATRIX(1 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(1+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(1+FTJMO,n+NUMTKONE)) ! 11
                            GMATRIX(1 + FTIMO, 2 + FTJMO) = GMATRIX(1 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(1+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(2+FTJMO,n+NUMTKONE)) ! 12
                            GMATRIX(1 + FTIMO, 3 + FTJMO) = GMATRIX(1 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(1+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(3+FTJMO,n+NUMTKONE)) ! 13
                            GMATRIX(1 + FTIMO, 4 + FTJMO) = GMATRIX(1 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(1+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(4+FTJMO,n+NUMTKONE)) ! 14

                            GMATRIX(2 + FTIMO, 1 + FTJMO) = GMATRIX(2 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(2+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(1+FTJMO,n+NUMTKONE)) ! 21
                            GMATRIX(2 + FTIMO, 2 + FTJMO) = GMATRIX(2 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(2+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(2+FTJMO,n+NUMTKONE)) ! 22
                            GMATRIX(2 + FTIMO, 3 + FTJMO) = GMATRIX(2 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(2+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(3+FTJMO,n+NUMTKONE)) ! 23
                            GMATRIX(2 + FTIMO, 4 + FTJMO) = GMATRIX(2 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(2+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(4+FTJMO,n+NUMTKONE)) ! 24

                            GMATRIX(3 + FTIMO, 1 + FTJMO) = GMATRIX(3 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(3+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(1+FTJMO,n+NUMTKONE)) ! 31
                            GMATRIX(3 + FTIMO, 2 + FTJMO) = GMATRIX(3 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(3+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(2+FTJMO,n+NUMTKONE)) ! 32
                            GMATRIX(3 + FTIMO, 3 + FTJMO) = GMATRIX(3 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(3+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(3+FTJMO,n+NUMTKONE)) ! 33
                            GMATRIX(3 + FTIMO, 4 + FTJMO) = GMATRIX(3 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(3+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(4+FTJMO,n+NUMTKONE)) ! 34

                            GMATRIX(4 + FTIMO, 1 + FTJMO) = GMATRIX(4 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(4+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(1+FTJMO,n+NUMTKONE)) ! 41
                            GMATRIX(4 + FTIMO, 2 + FTJMO) = GMATRIX(4 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(4+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(2+FTJMO,n+NUMTKONE)) ! 42
                            GMATRIX(4 + FTIMO, 3 + FTJMO) = GMATRIX(4 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(4+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(3+FTJMO,n+NUMTKONE)) ! 43
                            GMATRIX(4 + FTIMO, 4 + FTJMO) = GMATRIX(4 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(4+FTIMO,n + NUMTKONE)*CONJG(EIGENVECTORS(4+FTJMO,n+NUMTKONE)) ! 44

                        end do

                    end do
                end do

                if (dosorno == 0) then
                    do i = 1, NUMT ! Calculation of density for each atom
                        FTIMO = 4*(i-1)
                        do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                            greendensityperatom(1+i,IE) = greendensityperatom(1+i,IE) -&
                            &(1.0/PI)*AIMAG(GMATRIX(j + FTIMO, j + FTIMO))
                        end do
                    end do
                endif

                ini = 1 + (k-1)*4*NUMT
                fin = k*4*NUMT
                GK(:,ini:fin) = GMATRIX ! This is a table that contains all G(α,α',k;E) per energy E
            
            end do ! ends k-sum

            ! Fourier transform of G(k) into G(r-r')
            ! This constructs a 4*NUMIMP x 4*NUMIMP G(r-r') matrix for each energy EZ

            ! Startup
            GREENR(:,:) = (0.0,0.0)
            
            do i = 1, NUMIMP
                FTIMO = 4*(i-1)
                do j = 1, NUMIMP
                    FTJMO = 4*(j-1)

                    a = IMPPTSVAR(4,i)
                    aprime = IMPPTSVAR(4,j)
                    RPOINT = IMPPTSVAR(1,i)*a_1 + IMPPTSVAR(2,i)*a_2 + IMPPTSVAR(3,i)*a_3
                    RPRIMEPOINT = IMPPTSVAR(1,j)*a_1 + IMPPTSVAR(2,j)*a_2 + IMPPTSVAR(3,j)*a_3

                    FOURIERVEC(1:3) = RPRIMEPOINT(1:3) - RPOINT(1:3) + TPTS(1:3,aprime) - TPTS(1:3,a)

                    do k = 1, NUMK

                        KPOINT = KPTS(1:3,k)
                        expon = exp(-CI*DOT_PRODUCT(KPOINT,FOURIERVEC))

                        do l = 1, 4
                            do m = 1, 4
                                GREENR(m + FTIMO, l + FTJMO) = GREENR(m + FTIMO, l + FTJMO) +&
                                &expon*GK(m + 4*(a-1) , l + 4*(aprime-1) + NUMTKONE)                    
                            end do
                        end do

                    end do

                end do
            end do

            ! Writes the Green impurity elements on greenimp.txt
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1, '(F17.8,F17.8)') GREENR(i,j)
                end do
            end do

            if (dosorno == 0) then
                do i = 1, NUMT ! Calculation of full density
                    greendensityperatom(i+1,IE) = greendensityperatom(i+1,IE)/NUMK ! Normalization
                    greendensity(2,IE) = greendensity(2,IE) + greendensityperatom(1+i,IE)
                end do
            endif

        end do ! ends energies sum
        close(1) ! Closes the greenimp.txt file

        if (dosorno == 0) then
            open(1, file = 'greendensityperatom.txt', action = 'write')
            do j = 1, NUME+1 ! Energies = Intervals + 1
                do i = 1, NUMT+1
                    write (1,'(F17.8)',advance='no') greendensityperatom(i,j)
                end do
                write (1,*)
            end do
            close(1)

            open(1, file = 'greendensity.txt', action = 'write')
            do j = 1, NUME+1 ! Energies = Intervals + 1
                write (1,'(2F17.8)') greendensity(1,j), greendensity(2,j)
            end do
            close(1)
        endif

    end subroutine GREEN


    subroutine BANDS(NUMT,W,WORK,LWORK,RWORK,INFO,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0, &
        &ULCN,nu,nuzero,BETA,DELTA,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,RCONNECT,ALAT)
        implicit none

        integer :: NUMDIR, i, io, j, m, n, NUMT, LWORK, INFO, intpointer, checker, MAXNEIGHB, NEIGHBNUM(NUMT), &
        &JTYPE(MAXNEIGHB,NUMT)
        integer, allocatable, dimension(:) :: KINTS
        real*8 :: DIR(3), KPOINT(3), HORINT, W(4*NUMT), RWORK(3*(4*NUMT) - 2), chempot, E0(NUMT), ULCN(NUMT), &
        &nu(NUMT), nuzero(NUMT), BETA(3,NUMT), ALAT, HOPPVALS(MAXNEIGHB,NUMT), RCONNECT(3,MAXNEIGHB,NUMT)
        real*8, allocatable, dimension(:,:) :: INPOINT, OUTPOINT
        complex*16 :: HAMILTONIAN(4*NUMT,4*NUMT), WORK(LWORK), zPauli(2,2), IdentityPauli(2,2), xPauli(2,2), &
        &yPauli(2,2), DELTA(NUMT), HAMILTONIANPREP(4*NUMT,4*NUMT)

        ! The following reads the number of basis vectors from a file named basisvectors.dat
        NUMDIR = 0
        open (1, file = 'bandpoints.dat', action = 'read')
        do
            read(1,*,iostat=io)
            if (io/=0) exit
            NUMDIR = NUMDIR + 1
        end do
        close(1)

        NUMDIR = NUMDIR - 2 ! To ignore the final 2 lines in basisvectors.dat which are the configuration settings.

        allocate(INPOINT(3,NUMDIR))
        allocate(OUTPOINT(3,NUMDIR))
        allocate(KINTS(NUMDIR))

        open(1, file = 'bandpoints.dat', action = 'read')
        do i = 1, NUMDIR
            read (1,*) INPOINT(1:3,i), OUTPOINT(1:3,i), KINTS(i)
        end do
        close(1)

        ! Inverse space points, instead of configuration space input
        do i = 1, NUMDIR
            do j = 1, 3
                INPOINT(j,i) = ((2.0*PI)/ALAT)*INPOINT(j,i)
                OUTPOINT(j,i) = ((2.0*PI)/ALAT)*OUTPOINT(j,i)
            end do
        end do

        print *, 'Press 0 to print only intervals and any other number to print k-values as well.'
        read *, intpointer
        HORINT = 0.0 ! This is a parameter that adjusts the widths for each k-dimension window at the band diagram
        
        open(2, file = 'bands.txt', action = 'write')
        do i = 1, NUMDIR
            
            DIR = OUTPOINT(1:3,i) - INPOINT(1:3,i) ! Sets the direction along which we calculate K points
            KPOINT = INPOINT(1:3,i) ! Startup of each k

            ! In order to avoid taking some points twice
            if (i /= NUMDIR) then
                if (norm2(OUTPOINT(1:3,i) - INPOINT(1:3,i+1)) < 0.0001) then
                    checker = KINTS(i)
                else
                    checker = KINTS(i)+1
                endif
            else
                checker = KINTS(i)+1
            endif

            do j = 1, checker

                call HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,DELTA,HAMILTONIANPREP)

                call RANDOMKHAM(KPOINT,HAMILTONIANPREP,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,RCONNECT,NUMT,HAMILTONIAN)

                ! N because we only want eigenvalues and not eigenvectors
                call zheev ('N', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO) ! Don't forget to reconfigure those whenever the dimensions change!
                
                if (intpointer == 0) then
                    write (2,'(F17.8)', advance='no') HORINT
                    do m = 1, 4*NUMT
                        write (2,'(F17.8)', advance='no') W(m)
                    end do
                    write (2,*)
                else 
                    write (2,'(4F17.8)', advance='no') KPOINT, HORINT
                    do m = 1, 4*NUMT
                        write (2,'(F17.8)', advance='no') W(m)
                    end do
                    write (2,*)
                endif

                if (j /= checker) then
                    KPOINT = KPOINT + (1.0/KINTS(i))*DIR
                    HORINT = HORINT + (1.0/KINTS(i))*norm2(DIR)
                endif

            end do

        end do
        close(2)

    end subroutine BANDS

    subroutine INT_NUM_DEN(uniquecounter,EIGENVALUES,NUMT,NUMK,SORTEDEIGVALS,multiplicity)
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

        open(1, file = 'intnumdensity.txt', action = 'write')
        do j = 1, uniquecounter
            write (1,'(2F17.8)',advance='no') intnumdensity(1,j), intnumdensity(2,j)
        end do
        close(1)

    end subroutine INT_NUM_DEN

    ! --------------------------------------------------------------------------------------------------------------------------------
    subroutine NUM_DEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,NUME,lorentzbroad,metalorno)
        implicit none

        integer :: i, j, IE, NUMT, NUMK, NUME, metalorno, FTIMO
        real*8, allocatable, dimension(:,:) :: numdensityperatom, numdensity
        real*8 :: PI, lorentzbroad, EIGENVALUES(4*NUMT*NUMK), energyintervals
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT*NUMK)

        allocate(numdensityperatom(1+NUMT,NUME+1))
        allocate(numdensity(2,NUME+1))

        energyintervals = (MAXVAL(EIGENVALUES) - MINVAL(EIGENVALUES))/NUME

        ! These are the energies E
        do i = 0, NUME
            numdensityperatom(1,i+1) = MINVAL(EIGENVALUES) + energyintervals*i
            numdensity(1,i+1) = MINVAL(EIGENVALUES) + energyintervals*i
            numdensity(2,i+1) = 0.0
        end do

        ! Calculation of the DoS PER ATOM
        do i = 1, NUMT
            FTIMO = 4*(i-1)
            do IE = 0, NUME
                numdensityperatom(1+i,IE+1) = 0.0
                do j = 1, 4*NUMT*NUMK
                    numdensityperatom(1+i,IE+1) = numdensityperatom(1+i,IE+1) + (lorentzbroad/pi)*&
                    &((abs(EIGENVECTORS(1+FTIMO,j))**2 + abs(EIGENVECTORS(2+FTIMO,j))**2)/((numdensityperatom(1,IE+1) -&
                    &EIGENVALUES(j))**2 + lorentzbroad**2)) ! Sum over all eigenvalues with |u| as weights
                end do
            end do
        end do

        ! Normalization
        do i = 1, NUMT
            do j = 1, NUME+1
                numdensityperatom(i+1,j) = numdensityperatom(i+1,j)/NUMK
            end do
        end do

        if (metalorno == 0) then
            open(1, file = 'numdensityperatom.txt', action = 'write')
            do j = 1, NUME+1
                do i = 1, NUMT+1
                    write (1,'(F17.8)',advance='no') numdensityperatom(i,j)
                end do
                write (1,*)
            end do
            close(1)

            ! Calculation of the full density of states 
            do i = 1, NUMT
                numdensity(2,:) = numdensity(2,:) + (1.0/NUMT)*numdensityperatom(1+i,:) ! 1/NUMT for normalization
            end do

            open(1, file = 'numdensity.txt', action = 'write')
            do j = 1, NUME+1
                write (1,'(2F17.8)') numdensity(1,j), numdensity(2,j)
            end do
            close(1)
        else if (metalorno == 1) then ! A simple workaround so that we get different files for metals and SCs
            open(1, file = 'metalnumdensityperatom.txt', action = 'write')
            do j = 1, NUME+1
                do i = 1, NUMT+1
                    write (1,'(F17.8)',advance='no') numdensityperatom(i,j)
                end do
                write (1,*)
            end do
            close(1)

            ! Calculation of the full density of states 
            do i = 1, NUMT
                numdensity(2,:) = numdensity(2,:) + (1.0/NUMT)*numdensityperatom(1+i,:) ! 1/NUMT for normalization
            end do

            open(1, file = 'metalnumdensity.txt', action = 'write')
            do j = 1, NUME+1
                write (1,'(2F17.8)') numdensity(1,j), numdensity(2,j)
            end do
            close(1)
        endif

    end subroutine NUM_DEN
    
    subroutine BZ(a,b,c,N_x,N_y,N_z,PI,KPTS,Ntot,b_1,b_2,b_3)
        implicit none

        integer :: N_x, N_y, N_z, Ntot, counter, c_1, c_2, c_3
        real*8 :: a(3), b(3), c(3), b_1(3), b_2(3), b_3(3)
        real*8, allocatable, dimension(:,:) :: KPTS
        real*8 :: PI, volume
        counter = 1
		
        volume = DOT_PRODUCT(a, CROSS_PRODUCT(b,c)) !The volume of the unit cell

		!At this point we calculate the reciprocal space vectors.
        b_1 = 2.0*PI*CROSS_PRODUCT(b,c)/volume
        b_2 = 2.0*PI*CROSS_PRODUCT(c,a)/volume
        b_3 = 2.0*PI*CROSS_PRODUCT(a,b)/volume

        Ntot = N_x*N_y*N_z !The total number of different wavevectors in reciprocal space.

        allocate(KPTS(3,Ntot))
                    
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

    subroutine IDENTIFIER(imppointer,b_1,b_2,b_3,PI,a_1,a_2,a_3,TPTS,NUMT,NUMIMP,IMPPTSVAR)
        implicit none

        integer :: NUMIMP, i, j, io, NUMT, counter, imppointer, k
        real*8 :: b_1(3), b_2(3), b_3(3), PI, IMPPOINT(3), TPOINT(3), RPOINT(3), a_1(3), a_2(3), a_3(3), TPTS(3,NUMT),&
        &BASIS(3), eps
        real*8, allocatable, dimension(:,:) :: IMPPTS
        integer, allocatable, dimension(:,:) :: IMPPTSVAR

        eps = 0.0001

        NUMIMP = 0
        if (imppointer == 0) then
            ! Counts the number of impurities
            open (1, file = 'impurities.dat', action = 'read')
            do
                read(1,*,iostat=io)
                if (io/=0) exit
                NUMIMP = NUMIMP + 1
            end do
            close(1)

            NUMIMP = NUMIMP - 2

            allocate(IMPPTS(3,NUMIMP))
            allocate(IMPPTSVAR(4,NUMIMP)) ! This is an array with elements n_1, n_2, n_3, basis index, for each impurity site

            ! Reads the impurities coordinates
            open (1, file = 'impurities.dat', action = 'read')
                do i = 1, NUMIMP
                    read(1,*) IMPPTS(1:3,i)
                end do
            close(1)

            ! Discerns the R + τ part for each impurity vector
            do i = 1, NUMIMP
                IMPPOINT = IMPPTS(1:3,i)

                ! These find the R part
                IMPPTSVAR(1,i) = FLOOR((1/((2.0*PI)))*DOT_PRODUCT(IMPPOINT,b_1)) ! n_1
                IMPPTSVAR(2,i) = FLOOR((1/((2.0*PI)))*DOT_PRODUCT(IMPPOINT,b_2)) ! n_2
                IMPPTSVAR(3,i) = FLOOR((1/((2.0*PI)))*DOT_PRODUCT(IMPPOINT,b_3)) ! n_3
                RPOINT = IMPPTSVAR(1,i)*a_1 + IMPPTSVAR(2,i)*a_2 + IMPPTSVAR(3,i)*a_3

                ! The remainder is the IMPPOINT-R = τ'
                TPOINT = IMPPOINT - RPOINT

                ! This tries to identify that remainder with a basis atom in the unit cell
                counter = 0
                do j = 1, NUMT
                    BASIS = TPTS(1:3,j)

                    if (abs(TPOINT(1)-BASIS(1)) <= eps .and. abs(TPOINT(2)-BASIS(2)) <= eps &
                    &.and. abs(TPOINT(3)-BASIS(3)) <= eps) then
                        counter = j
                    endif
                end do

                if (counter == 0) then
                    print *, 'The impurity No.', i, 'does not correspond to a lattice atom.'
                    call exit(123)
                else
                    IMPPTSVAR(4,i) = counter
                endif
            end do
        else
            ! Counts the number of impurities
            open (1, file = 'impvals.txt', action = 'read')
            do
                read(1,*,iostat=io)
                if (io/=0) exit
                NUMIMP = NUMIMP + 1
            end do
            close(1)

            NUMIMP = NUMIMP - 2

            allocate(IMPPTSVAR(4,NUMIMP)) ! This is an array with elements n_1, n_2, n_3, basis index, for each impurity site

            ! Reads the impurities coordinates
            open (1, file = 'impvals.txt', action = 'read')
                do i = 1, NUMIMP
                    read(1,*) IMPPTSVAR(1:4,i)
                end do
            close(1)

            do i = 1, NUMIMP
                if (IMPPTSVAR(4,i) > NUMT .or. IMPPTSVAR(4,i) <= 0) then
                    print *, 'The impurity No.', i, 'does not correspond to a lattice atom.'
                    call exit(123)
                endif
            end do
        endif
    
    end subroutine IDENTIFIER

    subroutine CONSTANTS(s_0,s_1,s_2,s_3,CI,PI,KB) ! Sets some global constants
        implicit none
        complex*16 :: s_0(2,2), s_1(2,2), s_2(2,2), s_3(2,2), CI
        real*8 :: PI, KB
        
        s_0(1,1) = (1.0,0.0)
        s_0(1,2) = (0.0,0.0)
        s_0(2,1) = (0.0,0.0)
        s_0(2,2) = (1.0,0.0)
        s_1(1,1) = (0.0,0.0)
        s_1(1,2) = (1.0,0.0)
        s_1(2,1) = (1.0,0.0)
        s_1(2,2) = (0.0,0.0)
        s_2(1,1) = (0.0,0.0)
        s_2(1,2) = (0.0,-1.0)
        s_2(2,1) = (0.0,1.0)
        s_2(2,2) = (0.0,0.0)
        s_3(1,1) = (1.0,0.0)
        s_3(1,2) = (0.0,0.0)
        s_3(2,1) = (0.0,0.0)
        s_3(2,2) = (-1.0,0.0)

        CI = (0.0,1.0) ! setting the imaginary unit
        PI = 4.D0*atan(1.D0) ! setting π
        KB = 1.0 ! The value in eVs is 8.617385D-5

    end subroutine CONSTANTS

    function FERMI(E,T,KB) result(fE)
        implicit none
        real*8, intent(in) :: E, T
        real*8 :: fE, KB, term

        if (T == 0.0) then
            if (E < 0.0) then
                fE = 1.0
            else if (E == 0.0) then
                fE = 0.5
            else
                fE = 0.0
            endif
        else
            term = E/(KB*T)
            ! 460.0 is a check to see if this exponential is
            ! 200ln(10), i.e. so large that fortran can't process it, since the limit is ~300ln(10)
            if (term > 460.0) then
                fE = 0.0
            else if (term < -460.0) then
                fE = 1.0
            else
                fE = 1.0/(exp(term)+1.0)
            endif
        endif
		
    end function FERMI

    function CROSS_PRODUCT(x,y) result(cross)
        implicit none
        real*8, dimension(3), intent(in) :: x, y
        real*8, dimension(3) :: cross

        cross(1) = x(2)*y(3) - x(3)*y(2)
        cross(2) = x(3)*y(1) - x(1)*y(3)
        cross(3) = x(1)*y(2) - x(2)*y(1)
		
    end function CROSS_PRODUCT

end program TB