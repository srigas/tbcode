program GREEN
    implicit none
    integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, LWORK, INFO, NCELLS, &
    &NUMIMP, NUMCHEMTYPES, NUME, JATOM, rcheck, MAXNEIGHB, dosorno
    integer, allocatable, dimension(:) :: CHEMTYPE, NEIGHBNUM
    integer, allocatable, dimension(:,:) :: IMPPTSVAR, JTYPE
    real*8 :: ALAT, a_1(3), a_2(3), a_3(3), RMAX, R0, KPOINT(3), RPOINT(3), TTPRIME(3),&
    &chempot, T, PI, KB, b_1(3), b_2(3), b_3(3), eta, lambda, TOL, tempvalre, tempvalim
    real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, nu, nuzero, magnet, VSUPCOND
    real*8, allocatable, dimension(:,:) :: KPTS, TPTS, EIGENVALUES, RLATT, BETA, LHOPS, PREFACTORS, HOPPVALS
    real*8, allocatable, dimension(:,:,:) :: RCONNECT
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2)
    complex*16, allocatable, dimension(:) :: WORK, DELTA
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, HAMILTONIANPREP
    complex*16, allocatable, dimension(:,:,:) :: EXPONS, EIGENVECTORS
    character(len = 1) :: vecorints

    TOL = 0.0001 ! The fault tolerance for lattice vectors' norms

    call CONSTANTS(IdentityPauli,xPauli,yPauli,zPauli,CI,PI,KB)

    open(1, file = 'config.dat', action = 'read')
    read(1,*) ALAT
    read(1,*) a_1
    read(1,*) a_2
    read(1,*) a_3
    read(1,*) NUMKX,NUMKY,NUMKZ
    read(1,*) RMAX
    read(1,*) R0 
    read(1,*) NCELLS 
    read(1,*) T
    close(1)

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
    allocate(BETA(3,NUMT)) ! NUMT x 3 matrix with the B_x, B_y, B_z for each basis atom
    allocate(DELTA(NUMT))
    allocate(magnet(NUMT))
    allocate(VSUPCOND(NUMT))
    allocate(CHEMTYPE(NUMT)) ! Disciminates between different chemical elements

    open (1, file = 'basisvectors.dat', action = 'read')
    do i = 1, NUMT
        read(1,*) TPTS(1:3,i), CHEMTYPE(i), E0(i), ULCN(i), nuzero(i), BETA(1,i), BETA(2,i), BETA(3,i), VSUPCOND(i)
    end do
    close(1)

    open (1, file = 'scresults.dat', action = 'read')
    read(1,*) chempot
    do i = 1, NUMT
        read(1,*) nu(i), magnet(i), tempvalre, tempvalim
        DELTA(i) = dcmplx(tempvalre, tempvalim)
    end do
    close(1)

    NUMCHEMTYPES = MAXVAL(CHEMTYPE)
    allocate(LHOPS(NUMCHEMTYPES,NUMCHEMTYPES))
    allocate(PREFACTORS(NUMCHEMTYPES,NUMCHEMTYPES))

    ! The *4 factors are now due to spin and particle-hole
    allocate(EIGENVECTORS(4*NUMT,4*NUMT,NUMK))
    allocate(EIGENVALUES(4*NUMT,NUMK))
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

    print *, 'Initiating the diagonalization of the Hamiltonian...'

    ! Run the Hamiltonian to evaluate the Eigenvalues and Eigenvectors to be used for the Green calculation.
    call HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,DELTA,HAMILTONIANPREP)

    do kcounter = 1, NUMK

        call FOURIERHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)

        call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO)

        EIGENVALUES(:,kcounter) = W
        EIGENVECTORS(:,:,kcounter) = HAMILTONIAN
    end do

    print *, 'Diagonalization finished and eigenvectors/eigenvalues stored.'

    print *, 'Please insert the number of energies for the Green function.'
    read *, NUME
    print *, 'Please insert the imaginary part of the energies.'
    read *, eta

    print *, 'Choose the input way for the impurity sites, v/i.'
    print *, 'v to input via vectors or i to input via integers.'
    150 read *, vecorints

    if (vecorints /= 'v' .and. vecorints /= 'i') then
        print *, 'Invalid input. Please choose v or i.'
        goto 150
    endif

    call IDENTIFIER(vecorints,b_1,b_2,b_3,PI,a_1,a_2,a_3,TPTS,NUMT,NUMIMP,IMPPTSVAR)

    print *, 'If you want to get the full G(k) function and calculate the DoS press 0.'
    print *, 'If you only want to study the impurity problem press any other number.'
    read *, dosorno

    print *, 'Initiating Green functions calculations.'
    if (dosorno == 0) then
        call FULLGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NUME,eta)
    else
        call PARTIALGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NUME,eta)
    endif

    ! Prepares two config files to be used by the impurity program
    open (1, file = 'impconfig.dat', action = 'write')
        write (1,*) NUMIMP, '! Number of impurities.'
        write (1,*) NUME, '! Number of energy values.'
    close(1)

    open (1, file = 'impatoms.dat', action = 'write')
    do i = 1, NUMIMP
        j = IMPPTSVAR(4,i)
        write (1,'(6F15.7, A, I7)') E0(j), BETA(1,j), BETA(2,j), BETA(3,j), REAL(DELTA(j)), &
        &AIMAG(DELTA(j)), '    ! Impurity No. ', i
        write (1,*)
    end do
    write (1,'(A)') '------------------------------------------------------------'
    write (1,'(A)') 'Format: E_0,           B_x,           B_y,           B_z,           Re(D),         Im(D)'
    write (1,'(A)') 'Please insert the corresponding impurity value under each element.'
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
        integer :: NUMT, i
        real*8 :: chempot, E0(NUMT), ULCN(NUMT), nu(NUMT), nuzero(NUMT), BETA(3,NUMT)
        complex*16 :: xPauli(2,2), yPauli(2,2), zPauli(2,2), IdentityPauli(2,2), helperham(2,2), DELTA(NUMT),&
        & HAMILTONIANPREP(4*NUMT,4*NUMT)

        HAMILTONIANPREP(:,:) = (0.0,0.0)

        do i = 1, NUMT

            helperham = (E0(i) - chempot + ULCN(i)*(nu(i) - nuzero(i)))*IdentityPauli -&
            &BETA(1,i)*xPauli - BETA(2,i)*yPauli - BETA(3,i)*zPauli

            HAMILTONIANPREP(i, i) = helperham(1,1)
            HAMILTONIANPREP(i, i + NUMT) = helperham(1,2)
            !HAMILTONIANPREP(i, i + 2*NUMT) = (0.0,0.0)
            HAMILTONIANPREP(i, i + 3*NUMT) = DELTA(i)

            HAMILTONIANPREP(i + NUMT, i) = helperham(2,1)
            HAMILTONIANPREP(i + NUMT, i + NUMT) = helperham(2,2)
            HAMILTONIANPREP(i + NUMT, i + 2*NUMT) = DELTA(i)
            !HAMILTONIANPREP(i + NUMT, i + 3*NUMT) = (0.0,0.0)

            !HAMILTONIANPREP(i + 2*NUMT, i) = (0.0,0.0)
            HAMILTONIANPREP(i + 2*NUMT, i + NUMT) = CONJG(DELTA(i))
            HAMILTONIANPREP(i + 2*NUMT, i + 2*NUMT) = -CONJG(helperham(1,1))
            HAMILTONIANPREP(i + 2*NUMT, i + 3*NUMT) = CONJG(helperham(1,2))

            HAMILTONIANPREP(i + 3*NUMT, i) = CONJG(DELTA(i))
            !HAMILTONIANPREP(i + 3*NUMT, i + NUMT) = (0.0,0.0)
            HAMILTONIANPREP(i + 3*NUMT, i + 2*NUMT) = CONJG(helperham(2,1))
            HAMILTONIANPREP(i + 3*NUMT, i + 3*NUMT) = -CONJG(helperham(2,2))

        end do

    end subroutine HAMPREP

    subroutine FOURIERHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)
        implicit none
        integer, intent(in) :: kcounter
        integer :: NUMT, NUMK, i, jneighb, NEIGHBNUM(NUMT), MAXNEIGHB, JATOM, JTYPE(MAXNEIGHB,NUMT)
        real*8 :: HOPPVALS(MAXNEIGHB,NUMT)
        complex*16 :: HAMILTONIAN(4*NUMT,4*NUMT), HAMILTONIANPREP(4*NUMT,4*NUMT), EXPONS(NUMK,MAXNEIGHB,NUMT), term

        HAMILTONIAN(:,:) = HAMILTONIANPREP(:,:)

        do i = 1, NUMT
            do jneighb = 1, NEIGHBNUM(i)

                JATOM = JTYPE(jneighb,i)

                term = EXPONS(kcounter,jneighb,i)*HOPPVALS(jneighb,i)

                HAMILTONIAN(i,JATOM) = HAMILTONIAN(i,JATOM) + term
                HAMILTONIAN(i + NUMT,JATOM + NUMT) = HAMILTONIAN(i + NUMT,JATOM + NUMT) + term
                HAMILTONIAN(i + 2*NUMT,JATOM + 2*NUMT) = HAMILTONIAN(i + 2*NUMT,JATOM + 2*NUMT) - term
                HAMILTONIAN(i + 3*NUMT,JATOM + 3*NUMT) = HAMILTONIAN(i + 3*NUMT,JATOM + 3*NUMT) - term

            end do
        end do

    end subroutine FOURIERHAM

    subroutine FULLGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NUME,eta)
        implicit none

        integer :: NUMT, NUMK, NUME, IE, i, j, k, n, NUMIMP, IMPPTSVAR(4,NUMIMP), &
        &l, m, a, aprime, FTIMO, FTJMO, checker
        real*8 :: PI, EIGENVALUES(4*NUMT,NUMK), energyintervals, eta, greendensityperatom(1+NUMT,NUME+1), &
        &a_1(3), a_2(3), a_3(3), RPOINT(3), RPRIMEPOINT(3), FOURIERVEC(3), TPTS(3,NUMT), KPTS(3,NUMK), KPOINT(3),&
        &greendensity(2,NUME+1), densityperimpurity(1+NUMIMP,NUME+1)
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT,NUMK), GMATRIX(4*NUMT,4*NUMT), EZ, ENFRAC, GFK(NUMK,4*NUMT,4*NUMT),&
        &GREENR(4*NUMIMP,4*NUMIMP), FOUREXPONS(NUMK,NUMIMP**2), energies(NUME+1)
        character(len = 1) :: impdosyesorno

        ! This part sets up the E points to be used in plots concerning the Green's function
        ! At the moment it simply gives NUME real energies.

        print *, 'Should the program calculate the initial DoS for each future impurity site? y/n'
        170 read *, impdosyesorno

        if (impdosyesorno /= 'n' .and. impdosyesorno /= 'y') then
            print *, 'Invalid input. Please choose y or n.'
            goto 170
        endif

        energyintervals = (MAXVAL(EIGENVALUES) - MINVAL(EIGENVALUES))/NUME

        ! These are the "complex" energies E
        do i = 1, NUME+1
            energies(i) = dcmplx(MINVAL(EIGENVALUES) + energyintervals*(i-1), eta)
        end do

        ! This constructs the exponentials to be used in the Fourier transform
        checker = 1
        do j = 1, NUMIMP
            do i = 1, NUMIMP

                a = IMPPTSVAR(4,i)
                aprime = IMPPTSVAR(4,j)
                RPOINT = IMPPTSVAR(1,i)*a_1 + IMPPTSVAR(2,i)*a_2 + IMPPTSVAR(3,i)*a_3
                RPRIMEPOINT = IMPPTSVAR(1,j)*a_1 + IMPPTSVAR(2,j)*a_2 + IMPPTSVAR(3,j)*a_3

                FOURIERVEC(1:3) = RPRIMEPOINT(1:3) - RPOINT(1:3) + TPTS(1:3,aprime) - TPTS(1:3,a)

                do k = 1, NUMK

                    KPOINT = KPTS(1:3,k)
                    FOUREXPONS(k,checker) = exp(-CI*DOT_PRODUCT(KPOINT,FOURIERVEC))

                end do
                checker = checker + 1
            end do
        end do

        ! This part writes all the Fouriered Green function elements per energy at this text file
        open(1, file = 'greenhost.txt', action = 'write')

        do IE = 1, NUME+1 ! Begins a loop over the energies, in order to find G(E) for each E

            EZ = energies(IE)

            ! Startup values for the densities
            greendensityperatom(1,IE) = REAL(EZ)
            greendensity(1,IE) = REAL(EZ)
            greendensity(2,IE) = 0.0
            do i = 1, NUMT
                greendensityperatom(1+i,IE) = 0.0
            end do

            if (impdosyesorno == 'y') then
                densityperimpurity(1,IE) = REAL(EZ)
                do i = 1, NUMIMP
                    densityperimpurity(1+i,IE) = 0.0
                end do
            endif

            ! This initiates the calculation of the Green's function matrix G(α,α';E) per k-point
            do k = 1, NUMK ! k

                ! Set all G-matrix values equal to zero, so that the following summation can work
                GMATRIX(:,:) = (0.0,0.0)

                do i = 1, NUMT ! α
                    FTIMO = 4*(i-1)
                    do j = 1, NUMT ! α'
                        FTJMO = 4*(j-1)

                        do n = 1, 4*NUMT ! This is the sum over all eigenenergies per k

                            ENFRAC = (1.0/(EZ-EIGENVALUES(n,k)))
                            
                            GMATRIX(1 + FTIMO, 1 + FTJMO) = GMATRIX(1 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i,n,k)*CONJG(EIGENVECTORS(j,n,k)) ! 11
                            GMATRIX(1 + FTIMO, 2 + FTJMO) = GMATRIX(1 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i,n,k)*CONJG(EIGENVECTORS(j+NUMT,n,k)) ! 12
                            GMATRIX(1 + FTIMO, 3 + FTJMO) = GMATRIX(1 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i,n,k)*CONJG(EIGENVECTORS(j+2*NUMT,n,k)) ! 13
                            GMATRIX(1 + FTIMO, 4 + FTJMO) = GMATRIX(1 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i,n,k)*CONJG(EIGENVECTORS(j+3*NUMT,n,k)) ! 14

                            GMATRIX(2 + FTIMO, 1 + FTJMO) = GMATRIX(2 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+NUMT,n,k)*CONJG(EIGENVECTORS(j,n,k)) ! 21
                            GMATRIX(2 + FTIMO, 2 + FTJMO) = GMATRIX(2 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+NUMT,n,k)*CONJG(EIGENVECTORS(j+NUMT,n,k)) ! 22
                            GMATRIX(2 + FTIMO, 3 + FTJMO) = GMATRIX(2 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+NUMT,n,k)*CONJG(EIGENVECTORS(j+2*NUMT,n,k)) ! 23
                            GMATRIX(2 + FTIMO, 4 + FTJMO) = GMATRIX(2 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+NUMT,n,k)*CONJG(EIGENVECTORS(j+3*NUMT,n,k)) ! 24

                            GMATRIX(3 + FTIMO, 1 + FTJMO) = GMATRIX(3 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+2*NUMT,n,k)*CONJG(EIGENVECTORS(j,n,k)) ! 31
                            GMATRIX(3 + FTIMO, 2 + FTJMO) = GMATRIX(3 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+2*NUMT,n,k)*CONJG(EIGENVECTORS(j+NUMT,n,k)) ! 32
                            GMATRIX(3 + FTIMO, 3 + FTJMO) = GMATRIX(3 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+2*NUMT,n,k)*CONJG(EIGENVECTORS(j+2*NUMT,n,k)) ! 33
                            GMATRIX(3 + FTIMO, 4 + FTJMO) = GMATRIX(3 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+2*NUMT,n,k)*CONJG(EIGENVECTORS(j+3*NUMT,n,k)) ! 34

                            GMATRIX(4 + FTIMO, 1 + FTJMO) = GMATRIX(4 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+3*NUMT,n,k)*CONJG(EIGENVECTORS(j,n,k)) ! 41
                            GMATRIX(4 + FTIMO, 2 + FTJMO) = GMATRIX(4 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+3*NUMT,n,k)*CONJG(EIGENVECTORS(j+NUMT,n,k)) ! 42
                            GMATRIX(4 + FTIMO, 3 + FTJMO) = GMATRIX(4 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+3*NUMT,n,k)*CONJG(EIGENVECTORS(j+2*NUMT,n,k)) ! 43
                            GMATRIX(4 + FTIMO, 4 + FTJMO) = GMATRIX(4 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(i+3*NUMT,n,k)*CONJG(EIGENVECTORS(j+3*NUMT,n,k)) ! 44

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
                
                GFK(k,:,:) = GMATRIX ! This is a table that contains all G(α,α',k;E) per energy E
            
            end do ! ends k-loop

            ! Fourier transform of G(k) into G(r-r')
            ! This constructs a 4*NUMIMP x 4*NUMIMP G(r-r') matrix for each energy EZ

            ! Startup
            GREENR(:,:) = (0.0,0.0)
            checker = 1

            do j = 1, NUMIMP
                FTJMO = 4*(j-1)
                aprime = IMPPTSVAR(4,j)

                do i = 1, NUMIMP
                    FTIMO = 4*(i-1)
                    a = IMPPTSVAR(4,i)

                    do k = 1, NUMK

                        do m = 1, 4
                            do l = 1, 4
                                GREENR(m + FTIMO, l + FTJMO) = GREENR(m + FTIMO, l + FTJMO) +&
                                &(1.0/NUMK)*FOUREXPONS(k,checker)*GFK(k, m + 4*(a-1) , l + 4*(aprime-1))                    
                            end do
                        end do

                    end do

                    checker = checker + 1

                end do
            end do

            if (impdosyesorno == 'y') then
                do i = 1, NUMIMP ! Calculation of density for each impurity atom
                    FTIMO = 4*(i-1)
                    do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                        densityperimpurity(1+i,IE) = densityperimpurity(1+i,IE) -&
                        &(1.0/PI)*AIMAG(GREENR(j + FTIMO, j + FTIMO))
                    end do
                end do
            endif

            ! Writes the Green impurity elements on greenhost.txt
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1, '(F17.8,F17.8)') GREENR(i,j)
                end do
            end do

            do i = 1, NUMT ! Calculation of full density
                greendensityperatom(i+1,IE) = greendensityperatom(i+1,IE)/NUMK ! Normalization
                greendensity(2,IE) = greendensity(2,IE) + (1.0/NUMT)*greendensityperatom(1+i,IE)
            end do

        end do ! ends energies sum
        close(1) ! Closes the greenimp.txt file

        open(1, file = 'greendensityperatom.txt', action = 'write')
        do j = 1, NUME+1 ! Energies = Intervals + 1
            do i = 1, NUMT+1
                write (1,'(F17.8)',advance='no') greendensityperatom(i,j)
            end do
            write (1,*)
        end do
        close(1)

        ! The density per atom at the future impurity sites.
        if (impdosyesorno == 'y') then 
            open(1, file = 'hostdensities.txt', action = 'write')
            do j = 1, NUME+1 ! Energies = Intervals + 1
                do i = 1, NUMIMP+1
                    if (i == NUMIMP+1) then
                        write (1,'(F17.8)',advance='no') densityperimpurity(i,j)
                    else
                        write (1,'(F17.8, A)',advance='no') densityperimpurity(i,j), ','
                    endif
                end do
                write (1,*)
            end do
            close(1)
        endif

        open(1, file = 'energies.dat', action = 'write')
        do i = 1, NUME+1
            write (1, '(2F17.8)') energies(i)
        end do
        close(1)

        open(1, file = 'greendensity.txt', action = 'write')
        do j = 1, NUME+1 ! Energies = Intervals + 1
            write (1,'(2F17.8)') greendensity(1,j), greendensity(2,j)
        end do
        close(1)

    end subroutine FULLGREEN

    subroutine PARTIALGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NUME,eta)
        implicit none

        integer :: NUMT, NUMK, NUME, IE, i, j, k, n, NUMIMP, IMPPTSVAR(4,NUMIMP), min_val, max_val,&
        &l, m, a, aprime, FTIMO, FTJMO, checker, IMPATOMTYPE(NUMIMP), uniquecounter, NUMATOMS, IATOM, JATOM
        integer, allocatable, dimension(:) :: IMPATOMVALS, UNIQUEIMPATOMS
        real*8 :: EIGENVALUES(4*NUMT,NUMK), energyintervals, eta, densityperimpurity(1+NUMIMP,NUME+1),&
        &a_1(3), a_2(3), a_3(3), RPOINT(3), RPRIMEPOINT(3), FOURIERVEC(3), TPTS(3,NUMT), KPTS(3,NUMK), KPOINT(3)
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT,NUMK), GMATRIX(4*NUMT,4*NUMT), EZ, ENFRAC, GFK(NUMK,4*NUMT,4*NUMT),&
        &GREENR(4*NUMIMP,4*NUMIMP), FOUREXPONS(NUMK,NUMIMP**2), energies(NUME+1)
        character(len = 1) :: impdosyesorno

        ! This part sets up the E points to be used in plots concerning the Green's function
        ! At the moment it simply gives NUME real energies.

        print *, 'Should the program calculate the initial DoS for each future impurity site? y/n'
        170 read *, impdosyesorno

        if (impdosyesorno /= 'n' .and. impdosyesorno /= 'y') then
            print *, 'Invalid input. Please choose y or n.'
            goto 170
        endif
        
        energyintervals = (MAXVAL(EIGENVALUES) - MINVAL(EIGENVALUES))/NUME

        ! These are the "complex" energies E
        do i = 1, NUME+1
            energies(i) = dcmplx(MINVAL(EIGENVALUES) + energyintervals*(i-1), eta)
        end do

        ! Here we catalogue the different atom types, as numbered by their line on the basisvectors.dat file,
        ! that correspond to the impurities of the problem.
        IMPATOMTYPE(:) = IMPPTSVAR(4,:)
        allocate(UNIQUEIMPATOMS(NUMIMP)) ! Helper array

        uniquecounter = 0
        min_val = MINVAL(IMPATOMTYPE) - 1 ! -1.0 Is inserted for a case of complete degeneracy
        max_val = MAXVAL(IMPATOMTYPE)
        do while (min_val < max_val)
            uniquecounter = uniquecounter + 1
            min_val = MINVAL(IMPATOMTYPE, mask = IMPATOMTYPE > min_val)
            UNIQUEIMPATOMS(uniquecounter) = min_val
        enddo
        NUMATOMS = uniquecounter ! The number of basisvectos.dat atoms that enter the impurity problem. NUMATOMS < NUMT
        allocate(IMPATOMVALS(NUMATOMS))
        IMPATOMVALS = UNIQUEIMPATOMS(1:NUMATOMS)
        ! IMPATOMVALS now contains only the unique atom types that enter the impurity problem, in ascending order.

        ! The point of this subroutine, unlike the FULLGREEN subroutine, is to build the Green function only for the
        ! atoms that enter the impurity problem and not all the atoms in basisvectors.dat
        
        ! This constructs the exponentials to be used in the Fourier transform
        checker = 1
        do j = 1, NUMIMP
            do i = 1, NUMIMP

                a = IMPPTSVAR(4,i)
                aprime = IMPPTSVAR(4,j)
                RPOINT = IMPPTSVAR(1,i)*a_1 + IMPPTSVAR(2,i)*a_2 + IMPPTSVAR(3,i)*a_3
                RPRIMEPOINT = IMPPTSVAR(1,j)*a_1 + IMPPTSVAR(2,j)*a_2 + IMPPTSVAR(3,j)*a_3

                FOURIERVEC(1:3) = RPRIMEPOINT(1:3) - RPOINT(1:3) + TPTS(1:3,aprime) - TPTS(1:3,a)

                do k = 1, NUMK

                    KPOINT = KPTS(1:3,k)
                    FOUREXPONS(k,checker) = exp(-CI*DOT_PRODUCT(KPOINT,FOURIERVEC))

                end do
                checker = checker + 1
            end do
        end do

        ! This part writes all the Fouriered Green function elements per energy at this text file
        open(1, file = 'greenhost.txt', action = 'write')

        do IE = 1, NUME+1 ! Begins a loop over the energies, in order to find G(E) for each E

            EZ = energies(IE)

            if (impdosyesorno == 'y') then
                densityperimpurity(1,IE) = REAL(EZ)
                do i = 1, NUMIMP
                    densityperimpurity(1+i,IE) = 0.0
                end do
            endif

            ! This initiates the calculation of the Green's function matrix G(α,α';E) per k-point
            do k = 1, NUMK ! k

                ! Set all G-matrix values equal to zero, so that the following summation can work
                GMATRIX(:,:) = (0.0,0.0)

                do i = 1, NUMATOMS ! α
                    IATOM = IMPATOMVALS(i)
                    FTIMO = 4*(IATOM-1)

                    do j = 1, NUMATOMS ! α'
                        JATOM = IMPATOMVALS(j)
                        FTJMO = 4*(JATOM-1)

                        do n = 1, 4*NUMT ! This is the sum over all eigenenergies per k

                            ENFRAC = (1.0/(EZ-EIGENVALUES(n,k)))
                            
                            GMATRIX(1 + FTIMO, 1 + FTJMO) = GMATRIX(1 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM,n,k)*CONJG(EIGENVECTORS(JATOM,n,k)) ! 11
                            GMATRIX(1 + FTIMO, 2 + FTJMO) = GMATRIX(1 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM,n,k)*CONJG(EIGENVECTORS(JATOM+NUMT,n,k)) ! 12
                            GMATRIX(1 + FTIMO, 3 + FTJMO) = GMATRIX(1 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM,n,k)*CONJG(EIGENVECTORS(JATOM+2*NUMT,n,k)) ! 13
                            GMATRIX(1 + FTIMO, 4 + FTJMO) = GMATRIX(1 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM,n,k)*CONJG(EIGENVECTORS(JATOM+3*NUMT,n,k)) ! 14

                            GMATRIX(2 + FTIMO, 1 + FTJMO) = GMATRIX(2 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+NUMT,n,k)*CONJG(EIGENVECTORS(JATOM,n,k)) ! 21
                            GMATRIX(2 + FTIMO, 2 + FTJMO) = GMATRIX(2 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+NUMT,n,k)) ! 22
                            GMATRIX(2 + FTIMO, 3 + FTJMO) = GMATRIX(2 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+2*NUMT,n,k)) ! 23
                            GMATRIX(2 + FTIMO, 4 + FTJMO) = GMATRIX(2 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+3*NUMT,n,k)) ! 24

                            GMATRIX(3 + FTIMO, 1 + FTJMO) = GMATRIX(3 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+2*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM,n,k)) ! 31
                            GMATRIX(3 + FTIMO, 2 + FTJMO) = GMATRIX(3 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+2*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+NUMT,n,k)) ! 32
                            GMATRIX(3 + FTIMO, 3 + FTJMO) = GMATRIX(3 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+2*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+2*NUMT,n,k)) ! 33
                            GMATRIX(3 + FTIMO, 4 + FTJMO) = GMATRIX(3 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+2*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+3*NUMT,n,k)) ! 34

                            GMATRIX(4 + FTIMO, 1 + FTJMO) = GMATRIX(4 + FTIMO, 1 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+3*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM,n,k)) ! 41
                            GMATRIX(4 + FTIMO, 2 + FTJMO) = GMATRIX(4 + FTIMO, 2 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+3*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+NUMT,n,k)) ! 42
                            GMATRIX(4 + FTIMO, 3 + FTJMO) = GMATRIX(4 + FTIMO, 3 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+3*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+2*NUMT,n,k)) ! 43
                            GMATRIX(4 + FTIMO, 4 + FTJMO) = GMATRIX(4 + FTIMO, 4 + FTJMO) +&
                            &ENFRAC*EIGENVECTORS(IATOM+3*NUMT,n,k)*CONJG(EIGENVECTORS(JATOM+3*NUMT,n,k)) ! 44

                        end do

                    end do

                end do
                
                GFK(k,:,:) = GMATRIX ! This is a table that contains all G(α,α',k;E) per energy E
            
            end do ! ends k-loop

            ! Fourier transform of G(k) into G(r-r')
            ! This constructs a 4*NUMIMP x 4*NUMIMP G(r-r') matrix for each energy EZ

            ! Startup
            GREENR(:,:) = (0.0,0.0)
            checker = 1

            do j = 1, NUMIMP
                FTJMO = 4*(j-1)
                aprime = IMPPTSVAR(4,j)

                do i = 1, NUMIMP
                    FTIMO = 4*(i-1)
                    a = IMPPTSVAR(4,i)                

                    do k = 1, NUMK

                        do m = 1, 4
                            do l = 1, 4
                                GREENR(m + FTIMO, l + FTJMO) = GREENR(m + FTIMO, l + FTJMO) +&
                                &(1.0/NUMK)*FOUREXPONS(k,checker)*GFK(k, m + 4*(a-1) , l + 4*(aprime-1))                    
                            end do
                        end do

                    end do

                    checker = checker + 1

                end do
            end do

            if (impdosyesorno == 'y') then
                do i = 1, NUMIMP ! Calculation of density for each impurity atom
                    FTIMO = 4*(i-1)
                    do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                        densityperimpurity(1+i,IE) = densityperimpurity(1+i,IE) -&
                        &(1.0/PI)*AIMAG(GREENR(j + FTIMO, j + FTIMO))
                    end do
                end do
            endif

            ! Writes the Green impurity elements on greenhost.txt
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1, '(F17.8,F17.8)') GREENR(i,j)
                end do
            end do

        end do ! ends energies sum
        close(1) ! Closes the greenimp.txt file

        ! The density per atom at the future impurity sites.
        if (impdosyesorno == 'y') then 
            open(1, file = 'hostdensities.txt', action = 'write')
            do j = 1, NUME+1 ! Energies = Intervals + 1
                do i = 1, NUMIMP+1
                    if (i == NUMIMP+1) then
                        write (1,'(F17.8)',advance='no') densityperimpurity(i,j)
                    else
                        write (1,'(F17.8, A)',advance='no') densityperimpurity(i,j), ','
                    endif
                end do
                write (1,*)
            end do
            close(1)
        endif

        open(1, file = 'energies.dat', action = 'write')
        do i = 1, NUME+1
            write (1, '(2F17.8)') energies(i)
        end do
        close(1)

    end subroutine PARTIALGREEN

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

    subroutine IDENTIFIER(vecorints,b_1,b_2,b_3,PI,a_1,a_2,a_3,TPTS,NUMT,NUMIMP,IMPPTSVAR)
        implicit none

        integer :: NUMIMP, i, j, io, NUMT, counter
        real*8 :: b_1(3), b_2(3), b_3(3), PI, IMPPOINT(3), TPOINT(3), RPOINT(3), a_1(3), a_2(3), a_3(3), TPTS(3,NUMT),&
        &BASIS(3), eps
        real*8, allocatable, dimension(:,:) :: IMPPTS
        integer, allocatable, dimension(:,:) :: IMPPTSVAR
        character(len = 1) :: vecorints

        eps = 0.0001

        NUMIMP = 0
        if (vecorints == 'v') then
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

    function CROSS_PRODUCT(x,y) result(cross)
        implicit none
        real*8, dimension(3), intent(in) :: x, y
        real*8, dimension(3) :: cross

        cross(1) = x(2)*y(3) - x(3)*y(2)
        cross(2) = x(3)*y(1) - x(1)*y(3)
        cross(3) = x(1)*y(2) - x(2)*y(1)
		
    end function CROSS_PRODUCT

end program GREEN