program BANDS
    implicit none

    integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, m, n, IRLATT, IRLATTMAX, LWORK, INFO, NCELLS, &
    &NUMCHEMTYPES, JATOM, rcheck, MAXNEIGHB, NUMDIR, intpointer, checker, p
    integer, allocatable, dimension(:) :: CHEMTYPE, NEIGHBNUM, KINTS
    integer, allocatable, dimension(:,:) :: JTYPE
    real*8 :: ALAT, a_1(3), a_2(3), a_3(3), RMAX, R0, KPOINT(3), RPOINT(3), TTPRIME(3),&
    &chempot, T, PI, KB, lambda, TOL, HORINT, DIR(3), tempvalre, tempvalim
    real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, nu, nuzero, EIGENVALUES, VSUPCOND, magnet
    real*8, allocatable, dimension(:,:) :: TPTS, RLATT, BETA, LHOPS, PREFACTORS, NNDIS, HOPPVALS, INPOINT, &
    &OUTPOINT, ATWEIGHTS
    real*8, allocatable, dimension(:,:,:) :: RCONNECT
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2)
    complex*16, allocatable, dimension(:) :: WORK, DELTA
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, EIGENVECTORS, HAMILTONIANPREP
    character(len = 1) :: weightspointer

    TOL = 0.0001 ! The fault tolerance for lattice vectors' norms

    call CONSTANTS(IdentityPauli,xPauli,yPauli,zPauli,CI,PI,KB)

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
    close(1)

    ! This is because the lattice vectors are inserted through Bravais coordinates
    a_1 = ALAT*a_1
    a_2 = ALAT*a_2
    a_3 = ALAT*a_3

    IRLATTMAX = (2*NCELLS+1)**3 ! Configures how many neighbouring cells are taken into account

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
    allocate(ATWEIGHTS(NUMT,4*NUMT))

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

    ! We finally calculate the hopping elements.
    allocate(HOPPVALS(MAXNEIGHB,NUMT))
    do i = 1, NUMT
        do j = 1, NEIGHBNUM(i)

            RPOINT = RCONNECT(1:3,j,i)

            JATOM = JTYPE(j,i)
            lambda = PREFACTORS(CHEMTYPE(i),CHEMTYPE(JATOM))
                
            HOPPVALS(j,i) = -lambda*exp(-norm2(RPOINT)/R0)

        end do
    end do

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

    17 print *, 'Calculate the contribution of each atom in the band diagram? y/n'
    read (*, '(A)') weightspointer

    if (weightspointer /= 'y' .and. weightspointer /= 'n') then
        print *, 'Invalid value. Please enter y or n.'
        goto 17
    endif

    HORINT = 0.0 ! This is a parameter that adjusts the widths for each k-dimension window at the band diagram
    
    open(2, file = 'bands.txt', action = 'write')
    do i = 1, NUMDIR
        
        DIR = OUTPOINT(1:3,i) - INPOINT(1:3,i) ! Sets the direction along which we calculate K points
        KPOINT = INPOINT(1:3,i) ! Startup of each k

        ! In order to avoid taking some points twice
        if (i /= NUMDIR) then
            if (norm2(OUTPOINT(1:3,i) - INPOINT(1:3,i+1)) < TOL) then
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

            if (weightspointer == 'n') then
                ! N because we only want eigenvalues and not eigenvectors
                call zheev ('N', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO)

                if (intpointer == 0) then
                    do n = 1, 4*NUMT
                        write (2,'(F17.8, A, F17.8)', advance='no') HORINT, ',', W(n)
                        write (2,*)
                    end do
                else 
                    do n = 1, 4*NUMT
                        write (2,'(5F17.8)', advance='no') KPOINT, HORINT, W(n)
                        write (2,*)
                    end do
                endif
            else

                ! V because we also want the eigenvectors
                call zheev ('V', 'U', 4*NUMT, HAMILTONIAN, 4*NUMT, W, WORK, LWORK, RWORK, INFO)

                ! That part calculates the contribution of each atom to the band structure by asigning the
                ! square of its wavefunction as a weight.
                ATWEIGHTS(:,:) = 0.0
                do m = 1, NUMT
                    do n = 1, 4*NUMT

                        do p = 1, 4 ! Sums over the 4 components of the Ψ_m(n) spinor
                            ATWEIGHTS(m,n) = ATWEIGHTS(m,n) + (abs(HAMILTONIAN(m + (p-1)*NUMT, n)))**2
                        end do

                    end do
                end do

                if (intpointer == 0) then
                    do n = 1, 4*NUMT
                        write (2,'(F17.8, A, F17.8)', advance='no') HORINT, ',', W(n)
                        do m = 1, NUMT
                            write (2,'(A, F17.8)', advance='no') ',', ATWEIGHTS(m,n)
                        end do
                        write (2,*)
                    end do
                else 
                    do n = 1, 4*NUMT
                        write (2,'(5F17.8)', advance='no') KPOINT, HORINT, W(n)
                        do m = 1, NUMT
                            write (2,'(F17.8)', advance='no') ATWEIGHTS(m,n)
                        end do
                        write (2,*)
                    end do
                endif

            endif

            if (j /= checker) then
                KPOINT = KPOINT + (1.0/KINTS(i))*DIR
                HORINT = HORINT + (1.0/KINTS(i))*norm2(DIR)
            endif

        end do

    end do
    close(2)

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

    subroutine RANDOMKHAM(KPOINT,HAMILTONIANPREP,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,RCONNECT,NUMT,HAMILTONIAN)
        implicit none
        real*8, intent(in) :: KPOINT(3)
        integer :: NUMT, i, jneighb, MAXNEIGHB, NEIGHBNUM(NUMT), JATOM, JTYPE(MAXNEIGHB,NUMT)
        real*8 :: hopping, HOPPVALS(MAXNEIGHB,NUMT), RPOINT(3), RCONNECT(3,MAXNEIGHB,NUMT)
        complex*16 :: expon, HAMILTONIAN(4*NUMT,4*NUMT), HAMILTONIANPREP(4*NUMT,4*NUMT)

        HAMILTONIAN(:,:) = HAMILTONIANPREP(:,:)

        do i = 1, NUMT
            do jneighb = 1, NEIGHBNUM(i)

                RPOINT = RCONNECT(1:3,jneighb,i)
                JATOM = JTYPE(jneighb,i)

                hopping = HOPPVALS(jneighb,i)

                expon = exp(-CI*DOT_PRODUCT(KPOINT,RPOINT))

                HAMILTONIAN(i,JATOM) = HAMILTONIAN(i,JATOM) + expon*hopping
                HAMILTONIAN(i + NUMT,JATOM + NUMT) = HAMILTONIAN(i + NUMT,JATOM + NUMT) + expon*hopping
                HAMILTONIAN(i + 2*NUMT,JATOM + 2*NUMT) = HAMILTONIAN(i + 2*NUMT,JATOM + 2*NUMT) - CONJG(expon*hopping)
                HAMILTONIAN(i + 3*NUMT,JATOM + 3*NUMT) = HAMILTONIAN(i + 3*NUMT,JATOM + 3*NUMT) - CONJG(expon*hopping)

            end do
        end do        

    end subroutine RANDOMKHAM

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

end program BANDS