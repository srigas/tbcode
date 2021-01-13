program PFAFFIAN
    implicit none
    integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, k, IRLATT, IRLATTMAX, kcounter, NCELLS, &
    &NUMCHEMTYPES, JATOM, rcheck, MAXNEIGHB, chemtypeval, XINDEX, YINDEX, ZINDEX, checker, HALF, SCLAYERS, &
    &num, denom, BSTEPS, BCOUNTER, INFO, LWORK
    integer, allocatable, dimension(:) :: CHEMTYPE, NEIGHBNUM, IWORK
    integer, allocatable, dimension(:,:) :: JTYPE
    real*8 :: ALAT, a_1(3), a_2(3), a_3(3), RMAX, R0, KPOINT(3), RPOINT(3), TTPRIME(3), &
    &chempot, PI, KB, b_1(3), b_2(3), b_3(3), lambda, TOL, tempvalre, tempvalim, ezeroval, ulcnval, &
    &nuzeroval, zparam, yparam, randparam, ROTAXIS(3), theta, magB, TOTSIGN, BINC, BFINAL
    real*8, allocatable, dimension(:) :: E0, ULCN, nu, nuzero, readnu
    real*8, allocatable, dimension(:,:) :: KPTS, TPTS, RLATT, BETA, LHOPS, PREFACTORS, HOPPVALS, &
    &MAGBETS, PFSIGNS
    real*8, allocatable, dimension(:,:,:) :: RCONNECT
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2)
    complex*16, allocatable, dimension(:) :: DELTA, readdelta
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, HAMILTONIANPREP
    complex*16, allocatable, dimension(:,:,:) :: EXPONS
    character(len = 1) :: manyruns

    real*8, allocatable, dimension(:) :: RWORK
    complex*16 :: PFAFZERO(2), PFAFPI(2)
    complex*16, allocatable, dimension(:) :: WORK
    complex*16, allocatable, dimension(:,:) :: UMATRIX, UHERMMATRIX, ALPHAZERO, ALPHAPI

    TOL = 0.0001 ! The fault tolerance for lattice vectors' norms

    call CONSTANTS(IdentityPauli,xPauli,yPauli,zPauli,CI,PI,KB)

    open(1, file = 'pfafconfig.dat', action = 'read')
    read(1,*) ALAT
    read(1,*) a_1
    read(1,*) a_2
    read(1,*) a_3
    read(1,*) NUMKX,NUMKY,NUMKZ
    read(1,*) RMAX
    read(1,*) R0 
    read(1,*) NCELLS
    read(1,*) XINDEX
    read(1,*) YINDEX
    read(1,*) ZINDEX
    read(1,*) chemtypeval
    read(1,*) ezeroval
    read(1,*) ulcnval
    read(1,*) nuzeroval
    read(1,*) ROTAXIS
    read(1,*) num, denom
    read(1,*) magB
    read(1,*) manyruns
    close(1)

    theta = (num*PI)/denom

    NUMT = XINDEX*(2*YINDEX+1)*ZINDEX

    ! This is because the lattice vectors are inserted through Bravais coordinates
    a_1 = ALAT*a_1
    a_2 = ALAT*a_2
    a_3 = ALAT*a_3

    IRLATTMAX = (2*NCELLS+1)**3 ! Configures how many neighbouring cells are taken into account

    call BZ(a_1,a_2,a_3,NUMKX,NUMKY,NUMKZ,PI,KPTS,NUMK,b_1,b_2,b_3) ! We create the Brillouin Zone's k-points to be used later on

    allocate(CHEMTYPE(NUMT))
    allocate(E0(NUMT))
    allocate(ULCN(NUMT))
    allocate(nuzero(NUMT))

    do i = 1, NUMT
        CHEMTYPE(i) = chemtypeval
        E0(i) = ezeroval
        ULCN(i) = ulcnval
        nuzero(i) = nuzeroval
    end do

    allocate(TPTS(3,NUMT))

    allocate(readnu(ZINDEX))
    allocate(readdelta(ZINDEX))

    ! The following reads the number of entries from scresults.dat
    SCLAYERS = 0
    open (2, file = 'scresults.dat', action = 'read')
    do
        read(2,*,iostat=io)
        if (io/=0) exit
        SCLAYERS = SCLAYERS + 1
    end do
    close(2)

    SCLAYERS = SCLAYERS - 1

    HALF = ZINDEX/2
    open (1, file = 'scresults.dat', action = 'read')
    read(1,*) chempot
    do i = 1, HALF
        read(1,*) readnu(i), randparam, tempvalre, tempvalim
        readdelta(i) = dcmplx(tempvalre, tempvalim)
    end do
    do i = HALF+1, SCLAYERS-HALF
        read(1,*)
    end do
    do i = HALF+1, ZINDEX
        read(1,*) readnu(i), randparam, tempvalre, tempvalim
        readdelta(i) = dcmplx(tempvalre, tempvalim)
    end do
    close(1)

    ! The coordinates of the imp atoms are (:,0,ZINDEX-1), where : = 0, XINDEX-1
    ! Their NUMT starts from (ZINDEX-1)*(2*YINDEX+1)*XINDEX + XINDEX*YINDEX + 1
    ! and ends at (ZINDEX-1)*(2*YINDEX+1)*XINDEX + XINDEX*YINDEX + XINDEX

    allocate(nu(NUMT))
    allocate(DELTA(NUMT))

    checker = 1

    do k = 1, ZINDEX ! This labels different layers

        zparam = (k-1)*1.D0

        do j = -YINDEX, YINDEX

            yparam = -1.D0*j

                do i = 1, XINDEX

                    nu(checker) = readnu(k)
                    DELTA(checker) = readdelta(k)

                    TPTS(1,checker) = 1.D0*(i-1)
                    TPTS(2,checker) = yparam
                    TPTS(3,checker) = zparam

                    checker = checker + 1

                end do
        end do
    end do

    NUMCHEMTYPES = MAXVAL(CHEMTYPE)
    allocate(LHOPS(NUMCHEMTYPES,NUMCHEMTYPES))
    allocate(PREFACTORS(NUMCHEMTYPES,NUMCHEMTYPES))

    ! The *4 factors are now due to spin and particle-hole
    allocate(HAMILTONIAN(4*NUMT,4*NUMT))
    allocate(HAMILTONIANPREP(4*NUMT,4*NUMT))

    allocate(UMATRIX(4*NUMT,4*NUMT))
    allocate(UHERMMATRIX(4*NUMT,4*NUMT))
    allocate(ALPHAZERO(4*NUMT,4*NUMT))
    allocate(ALPHAPI(4*NUMT,4*NUMT))

    call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) ! Here we construct the RLATT Matrix consisting of the lattice sites

    call HOPPS(RLATT,IRLATTMAX,R0,NUMCHEMTYPES,LHOPS,NUMT,CHEMTYPE,TPTS,PREFACTORS)

    ! Setup of the unitary transformation matrices, U, U^\dagger
    do i = 1, NUMT

        UMATRIX(i, i) = 0.5*(1.0,0.0)
        UMATRIX(i, i + NUMT) = 0.5*(0.0,1.0)
        UMATRIX(i, i + 2*NUMT) = 0.5*(-1.0,0.0)
        UMATRIX(i, i + 3*NUMT) = 0.5*(0.0,-1.0)

        UMATRIX(i + NUMT, i) = 0.5*(0.0,1.0)
        UMATRIX(i + NUMT, i + NUMT) = 0.5*(1.0,0.0)
        UMATRIX(i + NUMT, i + 2*NUMT) = 0.5*(0.0,1.0)
        UMATRIX(i + NUMT, i + 3*NUMT) = 0.5*(1.0,0.0)

        UMATRIX(i + 2*NUMT, i) = 0.5*(0.0,1.0)
        UMATRIX(i + 2*NUMT, i + NUMT) = 0.5*(-1.0,0.0)
        UMATRIX(i + 2*NUMT, i + 2*NUMT) = 0.5*(0.0,1.0)
        UMATRIX(i + 2*NUMT, i + 3*NUMT) = 0.5*(-1.0,0.0)

        UMATRIX(i + 3*NUMT, i) = 0.5*(1.0,0.0)
        UMATRIX(i + 3*NUMT, i + NUMT) = 0.5*(0.0,-1.0)
        UMATRIX(i + 3*NUMT, i + 2*NUMT) = 0.5*(-1.0,0.0)
        UMATRIX(i + 3*NUMT, i + 3*NUMT) = 0.5*(0.0,1.0)

        UHERMMATRIX(i, i) = 0.5*(1.0,0.0)
        UHERMMATRIX(i, i + NUMT) = 0.5*(0.0,-1.0)
        UHERMMATRIX(i, i + 2*NUMT) = 0.5*(0.0,-1.0)
        UHERMMATRIX(i, i + 3*NUMT) = 0.5*(1.0,0.0)

        UHERMMATRIX(i + NUMT, i) = 0.5*(0.0,-1.0)
        UHERMMATRIX(i + NUMT, i + NUMT) = 0.5*(1.0,0.0)
        UHERMMATRIX(i + NUMT, i + 2*NUMT) = 0.5*(-1.0,0.0)
        UHERMMATRIX(i + NUMT, i + 3*NUMT) = 0.5*(0.0,1.0)

        UHERMMATRIX(i + 2*NUMT, i) = 0.5*(-1.0,0.0)
        UHERMMATRIX(i + 2*NUMT, i + NUMT) = 0.5*(0.0,-1.0)
        UHERMMATRIX(i + 2*NUMT, i + 2*NUMT) = 0.5*(0.0,-1.0)
        UHERMMATRIX(i + 2*NUMT, i + 3*NUMT) = 0.5*(-1.0,0.0)

        UHERMMATRIX(i + 3*NUMT, i) = 0.5*(0.0,1.0)
        UHERMMATRIX(i + 3*NUMT, i + NUMT) = 0.5*(1.0,0.0)
        UHERMMATRIX(i + 3*NUMT, i + 2*NUMT) = 0.5*(-1.0,0.0)
        UHERMMATRIX(i + 3*NUMT, i + 3*NUMT) = 0.5*(0.0,-1.0)

    end do

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

    allocate(BETA(3,NUMT))
    BETA(:,:) = 0.d0

    allocate(MAGBETS(3,XINDEX))

    if (manyruns == 'y') then
        print *, 'Enter the final magnitude for B.'
        read *, BFINAL

        print *, 'Enter the number of steps.'
        read *, BSTEPS

        BCOUNTER = 1
        BINC = (BFINAL-magB)/(BSTEPS-1)
        allocate(PFSIGNS(BSTEPS,2))
        PFSIGNS = 0.D0
    endif

    ! Configuration of ALTPFAF
    LWORK = 2*4*NUMT-1
    allocate(IWORK(4*NUMT))
    allocate(WORK(LWORK))
    allocate(RWORK(4*NUMT-1))
    
    180 call MAGCHAIN(XINDEX,ROTAXIS,theta,magB,MAGBETS)

    checker = (ZINDEX-1)*(2*YINDEX+1)*XINDEX + XINDEX*YINDEX + 1
    do i = 1, XINDEX
        DELTA(checker) = (0.D0,0.D0) ! <-- This needs to be checked
        BETA(:,checker) = MAGBETS(:,i)
        checker = checker + 1
    end do

    ! This prints the data we have gathered, in order to check if the atoms are okay.
    !open (1, file = 'basisvectors.dat', action = 'write')
    !do i = 1, NUMT
    !    write(1,'(3F7.1, I4, 3F7.1, 6F13.7)') TPTS(1:3,i), CHEMTYPE(i), E0(i), ULCN(i), nuzero(i), BETA(1,i), BETA(2,i), BETA(3,i), nu(i), DELTA(i)
    !end do
    !close(1)

    call HAMPREP(NUMT,xPauli,yPauli,zPauli,IdentityPauli,chempot,E0,ULCN,nu,nuzero,BETA,DELTA,HAMILTONIANPREP)
    
    ! ------------ k = 0 CALCULATIONS -----------------------------------------------!
    kcounter = NUMK/2
    call FINALHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)

    call MAKEALPHA(NUMT,UMATRIX,UHERMMATRIX,HAMILTONIAN,ALPHAZERO)
    
    call ZSKPF10('U', 'P', 4*NUMT, ALPHAZERO, 4*NUMT, PFAFZERO, IWORK, WORK, LWORK, RWORK, INFO)

    ! --------------------------------------------------------------------------------!

    ! ------------ k = pi CALCULATIONS -----------------------------------------------!
    kcounter = NUMK
    call FINALHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)

    call MAKEALPHA(NUMT,UMATRIX,UHERMMATRIX,HAMILTONIAN,ALPHAPI)

    call ZSKPF10('U', 'P', 4*NUMT, ALPHAPI, 4*NUMT, PFAFPI, IWORK, WORK, LWORK, RWORK, INFO)
    
    ! --------------------------------------------------------------------------------!

    ! ------------------ Majorana Number Calculations --------------------------------!

    TOTSIGN = REAL(PFAFPI(1))*REAL(PFAFZERO(1))

    if (TOTSIGN > 0.0) then
        TOTSIGN = 1.0
    else if (TOTSIGN < 0.0) then
        TOTSIGN = -1.0
    endif

    print *, 'Pfaffian for k = 0 is', PFAFZERO(1), 'times 10^', PFAFZERO(2)
    print *, 'Pfaffian for k = pi is', PFAFPI(1), 'times 10^', PFAFPI(2)

    if (manyruns == 'y') then

        if (BCOUNTER < BSTEPS) then

            PFSIGNS(BCOUNTER,1) = magB
            PFSIGNS(BCOUNTER,2) = TOTSIGN

            magB = magB + BINC
            BCOUNTER = BCOUNTER + 1

            print *, 'The Majorana number is equal to', TOTSIGN
            print *, 'Finished iteration No.', BCOUNTER-1
            print *, '-----------------------------------------------------------------'

            goto 180
        else
            open(1, file = 'Pfaffians.txt', action = 'write')
                do i = 1, BSTEPS
                    write (1,*) PFSIGNS(i,1), ',', PFSIGNS(i,2)
                end do
            close(1)
        endif

    else

        print *, 'The Majorana number is equal to ', TOTSIGN

    endif

    contains

    subroutine MAGCHAIN(NUM,ROTAXIS,theta,magB,NEWBETA)
        implicit none

        integer :: NUM, i
        real*8 :: theta, ROTAXIS(3), NEWBETA(3,NUM), magB, initheta, TOL

        TOL = 0.000001 ! The tolerance for the "if" checks.

        NEWBETA(:,:) = 0.d0

        if (abs(ROTAXIS(1) - 1.0) < TOL .and. abs(ROTAXIS(2)) < TOL .and. abs(ROTAXIS(3)) < TOL) then

            initheta = 0.d0

            do i = 1, NUM
                if (i == 1) then
                    NEWBETA(1,i) = 0.0
                    NEWBETA(2,i) = magB*SIN(initheta)
                    NEWBETA(3,i) = magB*COS(initheta)
                else
                    NEWBETA(1,i) = 0.0
                    NEWBETA(2,i) = magB*SIN(initheta + (i-1)*theta)
                    NEWBETA(3,i) = magB*COS(initheta + (i-1)*theta)
                endif
            end do
        else if (abs(ROTAXIS(1)) < TOL .and. abs(ROTAXIS(2) - 1.0) < TOL .and. abs(ROTAXIS(3)) < TOL) then

            initheta = 0.d0

            do i = 1, NUM
                if (i == 1) then
                    NEWBETA(1,i) = magB*COS(initheta)
                    NEWBETA(2,i) = 0.0 
                    NEWBETA(3,i) = magB*SIN(initheta)
                else
                    NEWBETA(1,i) = magB*COS(initheta + (i-1)*theta)
                    NEWBETA(2,i) = 0.0
                    NEWBETA(3,i) = magB*SIN(initheta + (i-1)*theta)
                endif
            end do
        else if (abs(ROTAXIS(1)) < TOL .and. abs(ROTAXIS(2)) < TOL .and. abs(ROTAXIS(3) - 1.0) < TOL) then

            initheta = 0.d0

            do i = 1, NUM
                if (i == 1) then
                    NEWBETA(1,i) = magB*COS(initheta)
                    NEWBETA(2,i) = magB*SIN(initheta)
                    NEWBETA(3,i) = 0.0
                else
                    NEWBETA(1,i) = magB*COS(initheta + (i-1)*theta)
                    NEWBETA(2,i) = magB*SIN(initheta + (i-1)*theta)
                    NEWBETA(3,i) = 0.0
                endif
            end do
        else
            ! This is the case where we want to rotate B along an arbitrary axis.
            print *, 'Invalid input for the magchain rotation axis.'

        endif

    end subroutine MAGCHAIN

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

    subroutine FINALHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)
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

    end subroutine FINALHAM

    subroutine MAKEALPHA(NUMT,UMATRIX,UHERMMATRIX,HAMILTONIAN,ALPHAMATRIX)
        implicit none
        
        integer :: NUMT, i, j
        complex*16 :: HAMILTONIAN(4*NUMT,4*NUMT), ALPHAMATRIX(4*NUMT,4*NUMT), UMATRIX(4*NUMT,4*NUMT), &
        &UHERMMATRIX(4*NUMT,4*NUMT), HELPERMAT(4*NUMT,4*NUMT), scalbeta, scalpha

        scalpha = (1.d0,0.d0)
        scalbeta = (0.d0,0.d0)

        call ZGEMM('N', 'N', 4*NUMT, 4*NUMT, 4*NUMT, scalpha, UMATRIX, 4*NUMT, HAMILTONIAN, 4*NUMT, scalbeta, &
        &HELPERMAT, 4*NUMT)
        HAMILTONIAN = (0.D0,0.D0)
        call ZGEMM('N', 'N', 4*NUMT, 4*NUMT, 4*NUMT, scalpha, HELPERMAT, 4*NUMT, UHERMMATRIX, 4*NUMT, scalbeta, &
        &HAMILTONIAN, 4*NUMT)
        ALPHAMATRIX = (0.D0,-1.D0)*HAMILTONIAN

    end subroutine MAKEALPHA

    subroutine BZ(a,b,c,N_x,N_y,N_z,PI,KPTS,Ntot,b_1,b_2,b_3)
        implicit none

        integer :: N_x, N_y, N_z, Ntot, counter, i, j, k
        real*8 :: a(3), b(3), c(3), b_1(3), b_2(3), b_3(3), flx, fly, flz, c_1, c_2, c_3
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
        
        flx = dfloat(N_x)
        fly = dfloat(N_y)
        flz = dfloat(N_z)
        
        c_1 = -flx/2.0+1.0
                    
        do i = 1, N_x
            c_2 = -fly/2.0
            do j = 1, N_y
                c_3 = -flz/2.0
                do k = 1, N_z
                    
                    KPTS(1, counter) = (b_1(1)*c_1)/flx
                    KPTS(2, counter) = 0.D0
                    KPTS(3, counter) = 0.D0
                    counter = counter + 1
                    c_3 = c_3 + 1.0
                end do
                c_2 = c_2 + 1.0
            end do
            c_1 = c_1 + 1.0
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
  
    subroutine ZMUL10(A, B)
        implicit none
        DOUBLE COMPLEX A(2)
        DOUBLE COMPLEX B
        DOUBLE PRECISION EXPONENT, DLAMCH
        INTEGER IEXPONENT
        EXTERNAL DLAMCH
  
        A(1) = A(1) * B
  
        IF( A(1).EQ.(0.0D0, 0.0D0) ) THEN
           A(1) = (0.0D0, 0.0D0)
           A(2) = (0.0D0, 0.0D0)
        ELSE
           EXPONENT = LOG10(ABS(A(1)))
  
           IF (EXPONENT .GE. 0D0) THEN
              IEXPONENT=INT(EXPONENT)
           ELSE
              IEXPONENT=INT(EXPONENT)-1
           END IF
  
        !     If 10**IEXPONENT is smaller than the safe minimum for inversion,
        !     the number is considered 0
           IF( DLAMCH( "S" ) > 10D0**IEXPONENT ) THEN
              A(1) = (0.0D0, 0.0D0)
              A(2) = (0.0D0, 0.0D0)
           ELSE
              A(1) = A(1)/10D0**IEXPONENT
              A(2) = A(2) + IEXPONENT
           END IF
        END IF
  
    end subroutine ZMUL10

    subroutine ZSKMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

        implicit none

        DOUBLE COMPLEX ALPHA,BETA
        INTEGER INCX,INCY,LDA,N
        CHARACTER UPLO
        
        DOUBLE COMPLEX A(LDA,*),X(*),Y(*)

        DOUBLE COMPLEX ONE
        PARAMETER (ONE= (1.0D+0,0.0D+0))
        DOUBLE COMPLEX ZERO
        PARAMETER (ZERO= (0.0D+0,0.0D+0))
        
        DOUBLE COMPLEX TEMP1,TEMP2
        INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
        
        LOGICAL LSAME
        EXTERNAL LSAME
        
        EXTERNAL XERBLA
        
        INTRINSIC MAX
        
        INFO = 0
        IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
            INFO = 1
        ELSE IF (N.LT.0) THEN
            INFO = 2
        ELSE IF (LDA.LT.MAX(1,N)) THEN
            INFO = 5
        ELSE IF (INCX.EQ.0) THEN
            INFO = 7
        ELSE IF (INCY.EQ.0) THEN
            INFO = 10
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('ZHEMV ',INFO)
            RETURN
        END IF
        
        IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
        
        IF (INCX.GT.0) THEN
            KX = 1
        ELSE
            KX = 1 - (N-1)*INCX
        END IF
        IF (INCY.GT.0) THEN
            KY = 1
        ELSE
            KY = 1 - (N-1)*INCY
        END IF
        
        IF (BETA.NE.ONE) THEN
            IF (INCY.EQ.1) THEN
                IF (BETA.EQ.ZERO) THEN
                    DO 10 I = 1,N
                        Y(I) = ZERO
        10             CONTINUE
                ELSE
                    DO 20 I = 1,N
                        Y(I) = BETA*Y(I)
        20             CONTINUE
                END IF
            ELSE
                IY = KY
                IF (BETA.EQ.ZERO) THEN
                    DO 30 I = 1,N
                        Y(IY) = ZERO
                        IY = IY + INCY
        30             CONTINUE
                ELSE
                    DO 40 I = 1,N
                        Y(IY) = BETA*Y(IY)
                        IY = IY + INCY
        40             CONTINUE
                END IF
            END IF
        END IF
        IF (ALPHA.EQ.ZERO) RETURN
        IF (LSAME(UPLO,'U')) THEN
            
            IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
                DO 60 J = 1,N
                    TEMP1 = ALPHA*X(J)
                    TEMP2 = ZERO
                    DO 50 I = 1,J - 1
                        Y(I) = Y(I) + TEMP1*A(I,J)
                        TEMP2 = TEMP2 - A(I,J)*X(I)
        50             CONTINUE
                    Y(J) = Y(J) + ALPHA*TEMP2
        60         CONTINUE
            ELSE
                JX = KX
                JY = KY
                DO 80 J = 1,N
                    TEMP1 = ALPHA*X(JX)
                    TEMP2 = ZERO
                    IX = KX
                    IY = KY
                    DO 70 I = 1,J - 1
                        Y(IY) = Y(IY) + TEMP1*A(I,J)
                        TEMP2 = TEMP2 - A(I,J)*X(IX)
                        IX = IX + INCX
                        IY = IY + INCY
        70             CONTINUE
                    Y(JY) = Y(JY) + ALPHA*TEMP2
                    JX = JX + INCX
                    JY = JY + INCY
        80         CONTINUE
            END IF
        ELSE
            
            IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
                DO 100 J = 1,N
                    TEMP1 = ALPHA*X(J)
                    TEMP2 = ZERO
                    DO 90 I = J + 1,N
                        Y(I) = Y(I) + TEMP1*A(I,J)
                        TEMP2 = TEMP2 - A(I,J)*X(I)
        90             CONTINUE
                    Y(J) = Y(J) + ALPHA*TEMP2
        100         CONTINUE
            ELSE
                JX = KX
                JY = KY
                DO 120 J = 1,N
                    TEMP1 = ALPHA*X(JX)
                    TEMP2 = ZERO
                    IX = JX
                    IY = JY
                    DO 110 I = J + 1,N
                        IX = IX + INCX
                        IY = IY + INCY
                        Y(IY) = Y(IY) + TEMP1*A(I,J)
                        TEMP2 = TEMP2 - A(I,J)*X(IX)
        110             CONTINUE
                    Y(JY) = Y(JY) + ALPHA*TEMP2
                    JX = JX + INCX
                    JY = JY + INCY
        120         CONTINUE
            END IF
        END IF
        
        RETURN
      
    end subroutine ZSKMV

    subroutine ZSKR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

        implicit none

        DOUBLE COMPLEX ALPHA
        DOUBLE COMPLEX BETA
        INTEGER K,LDA,LDB,LDC,N
        CHARACTER TRANS,UPLO
        
        DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)

        LOGICAL LSAME
        EXTERNAL LSAME
        
        EXTERNAL XERBLA
        
        INTRINSIC MAX
        
        DOUBLE COMPLEX TEMP1,TEMP2
        INTEGER I,INFO,J,L,NROWA
        LOGICAL UPPER
        
        DOUBLE COMPLEX ONE
        PARAMETER (ONE= (1.0D+0,0.0D+0))
        DOUBLE COMPLEX ZERO
        PARAMETER (ZERO= (0.0D+0,0.0D+0))
        
        IF (LSAME(TRANS,'N')) THEN
            NROWA = N
        ELSE
            NROWA = K
        END IF
        UPPER = LSAME(UPLO,'U')
        
        INFO = 0
        IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
            INFO = 1
        ELSE IF ((.NOT.LSAME(TRANS,'N')) .AND. (.NOT.LSAME(TRANS,'T'))) THEN
            INFO = 2
        ELSE IF (N.LT.0) THEN
            INFO = 3
        ELSE IF (K.LT.0) THEN
            INFO = 4
        ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
            INFO = 7
        ELSE IF (LDB.LT.MAX(1,NROWA)) THEN
            INFO = 9
        ELSE IF (LDC.LT.MAX(1,N)) THEN
            INFO = 12
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('ZSKR2K',INFO)
            RETURN
        END IF
        
        IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
        
        IF (ALPHA.EQ.ZERO) THEN
            IF (UPPER) THEN
                IF (BETA.EQ.ZERO) THEN
                    DO 20 J = 1,N
                        DO 10 I = 1,J
                            C(I,J) = ZERO
        10                 CONTINUE
        20             CONTINUE
                ELSE
                    DO 40 J = 1,N
                        DO 30 I = 1,J - 1
                            C(I,J) = BETA*C(I,J)
        30                 CONTINUE
                        C(J,J) = ZERO
        40             CONTINUE
                END IF
            ELSE
                IF (BETA.EQ.ZERO) THEN
                    DO 60 J = 1,N
                        DO 50 I = J,N
                            C(I,J) = ZERO
        50                 CONTINUE
        60             CONTINUE
                ELSE
                    DO 80 J = 1,N
                        C(J,J) = ZERO
                        DO 70 I = J + 1,N
                            C(I,J) = BETA*C(I,J)
        70                 CONTINUE
        80             CONTINUE
                END IF
            END IF
            RETURN
        END IF
        
        IF (LSAME(TRANS,'N')) THEN
            
            IF (UPPER) THEN
                DO 130 J = 1,N
                    IF (BETA.EQ.ZERO) THEN
                        DO 90 I = 1,J
                            C(I,J) = ZERO
        90                 CONTINUE
                    ELSE IF (BETA.NE.ONE) THEN
                        DO 100 I = 1,J - 1
                            C(I,J) = BETA*C(I,J)
        100                 CONTINUE
                        C(J,J) = ZERO
                    ELSE
                        C(J,J) = ZERO
                    END IF
                    DO 120 L = 1,K
                        IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                            TEMP1 = ALPHA*B(J,L)
                            TEMP2 = ALPHA*A(J,L)
                            DO 110 I = 1,J - 1
                                C(I,J) = C(I,J) + A(I,L)*TEMP1 - B(I,L)*TEMP2
        110                     CONTINUE
                            C(J,J) = ZERO
                        END IF
        120             CONTINUE
        130         CONTINUE
            ELSE
                DO 180 J = 1,N
                    IF (BETA.EQ.ZERO) THEN
                        DO 140 I = J,N
                            C(I,J) = ZERO
        140                 CONTINUE
                    ELSE IF (BETA.NE.ONE) THEN
                        DO 150 I = J + 1,N
                            C(I,J) = BETA*C(I,J)
        150                 CONTINUE
                        C(J,J) = ZERO
                    ELSE
                        C(J,J) = ZERO
                    END IF
                    DO 170 L = 1,K
                        IF ((A(J,L).NE.ZERO) .OR. (B(J,L).NE.ZERO)) THEN
                            TEMP1 = ALPHA*B(J,L)
                            TEMP2 = ALPHA*A(J,L)
                            DO 160 I = J + 1,N
                                C(I,J) = C(I,J) + A(I,L)*TEMP1 - B(I,L)*TEMP2
        160                     CONTINUE
                            C(J,J) = ZERO
                        END IF
        170             CONTINUE
        180         CONTINUE
            END IF
        ELSE
            
            IF (UPPER) THEN
                DO 210 J = 1,N
                    DO 200 I = 1,J
                        TEMP1 = ZERO
                        TEMP2 = ZERO
                        DO 190 L = 1,K
                            TEMP1 = TEMP1 + A(L,I)*B(L,J)
                            TEMP2 = TEMP2 + B(L,I)*A(L,J)
        190                 CONTINUE
                        IF (I.EQ.J) THEN
                            C(J,J) = ZERO
                        ELSE
                            IF (BETA.EQ.ZERO) THEN
                                C(I,J) = ALPHA*TEMP1 - ALPHA*TEMP2
                            ELSE
                                C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 - ALPHA*TEMP2
                            END IF
                        END IF
        200             CONTINUE
        210         CONTINUE
            ELSE
                DO 240 J = 1,N
                    DO 230 I = J,N
                        TEMP1 = ZERO
                        TEMP2 = ZERO
                        DO 220 L = 1,K
                            TEMP1 = TEMP1 + A(L,I)*B(L,J)
                            TEMP2 = TEMP2 + B(L,I)*A(L,J)
        220                 CONTINUE
                        IF (I.EQ.J) THEN
                            C(J,J) = ZERO
                        ELSE
                            IF (BETA.EQ.ZERO) THEN
                                C(I,J) = ALPHA*TEMP1 - ALPHA*TEMP2
                            ELSE
                                C(I,J) = BETA*C(I,J) + ALPHA*TEMP1 - ALPHA*TEMP2
                            END IF
                        END IF
        230             CONTINUE
        240         CONTINUE
            END IF
        END IF
        
        RETURN


    end subroutine ZSKR2K

    subroutine ZSKR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)

        implicit none

        DOUBLE COMPLEX ALPHA
        INTEGER INCX,INCY,LDA,N
        CHARACTER UPLO
        
        DOUBLE COMPLEX A(LDA,*),X(*),Y(*)

        DOUBLE COMPLEX ZERO
        PARAMETER (ZERO= (0.0D+0,0.0D+0))
        
        DOUBLE COMPLEX TEMP1,TEMP2
        INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
        
        LOGICAL LSAME
        EXTERNAL LSAME
        
        EXTERNAL XERBLA
        
        INTRINSIC MAX
        
        INFO = 0
        IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
            INFO = 1
        ELSE IF (N.LT.0) THEN
            INFO = 2
        ELSE IF (INCX.EQ.0) THEN
            INFO = 5
        ELSE IF (INCY.EQ.0) THEN
            INFO = 7
        ELSE IF (LDA.LT.MAX(1,N)) THEN
            INFO = 9
        END IF
        IF (INFO.NE.0) THEN
            CALL XERBLA('ZSKR2 ',INFO)
            RETURN
        END IF
        
        IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
        
        IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
            IF (INCX.GT.0) THEN
                KX = 1
            ELSE
                KX = 1 - (N-1)*INCX
            END IF
            IF (INCY.GT.0) THEN
                KY = 1
            ELSE
                KY = 1 - (N-1)*INCY
            END IF
            JX = KX
            JY = KY
        END IF
        
        IF (LSAME(UPLO,'U')) THEN
            
            IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
                DO 20 J = 1,N
                    IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                        TEMP1 = ALPHA*Y(J)
                        TEMP2 = ALPHA*X(J)
                        DO 10 I = 1,J - 1
                            A(I,J) = A(I,J) + X(I)*TEMP1 - Y(I)*TEMP2
        10                 CONTINUE
                        A(J,J) = ZERO
                    ELSE
                        A(J,J) = ZERO
                    END IF
        20         CONTINUE
            ELSE
                DO 40 J = 1,N
                    IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                        TEMP1 = ALPHA*Y(JY)
                        TEMP2 = ALPHA*X(JX)
                        IX = KX
                        IY = KY
                        DO 30 I = 1,J - 1
                            A(I,J) = A(I,J) + X(IX)*TEMP1 - Y(IY)*TEMP2
                            IX = IX + INCX
                            IY = IY + INCY
        30                 CONTINUE
                        A(J,J) = ZERO
                    ELSE
                        A(J,J) = ZERO
                    END IF
                    JX = JX + INCX
                    JY = JY + INCY
        40         CONTINUE
            END IF
        ELSE
            
            IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
                DO 60 J = 1,N
                    IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                        TEMP1 = ALPHA*Y(J)
                        TEMP2 = ALPHA*X(J)
                        A(J,J) = ZERO
                        DO 50 I = J + 1,N
                            A(I,J) = A(I,J) + X(I)*TEMP1 - Y(I)*TEMP2
        50                 CONTINUE
                    ELSE
                        A(J,J) = ZERO
                    END IF
        60         CONTINUE
            ELSE
                DO 80 J = 1,N
                    IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                        TEMP1 = ALPHA*Y(JY)
                        TEMP2 = ALPHA*X(JX)
                        A(J,J) = ZERO
                        IX = JX
                        IY = JY
                        DO 70 I = J + 1,N
                            IX = IX + INCX
                            IY = IY + INCY
                            A(I,J) = A(I,J) + X(IX)*TEMP1 - Y(IY)*TEMP2
        70                 CONTINUE
                    ELSE
                        A(J,J) = ZERO
                    END IF
                    JX = JX + INCX
                    JY = JY + INCY
        80         CONTINUE
            END IF
        END IF
        
        RETURN        

    end subroutine ZSKR2

    subroutine ZLASKTRF(UPLO, MODE, N, NB, A, LDA, IPIV, W, LDW, INFO)

        implicit none

        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, LDW, N, NB, STEP
        
        INTEGER            IPIV( * )
        DOUBLE COMPLEX     A( LDA, * )
        DOUBLE COMPLEX     W( LDW, * )

        DOUBLE COMPLEX     ZERO, ONE
        PARAMETER          ( ZERO = (0.0D+0, 0.0D+0), ONE = (1.0D+0,0.0D+0) )
        
        INTEGER            K, KK, KP, NPANEL, WK
        DOUBLE PRECISION   COLMAX
        DOUBLE COMPLEX     T
        
        LOGICAL            LSAME
        INTEGER            IZAMAX
        EXTERNAL           LSAME, IZAMAX
        
        EXTERNAL           ZSCAL, ZSWAP, ZCOPY, ZGEMV, XERBLA
        
        INTRINSIC          ABS, MAX
        

        INFO = 0

        IF( LSAME( MODE, 'P' ) ) THEN
            STEP = 2
        ELSE
            STEP = 1
        END IF

        
        NPANEL = NB * STEP

        IF( LSAME( UPLO, 'U' ) ) THEN
            

            WK = 0
            DO 10 K=N, MAX(N-NPANEL+1, 2), -1
                
                
                KK = K-1

                IF( K .LT. N) THEN

                IF( WK .GT. 0 ) THEN
                    A( K, K ) = ZERO
                    CALL ZGEMV( 'N', K, WK, +ONE, A( 1, N-(WK-1)*STEP ), LDA*STEP, W( K, NB-WK+1 ), &
                    &LDW, ONE, A( 1, K ), 1 )
                    CALL ZGEMV( 'N', K, WK, -ONE, W( 1, NB-WK+1 ), LDW, A( K, N-(WK-1)*STEP ), &
                    &LDA*STEP, ONE, A( 1, K ), 1 )
                    A( K, K ) = ZERO
                END IF

                
                IF( MOD(K, STEP).EQ.1 .OR. STEP.EQ.1 ) THEN
                    WK = WK + 1
                    CALL ZCOPY(K, A(1, K), 1, W(1, NB-WK+1 ), 1)
                END IF
                END IF

                IF( MOD(K, STEP) .EQ. 0) THEN
                    
                KP = IZAMAX(K-1, A( 1, K ), 1)
                COLMAX = ABS( A( KP, K ) )

                IF( COLMAX.EQ.ZERO ) THEN
                    
                    IF( INFO.EQ.0 ) THEN
                        INFO = K - 1
                    END IF
                    KP = KK
                END IF
                

                IF( KP .NE. KK ) THEN
                    CALL ZSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ),1)
                    CALL ZSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )

                    CALL ZSWAP( N-K+1, A( KK, K), LDA, A( KP, K), LDA)

                    CALL ZSCAL(KK-KP, -ONE, A(KP, KK), 1)
                    CALL ZSCAL(KK-KP-1, -ONE, A(KP, KP+1), LDA)

                    
                    IF( WK .GT. 0 ) THEN
                        CALL ZSWAP( WK, W( KK, NB-WK+1 ), LDW, W( KP, NB-WK+1 ), LDW)
                    END IF
                END IF

                
                IF( COLMAX .NE. ZERO ) THEN
                    
                    CALL ZSCAL(K-2, ONE/A( K-1, K ), A(1, K), 1)
                END IF
                
                IPIV( K-1 ) = KP

                ELSE
                    
                IPIV( K-1 ) = K-1
                END IF
        10      CONTINUE

    
            IF( N-NPANEL+1 .GT. 2 ) THEN

                
                T = A( N-NPANEL, N-NPANEL+1 )
                A( N-NPANEL, N-NPANEL+1 ) = ZERO

                IF( WK .LT. NB) THEN
                    
                CALL ZCOPY(N-NPANEL, A(1, N-NPANEL), 1, W(1, 1), 1)

                W( N-NPANEL, 1 ) = ZERO
                CALL ZGEMV( 'N', N-NPANEL, WK, +ONE, A( 1, N-(WK-1)*STEP ), LDA*STEP, &
                &W( N-NPANEL, NB-WK+1 ), LDW, ONE, W( 1, 1 ), 1 )
                CALL ZGEMV( 'N', N-NPANEL, WK, -ONE, W( 1, NB-WK+1 ), LDW, A( N-NPANEL, &
                &N-(WK-1)*STEP ), LDA*STEP, ONE, W( 1, 1 ), 1 )
                W( N-NPANEL, 1 ) = ZERO

                WK = WK + 1
                END IF

                
                CALL ZSKR2K( UPLO, "N", N-NPANEL, NB, ONE, A(1, N-(WK-1)*STEP), LDA*STEP, &
                &W(1,1), LDW, ONE, A(1, 1), LDA)

                A( N-NPANEL, N-NPANEL+1 ) = T
            END IF

        ELSE

            WK = 0
            DO 30 K=1, MIN(NPANEL, N-1)

                
                KK = K+1

                IF( K .GT. 1) THEN

                IF( WK .GT. 0) THEN
                    A( K, K ) = ZERO
                    CALL ZGEMV( 'N', N-K+1, WK, +ONE, A( K, 1 ), LDA*STEP, W( K, 1 ), LDW, ONE, A( K, K ), 1 )
                    CALL ZGEMV( 'N', N-K+1, WK, -ONE, W( K, 1 ), LDW, A( K, 1 ), LDA*STEP, ONE, A( K, K ), 1 )
                    A( K, K ) = ZERO
                END IF

                
                IF( MOD(K, STEP) .EQ. 0 ) THEN
                    WK = WK + 1
                    CALL ZCOPY(N-K+1, A(K, K), 1, W(K, WK), 1)
                END IF
                END IF

                IF( MOD(K, STEP) .EQ. 1 .OR. STEP .EQ. 1) THEN
                    
                KP = K + IZAMAX(N-K, A( K+1, K ), 1)
                COLMAX = ABS( A( KP, K ) )

                IF( COLMAX.EQ.ZERO ) THEN
                    
                    IF( INFO.EQ.0 ) THEN
                        INFO = K
                    END IF
                    KP = KK
                END IF
                

                IF( KP .NE. KK ) THEN
                    IF( KP.LT.N ) THEN
                        CALL ZSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ),1 )
                    END IF

                    CALL ZSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )

                    CALL ZSWAP( K, A( KK, 1), LDA, A( KP, 1), LDA)

                    CALL ZSCAL( KP-KK, -ONE, A(KK+1, KK), 1 )
                    CALL ZSCAL( KP-KK-1, -ONE, A(KP, KK+1), LDA )

                    
                    CALL ZSWAP( WK, W( KK, 1 ), LDW, W( KP, 1 ), LDW)
                END IF

                
                IF( COLMAX .NE. ZERO .AND. K .LE. N-2) THEN
                    
                    CALL ZSCAL( N-K-1, ONE/A( K+1, K ), A(K+2, K), 1 )
                END IF

                
                IPIV( K+1 ) = KP

                ELSE
                    
                IPIV(K+1) = K+1
                END IF
        30      CONTINUE
    

            IF( NPANEL .LT. N-1) THEN
                
                T = A( NPANEL+1, NPANEL )
                A( NPANEL+1, NPANEL ) = ZERO

                IF( WK .LT. NB) THEN
                    
                CALL ZCOPY(N-NPANEL, A(NPANEL+1, NPANEL+1), 1, W(NPANEL+1, NB), 1)

                W( NPANEL+1, NB ) = ZERO
                CALL ZGEMV( 'N', N-NPANEL, NB-1, +ONE, A( NPANEL+1, 1 ), LDA*STEP, &
                &W( NPANEL+1, 1 ), LDW, ONE, W( NPANEL+1, NB ), 1 )
                CALL ZGEMV( 'N', N-NPANEL, NB-1, -ONE, W( NPANEL+1, 1 ), LDW, A( NPANEL+1, 1 ), &
                &LDA*STEP, ONE, W( NPANEL+1, NB ), 1 )
                W( NPANEL+1, NB ) = ZERO
                END IF

                CALL ZSKR2K( UPLO, "N", N-NPANEL, NB, ONE, A(NPANEL+1,1), LDA*STEP, W(NPANEL+1,1), &
                &LDW, ONE, A(NPANEL+1, NPANEL+1), LDA)

                A(NPANEL+1, NPANEL)=T
            END IF

        END IF

    end subroutine ZLASKTRF

    subroutine ZSKTD2(UPLO, MODE, N, A, LDA, E, TAU, INFO )

        implicit none

        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, N
        
        DOUBLE PRECISION   E( * )
        DOUBLE COMPLEX     A( LDA, * ), TAU( * )

        DOUBLE COMPLEX     ONE, ZERO
        PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ))
        
        LOGICAL            UPPER, NORMAL
        INTEGER            I, K, STEP
        DOUBLE COMPLEX     ALPHA, TAUI
        
        EXTERNAL           XERBLA, ZLARFG
        
        LOGICAL            LSAME
        EXTERNAL           LSAME
        
        INTRINSIC          MAX, MIN, CONJG
        
        INFO = 0
        UPPER = LSAME( UPLO, 'U' )
        NORMAL = LSAME( MODE, 'N' )

        IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.NORMAL .AND. .NOT.LSAME( MODE, 'P' ) ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( .NOT.NORMAL .AND. MOD(N,2).NE.0 ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        END IF

        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'ZSKTD2', -INFO )
            RETURN
        END IF
        
        IF( N.LE.0 ) RETURN
        
        IF( .NOT. NORMAL ) THEN
            STEP = 2
            
            DO 5 I = 2, N-2, 2
                TAU( I ) = ZERO
        5       CONTINUE
        ELSE
            STEP = 1
        END IF

        IF( UPPER ) THEN
            
            A( N, N ) = ZERO
            DO 10 I = N - 1, 1, -STEP
                
                ALPHA = A( I, I+1 )
                CALL ZLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
                E( I ) = ALPHA
                
                IF( TAUI.NE.ZERO ) THEN
                    
                A( I, I+1 ) = ONE
                
                DO 12 K = 1, I
                    A( K, I+1 ) = CONJG(A( K, I+1 ))
        12            CONTINUE
    
                CALL ZSKMV( UPLO, I, CONJG(TAUI), A, LDA, A( 1, I+1 ),1, ZERO, TAU, 1 )

                DO 15 K = 1, I
                    A( K, I+1 ) = CONJG(A( K, I+1 ))
        15            CONTINUE

    
                CALL ZSKR2( UPLO, I-STEP+1, ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA )
        
                ELSE
                A( I, I ) = ZERO
                END IF
                A( I, I+1 ) = E( I )
                TAU( I ) = TAUI
        10    CONTINUE
        ELSE
            
            A( 1, 1 ) = ZERO
            DO 20 I = 1, N - 1, STEP
                
                ALPHA = A( I+1, I )
                CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
                E( I ) = ALPHA
                
                IF( TAUI.NE.ZERO ) THEN
                    
                A( I+1, I ) = ONE
                
                DO 30 K = I+2, N
                    A( K, I ) = CONJG(A( K, I ))
        30          CONTINUE

                CALL ZSKMV( UPLO, N-I, CONJG(TAUI), A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 )
        
                DO 40 K = I+2, N
                    A( K, I ) = CONJG(A( K, I ))
        40          CONTINUE

    
                IF( I.LT.N-1 ) THEN
                    CALL ZSKR2( UPLO, N-I-STEP+1, ONE, A( I+STEP, I ), 1, TAU( I+STEP-1 ), 1, &
                    &A( I+STEP, I+STEP ), LDA )
                END IF
                
                ELSE
                A( I+1, I+1 ) = ZERO
                END IF
                A( I+1, I ) = E( I )
                TAU( I ) = TAUI
        20    CONTINUE
        END IF
        
        RETURN

    end subroutine ZSKTD2

    subroutine ZLASKTRD(UPLO, MODE, N, NB, A, LDA, E, TAU, W, LDW)

        implicit none

        CHARACTER          UPLO, MODE
        INTEGER            LDA, LDW, N, NB
        
        DOUBLE PRECISION   E( * )
        DOUBLE COMPLEX     A( LDA, * ), TAU( * ), W( LDW, * )

        DOUBLE COMPLEX     ZERO, ONE
        PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ))
        
        INTEGER            I, NW, NW2, STEP, NPANEL
        DOUBLE COMPLEX     ALPHA
        
        EXTERNAL           ZGEMV, ZLACGV, ZLARFG
        
        LOGICAL            LSAME
        EXTERNAL           LSAME
        
        INTRINSIC          MIN, CONJG
        
        IF( N.LE.0 ) RETURN
        

        IF( LSAME( MODE, 'P' ) ) THEN
            STEP = 2
        ELSE
            STEP = 1
        END IF
        NPANEL = NB * STEP

        IF( LSAME( UPLO, 'U' ) ) THEN
            
            NW=0

            DO 10 I = N, MAX(N - NPANEL + 1, 2), -1

                NW2 = NW - MOD(I,STEP)
                IF( NW2 .GT. 0 ) THEN
                    
                A( I, I ) = ZERO
                CALL ZGEMV( 'No transpose', I, NW2, +ONE, A( 1, N-(NW2-1)*STEP ), &
                &LDA*STEP, W( I, NB-NW2+1 ), LDW, ONE, A( 1, I ), 1 )
                CALL ZGEMV( 'No transpose', I, NW2, -ONE, W( 1, NB-NW2+1 ), LDW, &
                &A( I, N-(NW2-1)*STEP ), LDA*STEP, ONE, A( 1, I ), 1 )
                A( I, I ) = ZERO
                END IF

                
                IF( STEP.EQ.2 .AND. MOD( I, STEP ).EQ.1 ) THEN
                TAU( I-1 ) = ZERO
                GOTO 10
                END IF

                IF( I.GT.1 ) THEN
                    
                ALPHA = A( I-1, I )
                CALL ZLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
                E( I-1 ) = ALPHA
                A( I-1, I ) = ONE
                
                CALL ZLACGV( I-1, A( 1, I ), 1)

                CALL ZSKMV( 'Upper', I-1, CONJG(TAU(I-1)), A, LDA, A( 1, I ), 1, ZERO, W( 1, NB-NW ), 1 )
                IF( NW .GT. 0 ) THEN
                    CALL ZGEMV( 'Transpose', I-1, NW, ONE, W( 1, NB-NW+1 ), LDW, &
                    &A( 1, I ), 1, ZERO, W( I+1, NB-NW ), 1 )
                    CALL ZGEMV( 'No transpose', I-1, NW, CONJG(TAU(I-1)), A( 1, N-(NW-1)*STEP ), &
                    &LDA*STEP, W( I+1, NB-NW ), 1, ONE, W( 1, NB-NW ), 1 )
                    CALL ZGEMV( 'Transpose', I-1, NW, ONE, A( 1, N-(NW-1)*STEP ), LDA*STEP, &
                    &A( 1, I ), 1, ZERO, W( I+1, NB-NW ), 1 )
                    CALL ZGEMV( 'No transpose', I-1, NW, -CONJG(TAU(I-1)), W( 1, NB-NW+1 ), &
                    &LDW, W( I+1, NB-NW ), 1, ONE, W( 1, NB-NW ), 1 )
                END IF

                
                CALL ZLACGV( I-1, A( 1, I ), 1)

                
                NW = NW + 1

                END IF

        10    CONTINUE
        ELSE
            
            NW = 0

            DO 20 I = 1, MIN(NPANEL, N-1)
                
                NW2 = NW - MOD(I+1,STEP)
                IF( NW2 .GT. 0 ) THEN
                A( I, I ) = ZERO
                CALL ZGEMV( 'No transpose', N-I+1, NW2, +ONE, A( I, 1 ), &
                &LDA*STEP, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
                CALL ZGEMV( 'No transpose', N-I+1, NW2, -ONE, W( I, 1 ), &
                &LDW, A( I, 1 ), LDA*STEP, ONE, A( I, I ), 1 )
                A( I, I ) = ZERO
                END IF

                
                IF( STEP.EQ.2 .AND. MOD( I, STEP ).EQ.0 ) THEN
                TAU( I ) = ZERO
                GOTO 20
                END IF

                IF( I.LT.N ) THEN
                    
                ALPHA = A( I+1, I )
                CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
                E( I ) = ALPHA
                A( I+1, I ) = ONE
                
                CALL ZLACGV( N-I, A( I+1, I ), 1 )

                CALL ZSKMV( 'Lower', N-I, CONJG(TAU( I )), A( I+1, I+1 ), LDA, A( I+1, I ), &
                &1, ZERO, W( I+1, NW+1 ), 1 )
                IF( NW .GT. 0 ) THEN
                    CALL ZGEMV( 'Transpose', N-I, NW, ONE, W( I+1, 1 ), LDW, A( I+1, I ), &
                    &1, ZERO, W( 1, NW+1 ), 1 )
                    CALL ZGEMV( 'No transpose', N-I, NW, CONJG(TAU( I )), A( I+1, 1 ), &
                    &LDA*STEP, W( 1, NW+1 ), 1, ONE, W( I+1, NW+1 ), 1 )
                    CALL ZGEMV( 'Transpose', N-I, NW, ONE, A( I+1, 1 ), LDA*STEP, &
                    &A( I+1, I ), 1, ZERO, W( 1, NW+1 ), 1 )
                    CALL ZGEMV( 'No transpose', N-I, NW, -CONJG(TAU( I )), W( I+1, 1 ), &
                    &LDW, W( 1, NW+1 ), 1, ONE, W( I+1, NW+1 ), 1 )
                END IF

                
                CALL ZLACGV( N-I, A( I+1, I ), 1 )

                
                NW = NW + 1
                END IF
                
        20    CONTINUE
        END IF
        
        RETURN

    end subroutine ZLASKTRD

    subroutine ZSKTF2(UPLO, MODE, N, A, LDA, IPIV, INFO )

        implicit none

        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, N
        
        INTEGER            IPIV( * )
        DOUBLE COMPLEX     A( LDA, * )

        DOUBLE COMPLEX     ZERO, ONE
        PARAMETER          ( ZERO = (0.0D+0, 0.0D+0), ONE = (1.0D+0, 0.0D+0) )

        LOGICAL            UPPER, NORMAL
        INTEGER            K, KK, KP, STEP
        DOUBLE PRECISION   COLMAX
        
        LOGICAL            LSAME
        INTEGER            IZAMAX
        EXTERNAL           LSAME, IZAMAX
        
        EXTERNAL           ZSCAL, ZSWAP, XERBLA
        
        INTRINSIC          ABS, MAX, SQRT
        
        INFO = 0
        UPPER = LSAME( UPLO, 'U' )
        NORMAL = LSAME( MODE, 'N' )

        IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.NORMAL .AND. .NOT.LSAME( MODE, 'P' ) ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( .NOT.NORMAL .AND. MOD(N,2).EQ.1 ) THEN
            
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'ZSKTF2', -INFO )
            RETURN
        END IF

        
        IF( N .EQ. 0 ) RETURN

        IF( NORMAL ) THEN
            STEP = 1
        ELSE
            STEP = 2
        END IF

        IF( UPPER ) THEN
            
            IPIV( N ) = N

            DO 10 K=N, 2, -1
                
                IF( MOD(K, STEP) .EQ. 0) THEN
                    
                KP = IZAMAX(K-1, A( 1, K ), 1)
                COLMAX = ABS( A( KP, K ) )

                IF( COLMAX.EQ.ZERO ) THEN
                    
                    IF( INFO.EQ.0 ) THEN
                        INFO = K - 1
                    END IF
                    KP = K-1
                END IF

                
                KK = K-1

                IF( KP .NE. KK ) THEN
                    CALL ZSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1)
                    CALL ZSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )

                    CALL ZSWAP( N-K+1, A( KK, K), LDA, A( KP, K), LDA)

                    CALL ZSCAL(KK-KP, -ONE, A(KP, KK), 1)
                    CALL ZSCAL(KK-KP-1, -ONE, A(KP, KP+1), LDA)
                END IF

                
                IF( COLMAX .NE. ZERO ) THEN
                    CALL ZSKR2( UPLO, K-2, ONE/A( K-1,K ), A( 1, K ), 1, A( 1, K-1 ), 1, A( 1, 1 ), LDA )

        
                    CALL ZSCAL(K-2, ONE/A( K-1, K ), A(1, K), 1)
                END IF
                
                IPIV( K-1 ) = KP
                ELSE
                IPIV( K-1 ) = K-1
                END IF
        10      CONTINUE

        ELSE
            

            IPIV( 1 ) = 1

            DO 20 K=1, N-1
                
                IF( MOD(K, STEP).EQ.1 .OR. STEP.EQ.1 ) THEN
                    
                KP = K + IZAMAX(N-K, A( K+1, K ), 1)
                COLMAX = ABS( A( KP, K ) )

                IF( COLMAX.EQ.ZERO ) THEN
                    
                    IF( INFO.EQ.0 ) THEN
                        INFO = K
                    END IF
                    KP = K+1
                END IF

                
                KK = K+1

                IF( KP .NE. KK ) THEN
                    IF( KP.LT.N ) THEN
                        CALL ZSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ),1 )
                    END IF

                    CALL ZSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )

                    CALL ZSWAP( K, A( KK, 1), LDA, A( KP, 1), LDA)

                    CALL ZSCAL(KP-KK, -ONE, A(KK+1, KK), 1)
                    CALL ZSCAL(KP-KK-1, -ONE, A(KP, KK+1), LDA)
                END IF

                
                IF( COLMAX .NE. ZERO .AND. K+2 .LE. N) THEN
                    CALL ZSKR2( UPLO, N-K-1, ONE/A( K+1,K ), A( K+2, K ), 1, &
                    &A( K+2, K+1 ), 1, A( K+2, K+2 ), LDA )

        
                    CALL ZSCAL(N-K-1, ONE/A( K+1, K ), A(K+2, K), 1)
                END IF

                
                IPIV( K+1 ) = KP
                ELSE
                IPIV( K+1 ) = K+1
                END IF
        20      CONTINUE

        END IF

    end subroutine ZSKTF2

    subroutine ZSKTRF(UPLO, MODE, N, A, LDA, IPIV, WORK, LWORK, INFO)

        implicit none

        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, LWORK, N
        
        INTEGER            IPIV( * )
        DOUBLE COMPLEX     A( LDA, * ), WORK( * )

        LOGICAL            LQUERY, UPPER, NORMAL
        INTEGER            IINFO, J, K, K2, PIV, LWKOPT, NB, NBMIN, NPANEL
        
        LOGICAL            LSAME
        INTEGER            ILAENV
        EXTERNAL           LSAME, ILAENV
        
        EXTERNAL           ZSWAP, XERBLA
        
        INTRINSIC          MAX
        
        INFO = 0
        UPPER = LSAME( UPLO, 'U' )
        NORMAL = LSAME( MODE, 'N' )
        LQUERY = ( LWORK.EQ.-1 )
        IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.NORMAL .AND. .NOT.LSAME( MODE, 'P' ) ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( .NOT.NORMAL .AND. MOD(N,2).EQ.1 ) THEN
            
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
            INFO = -8
        END IF
        

        IF( INFO.EQ.0 ) THEN
            
            NB = ILAENV( 1, 'ZHETRF', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
            WORK( 1 ) = LWKOPT
        END IF
        
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'ZSKTRF', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        

        NBMIN=NB
        IF( NB.GT.1 .AND. NB.LT.N ) THEN
            IF( LWORK .LT. N*NB ) THEN
                NB = MAX( LWORK / N, 1 )
                NBMIN = MAX( 2, ILAENV( 2, 'ZHETRF', UPLO, N, -1, -1, -1 ) )
            END IF
        ELSE
            NB = N
        END IF

        IF( NB.LT.NBMIN ) NB = N
        
        IF( N .EQ. 0 ) RETURN

        IF( LSAME( MODE, 'N' ) ) THEN
            NPANEL = NB
        ELSE
            NPANEL = MIN(NB*2, N)
        END IF

        IF( UPPER ) THEN
            

            IPIV( N ) = N
            
            DO 10 K = N, MAX(NPANEL, 1), -NPANEL
                
                IF( K.GE.NPANEL*2 ) THEN
                    
                CALL ZLASKTRF( UPLO, MODE, K, NB, A, LDA, IPIV, WORK, N, IINFO )

                K2 = K-NPANEL
                ELSE
                    
                PIV = IPIV( K )

                CALL ZSKTF2( UPLO, MODE, K, A, LDA, IPIV, IINFO )

                IPIV( K ) = PIV

                K2 = 1
                END IF
                

                IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO

        
                IF( K .LT. N ) THEN
                DO 20 J=K-1, K2, -1
                    CALL ZSWAP( N-K, A( J, K+1 ), LDA, A( IPIV( J ), K+1 ), LDA )
        20            CONTINUE
                END IF

        10   CONTINUE
    
        ELSE
            

            IPIV( 1 ) = 1
            
            DO 30 K = 1, MIN(N-NPANEL+1, N-1), NPANEL
                
                IF( K.LE.N-NPANEL*2+1 ) THEN
                    
                CALL ZLASKTRF( UPLO, MODE, N-K+1, NB, A( K, K ), LDA, IPIV( K ), WORK, N, IINFO )

                K2 = K + NPANEL
                ELSE
                    
                PIV = IPIV( K )

                CALL ZSKTF2( UPLO, MODE, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )

                IPIV( K ) = PIV

                K2 = N
                END IF
                
                IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + K - 1
        
                DO 40 J = K+1, K2
                IPIV( J ) = IPIV( J ) + K - 1
        40         CONTINUE

    
                IF( K .GT. 1 ) THEN
                DO 50 J=K+1, K2
                    CALL ZSWAP( K-1, A( J, 1 ), LDA, A( IPIV( J ), 1 ), LDA )
        50            CONTINUE
                END IF


        30   CONTINUE
    
        END IF
        
        WORK( 1 ) = LWKOPT
        RETURN

    end subroutine ZSKTRF

    subroutine ZSKTRD(UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK, INFO )

        implicit none

        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, LWORK, N
        
        DOUBLE PRECISION   E( * )
        DOUBLE COMPLEX     A( LDA, * ), TAU( * ), WORK( * )

        DOUBLE COMPLEX     CONE
        PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
        
        LOGICAL            LQUERY, UPPER, NORMAL
        INTEGER            I, IINFO, IWS, J, LDWORK, LWKOPT, NB, NBMIN, NX, STEP, NPANEL, NXPANEL
        
        EXTERNAL           XERBLA
        
        INTRINSIC          MAX
        
        LOGICAL            LSAME
        INTEGER            ILAENV
        EXTERNAL           LSAME, ILAENV
        
        INFO = 0
        UPPER = LSAME( UPLO, 'U' )
        NORMAL = LSAME( MODE, 'N' )

        LQUERY = ( LWORK.EQ.-1 )
        IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.NORMAL .AND. .NOT.LSAME( MODE, 'P' ) ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( .NOT.NORMAL .AND. MOD(N,2).NE.0 ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
            INFO = -9
        END IF
        
        IF( INFO.EQ.0 ) THEN
            
            NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
            WORK( 1 ) = LWKOPT
        END IF
        
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'ZSKTRD', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        
        IF( N.EQ.0 ) THEN
            WORK( 1 ) = 1
            RETURN
        END IF
        
        NX = N
        IWS = 1
        IF( NB.GT.1 .AND. NB.LT.N ) THEN
            
            NX = MAX( NB, ILAENV( 3, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
            IF( NX.LT.N ) THEN
                
                LDWORK = N
                IWS = LDWORK*NB
                IF( LWORK.LT.IWS ) THEN
                    
                NB = MAX( LWORK / LDWORK, 1 )
                NBMIN = ILAENV( 2, 'ZHETRD', UPLO, N, -1, -1, -1 )
                IF( NB.LT.NBMIN .OR. NB.LE.1 ) NX = N
                END IF
            ELSE
                NX = N
            END IF
        ELSE
            NB = 1
        END IF
        

        IF( .NOT.NORMAL ) THEN
            STEP = 2
        ELSE
            STEP = 1
        END IF

        NPANEL = NB * STEP
        NXPANEL = NX * STEP

        IF( UPPER ) THEN
            
            DO 20 I = N, NXPANEL + NPANEL, -NPANEL
                
                CALL ZLASKTRD( UPLO, MODE, I, NB, A, LDA, E, TAU, WORK, LDWORK )
        
                CALL ZSKR2K( UPLO, 'No transpose', I-NPANEL, NB, CONE, A( 1, I-NPANEL+STEP ), &
                &LDA*STEP, WORK, LDWORK, CONE, A, LDA )
        
                DO 10 J = I-NPANEL+1+STEP-1, I, STEP
                A( J-1, J ) = E( J-1 )
        10       CONTINUE
        20    CONTINUE
    
            CALL ZSKTD2( UPLO, MODE, I, A, LDA, E, TAU, IINFO )
        ELSE
            
            DO 40 I = 1, N - NXPANEL, NPANEL
                
                CALL ZLASKTRD( UPLO, MODE, N-I+1, NB, A( I, I ), LDA, E( I ), TAU( I ), WORK, LDWORK )
        
                CALL ZSKR2K( UPLO, 'No transpose', N-I-NPANEL+1, NB, CONE, A( I+NPANEL, I ), &
                &LDA*STEP, WORK( NPANEL+1 ), LDWORK, CONE, A( I+NPANEL, I+NPANEL ), LDA )
        
                DO 30 J = I, I + NPANEL - 1, STEP
                A( J+1, J ) = E( J )
        30       CONTINUE
        40    CONTINUE
    
            CALL ZSKTD2( UPLO, MODE, N-I+1, A( I, I ), LDA, E( I ), TAU( I ), IINFO )
        END IF
        
        WORK( 1 ) = LWKOPT
        RETURN

    end subroutine ZSKTRD

    subroutine ZSKPF10( UPLO, MTHD, N, A, LDA, PFAFF, IWORK, WORK, LWORK, RWORK, INFO)

        implicit none

        CHARACTER          UPLO, MTHD
        INTEGER            INFO, LDA, LWORK, N
        
        INTEGER            IWORK( * )
        DOUBLE PRECISION   RWORK( * )
        DOUBLE COMPLEX     PFAFF( 2 )
        DOUBLE COMPLEX     A( LDA, * ), WORK( * )

        DOUBLE COMPLEX       ONE, ZERO
        PARAMETER          ( ONE = (1.0D+0, 0.0D+0) )
        PARAMETER          ( ZERO = (0.0D+0, 0.0D+0) )

        DOUBLE PRECISION   RONE
        PARAMETER          ( RONE = 1.0D+0 )

        INTEGER            I,K
        DOUBLE PRECISION   TEMP
        
        LOGICAL            LQUERY, UPPER, LTL

        EXTERNAL           XERBLA
        
        LOGICAL            LSAME
        EXTERNAL           LSAME
        INTRINSIC          CONJG, CMPLX

        INFO = 0
        UPPER = LSAME( UPLO, 'U' )
        LTL = LSAME( MTHD, 'P' )
        LQUERY = ( LWORK.EQ.-1 )

        IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
            INFO = -1
        ELSE IF( .NOT.LTL .AND. .NOT.LSAME( MTHD, 'H' ) ) THEN
            INFO = -2
        ELSE IF( N.LT.0 ) THEN
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
            INFO = -9
        ELSE IF( MOD(N,2).NE.1 .AND. .NOT.LTL .AND. LWORK.LT.N .AND. .NOT.LQUERY ) THEN
            INFO = -9
        END IF

        IF( INFO.EQ.0 .AND. LQUERY) THEN
            IF( MOD(N,2).EQ.1 ) THEN
                WORK(1) = 1
            ELSE IF( LTL ) THEN
                
                CALL ZSKTRF( UPLO, "P", N, A, LDA, IWORK, WORK, LWORK, INFO )
            ELSE
                
                CALL ZSKTRD( UPLO, "P", N, A, LDA, RWORK, WORK, WORK, LWORK, INFO)
                WORK(1) = WORK(1) + N - 1
            END IF
        END IF
        
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'ZSKPF10', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF

        PFAFF( 1 ) = ONE
        PFAFF( 2 ) = ZERO

        
        IF( N.EQ.0 ) THEN
            RETURN
        ELSE IF( MOD(N,2).EQ.1 ) THEN
            PFAFF( 1 ) = ZERO
            RETURN
        END IF

        IF( LTL ) THEN
            
            CALL ZSKTRF( UPLO, "P", N, A, LDA, IWORK, WORK, LWORK, INFO )

            
            IF( INFO .GT. 0 ) THEN
                PFAFF( 1 ) = ZERO
                PFAFF( 2 ) = ZERO
                INFO = 0
            ELSE
                IF( UPPER ) THEN

                DO 10 I = 1, N-1, 2
                    CALL ZMUL10( PFAFF, A( I, I+1 ) )

                    
                    IF( IWORK( I ) .NE. I ) PFAFF( 1 ) = -PFAFF( 1 )
        10            CONTINUE

                ELSE

                DO 20 I = 1, N-1, 2
                    CALL ZMUL10( PFAFF, -A( I+1, I ) )

                    
                    IF( IWORK( I+1 ) .NE. I+1 ) PFAFF( 1 ) = -PFAFF( 1 )
        20            CONTINUE

                END IF
            END IF
        ELSE

            
            CALL ZSKTRD(UPLO, "P", N, A, LDA, RWORK, WORK, WORK( N ), LWORK-N+1, INFO)

            IF( UPPER ) THEN
                
                DO 30 I = 1, N-1, 2
                CALL ZMUL10( PFAFF, CMPLX(RWORK( I ),KIND=KIND(RWORK)) )

                
                TEMP = RONE
                DO 40 K=1, I-1
                    TEMP = TEMP + CONJG(A(K,I+1))*A(K,I+1)
        40            CONTINUE

                PFAFF( 1 ) = PFAFF( 1 ) * (ONE - WORK( I ) * TEMP)
        30         CONTINUE

            ELSE

                
                DO 50 I = 1, N-1, 2
                CALL ZMUL10( PFAFF,-CMPLX(RWORK( I ),KIND=KIND(RWORK)))

                
                TEMP = RONE
                DO 60 K=I+2, N
                    TEMP = TEMP + CONJG(A(K,I))*A(K,I)
        60            CONTINUE

                PFAFF( 1 ) = PFAFF( 1 ) * (ONE - WORK( I ) * TEMP)
        50         CONTINUE

            END IF

            
            WORK( 1 ) = WORK( N ) + N-1
        END IF

        RETURN

    end subroutine ZSKPF10

end program PFAFFIAN