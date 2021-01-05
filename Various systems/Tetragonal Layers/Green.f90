program GREEN
    implicit none
    integer :: NUMKX, NUMKY, NUMKZ, NUMK, NUMT, io, i, j, IRLATT, IRLATTMAX, kcounter, LWORK, INFO, NCELLS, &
    &NUMIMP, NUMCHEMTYPES, JATOM, rcheck, MAXNEIGHB, NPNT, NPNT1, NPNT2, NPNT3, &
    &NPOL, NUMEDOS, ZEROPOL, IEMXD, IEMXDZ, FTIMO
    integer, allocatable, dimension(:) :: CHEMTYPE, NEIGHBNUM
    integer, allocatable, dimension(:,:) :: IMPPTSVAR, JTYPE
    real*8 :: ALAT, a_1(3), a_2(3), a_3(3), RMAX, R0, KPOINT(3), RPOINT(3), TTPRIME(3), ENERGY, EMIN,&
    &chempot, T, TBROD, PI, KB, b_1(3), b_2(3), b_3(3), etados, lambda, TOL, tempvalre, tempvalim, EBOT, EMU
    real*8, allocatable, dimension(:) :: W, RWORK, E0, ULCN, nu, nuzero, magnet, VSUPCOND
    real*8, allocatable, dimension(:,:) :: KPTS, TPTS, EIGENVALUES, RLATT, BETA, LHOPS, PREFACTORS, HOPPVALS
    real*8, allocatable, dimension(:,:,:) :: RCONNECT
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), tempparam
    complex*16, allocatable, dimension(:) :: WORK, DELTA, DF, EZ, DFDOS, EZDOS
    complex*16, allocatable, dimension(:,:) :: HAMILTONIAN, HAMILTONIANPREP
    complex*16, allocatable, dimension(:,:,:) :: EXPONS, EIGENVECTORS
    character(len = 1) :: yesorno, vecorints, dosorno, magorno
    character(len = 4) :: fullorp

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
    read(1,*) ! chempot 
    read(1,*) ! inicharge
    read(1,*) ! inidelta
    read(1,*) ! mixfactorN
    read(1,*) ! mixfactorD
    read(1,*) NUMEDOS
    read(1,*) etados
    read(1,*) ! metalepsilon
    read(1,*) ! scepsilon
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

    do i = 1, NUMT
        nu(i) = 0.d0
        DELTA(i) = (0.d0,0.d0)

        FTIMO = 4*(i-1)

        do kcounter = 1, NUMK
            do j = 1, 4*NUMT

                ENERGY = EIGENVALUES(j,kcounter)

                nu(i) = nu(i) + FERMI(ENERGY,T,KB)*abs(EIGENVECTORS(i,j,kcounter))**2 +&
                & FERMI(ENERGY,T,KB)*abs(EIGENVECTORS(i+NUMT,j,kcounter))**2

                DELTA(i) = DELTA(i) -FERMI(ENERGY,T,KB)*VSUPCOND(i)*&
                &EIGENVECTORS(i,j,kcounter)*CONJG(EIGENVECTORS(i+3*NUMT,j,kcounter))
            end do
        end do

        nu(i) = nu(i)/NUMK
        DELTA(i) = DELTA(i)/NUMK
    end do

    print *, 'Diagonalization finished and eigenvectors/eigenvalues stored.'

    print *, 'We proceed with the setup of the energies to be used for the calculation of G.'

    open(1, file = 'greenconfig.dat', action = 'read')
    read(1,*) NPOL
    read(1,*) dosorno
    read(1,*) NPNT1
    read(1,*) NPNT2
    read(1,*) NPNT3
    read(1,*) EBOT
    read(1,*) vecorints
    read(1,*) fullorp
    read(1,*) magorno
    close(1)

    call IDENTIFIER(vecorints,b_1,b_2,b_3,PI,a_1,a_2,a_3,TPTS,NUMT,NUMIMP,IMPPTSVAR)

    if (EBOT == -100.0) then
        EBOT = MINVAL(EIGENVALUES) - 0.5 ! The - 0.5 factor is inserted to avoid broadening cutoffs
    endif
    
    IEMXD = (NPNT1+NPNT2+NPNT3+NPOL)*2 ! In order to setup the EZ and DF arrays

    allocate(DF(IEMXD))
    allocate(EZ(IEMXD))

    ! If NPOL == 0 we check if the user wants to perform a DoS calculation apart from the other calculations.
    if (NPOL /= 0 .and. dosorno == 'y') then
        ZEROPOL = 0
        TBROD = etados/(KB*PI)

        IEMXDZ = NUMEDOS*2

        allocate(DFDOS(IEMXDZ))
        allocate(EZDOS(IEMXDZ))

        print *, 'It appears that you want to perform a DoS calculation apart from the other calculations.'
        print *, 'Should the minimum energy for the DoS be the same as EBOT ? y/n'
        171 read *, yesorno
        if (yesorno == 'y') then
            EMIN = EBOT
        else if (yesorno == 'n') then
            print *, 'Please enter the minimum energy to be used for the DoS.'
            read *, EMIN
        else
            print *, 'Invalid input. Please enter y or n.'
            goto 171
        endif
        print *, 'Should the maximum eigenvalue of the spectrum be used as the maximum energy? y/n'
        172 read *, yesorno
        if (yesorno == 'y') then
            EMU = MAXVAL(EIGENVALUES) + 0.5 ! The + 0.5 factor is inserted to avoid broadening cutoffs
        else if (yesorno == 'n') then
            print *, 'Please enter the maximum energy to be used for the DoS.'
            read *, EMU
        else
            print *, 'Invalid input. Please enter y or n.'
            goto 172
        endif

        call EMESHT(EZDOS,DFDOS,NPNT,EMIN,EMU,TBROD,ZEROPOL,NPNT1,NUMEDOS,NPNT3,PI,KB,IEMXDZ)

        NUMEDOS = NPNT ! To be written later on the file impconfig.dat

        open(1, file = 'energiesfordos.dat', action = 'write')
            do i = 1, NPNT
                tempparam = EZDOS(i)
                write (1, '(2F17.8)') REAL(tempparam), AIMAG(tempparam)
            end do
        close(1)

        print *, 'Energies and weights printed. Initiating Green functions DoS calculations.'

        if (fullorp == 'full') then
            call FULLGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NPNT,ZEROPOL)
        else
            call PARTIALGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NPNT,ZEROPOL)
        endif

        print *, 'DoS calculated.'
        print *, 'We proceed with the setup of the energies to be used for the contours.'

    endif
    
    EMU = 0.0 ! Because we are studying superconductivity
    call EMESHT(EZ,DF,NPNT,EBOT,EMU,T,NPOL,NPNT1,NPNT2,NPNT3,PI,KB,IEMXD)

    open(1, file = 'energies.dat', action = 'write')
        do i = 1, NPNT
            tempparam = EZ(i)
            write (1, '(2F17.8)') REAL(tempparam), AIMAG(tempparam)
        end do
    close(1)

    open(1, file = 'fermiweights.dat', action = 'write')
        do i = 1, NPNT
            tempparam = DF(i)
            write (1, '(2F17.8)') REAL(tempparam), AIMAG(tempparam)
        end do
    close(1)

    print *, 'Energies and weights printed.'

    print *, 'Initiating Green functions calculations.'
    if (fullorp == 'full') then
        call FULLGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NPNT,NPOL)
    else
        call PARTIALGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NPNT,NPOL)
    endif

    ! Prepares two config files to be used by the impurity program
    open (1, file = 'impconfig.dat', action = 'write')
        write (1,*) NUMIMP, '! Number of impurities.'
        write (1,*) NPNT, '! Number of energy values.'
        write (1,*) NUMEDOS, '! Number of energy values for the DoS.'
        write (1,*) '          n ! y to print G_imp elements on a file and n to ignore.'
        write (1,*) '   0.00001  ! epsilon for impurity self-consistency.'
        do i = 1, NUMIMP
            j = IMPPTSVAR(4,i)
            write (1,'(3F15.9, A, I7)') nu(j), REAL(DELTA(j)), AIMAG(DELTA(j)), ' ! n and D for atom No. ', i
        end do
    close(1)

    if (magorno == 'n') then
        open (1, file = 'impatoms.dat', action = 'write')
        do i = 1, NUMIMP
            j = IMPPTSVAR(4,i)
            write (1,'(7F15.9, A, I7)') E0(j), BETA(1,j), BETA(2,j), BETA(3,j), VSUPCOND(j), ULCN(j), nuzero(j),&
            & '    ! Host No. ', i
            write (1,'(7F15.9, A, I7)') E0(j), BETA(1,j), BETA(2,j), BETA(3,j), VSUPCOND(j), ULCN(j), nuzero(j),&
            & '    ! Impurity No. ', i
        end do
        write (1,'(A)') '---------------------------------------------------------------------------------------&
        &-------------------'
        write (1,'(A)') 'Format: E_0,           B_x,           B_y,           B_z,            V,&
        &             U,            n_0'
        close(1)
    else
        call MAGCHAIN(NUMT,NUMIMP,PI,E0,BETA,VSUPCOND,ULCN,nuzero,IMPPTSVAR)
    endif

    contains

    subroutine MAGCHAIN(NUMT,NUMIMP,PI,E0,BETA,VSUPCOND,ULCN,nuzero,IMPPTSVAR)
        implicit none

        integer :: NUMT, NUMIMP, num, denom, i, IMPPTSVAR(4,NUMIMP)
        real*8 :: ROTAXIS(3), theta, E0(NUMT), BETA(3,NUMT), NEWBETA(3,NUMIMP), PI, magB, initheta, TOL, VSUPCONDIMP(NUMT),&
        &ROTMATRIX(3,3), INBETA(3), OUTBETA(3), one, zero, sine, cosine, ULCN(NUMT), nuzero(NUMT), VSUPCOND(NUMT)

        TOL = 0.000001 ! The tolerance for the "if" checks.

        print *, 'Preparing the magnetic chain inputs.'

        open(1, file = 'greenconfig.dat', action = 'read')
            do i = 1, 9
                read(1,*)
            end do ! Skips the file's first 9 lines
            read(1,*) ROTAXIS
            read(1,*) num, denom
        close(1)

        theta = (num*PI)/denom

        if (abs(ROTAXIS(1) - 1.0) < TOL .and. abs(ROTAXIS(2)) < TOL .and. abs(ROTAXIS(3)) < TOL) then
            ! This is a rotation equivalent to a rotation by the x-axis (i.e. in the y-z plane)
            print *, 'Please enter the starting angle with respect to the z-axis.'
            print *, 'Format a*PI/b, with a, b: integers. Please enter a.'
            read *, num
            print *, 'Please enter b.'
            read *, denom

            initheta = (num*PI)/denom

            print *, 'Please enter the magnitude of B.'
            read *, magB

            do i = 1, NUMIMP
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
            ! This is a rotation equivalent to a rotation by the y-axis (i.e. in the x-z plane)
            print *, 'Please enter the starting angle with respect to the x-axis.'
            print *, 'Format a*PI/b, with a, b: integers. Please enter a.'
            read *, num
            print *, 'Please enter b.'
            read *, denom

            initheta = (num*PI)/denom

            print *, 'Please enter the magnitude of B.'
            read *, magB

            do i = 1, NUMIMP
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
            ! This is a rotation equivalent to a rotation by the z-axis (i.e. in the x-y plane)
            print *, 'Please enter the starting angle with respect to the x-axis.'
            print *, 'Format a*PI/b, with a, b: integers. Please enter a.'
            read *, num
            print *, 'Please enter b.'
            read *, denom

            initheta = (num*PI)/denom

            print *, 'Please enter the magnitude of B.'
            read *, magB

            do i = 1, NUMIMP
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
            print *, 'Please enter the value of B for the first atom of the chain. Format Bx,By,Bz.'
            read *, NEWBETA(:,1)

            cosine = COS(theta)
            sine = SIN(theta)

            ! This sets up the rotation matrix.
            ROTMATRIX(1,1) = cosine + ROTAXIS(1)**2*(1.0-cosine)
            ROTMATRIX(1,2) = ROTAXIS(1)*ROTAXIS(2)*(1.0-cosine) - ROTAXIS(3)*sine
            ROTMATRIX(1,3) = ROTAXIS(1)*ROTAXIS(3)*(1.0-cosine) + ROTAXIS(2)*sine

            ROTMATRIX(2,1) = ROTAXIS(1)*ROTAXIS(2)*(1.0-cosine) + ROTAXIS(3)*sine
            ROTMATRIX(2,2) = cosine + ROTAXIS(2)**2*(1.0-cosine)
            ROTMATRIX(2,3) = ROTAXIS(2)*ROTAXIS(3)*(1.0-cosine) - ROTAXIS(1)*sine

            ROTMATRIX(3,1) = ROTAXIS(1)*ROTAXIS(3)*(1.0-cosine) - ROTAXIS(2)*sine
            ROTMATRIX(3,2) = ROTAXIS(3)*ROTAXIS(2)*(1.0-cosine) + ROTAXIS(1)*sine
            ROTMATRIX(3,3) = cosine + ROTAXIS(3)**2*(1.0-cosine)

            ! Each Beta is multiplied by the rotation matrix.
            if (NUMIMP > 1) then
                one = 1.0
                zero = 0.0

                do i = 2, NUMIMP
                    INBETA(:) = NEWBETA(:,i-1)
                    call DGEMV ('N', 3, 3, one, ROTMATRIX, 3, INBETA, 1, zero, OUTBETA, 1)
                    NEWBETA(:,i) = OUTBETA(:)
                enddo
            endif

        endif

        open (1, file = 'impatoms.dat', action = 'write')
        do i = 1, NUMIMP
            j = IMPPTSVAR(4,i)
            VSUPCONDIMP(j) = 0.0
            write (1,'(7F15.9, A, I7)') E0(j), BETA(1,j), BETA(2,j), BETA(3,j), VSUPCOND(j), ULCN(j), nuzero(j),&
            &'    ! Host No. ', i
            write (1,'(7F15.9, A, I7)') E0(j), NEWBETA(1,i), NEWBETA(2,i), NEWBETA(3,i), VSUPCONDIMP(j), ULCN(j), nuzero(j),&
            &'    ! Impurity No. ', i
        end do
        write (1,'(A)') '---------------------------------------------------------------------------------------&
        &-------------------'
        write (1,'(A)') 'Format: E_0,           B_x,           B_y,           B_z,            V,&
        &             U,            n_0'
        close(1)

    endsubroutine MAGCHAIN

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

    subroutine FULLGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,PI,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NPNT,NPOL)
        implicit none

        integer :: NUMT, NUMK, IE, i, j, k, n, NUMIMP, IMPPTSVAR(4,NUMIMP), &
        &l, m, a, aprime, FTIMO, FTJMO, checker, NPNT, NPOL
        real*8 :: PI, EIGENVALUES(4*NUMT,NUMK), greendensityperatom(1+NUMT,NPNT), &
        &a_1(3), a_2(3), a_3(3), RPOINT(3), RPRIMEPOINT(3), FOURIERVEC(3), TPTS(3,NUMT), KPTS(3,NUMK), KPOINT(3),&
        &greendensity(2,NPNT), densityperimpurity(1+NUMIMP,NPNT), tempval1, tempval2
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT,NUMK), GMATRIX(4*NUMT,4*NUMT), EZ, ENFRAC, GFK(NUMK,4*NUMT,4*NUMT),&
        &GREENR(4*NUMIMP,4*NUMIMP), FOUREXPONS(NUMK,NUMIMP**2), energies(NPNT), WEIGHTS(NPNT)

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

        if (NPOL == 0) then

            open(1, file = 'energiesfordos.dat', action = 'read')
            do i = 1, NPNT
                read(1,*) tempval1, tempval2
                energies(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

        else

            open(1, file = 'energies.dat', action = 'read')
            do i = 1, NPNT
                read(1,*) tempval1, tempval2
                energies(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

            open(1, file = 'fermiweights.dat', action = 'read')
            do i = 1, NPNT
                read(1,*) tempval1, tempval2
                WEIGHTS(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

        endif

        ! This part writes all the Fouriered Green function elements per energy at this text file
        if (NPOL == 0) then
            open(1, file = 'greenhostfordos.txt', action = 'write')
        else
            open(1, file = 'greenhost.txt', action = 'write')
        endif

        do IE = 1, NPNT ! Begins a loop over the energies, in order to find G(E) for each E

            EZ = energies(IE)

            if (NPOL == 0) then
                greendensityperatom(1,IE) = REAL(EZ)
                greendensity(1,IE) = REAL(EZ)
                greendensity(2,IE) = 0.0
                do i = 1, NUMT
                    greendensityperatom(1+i,IE) = 0.0
                end do
                
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

                if (NPOL == 0) then
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

            if (NPOL == 0) then
                do i = 1, NUMIMP ! Calculation of density for each impurity atom
                    FTIMO = 4*(i-1)
                    do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                        densityperimpurity(1+i,IE) = densityperimpurity(1+i,IE) -&
                        &(1.0/PI)*AIMAG(GREENR(j + FTIMO, j + FTIMO))
                    end do
                end do
            endif

            ! Writes the Green impurity elements on the corresponding .txt file
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1, '(F17.8, F17.8)') GREENR(i,j)
                end do
            end do

            if (NPOL == 0) then
                do i = 1, NUMT ! Calculation of full density
                    greendensityperatom(i+1,IE) = greendensityperatom(i+1,IE)/NUMK ! Normalization
                    greendensity(2,IE) = greendensity(2,IE) + (1.0/NUMT)*greendensityperatom(1+i,IE)
                end do
            endif

        end do ! ends energies sum
        close(1) ! Closes the .txt file

        if (NPOL == 0) then
            open(1, file = 'greendensityperatom.txt', action = 'write')
            do j = 1, NPNT
                do i = 1, NUMT+1
                    if (i == NUMT+1) then
                        write (1,'(F17.8)',advance='no') greendensityperatom(i,j)
                    else
                        write (1,'(F17.8, A)',advance='no') greendensityperatom(i,j), ','
                    endif
                end do
                write (1,*)
            end do
            close(1)

            ! The density per atom at the future impurity sites.
            open(1, file = 'hostdensities.txt', action = 'write')
            do j = 1, NPNT
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

            open(1, file = 'greendensity.txt', action = 'write')
            do j = 1, NPNT ! Energies = Intervals + 1
                write (1,'(F17.8, A, F17.8)') greendensity(1,j), ',', greendensity(2,j)
            end do
            close(1)
        endif

    end subroutine FULLGREEN

    subroutine PARTIALGREEN(EIGENVALUES,EIGENVECTORS,NUMT,NUMK,TPTS,a_1,a_2,a_3,KPTS,NUMIMP,IMPPTSVAR,NPNT,NPOL)
        implicit none

        integer :: NUMT, NUMK, IE, i, j, k, n, NUMIMP, IMPPTSVAR(4,NUMIMP), min_val, max_val, l, m, a, aprime,&
        &FTIMO, FTJMO, checker, IMPATOMTYPE(NUMIMP), uniquecounter, NUMATOMS, IATOM, JATOM, NPNT, NPOL
        integer, allocatable, dimension(:) :: IMPATOMVALS, UNIQUEIMPATOMS
        real*8 :: EIGENVALUES(4*NUMT,NUMK), densityperimpurity(1+NUMIMP,NPNT), tempval1, tempval2,&
        &a_1(3), a_2(3), a_3(3), RPOINT(3), RPRIMEPOINT(3), FOURIERVEC(3), TPTS(3,NUMT), KPTS(3,NUMK), KPOINT(3)
        complex*16 :: EIGENVECTORS(4*NUMT,4*NUMT,NUMK), GMATRIX(4*NUMT,4*NUMT), EZ, ENFRAC, GFK(NUMK,4*NUMT,4*NUMT),&
        &GREENR(4*NUMIMP,4*NUMIMP), FOUREXPONS(NUMK,NUMIMP**2), energies(NPNT), WEIGHTS(NPNT)

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

        if (NPOL == 0) then

            open(1, file = 'energiesfordos.dat', action = 'read')
            do i = 1, NPNT
                read(1,*) tempval1, tempval2
                energies(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

        else

            open(1, file = 'energies.dat', action = 'read')
            do i = 1, NPNT
                read(1,*) tempval1, tempval2
                energies(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

            open(1, file = 'fermiweights.dat', action = 'read')
            do i = 1, NPNT
                read(1,*) tempval1, tempval2
                WEIGHTS(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

        endif

        ! This part writes all the Fouriered Green function elements per energy at this text file
        if (NPOL == 0) then
            open(1, file = 'greenhostfordos.txt', action = 'write')
        else
            open(1, file = 'greenhost.txt', action = 'write')
        endif

        do IE = 1, NPNT ! Begins a loop over the energies, in order to find G(E) for each E

            EZ = energies(IE)

            if (NPOL == 0) then
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

            if (NPOL == 0) then
                do i = 1, NUMIMP ! Calculation of density for each impurity atom
                    FTIMO = 4*(i-1)
                    do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                        densityperimpurity(1+i,IE) = densityperimpurity(1+i,IE) -&
                        &(1.0/PI)*AIMAG(GREENR(j + FTIMO, j + FTIMO))
                    end do
                end do
            endif

            ! Writes the Green impurity elements on the corresponding .txt file
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1, '(F17.8, F17.8)') GREENR(i,j)
                end do
            end do

        end do ! ends energies sum
        close(1) ! Closes the .txt file

        ! The density per atom at the future impurity sites.
        if (NPOL == 0) then 
            open(1, file = 'hostdensities.txt', action = 'write')
            do j = 1, NPNT ! Energies = Intervals + 1
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

    end subroutine PARTIALGREEN

    subroutine GAUFD(XI,WI,N)

        integer :: N
        real*8 :: WI(*), XI(*)
        
        if (N == 1) then
            XI(1) = -49817229548128141768.D-20
            WI(1) = 10000000000000031192.D-19
        else if (N == 2) then
            XI(1) = -78465071850839016234.D-20
            XI(2) = -20091536266094051757.D-20
            WI(1) = 50923235990870048433.D-20
            WI(2) = 49076764009130263488.D-20
        else if (N == 3) then
            XI(1) = -88288518955458358024.D-20
            XI(2) = -48117621892777473749.D-20
            XI(3) = -88198184413497647625.D-21
            WI(1) = 28858444436509900908.D-20
            WI(2) = 45966895698954759346.D-20
            WI(3) = 25174659864535651667.D-20
        else if (N == 4) then
            XI(1) = -92613063531202843773.D-20
            XI(2) = -64918327008663578157.D-20
            XI(3) = -28982568853420020298.D-20
            XI(4) = -24595209663255169680.D-21
            WI(1) = 18501429405165520392.D-20
            WI(2) = 34614391006511784214.D-20
            WI(3) = 34152482191988153127.D-20
            WI(4) = 12731697396334854188.D-20
        else if (N == 5) then
            XI(1) = -94875333872503463082.D-20
            XI(2) = -74805843506753178608.D-20
            XI(3) = -45504655263391074765.D-20
            XI(4) = -16657582360358973599.D-20
            XI(5) = 27402283545708211900.D-21
            WI(1) = 12939804504572789754.D-20
            WI(2) = 26102400189213290231.D-20
            WI(3) = 30851911091450589451.D-20
            WI(4) = 24746815229701880449.D-20
            WI(5) = 53590689850617620359.D-21
        else if (N == 6) then
            XI(1) = -96204950250095729781.D-20
            XI(2) = -80971428101130972258.D-20
            XI(3) = -57293627456482418171.D-20
            XI(4) = -30755197635518367504.D-20
            XI(5) = -82123839469384988331.D-21
            XI(6) = 83748358371240941581.D-21
            WI(1) = 96268650841705383829.D-21
            WI(2) = 20246201047059595265.D-20
            WI(3) = 26160719441051813381.D-20
            WI(4) = 25781980698475975536.D-20
            WI(5) = 16683001513553609336.D-20
            WI(6) = 15012322156887800205.D-21
        else if (N == 7) then
            XI(1) = -97053934379083423143.D-20
            XI(2) = -85045695849615413757.D-20
            XI(3) = -65665104053460540522.D-20
            XI(4) = -42357896269371657364.D-20
            XI(5) = -19472732441816555564.D-20
            XI(6) = -19669621223691542539.D-21
            XI(7) = 15142830586888806919.D-20
            WI(1) = 74948008822570509041.D-21
            WI(2) = 16170863905729061704.D-20
            WI(3) = 22007120289205973485.D-20
            WI(4) = 23880411919774885893.D-20
            WI(5) = 20952460047488907594.D-20
            WI(6) = 92465405554445737538.D-21
            WI(7) = 24780240009985858690.D-22
        else if (N == 8) then
            XI(1) = -97630544447925725992.D-20
            XI(2) = -87873822716479965943.D-20
            XI(3) = -71736329217593360204.D-20
            XI(4) = -51463306578144813387.D-20
            XI(5) = -29967081434747298359.D-20
            XI(6) = -10763455942936048359.D-20
            XI(7) = 35963113675701677498.D-21
            XI(8) = 23003149140664609750.D-20
            WI(1) = 60394634019629989770.D-21
            WI(2) = 13252509350880929004.D-20
            WI(3) = 18643612522057003210.D-20
            WI(4) = 21413715867867937533.D-20
            WI(5) = 21005092708864293339.D-20
            WI(6) = 16003068683842947897.D-20
            WI(7) = 36159126989806650464.D-21
            WI(8) = 26624765543536915040.D-23
        else if (N == 9) then
            XI(1) = -98041275487012188695.D-20
            XI(2) = -89918326179154863440.D-20
            XI(3) = -76254129548477842110.D-20
            XI(4) = -58579104527384144901.D-20
            XI(5) = -38924212142470946276.D-20
            XI(6) = -19724340764961096691.D-20
            XI(7) = -40039281758884590381.D-21
            XI(8) = 97228170103579374416.D-21
            XI(9) = 31678885353558278864.D-20
            WI(1) = 49992516372028853833.D-21
            WI(2) = 11099301824870447793.D-20
            WI(3) = 15971411690431220541.D-20
            WI(4) = 19037877203046567198.D-20
            WI(5) = 19869087157813151863.D-20
            WI(6) = 17972334325952047726.D-20
            WI(7) = 10203571121909080322.D-20
            WI(8) = 84501828581921130722.D-22
            WI(9) = 21467529556997868476.D-24
        else if (N == 10) then
            XI(1) = -98345122025502045873.D-20
            XI(2) = -91446749996879318119.D-20
            XI(3) = -79700500547314513626.D-20
            XI(4) = -64189534981349313375.D-20
            XI(5) = -46376588343242516012.D-20
            XI(6) = -28030431525349494354.D-20
            XI(7) = -11327091328726333942.D-20
            XI(8) = 17437648086722052805.D-21
            XI(9) = 16877498338102917782.D-20
            XI(10) = 40960465258252015313.D-20
            WI(1) = 42278597323639457484.D-21
            WI(2) = 94666349251635366832.D-21
            WI(3) = 13843777024241956101.D-20
            WI(4) = 16932936699837666261.D-20
            WI(5) = 18398357022114735352.D-20
            WI(6) = 17939886390638648260.D-20
            WI(7) = 14468854182396060463.D-20
            WI(8) = 46026485095922891703.D-21
            WI(9) = 11890402956686871419.D-22
            WI(10) = 14148408460516817666.D-25
        else if (N == 11) then
            XI(1) = -98576901837451635280.D-20
            XI(2) = -92621727156102677473.D-20
            XI(3) = -82389243156123939088.D-20
            XI(4) = -68670708816882492198.D-20
            XI(5) = -52549052940365991088.D-20
            XI(6) = -35349156561982307316.D-20
            XI(7) = -18652071146560858606.D-20
            XI(8) = -45389164233559550280.D-21
            XI(9) = 76984180593432347734.D-21
            XI(10) = 24899533750455431614.D-20
            XI(11) = 50711636785486806957.D-20
            WI(1) = 36383684790132198923.D-21
            WI(2) = 81985364434128201418.D-21
            WI(3) = 12133566247788805356.D-20
            WI(4) = 15122112006362489825.D-20
            WI(5) = 16900090791849557413.D-20
            WI(6) = 17240157268363508589.D-20
            WI(7) = 15745585899461757802.D-20
            WI(8) = 97600157144810676257.D-21
            WI(9) = 12496828256639735424.D-21
            WI(10) = 11876318920871395759.D-23
            WI(11) = 80046822403386311030.D-27
        else if (N == 12) then
            XI(1) = -98758247347129831371.D-20
            XI(2) = -93546465146779806654.D-20
            XI(3) = -84528996754470930223.D-20
            XI(4) = -72299594230844519839.D-20
            XI(5) = -57679398168141327066.D-20
            XI(6) = -41683730779892996801.D-20
            XI(7) = -25514627335790291149.D-20
            XI(8) = -10710838211747769681.D-20
            XI(9) = 12720145729326415607.D-21
            XI(10) = 14540842218988328389.D-20
            XI(11) = 33552500235752414908.D-20
            XI(12) = 60838109964484063119.D-20
            WI(1) = 31765161579790701148.D-21
            WI(2) = 71927618746964313778.D-21
            WI(3) = 10742555378156694842.D-20
            WI(4) = 13578811351554214795.D-20
            WI(5) = 15492042553417744038.D-20
            WI(6) = 16300300254834219520.D-20
            WI(7) = 15784577013790806216.D-20
            WI(8) = 12921482926208917372.D-20
            WI(9) = 46096943233133302568.D-21
            WI(10) = 20030610755774790850.D-22
            WI(11) = 95165705752725893549.D-25
            WI(12) = 40143360822128708729.D-28
        else if (N == 13) then
            XI(1) = -98903182721370020265.D-20
            XI(2) = -94288936524363459773.D-20
            XI(3) = -86261843870640242196.D-20
            XI(4) = -75277808759167753869.D-20
            XI(5) = -61972590294795871779.D-20
            XI(6) = -47139332563986024748.D-20
            XI(7) = -31718188942187627557.D-20
            XI(8) = -16854863011308355787.D-20
            XI(9) = -41195843159851553906.D-21
            XI(10) = 71957380142115164738.D-21
            XI(11) = 22223926926874000328.D-20
            XI(12) = 42682885634093164862.D-20
            XI(13) = 71270930856714354732.D-20
            WI(1) = 28069991026027589482.D-21
            WI(2) = 63803895087070663653.D-21
            WI(3) = 95973484361405430270.D-21
            WI(4) = 12264378189747678145.D-20
            WI(5) = 14213612346123977130.D-20
            WI(6) = 15296686007570952707.D-20
            WI(7) = 15358437552921000921.D-20
            WI(8) = 14007635729175637795.D-20
            WI(9) = 87531230524252970103.D-21
            WI(10) = 12989730151883234012.D-21
            WI(11) = 22351943999969127535.D-23
            WI(12) = 65097139765619073344.D-26
            WI(13) = 18257341724040876662.D-29
        else if (N == 14) then
            XI(1) = -99021130855943209687.D-20
            XI(2) = -94895368426058288869.D-20
            XI(3) = -87686856465753704289.D-20
            XI(4) = -77752669471002194917.D-20
            XI(5) = -65594116901081876554.D-20
            XI(6) = -51841232227159879604.D-20
            XI(7) = -37243750660439082187.D-20
            XI(8) = -22693429290756856295.D-20
            XI(9) = -93940943648510570987.D-21
            XI(10) = 16521198218716065629.D-21
            XI(11) = 13919799114797561344.D-20
            XI(12) = 30521886852802066309.D-20
            XI(13) = 52192337126752562221.D-20
            XI(14) = 81957965081548293179.D-20
            WI(1) = 25060310888021301605.D-21
            WI(2) = 57137272611562033779.D-21
            WI(3) = 86434450014324433897.D-21
            WI(4) = 11141118228632175288.D-20
            WI(5) = 13070790263291078499.D-20
            WI(6) = 14310195071194851995.D-20
            WI(7) = 14737968606274298328.D-20
            WI(8) = 14154903694980505066.D-20
            WI(9) = 11456160782223814050.D-20
            WI(10) = 40466499493397342820.D-21
            WI(11) = 21701008894932486895.D-22
            WI(12) = 19960253076851250807.D-24
            WI(13) = 39376501060604877095.D-27
            WI(14) = 76596142918862399780.D-31
        else if (N == 15) then
            XI(1) = -99118619138431485634.D-20
            XI(2) = -95398089203095832045.D-20
            XI(3) = -88874665207045485764.D-20
            XI(4) = -79832886799647722652.D-20
            XI(5) = -68674462209286747178.D-20
            XI(6) = -55907326778454372362.D-20
            XI(7) = -42138595122137487519.D-20
            XI(8) = -28083407355763995168.D-20
            XI(9) = -14649293944496725019.D-20
            XI(10) = -30865949117072113052.D-21
            XI(11) = 75989566859912966734.D-21
            XI(12) = 21425891814116860148.D-20
            XI(13) = 39280262275215780450.D-20
            XI(14) = 62012182191671475949.D-20
            XI(15) = 92858877219218103945.D-20
            WI(1) = 22570991165870390473.D-21
            WI(2) = 51589746641923392000.D-21
            WI(3) = 78401918844466166239.D-21
            WI(4) = 10176234626640128024.D-20
            WI(5) = 12055819130110177262.D-20
            WI(6) = 13377324647273569326.D-20
            WI(7) = 14041818603247829422.D-20
            WI(8) = 13919569003129657925.D-20
            WI(9) = 12562361445602688222.D-20
            WI(10) = 74852662340708470150.D-21
            WI(11) = 10996744175647251144.D-21
            WI(12) = 25513307315040157893.D-23
            WI(13) = 15270418102934789627.D-25
            WI(14) = 21560859319293022163.D-28
            WI(15) = 30032040385910287756.D-32
        else if (N == 16) then
            XI(1) = -99200289748411473927.D-20
            XI(2) = -95820266378296634182.D-20
            XI(3) = -89876661129475763142.D-20
            XI(4) = -81599671254865832971.D-20
            XI(5) = -71315812647978444249.D-20
            XI(6) = -59440032425488487666.D-20
            XI(7) = -46470396871945791541.D-20
            XI(8) = -32991653294098863600.D-20
            XI(9) = -19716091326862980561.D-20
            XI(10) = -76605243443508959615.D-21
            XI(11) = 26155046503992069925.D-21
            XI(12) = 14307776307824938682.D-20
            XI(13) = 29506185654032182160.D-20
            XI(14) = 48403577800553841578.D-20
            XI(15) = 72091584865612160132.D-20
            XI(16) = 10394188783181811718.D-19
            WI(1) = 20484388078614008045.D-21
            WI(2) = 46916532350372347409.D-21
            WI(3) = 71569877291069983495.D-21
            WI(4) = 93424466379672137196.D-21
            WI(5) = 11156011364306951297.D-20
            WI(6) = 12512553084306063601.D-20
            WI(7) = 13329704953113185969.D-20
            WI(8) = 13510959073859290681.D-20
            WI(9) = 12840858805365359846.D-20
            WI(10) = 10016528657871746742.D-20
            WI(11) = 32102655847303900301.D-21
            WI(12) = 18115418480524121495.D-22
            WI(13) = 24274994772381143993.D-24
            WI(14) = 10371321943363515335.D-26
            WI(15) = 10868941709467004901.D-29
            WI(16) = 11117372791599461059.D-33
        endif

    end subroutine GAUFD

    subroutine GAULEG(XI,WI,N)

        integer :: N, I
        real*8 :: WI(*), XI(*)
        
        if (N > 32) N = ((N-1)/4+1)*4
        if (N > 64) N = ((N-1)/8+1)*8

        if (N == 1) then
            XI(1) = 0.D0
            WI(1) = 2.D0
        else if (N == 2) then
            XI(1) = -57735026918962576451.D-20
            WI(1) = 10000000000000000000.D-19
        else if (N == 3) then
            XI(1) = -77459666924148337704.D-20
            XI(2) = 0.D0
            WI(1) = 55555555555555555556.D-20
            WI(2) = 88888888888888888889.D-20
        else if (N == 4) then
            XI(1) = -86113631159405257522.D-20
            XI(2) = -33998104358485626480.D-20
            WI(1) = 34785484513745385737.D-20
            WI(2) = 65214515486254614263.D-20
        else if (N == 5) then
            XI(1) = -90617984593866399280.D-20
            XI(2) = -53846931010568309104.D-20
            XI(3) = 0.D0
            WI(1) = 23692688505618908751.D-20
            WI(2) = 47862867049936646804.D-20
            WI(3) = 56888888888888888889.D-20
        else if (N == 6) then
            XI(1) = -93246951420315202781.D-20
            XI(2) = -66120938646626451366.D-20
            XI(3) = -23861918608319690863.D-20
            WI(1) = 17132449237917034504.D-20
            WI(2) = 36076157304813860757.D-20
            WI(3) = 46791393457269104739.D-20
        else if (N == 7) then
            XI(1) = -94910791234275852453.D-20
            XI(2) = -74153118559939443986.D-20
            XI(3) = -40584515137739716691.D-20
            XI(4) = 0.D0
            WI(1) = 12948496616886969327.D-20
            WI(2) = 27970539148927666790.D-20
            WI(3) = 38183005050511894495.D-20
            WI(4) = 41795918367346938776.D-20
        else if (N == 8) then
            XI(1) = -96028985649753623168.D-20
            XI(2) = -79666647741362673959.D-20
            XI(3) = -52553240991632898582.D-20
            XI(4) = -18343464249564980494.D-20
            WI(1) = 10122853629037625915.D-20
            WI(2) = 22238103445337447054.D-20
            WI(3) = 31370664587788728734.D-20
            WI(4) = 36268378337836198297.D-20
        else if (N == 9) then
            XI(1) = -96816023950762608984.D-20
            XI(2) = -83603110732663579430.D-20
            XI(3) = -61337143270059039731.D-20
            XI(4) = -32425342340380892904.D-20
            XI(5) = 0.D0
            WI(1) = 81274388361574411972.D-21
            WI(2) = 18064816069485740406.D-20
            WI(3) = 26061069640293546232.D-20
            WI(4) = 31234707704000284007.D-20
            WI(5) = 33023935500125976316.D-20
        else if (N == 10) then
            XI(1) = -97390652851717172008.D-20
            XI(2) = -86506336668898451073.D-20
            XI(3) = -67940956829902440623.D-20
            XI(4) = -43339539412924719080.D-20
            XI(5) = -14887433898163121088.D-20
            WI(1) = 66671344308688137594.D-21
            WI(2) = 14945134915058059315.D-20
            WI(3) = 21908636251598204400.D-20
            WI(4) = 26926671930999635509.D-20
            WI(5) = 29552422471475287017.D-20
        else if (N == 11) then
            XI(1) = -97822865814605699280.D-20
            XI(2) = -88706259976809529908.D-20
            XI(3) = -73015200557404932409.D-20
            XI(4) = -51909612920681181593.D-20
            XI(5) = -26954315595234497233.D-20
            XI(6) = 0.D0
            WI(1) = 55668567116173666483.D-21
            WI(2) = 12558036946490462463.D-20
            WI(3) = 18629021092773425143.D-20
            WI(4) = 23319376459199047992.D-20
            WI(5) = 26280454451024666218.D-20
            WI(6) = 27292508677790063071.D-20
        else if (N == 12) then
            XI(1) = -98156063424671925069.D-20
            XI(2) = -90411725637047485668.D-20
            XI(3) = -76990267419430468704.D-20
            XI(4) = -58731795428661744730.D-20
            XI(5) = -36783149899818019375.D-20
            XI(6) = -12523340851146891547.D-20
            WI(1) = 47175336386511827195.D-21
            WI(2) = 10693932599531843096.D-20
            WI(3) = 16007832854334622633.D-20
            WI(4) = 20316742672306592175.D-20
            WI(5) = 23349253653835480876.D-20
            WI(6) = 24914704581340278500.D-20
        else if (N == 13) then
            XI(1) = -98418305471858814947.D-20
            XI(2) = -91759839922297796521.D-20
            XI(3) = -80157809073330991279.D-20
            XI(4) = -64234933944034022064.D-20
            XI(5) = -44849275103644685288.D-20
            XI(6) = -23045831595513479407.D-20
            XI(7) = 0.D0
            WI(1) = 40484004765315879520.D-21
            WI(2) = 92121499837728447914.D-21
            WI(3) = 13887351021978723846.D-20
            WI(4) = 17814598076194573828.D-20
            WI(5) = 20781604753688850231.D-20
            WI(6) = 22628318026289723841.D-20
            WI(7) = 23255155323087391019.D-20
        else if (N == 14) then
            XI(1) = -98628380869681233884.D-20
            XI(2) = -92843488366357351734.D-20
            XI(3) = -82720131506976499319.D-20
            XI(4) = -68729290481168547015.D-20
            XI(5) = -51524863635815409197.D-20
            XI(6) = -31911236892788976044.D-20
            XI(7) = -10805494870734366207.D-20
            WI(1) = 35119460331751863032.D-21
            WI(2) = 80158087159760209806.D-21
            WI(3) = 12151857068790318469.D-20
            WI(4) = 15720316715819353457.D-20
            WI(5) = 18553839747793781374.D-20
            WI(6) = 20519846372129560397.D-20
            WI(7) = 21526385346315779020.D-20
        else if (N == 15) then
            XI(1) = -98799251802048542849.D-20
            XI(2) = -93727339240070590431.D-20
            XI(3) = -84820658341042721620.D-20
            XI(4) = -72441773136017004742.D-20
            XI(5) = -57097217260853884754.D-20
            XI(6) = -39415134707756336990.D-20
            XI(7) = -20119409399743452230.D-20
            XI(8) = 0.D0
            WI(1) = 30753241996117268355.D-21
            WI(2) = 70366047488108124709.D-21
            WI(3) = 10715922046717193501.D-20
            WI(4) = 13957067792615431445.D-20
            WI(5) = 16626920581699393355.D-20
            WI(6) = 18616100001556221103.D-20
            WI(7) = 19843148532711157646.D-20
            WI(8) = 20257824192556127288.D-20
        else if (N == 16) then
            XI(1) = -98940093499164993260.D-20
            XI(2) = -94457502307323257608.D-20
            XI(3) = -86563120238783174388.D-20
            XI(4) = -75540440835500303390.D-20
            XI(5) = -61787624440264374845.D-20
            XI(6) = -45801677765722738634.D-20
            XI(7) = -28160355077925891323.D-20
            XI(8) = -95012509837637440185.D-21
            WI(1) = 27152459411754094852.D-21
            WI(2) = 62253523938647892863.D-21
            WI(3) = 95158511682492784810.D-21
            WI(4) = 12462897125553387205.D-20
            WI(5) = 14959598881657673208.D-20
            WI(6) = 16915651939500253819.D-20
            WI(7) = 18260341504492358887.D-20
            WI(8) = 18945061045506849629.D-20
        else if (N == 17) then
            XI(1) = -99057547531441733568.D-20
            XI(2) = -95067552176876776122.D-20
            XI(3) = -88023915372698590212.D-20
            XI(4) = -78151400389680140693.D-20
            XI(5) = -65767115921669076585.D-20
            XI(6) = -51269053708647696789.D-20
            XI(7) = -35123176345387631530.D-20
            XI(8) = -17848418149584785585.D-20
            XI(9) = 0.D0
            WI(1) = 24148302868547931960.D-21
            WI(2) = 55459529373987201129.D-21
            WI(3) = 85036148317179180884.D-21
            WI(4) = 11188384719340397109.D-20
            WI(5) = 13513636846852547329.D-20
            WI(6) = 15404576107681028808.D-20
            WI(7) = 16800410215645004451.D-20
            WI(8) = 17656270536699264633.D-20
            WI(9) = 17944647035620652546.D-20
        else if (N == 18) then
            XI(1) = -99156516842093094673.D-20
            XI(2) = -95582394957139775518.D-20
            XI(3) = -89260246649755573921.D-20
            XI(4) = -80370495897252311568.D-20
            XI(5) = -69168704306035320787.D-20
            XI(6) = -55977083107394753461.D-20
            XI(7) = -41175116146284264604.D-20
            XI(8) = -25188622569150550959.D-20
            XI(9) = -84775013041735301242.D-21
            WI(1) = 21616013526483310313.D-21
            WI(2) = 49714548894969796453.D-21
            WI(3) = 76425730254889056529.D-21
            WI(4) = 10094204410628716556.D-20
            WI(5) = 12255520671147846018.D-20
            WI(6) = 14064291467065065120.D-20
            WI(7) = 15468467512626524493.D-20
            WI(8) = 16427648374583272299.D-20
            WI(9) = 16914238296314359184.D-20
        else if (N == 19) then
            XI(1) = -99240684384358440319.D-20
            XI(2) = -96020815213483003085.D-20
            XI(3) = -90315590361481790164.D-20
            XI(4) = -82271465653714282498.D-20
            XI(5) = -72096617733522937862.D-20
            XI(6) = -60054530466168102347.D-20
            XI(7) = -46457074137596094572.D-20
            XI(8) = -31656409996362983199.D-20
            XI(9) = -16035864564022537587.D-20
            XI(10) = 0.D0
            WI(1) = 19461788229726477036.D-21
            WI(2) = 44814226765699600333.D-21
            WI(3) = 69044542737641226581.D-21
            WI(4) = 91490021622449999464.D-21
            WI(5) = 11156664554733399472.D-20
            WI(6) = 12875396253933622768.D-20
            WI(7) = 14260670217360661178.D-20
            WI(8) = 15276604206585966678.D-20
            WI(9) = 15896884339395434765.D-20
            WI(10) = 16105444984878369598.D-20
        else if (N == 20) then
            XI(1) = -99312859918509492479.D-20
            XI(2) = -96397192727791379127.D-20
            XI(3) = -91223442825132590587.D-20
            XI(4) = -83911697182221882339.D-20
            XI(5) = -74633190646015079261.D-20
            XI(6) = -63605368072651502545.D-20
            XI(7) = -51086700195082709800.D-20
            XI(8) = -37370608871541956067.D-20
            XI(9) = -22778585114164507808.D-20
            XI(10) = -76526521133497333755.D-21
            WI(1) = 17614007139152118312.D-21
            WI(2) = 40601429800386941331.D-21
            WI(3) = 62672048334109063570.D-21
            WI(4) = 83276741576704748725.D-21
            WI(5) = 10193011981724043504.D-20
            WI(6) = 11819453196151841731.D-20
            WI(7) = 13168863844917662690.D-20
            WI(8) = 14209610931838205133.D-20
            WI(9) = 14917298647260374679.D-20
            WI(10) = 15275338713072585070.D-20
        else if (N == 21) then
            XI(1) = -99375217062038950026.D-20
            XI(2) = -96722683856630629432.D-20
            XI(3) = -92009933415040082879.D-20
            XI(4) = -85336336458331728365.D-20
            XI(5) = -76843996347567790862.D-20
            XI(6) = -66713880419741231931.D-20
            XI(7) = -55161883588721980706.D-20
            XI(8) = -42434212020743878357.D-20
            XI(9) = -28802131680240109660.D-20
            XI(10) = -14556185416089509094.D-20
            XI(11) = 0.D0
            WI(1) = 16017228257774333324.D-21
            WI(2) = 36953789770852493800.D-21
            WI(3) = 57134425426857208284.D-21
            WI(4) = 76100113628379302017.D-21
            WI(5) = 93444423456033861553.D-21
            WI(6) = 10879729916714837766.D-20
            WI(7) = 12183141605372853420.D-20
            WI(8) = 13226893863333746178.D-20
            WI(9) = 13988739479107315472.D-20
            WI(10) = 14452440398997005906.D-20
            WI(11) = 14608113364969042719.D-20
        else if (N == 22) then
            XI(1) = -99429458548239929207.D-20
            XI(2) = -97006049783542872712.D-20
            XI(3) = -92695677218717400052.D-20
            XI(4) = -86581257772030013654.D-20
            XI(5) = -78781680597920816200.D-20
            XI(6) = -69448726318668278005.D-20
            XI(7) = -58764040350691159296.D-20
            XI(8) = -46935583798675702641.D-20
            XI(9) = -34193582089208422516.D-20
            XI(10) = -20786042668822128548.D-20
            XI(11) = -69739273319722221214.D-21
            WI(1) = 14627995298272200685.D-21
            WI(2) = 33774901584814154793.D-21
            WI(3) = 52293335152683285940.D-21
            WI(4) = 69796468424520488095.D-21
            WI(5) = 85941606217067727414.D-21
            WI(6) = 10041414444288096493.D-20
            WI(7) = 11293229608053921839.D-20
            WI(8) = 12325237681051242429.D-20
            WI(9) = 13117350478706237073.D-20
            WI(10) = 13654149834601517135.D-20
            WI(11) = 13925187285563199338.D-20
        else if (N == 23) then
            XI(1) = -99476933499755212352.D-20
            XI(2) = -97254247121811523196.D-20
            XI(3) = -93297108682601610235.D-20
            XI(4) = -87675235827044166738.D-20
            XI(5) = -80488840161883989215.D-20
            XI(6) = -71866136313195019446.D-20
            XI(7) = -61960987576364615639.D-20
            XI(8) = -50950147784600754969.D-20
            XI(9) = -39030103803029083142.D-20
            XI(10) = -26413568097034493053.D-20
            XI(11) = -13325682429846611093.D-20
            XI(12) = 0.D0
            WI(1) = 13411859487141772081.D-21
            WI(2) = 30988005856979444311.D-21
            WI(3) = 48037671731084668572.D-21
            WI(4) = 64232421408525852127.D-21
            WI(5) = 79281411776718954923.D-21
            WI(6) = 92915766060035147477.D-21
            WI(7) = 10489209146454141007.D-20
            WI(8) = 11499664022241136494.D-20
            WI(9) = 12304908430672953047.D-20
            WI(10) = 12890572218808214998.D-20
            WI(11) = 13246203940469661737.D-20
            WI(12) = 13365457218610617535.D-20
        else if (N == 24) then
            XI(1) = -99518721999702136018.D-20
            XI(2) = -97472855597130949820.D-20
            XI(3) = -93827455200273275852.D-20
            XI(4) = -88641552700440103421.D-20
            XI(5) = -82000198597390292195.D-20
            XI(6) = -74012419157855436424.D-20
            XI(7) = -64809365193697556925.D-20
            XI(8) = -54542147138883953566.D-20
            XI(9) = -43379350762604513849.D-20
            XI(10) = -31504267969616337439.D-20
            XI(11) = -19111886747361630916.D-20
            XI(12) = -64056892862605626085.D-21
            WI(1) = 12341229799987199547.D-21
            WI(2) = 28531388628933663181.D-21
            WI(3) = 44277438817419806169.D-21
            WI(4) = 59298584915436780746.D-21
            WI(5) = 73346481411080305734.D-21
            WI(6) = 86190161531953275917.D-21
            WI(7) = 97618652104113888270.D-21
            WI(8) = 10744427011596563478.D-20
            WI(9) = 11550566805372560135.D-20
            WI(10) = 12167047292780339120.D-20
            WI(11) = 12583745634682829612.D-20
            WI(12) = 12793819534675215697.D-20
        else if (N == 25) then
            XI(1) = -99555696979049809791.D-20
            XI(2) = -97666392145951751150.D-20
            XI(3) = -94297457122897433941.D-20
            XI(4) = -89499199787827536885.D-20
            XI(5) = -83344262876083400142.D-20
            XI(6) = -75925926303735763058.D-20
            XI(7) = -67356636847346836449.D-20
            XI(8) = -57766293024122296772.D-20
            XI(9) = -47300273144571496052.D-20
            XI(10) = -36117230580938783774.D-20
            XI(11) = -24386688372098843205.D-20
            XI(12) = -12286469261071039639.D-20
            XI(13) = 0.D0
            WI(1) = 11393798501026287948.D-21
            WI(2) = 26354986615032137262.D-21
            WI(3) = 40939156701306312656.D-21
            WI(4) = 54904695975835191926.D-21
            WI(5) = 68038333812356917207.D-21
            WI(6) = 80140700335001018013.D-21
            WI(7) = 91028261982963649811.D-21
            WI(8) = 10053594906705064420.D-20
            WI(9) = 10851962447426365312.D-20
            WI(10) = 11485825914571164834.D-20
            WI(11) = 11945576353578477223.D-20
            WI(12) = 12224244299031004169.D-20
            WI(13) = 12317605372671545120.D-20
        else if (N == 26) then
            XI(1) = -99588570114561692900.D-20
            XI(2) = -97838544595647099110.D-20
            XI(3) = -94715906666171425014.D-20
            XI(4) = -90263786198430707422.D-20
            XI(5) = -84544594278849801880.D-20
            XI(6) = -77638594882067885619.D-20
            XI(7) = -69642726041995726486.D-20
            XI(8) = -60669229301761806323.D-20
            XI(9) = -50844071482450571770.D-20
            XI(10) = -40305175512348630648.D-20
            XI(11) = -29200483948595689514.D-20
            XI(12) = -17685882035689018397.D-20
            XI(13) = -59230093429313207094.D-21
            WI(1) = 10551372617343007156.D-21
            WI(2) = 24417851092631908790.D-21
            WI(3) = 37962383294362763950.D-21
            WI(4) = 50975825297147811998.D-21
            WI(5) = 63274046329574835539.D-21
            WI(6) = 74684149765659745887.D-21
            WI(7) = 85045894313485239210.D-21
            WI(8) = 94213800355914148464.D-21
            WI(9) = 10205916109442542324.D-20
            WI(10) = 10847184052857659066.D-20
            WI(11) = 11336181654631966655.D-20
            WI(12) = 11666044348529658204.D-20
            WI(13) = 11832141527926227652.D-20
        else if (N == 27) then
            XI(1) = -99617926288898856694.D-20
            XI(2) = -97992347596150122286.D-20
            XI(3) = -95090055781470500685.D-20
            XI(4) = -90948232067749110430.D-20
            XI(5) = -85620790801829449030.D-20
            XI(6) = -79177163907050822714.D-20
            XI(7) = -71701347373942369929.D-20
            XI(8) = -63290797194649514093.D-20
            XI(9) = -54055156457945689490.D-20
            XI(10) = -44114825175002688059.D-20
            XI(11) = -33599390363850889973.D-20
            XI(12) = -22645936543953685886.D-20
            XI(13) = -11397258560952996693.D-20
            XI(14) = 0.D0
            WI(1) = 97989960512943602612.D-22
            WI(2) = 22686231596180623196.D-21
            WI(3) = 35297053757419711023.D-21
            WI(4) = 47449412520615062704.D-21
            WI(5) = 58983536859833599110.D-21
            WI(6) = 69748823766245592984.D-21
            WI(7) = 79604867773057771263.D-21
            WI(8) = 88423158543756950194.D-21
            WI(9) = 96088727370028507566.D-21
            WI(10) = 10250163781774579867.D-20
            WI(11) = 10757828578853318721.D-20
            WI(12) = 11125248835684519267.D-20
            WI(13) = 11347634610896514862.D-20
            WI(14) = 11422086737895698905.D-20
        else if (N == 28) then
            XI(1) = -99644249757395444995.D-20
            XI(2) = -98130316537087275369.D-20
            XI(3) = -95425928062893819725.D-20
            XI(4) = -91563302639213207387.D-20
            XI(5) = -86589252257439504894.D-20
            XI(6) = -80564137091717917145.D-20
            XI(7) = -73561087801363177203.D-20
            XI(8) = -65665109403886496122.D-20
            XI(9) = -56972047181140171931.D-20
            XI(10) = -47587422495511826103.D-20
            XI(11) = -37625151608907871022.D-20
            XI(12) = -27206162763517807768.D-20
            XI(13) = -16456928213338077128.D-20
            XI(14) = -55079289884034270427.D-21
            WI(1) = 91242825930945177388.D-22
            WI(2) = 21132112592771259752.D-21
            WI(3) = 32901427782304379978.D-21
            WI(4) = 44272934759004227840.D-21
            WI(5) = 55107345675716745431.D-21
            WI(6) = 65272923966999595793.D-21
            WI(7) = 74646214234568779024.D-21
            WI(8) = 83113417228901218390.D-21
            WI(9) = 90571744393032840942.D-21
            WI(10) = 96930657997929915850.D-21
            WI(11) = 10211296757806076981.D-20
            WI(12) = 10605576592284641791.D-20
            WI(13) = 10871119225829413525.D-20
            WI(14) = 11004701301647519628.D-20
        else if (N == 29) then
            XI(1) = -99667944226059658616.D-20
            XI(2) = -98254550526141317487.D-20
            XI(3) = -95728559577808772580.D-20
            XI(4) = -92118023295305878509.D-20
            XI(5) = -87463780492010279042.D-20
            XI(6) = -81818548761525244499.D-20
            XI(7) = -75246285173447713391.D-20
            XI(8) = -67821453760268651516.D-20
            XI(9) = -59628179713822782038.D-20
            XI(10) = -50759295512422764210.D-20
            XI(11) = -41315288817400866389.D-20
            XI(12) = -31403163786763993495.D-20
            XI(13) = -21135228616600107451.D-20
            XI(14) = -10627823013267923017.D-20
            XI(15) = 0.D0
            WI(1) = 85169038787464096543.D-22
            WI(2) = 19732085056122705984.D-21
            WI(3) = 30740492202093622644.D-21
            WI(4) = 41402062518682836105.D-21
            WI(5) = 51594826902497923913.D-21
            WI(6) = 61203090657079138542.D-21
            WI(7) = 70117933255051278570.D-21
            WI(8) = 78238327135763783828.D-21
            WI(9) = 85472257366172527545.D-21
            WI(10) = 91737757139258763348.D-21
            WI(11) = 96963834094408606302.D-21
            WI(12) = 10109127375991496612.D-20
            WI(13) = 10407331007772937391.D-20
            WI(14) = 10587615509732094141.D-20
            WI(15) = 10647938171831424425.D-20
        else if (N == 30) then
            XI(1) = -99689348407464954027.D-20
            XI(2) = -98366812327974720997.D-20
            XI(3) = -96002186496830751222.D-20
            XI(4) = -92620004742927432588.D-20
            XI(5) = -88256053579205268154.D-20
            XI(6) = -82956576238276839744.D-20
            XI(7) = -76777743210482619492.D-20
            XI(8) = -69785049479331579693.D-20
            XI(9) = -62052618298924286114.D-20
            XI(10) = -53662414814201989926.D-20
            XI(11) = -44703376953808917678.D-20
            XI(12) = -35270472553087811347.D-20
            XI(13) = -25463692616788984644.D-20
            XI(14) = -15386991360858354696.D-20
            XI(15) = -51471842555317695833.D-21
            WI(1) = 79681924961666056155.D-22
            WI(2) = 18466468311090959142.D-21
            WI(3) = 28784707883323369350.D-21
            WI(4) = 38799192569627049597.D-21
            WI(5) = 48402672830594052903.D-21
            WI(6) = 57493156217619066482.D-21
            WI(7) = 65974229882180495128.D-21
            WI(8) = 73755974737705206268.D-21
            WI(9) = 80755895229420215355.D-21
            WI(10) = 86899787201082979802.D-21
            WI(11) = 92122522237786128718.D-21
            WI(12) = 96368737174644259639.D-21
            WI(13) = 99593420586795267063.D-21
            WI(14) = 10176238974840550460.D-20
            WI(15) = 10285265289355884034.D-20
        else if (N == 31) then
            XI(1) = -99708748181947707406.D-20
            XI(2) = -98468590966515248400.D-20
            XI(3) = -96250392509294966179.D-20
            XI(4) = -93075699789664816496.D-20
            XI(5) = -88976002994827104337.D-20
            XI(6) = -83992032014626734009.D-20
            XI(7) = -78173314841662494041.D-20
            XI(8) = -71577678458685328391.D-20
            XI(9) = -64270672292426034618.D-20
            XI(10) = -56324916140714926272.D-20
            XI(11) = -47819378204490248044.D-20
            XI(12) = -38838590160823294306.D-20
            XI(13) = -29471806998170161662.D-20
            XI(14) = -19812119933557062877.D-20
            XI(15) = -99555312152341520325.D-21
            XI(16) = 0.D0
            WI(1) = 74708315792487758587.D-22
            WI(2) = 17318620790310582463.D-21
            WI(3) = 27009019184979421801.D-21
            WI(4) = 36432273912385464024.D-21
            WI(5) = 45493707527201102902.D-21
            WI(6) = 54103082424916853712.D-21
            WI(7) = 62174786561028426910.D-21
            WI(8) = 69628583235410366168.D-21
            WI(9) = 76390386598776616426.D-21
            WI(10) = 82392991761589263904.D-21
            WI(11) = 87576740608477876126.D-21
            WI(12) = 91890113893641478215.D-21
            WI(13) = 95290242912319512807.D-21
            WI(14) = 97743335386328725093.D-21
            WI(15) = 99225011226672307875.D-21
            WI(16) = 99720544793426451428.D-21
        else if (N == 32) then
            XI(1) = -99726386184948156354.D-20
            XI(2) = -98561151154526833540.D-20
            XI(3) = -96476225558750643077.D-20
            XI(4) = -93490607593773968917.D-20
            XI(5) = -89632115576605212397.D-20
            XI(6) = -84936761373256997013.D-20
            XI(7) = -79448379596794240696.D-20
            XI(8) = -73218211874028968039.D-20
            XI(9) = -66304426693021520098.D-20
            XI(10) = -58771575724076232904.D-20
            XI(11) = -50689990893222939002.D-20
            XI(12) = -42135127613063534536.D-20
            XI(13) = -33186860228212764978.D-20
            XI(14) = -23928736225213707454.D-20
            XI(15) = -14447196158279649349.D-20
            XI(16) = -48307665687738316235.D-21
            WI(1) = 70186100094700966004.D-22
            WI(2) = 16274394730905670605.D-21
            WI(3) = 25392065309262059456.D-21
            WI(4) = 34273862913021433103.D-21
            WI(5) = 42835898022226680657.D-21
            WI(6) = 50998059262376176196.D-21
            WI(7) = 58684093478535547145.D-21
            WI(8) = 65822222776361846838.D-21
            WI(9) = 72345794108848506225.D-21
            WI(10) = 78193895787070306472.D-21
            WI(11) = 83311924226946755222.D-21
            WI(12) = 87652093004403811143.D-21
            WI(13) = 91173878695763884713.D-21
            WI(14) = 93844399080804565639.D-21
            WI(15) = 95638720079274859419.D-21
            WI(16) = 96540088514727800567.D-21
        else if (N == 36) then
            XI(1) = -99783046248408583620.D-20
            XI(2) = -98858647890221223807.D-20
            XI(3) = -97202769104969794934.D-20
            XI(4) = -94827298439950754520.D-20
            XI(5) = -91749777451565906608.D-20
            XI(6) = -87992980089039713198.D-20
            XI(7) = -83584716699247530642.D-20
            XI(8) = -78557623013220651283.D-20
            XI(9) = -72948917159355658209.D-20
            XI(10) = -66800123658552106210.D-20
            XI(11) = -60156765813598053508.D-20
            XI(12) = -53068028592624516164.D-20
            XI(13) = -45586394443342026721.D-20
            XI(14) = -37767254711968921632.D-20
            XI(15) = -29668499534402827050.D-20
            XI(16) = -21350089231686557894.D-20
            XI(17) = -12873610380938478865.D-20
            XI(18) = -43018198473708607227.D-21
            WI(1) = 55657196642450453613.D-22
            WI(2) = 12915947284065574405.D-21
            WI(3) = 20181515297735471532.D-21
            WI(4) = 27298621498568779094.D-21
            WI(5) = 34213810770307229921.D-21
            WI(6) = 40875750923644895474.D-21
            WI(7) = 47235083490265978417.D-21
            WI(8) = 53244713977759919092.D-21
            WI(9) = 58860144245324817310.D-21
            WI(10) = 64039797355015489556.D-21
            WI(11) = 68745323835736442614.D-21
            WI(12) = 72941885005653061354.D-21
            WI(13) = 76598410645870674529.D-21
            WI(14) = 79687828912071601909.D-21
            WI(15) = 82187266704339709517.D-21
            WI(16) = 84078218979661934933.D-21
            WI(17) = 85346685739338627492.D-21
            WI(18) = 85983275670394747490.D-21
        else if (N == 40) then
            XI(1) = -99823770971055920035.D-20
            XI(2) = -99072623869945700645.D-20
            XI(3) = -97725994998377426266.D-20
            XI(4) = -95791681921379165580.D-20
            XI(5) = -93281280827867653336.D-20
            XI(6) = -90209880696887429673.D-20
            XI(7) = -86595950321225950382.D-20
            XI(8) = -82461223083331166320.D-20
            XI(9) = -77830565142651938769.D-20
            XI(10) = -72731825518992710328.D-20
            XI(11) = -67195668461417954838.D-20
            XI(12) = -61255388966798023795.D-20
            XI(13) = -54946712509512820208.D-20
            XI(14) = -48307580168617871291.D-20
            XI(15) = -41377920437160500152.D-20
            XI(16) = -34199409082575847301.D-20
            XI(17) = -26815218500725368114.D-20
            XI(18) = -19269758070137109972.D-20
            XI(19) = -11608407067525520848.D-20
            XI(20) = -38772417506050821933.D-21
            WI(1) = 45212770985331912585.D-22
            WI(2) = 10498284531152813615.D-21
            WI(3) = 16421058381907888713.D-21
            WI(4) = 22245849194166957262.D-21
            WI(5) = 27937006980023401098.D-21
            WI(6) = 33460195282547847393.D-21
            WI(7) = 38782167974472017640.D-21
            WI(8) = 43870908185673271992.D-21
            WI(9) = 48695807635072232061.D-21
            WI(10) = 53227846983936824355.D-21
            WI(11) = 57439769099391551367.D-21
            WI(12) = 61306242492928939167.D-21
            WI(13) = 64804013456601038075.D-21
            WI(14) = 67912045815233903826.D-21
            WI(15) = 70611647391286779695.D-21
            WI(16) = 72886582395804059061.D-21
            WI(17) = 74723169057968264200.D-21
            WI(18) = 76110361900626242372.D-21
            WI(19) = 77039818164247965588.D-21
            WI(20) = 77505947978424811264.D-21
        else if (N == 44) then
            XI(1) = -99854020063677422494.D-20
            XI(2) = -99231639213851580848.D-20
            XI(3) = -98115183307791396666.D-20
            XI(4) = -96509965042249313939.D-20
            XI(5) = -94423950911819409920.D-20
            XI(6) = -91867525998417577432.D-20
            XI(7) = -88853423828604320234.D-20
            XI(8) = -85396659500471037873.D-20
            XI(9) = -81514453964513501049.D-20
            XI(10) = -77226147924875589902.D-20
            XI(11) = -72553105366071700261.D-20
            XI(12) = -67518607066612236533.D-20
            XI(13) = -62147734590357584780.D-20
            XI(14) = -56467245318547076842.D-20
            XI(15) = -50505439138820231798.D-20
            XI(16) = -44292017452541148383.D-20
            XI(17) = -37857935201470713251.D-20
            XI(18) = -31235246650278581224.D-20
            XI(19) = -24456945692820125151.D-20
            XI(20) = -17556801477551678575.D-20
            XI(21) = -10569190170865324712.D-20
            XI(22) = -35289236964135359058.D-21
            WI(1) = 37454048031127775152.D-22
            WI(2) = 87004813675248441226.D-22
            WI(3) = 13619586755579985520.D-21
            WI(4) = 18471481736814749172.D-21
            WI(5) = 23231481902019210629.D-21
            WI(6) = 27875782821281010081.D-21
            WI(7) = 32381222812069820881.D-21
            WI(8) = 36725347813808873643.D-21
            WI(9) = 40886512310346218908.D-21
            WI(10) = 44843984081970031446.D-21
            WI(11) = 48578046448352037528.D-21
            WI(12) = 52070096091704461881.D-21
            WI(13) = 55302735563728052549.D-21
            WI(14) = 58259859877595495334.D-21
            WI(15) = 60926736701561968039.D-21
            WI(16) = 63290079733203854950.D-21
            WI(17) = 65338114879181434984.D-21
            WI(18) = 67060638906293652396.D-21
            WI(19) = 68449070269366660985.D-21
            WI(20) = 69496491861572578037.D-21
            WI(21) = 70197685473558212587.D-21
            WI(22) = 70549157789354068811.D-21
        else if (N == 48) then
            XI(1) = -99877100725242611860.D-20
            XI(2) = -99353017226635075755.D-20
            XI(3) = -98412458372282685774.D-20
            XI(4) = -97059159254624725046.D-20
            XI(5) = -95298770316043086072.D-20
            XI(6) = -93138669070655433311.D-20
            XI(7) = -90587913671556967282.D-20
            XI(8) = -87657202027424788591.D-20
            XI(9) = -84358826162439353071.D-20
            XI(10) = -80706620402944262708.D-20
            XI(11) = -76715903251574033925.D-20
            XI(12) = -72403413092381465467.D-20
            XI(13) = -67787237963266390521.D-20
            XI(14) = -62886739677651362400.D-20
            XI(15) = -57722472608397270382.D-20
            XI(16) = -52316097472223303368.D-20
            XI(17) = -46690290475095840454.D-20
            XI(18) = -40868648199071672992.D-20
            XI(19) = -34875588629216073816.D-20
            XI(20) = -28736248735545557674.D-20
            XI(21) = -22476379039468906122.D-20
            XI(22) = -16122235606889171806.D-20
            XI(23) = -97004699209462698930.D-21
            XI(24) = -32380170962869362033.D-21
            WI(1) = 31533460523058386327.D-22
            WI(2) = 73275539012762621024.D-22
            WI(3) = 11477234579234539490.D-21
            WI(4) = 15579315722943848728.D-21
            WI(5) = 19616160457355527814.D-21
            WI(6) = 23570760839324379141.D-21
            WI(7) = 27426509708356948200.D-21
            WI(8) = 31167227832798088902.D-21
            WI(9) = 34777222564770438893.D-21
            WI(10) = 38241351065830706317.D-21
            WI(11) = 41545082943464749214.D-21
            WI(12) = 44674560856694280419.D-21
            WI(13) = 47616658492490474826.D-21
            WI(14) = 50359035553854474958.D-21
            WI(15) = 52890189485193667096.D-21
            WI(16) = 55199503699984162868.D-21
            WI(17) = 57277292100403215705.D-21
            WI(18) = 59114839698395635746.D-21
            WI(19) = 60704439165893880053.D-21
            WI(20) = 62039423159892663904.D-21
            WI(21) = 63114192286254025657.D-21
            WI(22) = 63924238584648186624.D-21
            WI(23) = 64466164435950082207.D-21
            WI(24) = 64737696812683922503.D-21
        else if (N == 52) then
            XI(1) = -99895111110395027809.D-20
            XI(2) = -99447759092921602925.D-20
            XI(3) = -98644619565154984065.D-20
            XI(4) = -97488388422174450314.D-20
            XI(5) = -95983182693308655253.D-20
            XI(6) = -94134385364135905684.D-20
            XI(7) = -91948612891642453989.D-20
            XI(8) = -89433689053449532252.D-20
            XI(9) = -86598616284606758524.D-20
            XI(10) = -83453543232673453496.D-20
            XI(11) = -80009728343046832433.D-20
            XI(12) = -76279499519374496028.D-20
            XI(13) = -72276209974998319368.D-20
            XI(14) = -68014190422716770209.D-20
            XI(15) = -63508697769524592430.D-20
            XI(16) = -58775860497957906990.D-20
            XI(17) = -53832620928582743838.D-20
            XI(18) = -48696674569809607778.D-20
            XI(19) = -43386406771876167031.D-20
            XI(20) = -37920826911609366925.D-20
            XI(21) = -32319500343480782550.D-20
            XI(22) = -26602478360500182747.D-20
            XI(23) = -20790226415636605969.D-20
            XI(24) = -14903550860694918049.D-20
            XI(25) = -89635244648900565489.D-21
            XI(26) = -29914109797338766044.D-21
            WI(1) = 26913169500471111189.D-22
            WI(2) = 62555239629732768999.D-22
            WI(3) = 98026345794627520617.D-22
            WI(4) = 13315114982340960657.D-21
            WI(5) = 16780023396300735678.D-21
            WI(6) = 20184891507980792203.D-21
            WI(7) = 23517513553984461590.D-21
            WI(8) = 26765953746504013449.D-21
            WI(9) = 29918581147143946641.D-21
            WI(10) = 32964109089718797915.D-21
            WI(11) = 35891634835097232942.D-21
            WI(12) = 38690678310423978985.D-21
            WI(13) = 41351219500560271679.D-21
            WI(14) = 43863734259000407995.D-21
            WI(15) = 46219228372784793508.D-21
            WI(16) = 48409269744074896854.D-21
            WI(17) = 50426018566342377218.D-21
            WI(18) = 52262255383906993034.D-21
            WI(19) = 53911406932757264751.D-21
            WI(20) = 55367569669302652549.D-21
            WI(21) = 56625530902368597191.D-21
            WI(22) = 57680787452526827654.D-21
            WI(23) = 58529561771813868550.D-21
            WI(24) = 59168815466042970369.D-21
            WI(25) = 59596260171248158258.D-21
            WI(26) = 59810365745291860248.D-21
        else if (N == 56) then
            XI(1) = -99909434380146558435.D-20
            XI(2) = -99523122608106974722.D-20
            XI(3) = -98829371554016151109.D-20
            XI(4) = -97830170914025638338.D-20
            XI(5) = -96528590190549018363.D-20
            XI(6) = -94928647956196263565.D-20
            XI(7) = -93035288024749630055.D-20
            XI(8) = -90854362042065549085.D-20
            XI(9) = -88392610832782754079.D-20
            XI(10) = -85657643376274863540.D-20
            XI(11) = -82657913214288165167.D-20
            XI(12) = -79402692289386649803.D-20
            XI(13) = -75902042270512890220.D-20
            XI(14) = -72166783445018808352.D-20
            XI(15) = -68208461269447045550.D-20
            XI(16) = -64039310680700689427.D-20
            XI(17) = -59672218277066332010.D-20
            XI(18) = -55120682485553461875.D-20
            XI(19) = -50398771838438171420.D-20
            XI(20) = -45521081487845957895.D-20
            XI(21) = -40502688092709127812.D-20
            XI(22) = -35359103217495452097.D-20
            XI(23) = -30106225386722066905.D-20
            XI(24) = -24760290943433720397.D-20
            XI(25) = -19337823863527525824.D-20
            XI(26) = -13855584681037624201.D-20
            XI(27) = -83305186822435374440.D-21
            XI(28) = -27797035287275437094.D-21
            WI(1) = 23238553757732154976.D-22
            WI(2) = 54025222460153377612.D-22
            WI(3) = 84690631633078876618.D-22
            WI(4) = 11509824340383382175.D-21
            WI(5) = 14515089278021471807.D-21
            WI(6) = 17475512911400946507.D-21
            WI(7) = 20381929882402572635.D-21
            WI(8) = 23225351562565316937.D-21
            WI(9) = 25996987058391952192.D-21
            WI(10) = 28688268473822741730.D-21
            WI(11) = 31290876747310447868.D-21
            WI(12) = 33796767115611761295.D-21
            WI(13) = 36198193872315186036.D-21
            WI(14) = 38487734259247662487.D-21
            WI(15) = 40658311384744517880.D-21
            WI(16) = 42703216084667086511.D-21
            WI(17) = 44616127652692283213.D-21
            WI(18) = 46391133373001896762.D-21
            WI(19) = 48022746793600258121.D-21
            WI(20) = 49505924683047578920.D-21
            WI(21) = 50836082617798480560.D-21
            WI(22) = 52009109151741399843.D-21
            WI(23) = 53021378524010763968.D-21
            WI(24) = 53869761865714485709.D-21
            WI(25) = 54551636870889421062.D-21
            WI(26) = 55064895901762425796.D-21
            WI(27) = 55407952503245123218.D-21
            WI(28) = 55579746306514395846.D-21
        else if (N == 60) then
            XI(1) = -99921012322743602204.D-20
            XI(2) = -99584052511883817387.D-20
            XI(3) = -98978789522222171737.D-20
            XI(4) = -98106720175259818562.D-20
            XI(5) = -96970178876505273373.D-20
            XI(6) = -95572225583999610740.D-20
            XI(7) = -93916627611642324950.D-20
            XI(8) = -92007847617762755286.D-20
            XI(9) = -89851031081004594194.D-20
            XI(10) = -87451992264689831513.D-20
            XI(11) = -84817198478592963249.D-20
            XI(12) = -81953752616214575937.D-20
            XI(13) = -78869373993226405457.D-20
            XI(14) = -75572377530658568687.D-20
            XI(15) = -72071651335573039944.D-20
            XI(16) = -68376632738135543722.D-20
            XI(17) = -64497282848947706781.D-20
            XI(18) = -60444059704851036344.D-20
            XI(19) = -56227890075394453918.D-20
            XI(20) = -51860140005856974742.D-20
            XI(21) = -47352584176170711111.D-20
            XI(22) = -42717374158307838931.D-20
            XI(23) = -37967005657679797715.D-20
            XI(24) = -33114284826844819425.D-20
            XI(25) = -28172293742326169169.D-20
            XI(26) = -23154355137602933801.D-20
            XI(27) = -18073996487342541724.D-20
            XI(28) = -12944913539694500315.D-20
            XI(29) = -77809333949536569419.D-21
            XI(30) = -25959772301247798589.D-21
            WI(1) = 20268119688737582862.D-22
            WI(2) = 47127299269535689099.D-22
            WI(3) = 73899311633454553761.D-22
            WI(4) = 10047557182287984515.D-21
            WI(5) = 12678166476815960077.D-21
            WI(6) = 15274618596784799378.D-21
            WI(7) = 17829901014207720191.D-21
            WI(8) = 20337120729457286785.D-21
            WI(9) = 22789516943997819862.D-21
            WI(10) = 25180477621521248381.D-21
            WI(11) = 27503556749924791635.D-21
            WI(12) = 29752491500788945241.D-21
            WI(13) = 31921219019296328949.D-21
            WI(14) = 34003892724946422835.D-21
            WI(15) = 35994898051084503067.D-21
            WI(16) = 37888867569243444031.D-21
            WI(17) = 39680695452380799470.D-21
            WI(18) = 41365551235584755613.D-21
            WI(19) = 42938892835935641954.D-21
            WI(20) = 44396478795787113328.D-21
            WI(21) = 45734379716114486647.D-21
            WI(22) = 46948988848912204847.D-21
            WI(23) = 48037031819971180964.D-21
            WI(24) = 48995575455756835389.D-21
            WI(25) = 49822035690550181011.D-21
            WI(26) = 50514184532509374598.D-21
            WI(27) = 51070156069855627405.D-21
            WI(28) = 51488451500980933995.D-21
            WI(29) = 51767943174910187544.D-21
            WI(30) = 51907877631220639733.D-21
        else if (N == 64) then
            XI(1) = -99930504173577213947.D-20
            XI(2) = -99634011677195527930.D-20
            XI(3) = -99101337147674432089.D-20
            XI(4) = -98333625388462595699.D-20
            XI(5) = -97332682778991096379.D-20
            XI(6) = -96100879965205371896.D-20
            XI(7) = -94641137485840281604.D-20
            XI(8) = -92956917213193957580.D-20
            XI(9) = -91052213707850280575.D-20
            XI(10) = -88931544599511410585.D-20
            XI(11) = -86599939815409281976.D-20
            XI(12) = -84062929625258036275.D-20
            XI(13) = -81326531512279755974.D-20
            XI(14) = -78397235894334140761.D-20
            XI(15) = -75281990726053189661.D-20
            XI(16) = -71988185017161082685.D-20
            XI(17) = -68523631305423324256.D-20
            XI(18) = -64896547125465733986.D-20
            XI(19) = -61115535517239325025.D-20
            XI(20) = -57189564620263403428.D-20
            XI(21) = -53127946401989454566.D-20
            XI(22) = -48940314570705295748.D-20
            XI(23) = -44636601725346408798.D-20
            XI(24) = -40227015796399160370.D-20
            XI(25) = -35722015833766811595.D-20
            XI(26) = -31132287199021095616.D-20
            XI(27) = -26468716220876741637.D-20
            XI(28) = -21742364374000708415.D-20
            XI(29) = -16964442042399281804.D-20
            XI(30) = -12146281929612055447.D-20
            XI(31) = -72993121787799039450.D-21
            XI(32) = -24350292663424432509.D-21
            WI(1) = 17832807216964305761.D-22
            WI(2) = 41470332605624720784.D-22
            WI(3) = 65044579689783580339.D-22
            WI(4) = 88467598263639498816.D-22
            WI(5) = 11168139460131126446.D-21
            WI(6) = 13463047896718642240.D-21
            WI(7) = 15726030476024719270.D-21
            WI(8) = 17951715775697343343.D-21
            WI(9) = 20134823153530209297.D-21
            WI(10) = 22270173808383254264.D-21
            WI(11) = 24352702568710873323.D-21
            WI(12) = 26377469715054658680.D-21
            WI(13) = 28339672614259483226.D-21
            WI(14) = 30234657072402478867.D-21
            WI(15) = 32057928354851553585.D-21
            WI(16) = 33805161837141609392.D-21
            WI(17) = 35472213256882383811.D-21
            WI(18) = 37055128540240046040.D-21
            WI(19) = 38550153178615629129.D-21
            WI(20) = 39953741132720341387.D-21
            WI(21) = 41262563242623528610.D-21
            WI(22) = 42473515123653589007.D-21
            WI(23) = 43583724529323453377.D-21
            WI(24) = 44590558163756563060.D-21
            WI(25) = 45491627927418144480.D-21
            WI(26) = 46284796581314417296.D-21
            WI(27) = 46968182816210017325.D-21
            WI(28) = 47540165714830308662.D-21
            WI(29) = 47999388596458307728.D-21
            WI(30) = 48344762234802957170.D-21
            WI(31) = 48575467441503426935.D-21
            WI(32) = 48690957009139720383.D-21
        else if (N == 72) then
            XI(1) = -99944993445296262425.D-20
            XI(2) = -99710287164272906797.D-20
            XI(3) = -99288495101680195780.D-20
            XI(4) = -98680315237583044393.D-20
            XI(5) = -97886877855723380415.D-20
            XI(6) = -96909669799878047771.D-20
            XI(7) = -95750524757769826743.D-20
            XI(8) = -94411618527253794342.D-20
            XI(9) = -92895464588091801987.D-20
            XI(10) = -91204909268867146184.D-20
            XI(11) = -89343126358809125431.D-20
            XI(12) = -87313611129877890577.D-20
            XI(13) = -85120173765443785094.D-20
            XI(14) = -82766932202275431476.D-20
            XI(15) = -80258304396929185201.D-20
            XI(16) = -77599000029998256214.D-20
            XI(17) = -74794011663283239067.D-20
            XI(18) = -71848605366223496136.D-20
            XI(19) = -68768310829046780749.D-20
            XI(20) = -65558910981120105615.D-20
            XI(21) = -62226431133946819101.D-20
            XI(22) = -58777127669165559920.D-20
            XI(23) = -55217476292771439740.D-20
            XI(24) = -51554159877600352386.D-20
            XI(25) = -47794055916894045718.D-20
            XI(26) = -43944223612496073208.D-20
            XI(27) = -40011890621916161270.D-20
            XI(28) = -36004439489141923754.D-20
            XI(29) = -31929393784671217133.D-20
            XI(30) = -27794403980784760698.D-20
            XI(31) = -23607233088575992506.D-20
            XI(32) = -19375742083702605980.D-20
            XI(33) = -15107875148221003033.D-20
            XI(34) = -10811644756210281481.D-20
            XI(35) = -64951166311857114075.D-21
            XI(36) = -21663946035424044670.D-21
            WI(1) = 14115163939753270721.D-22
            WI(2) = 32831697746674940448.D-22
            WI(3) = 51514360187894788637.D-22
            WI(4) = 70102723218632486081.D-22
            WI(5) = 88559960737053994464.D-22
            WI(6) = 10685108165352501353.D-21
            WI(7) = 12494165619873090060.D-21
            WI(8) = 14279769054554032189.D-21
            WI(9) = 16038564950285132556.D-21
            WI(10) = 17767250789200705073.D-21
            WI(11) = 19462580863294276956.D-21
            WI(12) = 21121372216440556069.D-21
            WI(13) = 22740510555035754038.D-21
            WI(14) = 24316956064419165043.D-21
            WI(15) = 25847749100655890089.D-21
            WI(16) = 27330015738950934497.D-21
            WI(17) = 28760973164701761061.D-21
            WI(18) = 30137934895375479293.D-21
            WI(19) = 31458315822561813978.D-21
            WI(20) = 32719637064293846704.D-21
            WI(21) = 33919530618286059497.D-21
            WI(22) = 35055743807217870434.D-21
            WI(23) = 36126143507637992986.D-21
            WI(24) = 37128720154502899461.D-21
            WI(25) = 38061591513802163834.D-21
            WI(26) = 38923006216169663800.D-21
            WI(27) = 39711347044834901782.D-21
            WI(28) = 40425133971733970043.D-21
            WI(29) = 41063026936075061102.D-21
            WI(30) = 41623828360138598208.D-21
            WI(31) = 42106485397586464147.D-21
            WI(32) = 42510091910057720078.D-21
            WI(33) = 42833890168338813667.D-21
            WI(34) = 43077272274913699745.D-21
            WI(35) = 43239781305222617485.D-21
            WI(36) = 43321112165486537076.D-21
        else if (N == 80) then
            XI(1) = -99955382265162826808.D-20
            XI(2) = -99764986439820796834.D-20
            XI(3) = -99422754096571174921.D-20
            XI(4) = -98929130249972665638.D-20
            XI(5) = -98284857273860426984.D-20
            XI(6) = -97490914058572041817.D-20
            XI(7) = -96548508904379327573.D-20
            XI(8) = -95459076634363795066.D-20
            XI(9) = -94224276130985166430.D-20
            XI(10) = -92845987717243869507.D-20
            XI(11) = -91326310257175780181.D-20
            XI(12) = -89667557943877142823.D-20
            XI(13) = -87872256767821427668.D-20
            XI(14) = -85943140666311134922.D-20
            XI(15) = -83883147358025523066.D-20
            XI(16) = -81695413868146347481.D-20
            XI(17) = -79383271750460544325.D-20
            XI(18) = -76950242013504137461.D-20
            XI(19) = -74400029758359727224.D-20
            XI(20) = -71736518536209988023.D-20
            XI(21) = -68963764434202760078.D-20
            XI(22) = -66085989898611980174.D-20
            XI(23) = -63107577304687196625.D-20
            XI(24) = -60033062282975174315.D-20
            XI(25) = -56867126812270978473.D-20
            XI(26) = -53614592089713193202.D-20
            XI(27) = -50280411188878498759.D-20
            XI(28) = -46869661517054447704.D-20
            XI(29) = -43387537083175609306.D-20
            XI(30) = -39839340588196922702.D-20
            XI(31) = -36230475349948731562.D-20
            XI(32) = -32566437074770191462.D-20
            XI(33) = -28852805488451185311.D-20
            XI(34) = -25095235839227212049.D-20
            XI(35) = -21299450285766613257.D-20
            XI(36) = -17471229183264681256.D-20
            XI(37) = -13616402280914388656.D-20
            XI(38) = -97408398441584599063.D-21
            XI(39) = -58504437152420668629.D-21
            XI(40) = -19511383256793997654.D-21
            WI(1) = 11449500037252453354.D-22
            WI(2) = 26635335911449182005.D-22
            WI(3) = 41803131248325406933.D-22
            WI(4) = 56909224518929595161.D-22
            WI(5) = 71929047688527908617.D-22
            WI(6) = 86839452691609565673.D-22
            WI(7) = 10161766041219293248.D-21
            WI(8) = 11624114120416127748.D-21
            WI(9) = 13068761592477329935.D-21
            WI(10) = 14493508040524581339.D-21
            WI(11) = 15896183583752925452.D-21
            WI(12) = 17274652056269915690.D-21
            WI(13) = 18626814208301626849.D-21
            WI(14) = 19950610878140102324.D-21
            WI(15) = 21244026115781505420.D-21
            WI(16) = 22505090246332526056.D-21
            WI(17) = 23731882865930076655.D-21
            WI(18) = 24922535764115501728.D-21
            WI(19) = 26075235767565117065.D-21
            WI(20) = 27188227500486381444.D-21
            WI(21) = 28259816057276862255.D-21
            WI(22) = 29288369583267847694.D-21
            WI(23) = 30272321759557980659.D-21
            WI(24) = 31210174188114701643.D-21
            WI(25) = 32100498673487773148.D-21
            WI(26) = 32941939397645401383.D-21
            WI(27) = 33733214984611522817.D-21
            WI(28) = 34473120451753928794.D-21
            WI(29) = 35160529044747593496.D-21
            WI(30) = 35794393953416054603.D-21
            WI(31) = 36373749905835978044.D-21
            WI(32) = 36897714638276008839.D-21
            WI(33) = 37365490238730490027.D-21
            WI(34) = 37776364362001397490.D-21
            WI(35) = 38129711314477638344.D-21
            WI(36) = 38424993006959423185.D-21
            WI(37) = 38661759774076463327.D-21
            WI(38) = 38839651059051968932.D-21
            WI(39) = 38958395962769531199.D-21
            WI(40) = 39017813656306654811.D-21
        else if (N == 88) then
            XI(1) = -99963083606662645367.D-20
            XI(2) = -99805540804996391115.D-20
            XI(3) = -99522317584752983826.D-20
            XI(4) = -99113707813467295361.D-20
            XI(5) = -98580218616336573201.D-20
            XI(6) = -97922520312124670724.D-20
            XI(7) = -97141441021462325102.D-20
            XI(8) = -96237964604708261699.D-20
            XI(9) = -95213229349011862112.D-20
            XI(10) = -94068526328823504199.D-20
            XI(11) = -92805297840849639853.D-20
            XI(12) = -91425135510152027922.D-20
            XI(13) = -89929778323076240540.D-20
            XI(14) = -88321110405889090804.D-20
            XI(15) = -86601158661395039865.D-20
            XI(16) = -84772090209816777675.D-20
            XI(17) = -82836209659047569960.D-20
            XI(18) = -80795956199923475240.D-20
            XI(19) = -78653900532828988584.D-20
            XI(20) = -76412741628420688124.D-20
            XI(21) = -74075303326860783021.D-20
            XI(22) = -71644530779725048283.D-20
            XI(23) = -69123486739088043741.D-20
            XI(24) = -66515347698445317191.D-20
            XI(25) = -63823399890332645027.D-20
            XI(26) = -61051035145681898624.D-20
            XI(27) = -58201746620128735758.D-20
            XI(28) = -55279124392655547280.D-20
            XI(29) = -52286850942114375373.D-20
            XI(30) = -49228696507328672289.D-20
            XI(31) = -46108514336619651061.D-20
            XI(32) = -42930235832742437021.D-20
            XI(33) = -39697865599349105025.D-20
            XI(34) = -36415476395219828782.D-20
            XI(35) = -33087204002619627511.D-20
            XI(36) = -29717242016246430559.D-20
            XI(37) = -26309836559336260023.D-20
            XI(38) = -22869280933583131390.D-20
            XI(39) = -19399910209614678858.D-20
            XI(40) = -15906095764839421594.D-20
            XI(41) = -12392239775547906165.D-20
            XI(42) = -88627696702076058805.D-21
            XI(43) = -53221325509403576522.D-21
            XI(44) = -17747895902112098289.D-21
            WI(1) = 94733513612529545899.D-23
            WI(2) = 22040586988692204378.D-22
            WI(3) = 34598684102906050894.D-22
            WI(4) = 47114803663497454871.D-22
            WI(5) = 59571849982300390278.D-22
            WI(6) = 71953989055383030738.D-22
            WI(7) = 84245466103046449651.D-22
            WI(8) = 96430834163440844925.D-22
            WI(9) = 10849469467224321456.D-21
            WI(10) = 12042186598227368180.D-21
            WI(11) = 13219730188112258690.D-21
            WI(12) = 14380617607801932514.D-21
            WI(13) = 15523385539474914063.D-21
            WI(14) = 16646594210848201219.D-21
            WI(15) = 17748828358570096473.D-21
            WI(16) = 18828699173754781504.D-21
            WI(17) = 19884846011190532102.D-21
            WI(18) = 20915938130610818211.D-21
            WI(19) = 21920676359974122618.D-21
            WI(20) = 22897794734780677291.D-21
            WI(21) = 23846062091859655017.D-21
            WI(22) = 24764283620768778908.D-21
            WI(23) = 25651302368961951718.D-21
            WI(24) = 26506000699434738764.D-21
            WI(25) = 27327301698855331284.D-21
            WI(26) = 28114170534408613493.D-21
            WI(27) = 28865615757635429249.D-21
            WI(28) = 29580690553619349115.D-21
            WI(29) = 30258493933943525335.D-21
            WI(30) = 30898171871912197634.D-21
            WI(31) = 31498918378604892320.D-21
            WI(32) = 32059976518406388069.D-21
            WI(33) = 32580639362732108686.D-21
            WI(34) = 33060250880746700145.D-21
            WI(35) = 33498206765953092528.D-21
            WI(36) = 33893955197610259240.D-21
            WI(37) = 34246997536020078737.D-21
            WI(38) = 34556888950807084135.D-21
            WI(39) = 34823238981399354993.D-21
            WI(40) = 35045712029004261397.D-21
            WI(41) = 35224027779459108533.D-21
            WI(42) = 35357961556423843794.D-21
            WI(43) = 35447344604470769706.D-21
            WI(44) = 35492064301714545296.D-21
        else if (N == 96) then
            XI(1) = -99968950458161677884.D-20
            XI(2) = -99836436150026625442.D-20
            XI(3) = -99598185401837168989.D-20
            XI(4) = -99254388318724897532.D-20
            XI(5) = -98805415112002250496.D-20
            XI(6) = -98251725228772058939.D-20
            XI(7) = -97593918962884178795.D-20
            XI(8) = -96832679623798884729.D-20
            XI(9) = -95968830381449102833.D-20
            XI(10) = -95003272972082660508.D-20
            XI(11) = -93937032620255700317.D-20
            XI(12) = -92771245540840037901.D-20
            XI(13) = -91507142091195428283.D-20
            XI(14) = -90146063448296036815.D-20
            XI(15) = -88689451757413437798.D-20
            XI(16) = -87138850590440317660.D-20
            XI(17) = -85495903346571022516.D-20
            XI(18) = -83762351122975346770.D-20
            XI(19) = -81940031074042640835.D-20
            XI(20) = -80030874413847743608.D-20
            XI(21) = -78036904386796125576.D-20
            XI(22) = -75960234117661479314.D-20
            XI(23) = -73803064374440921909.D-20
            XI(24) = -71567681234896650925.D-20
            XI(25) = -69256453664217173081.D-20
            XI(26) = -66871831004391619210.D-20
            XI(27) = -64416340378496711088.D-20
            XI(28) = -61892584012546857021.D-20
            XI(29) = -59303236477757208075.D-20
            XI(30) = -56651041856139716841.D-20
            XI(31) = -53938810832435743623.D-20
            XI(32) = -51169417715466767359.D-20
            XI(33) = -48345797392059635977.D-20
            XI(34) = -45470942216774300864.D-20
            XI(35) = -42547898840730054536.D-20
            XI(36) = -39579764982890860329.D-20
            XI(37) = -36569686147231363503.D-20
            XI(38) = -33520852289262542262.D-20
            XI(39) = -30436494435449635302.D-20
            XI(40) = -27319881259104914149.D-20
            XI(41) = -24174315616384001233.D-20
            XI(42) = -21003131046056720360.D-20
            XI(43) = -17809688236761860276.D-20
            XI(44) = -14597371465489694199.D-20
            XI(45) = -11369585011066592091.D-20
            XI(46) = -81297495464425558994.D-21
            XI(47) = -48812985136049731112.D-21
            XI(48) = -16276744849602969579.D-21
            WI(1) = 79647473413013308824.D-23
            WI(2) = 18545712028943610772.D-22
            WI(3) = 29096372771815468255.D-22
            WI(4) = 39656312639089457047.D-22
            WI(5) = 50134769158268190741.D-22
            WI(6) = 60590278611871775766.D-22
            WI(7) = 70960176225543354086.D-22
            WI(8) = 81274503287260803750.D-22
            WI(9) = 91484770246534778594.D-22
            WI(10) = 10160502506419551876.D-21
            WI(11) = 11162229814602627176.D-21
            WI(12) = 12151660787661414925.D-21
            WI(13) = 13128256444197094780.D-21
            WI(14) = 14090951969013093016.D-21
            WI(15) = 15038721680197214168.D-21
            WI(16) = 15970562935363605962.D-21
            WI(17) = 16885479435374326053.D-21
            WI(18) = 17782502280253696814.D-21
            WI(19) = 18660679624254468586.D-21
            WI(20) = 19519081141356853600.D-21
            WI(21) = 20356797152792231598.D-21
            WI(22) = 21172939892395877365.D-21
            WI(23) = 21966644438730298358.D-21
            WI(24) = 22737069658335806966.D-21
            WI(25) = 23483399085928249636.D-21
            WI(26) = 24204841792364550405.D-21
            WI(27) = 24900633222483566277.D-21
            WI(28) = 25570036005349360697.D-21
            WI(29) = 26212340735672413804.D-21
            WI(30) = 26826866725591762076.D-21
            WI(31) = 27412962726029242828.D-21
            WI(32) = 27970007616848334438.D-21
            WI(33) = 28497411065085385645.D-21
            WI(34) = 28994614150555236543.D-21
            WI(35) = 29461089958167905970.D-21
            WI(36) = 29896344136328385984.D-21
            WI(37) = 30299915420827593794.D-21
            WI(38) = 30671376123669149014.D-21
            WI(39) = 31010332586313837423.D-21
            WI(40) = 31316425596861355813.D-21
            WI(41) = 31589330770727168558.D-21
            WI(42) = 31828758894411006535.D-21
            WI(43) = 32034456231992663218.D-21
            WI(44) = 32206204794030250669.D-21
            WI(45) = 32343822568575928429.D-21
            WI(46) = 32447163714064269364.D-21
            WI(47) = 32516118713868835987.D-21
            WI(48) = 32550614492363166242.D-21
        else if (N == 104) then
            XI(1) = -99974025462451620190.D-20
            XI(2) = -99861349070638500381.D-20
            XI(3) = -99660327736123583184.D-20
            XI(4) = -99361508866476523597.D-20
            XI(5) = -98976332950492964529.D-20
            XI(6) = -98504935190347884894.D-20
            XI(7) = -97950469321598158611.D-20
            XI(8) = -97301711615545710014.D-20
            XI(9) = -96557677384453027357.D-20
            XI(10) = -95730762749516840775.D-20
            XI(11) = -94823341083979390414.D-20
            XI(12) = -93824333288826003337.D-20
            XI(13) = -92742344950081829375.D-20
            XI(14) = -91576253747884907173.D-20
            XI(15) = -90327594749081629450.D-20
            XI(16) = -88997023646288045179.D-20
            XI(17) = -87586186756351095122.D-20
            XI(18) = -86096161127797694571.D-20
            XI(19) = -84528340899234132869.D-20
            XI(20) = -82884124840218026839.D-20
            XI(21) = -81165006660051550143.D-20
            XI(22) = -79372537644043611522.D-20
            XI(23) = -77508338027046977385.D-20
            XI(24) = -75574092471832298179.D-20
            XI(25) = -73571549017801969428.D-20
            XI(26) = -71502517397363087234.D-20
            XI(27) = -69368867439397065602.D-20
            XI(28) = -67172527368333210088.D-20
            XI(29) = -64915482063111896041.D-20
            XI(30) = -62599771263251515327.D-20
            XI(31) = -60227487725540047159.D-20
            XI(32) = -57800775332757748445.D-20
            XI(33) = -55321827156203442307.D-20
            XI(34) = -52792883473767725346.D-20
            XI(35) = -50216229745345025782.D-20
            XI(36) = -47594194547413982495.D-20
            XI(37) = -44929147468652664376.D-20
            XI(38) = -42223496968490362523.D-20
            XI(39) = -39479688200531193751.D-20
            XI(40) = -36700200802816507262.D-20
            XI(41) = -33887546656923060512.D-20
            XI(42) = -31044267617922098404.D-20
            XI(43) = -28172933217250806882.D-20
            XI(44) = -25276138340572094235.D-20
            XI(45) = -22356500882721258998.D-20
            XI(46) = -19416659381858811954.D-20
            XI(47) = -16459270634967512819.D-20
            XI(48) = -13487007296848542737.D-20
            XI(49) = -10502555464786646703.D-20
            XI(50) = -75086122510670317989.D-21
            XI(51) = -45078833455377862647.D-21
            XI(52) = -15030805704205808070.D-21
            WI(1) = 80537301359223283137.D-23
            WI(2) = 88612279172572261915.D-23
            WI(3) = 29505344364932164953.D-22
            WI(4) = 29860447460835539614.D-22
            WI(5) = 44717827658964106576.D-22
            WI(6) = 54373252304872340492.D-22
            WI(7) = 63840289565270698417.D-22
            WI(8) = 64909501207992000200.D-22
            WI(9) = 83745523985271410846.D-22
            WI(10) = 85549571290963629793.D-22
            WI(11) = 94157585502386417660.D-22
            WI(12) = 10364141417799710457.D-21
            WI(13) = 11288156174522055729.D-21
            WI(14) = 12085767303575836382.D-21
            WI(15) = 12898236413125625150.D-21
            WI(16) = 13709097817130044944.D-21
            WI(17) = 14507175337544532779.D-21
            WI(18) = 15291431851904123958.D-21
            WI(19) = 16062483287466054902.D-21
            WI(20) = 16819201771329851973.D-21
            WI(21) = 17560584847393953123.D-21
            WI(22) = 18286096131447045103.D-21
            WI(23) = 18995087062019415460.D-21
            WI(24) = 19686910416592402212.D-21
            WI(25) = 20360942228950003547.D-21
            WI(26) = 21016573507412752623.D-21
            WI(27) = 21653211654050478214.D-21
            WI(28) = 22270281336456832151.D-21
            WI(29) = 22867224894903991760.D-21
            WI(30) = 23443502859080308868.D-21
            WI(31) = 23998594434228088202.D-21
            WI(32) = 24531997972222651431.D-21
            WI(33) = 25043231424904261260.D-21
            WI(34) = 25531832779704979249.D-21
            WI(35) = 25997360477175197755.D-21
            WI(36) = 26439393810028101126.D-21
            WI(37) = 26857533303339869813.D-21
            WI(38) = 27251401075562153020.D-21
            WI(39) = 27620641180020444979.D-21
            WI(40) = 27964919926589683664.D-21
            WI(41) = 28283926183256317188.D-21
            WI(42) = 28577371657294272064.D-21
            WI(43) = 28844991155800689526.D-21
            WI(44) = 29086542825355955653.D-21
            WI(45) = 29301808370591421861.D-21
            WI(46) = 29490593251467277733.D-21
            WI(47) = 29652726859082281213.D-21
            WI(48) = 29788062669856454763.D-21
            WI(49) = 29896478377947402691.D-21
            WI(50) = 29977876005780577121.D-21
            WI(51) = 30032181992593600121.D-21
            WI(52) = 30059347260914619701.D-21
        else if (N == 112) then
            XI(1) = -80488595332676199779.D-20
            XI(2) = -10160478627668805531.D-19
            XI(3) = -95439318547860299986.D-20
            XI(4) = -91546298836742760945.D-20
            XI(5) = -10340513562261696462.D-19
            XI(6) = -92662585778661689507.D-20
            XI(7) = -90489298887912825445.D-20
            XI(8) = -10131949463069062515.D-19
            XI(9) = -90464129755666057124.D-20
            XI(10) = -86603614404095298415.D-20
            XI(11) = -91636186250656082289.D-20
            XI(12) = -10121834958204258714.D-19
            XI(13) = -10078277129241111043.D-19
            XI(14) = -87962237189807963118.D-20
            XI(15) = -91677434744694723682.D-20
            XI(16) = -90487658323691933037.D-20
            XI(17) = -89271920604983483215.D-20
            XI(18) = -87962151397116941034.D-20
            XI(19) = -86599374258026102045.D-20
            XI(20) = -85171753748871651951.D-20
            XI(21) = -83675665375815917519.D-20
            XI(22) = -82114054838778695246.D-20
            XI(23) = -80488575638784496812.D-20
            XI(24) = -78800297986553099992.D-20
            XI(25) = -77050563606913539332.D-20
            XI(26) = -75240747787024521741.D-20
            XI(27) = -73372261848570062821.D-20
            XI(28) = -71446562385808640412.D-20
            XI(29) = -69465151030101343512.D-20
            XI(30) = -67429572842563396177.D-20
            XI(31) = -65341415102530449682.D-20
            XI(32) = -63202306093043949251.D-20
            XI(33) = -61013913826811007054.D-20
            XI(34) = -58777944746249047206.D-20
            XI(35) = -56496142392748458796.D-20
            XI(36) = -54170286047113786527.D-20
            XI(37) = -51802189342131382164.D-20
            XI(38) = -49393698848352747124.D-20
            XI(39) = -46946692634193765874.D-20
            XI(40) = -44463078801473008955.D-20
            XI(41) = -41944793997530961580.D-20
            XI(42) = -39393801905090372790.D-20
            XI(43) = -36812091711035275365.D-20
            XI(44) = -34201676555302667702.D-20
            XI(45) = -31564591961096358173.D-20
            XI(46) = -28902894247647038221.D-20
            XI(47) = -26218658926756261480.D-20
            XI(48) = -23513979084374651836.D-20
            XI(49) = -20790963748476333913.D-20
            XI(50) = -18051736244502265770.D-20
            XI(51) = -15298432539654847416.D-20
            XI(52) = -12533199577334872554.D-20
            XI(53) = -97581936030195779185.D-21
            XI(54) = -69755784828872187578.D-21
            XI(55) = -41875240164992552336.D-21
            XI(56) = -13962042448558683275.D-21
            WI(1) = 16569636576973647761.D-21
            WI(2) = 14402201696983724428.D-39
            WI(3) = 34709520255747461200.D-22
            WI(4) = 56778376879010079815.D-21
            WI(5) = 37731161063754411165.D-47
            WI(6) = 45632271760860013289.D-22
            WI(7) = 11520495792117312626.D-21
            WI(8) = 49246746812766292842.D-38
            WI(9) = 14372045777313838419.D-21
            WI(10) = 13870830688358627358.D-21
            WI(11) = -37371950806309312086.D-21
            WI(12) = 18775853756740260905.D-37
            WI(13) = 12261054337390208321.D-34
            WI(14) = 13637759622784054281.D-21
            WI(15) = 20686815045027872828.D-21
            WI(16) = 84327272746907118565.D-22
            WI(17) = 11759939136603779624.D-21
            WI(18) = 13945650387627380053.D-21
            WI(19) = 14124057734648147892.D-21
            WI(20) = 14634623094715024217.D-21
            WI(21) = 15300784134556015037.D-21
            WI(22) = 15940621199485025417.D-21
            WI(23) = 16571845899450453236.D-21
            WI(24) = 17192076898900235720.D-21
            WI(25) = 17800021943506593268.D-21
            WI(26) = 18393902644953860296.D-21
            WI(27) = 18973392242557999400.D-21
            WI(28) = 19538093130310137860.D-21
            WI(29) = 20087558367323984780.D-21
            WI(30) = 20621359698898228923.D-21
            WI(31) = 21139081089410710408.D-21
            WI(32) = 21640318863160403311.D-21
            WI(33) = 22124682168949492325.D-21
            WI(34) = 22591793313385707526.D-21
            WI(35) = 23041288057354547180.D-21
            WI(36) = 23472815898299421860.D-21
            WI(37) = 23886040343752650004.D-21
            WI(38) = 24280639173689438499.D-21
            WI(39) = 24656304691784812879.D-21
            WI(40) = 25012743965344867773.D-21
            WI(41) = 25349679053726235917.D-21
            WI(42) = 25666847225065297520.D-21
            WI(43) = 25964001161148146161.D-21
            WI(44) = 26240909150261532990.D-21
            WI(45) = 26497355267874394909.D-21
            WI(46) = 26733139545009064708.D-21
            WI(47) = 26948078124170862883.D-21
            WI(48) = 27142003402714474402.D-21
            WI(49) = 27314764163535311670.D-21
            WI(50) = 27466225692983949724.D-21
            WI(51) = 27596269885911683746.D-21
            WI(52) = 27704795337765294488.D-21
            WI(53) = 27791717423659206498.D-21
            WI(54) = 27856968364363379146.D-21
            WI(55) = 27900497279155473642.D-21
            WI(56) = 27922270225496082396.D-21
        else
            write (6,FMT=*) 'CASE N=', N, ' IS NOT PROVIDED IN GAULEG'
            STOP 'GAULEG'
        endif

        do I = N, (N+3)/2,-1
            XI(I) = -XI(N+1-I)
            WI(I) = WI(N+1-I)
        end do

    end subroutine GAULEG
    
    subroutine EMESHT(EZ,DF,NPNT,EBOT,EMU,TK,NPOL,NPNT1,NPNT2,NPNT3,PI,KB,IEMXD)
        ! **********************************************************************
        ! *                                                                    *
        ! * This subroutine provides the energy mesh in array EZ and the       *
        ! * appropriate integration weights in array DF.                       *
        ! *                                                                    *
        ! * Poles of the Fermi function C (Matsubara frequencies) and          *
        ! * a contour in the complex energy are used as described in (????).   *
        ! *                                                                    *
        ! * The contour consists of three straight lines with                  *
        ! * NPNT1, NPNT2, and NPNT3 integration points and is determined by    *
        ! * the input arguments: EBOT, EMU, TK, and NPOL.                      *
        ! *                                                                    *
        ! *            TK   = temperature in K                                 *
        ! *            EMU  = chemical potential in Ry                         *
        ! *            EBOT = bottom of contour in Ry                          *
        ! *            NPOL = number of Matsubara frequencies                  *
        ! *                                                                    *
        ! * The three lines are defined by:                                    *
        ! *                                                                    *
        ! *  1. the line from EBOT to EBOT+2*NPOL*pi*i*k*TK                    *
        ! *              with NPNT1 integration points (Gauss-Legendre rule)   *
        ! *                                                                    *
        ! *  2. the line from EBOT+2*NPOL*pi*i*k*TK to                         *
        ! *                   EMU+(2*NPOL*pi*i-30)*k*TK                        *
        ! *              with NPNT2 integration points (Gauss-Legendre rule)   *
        ! *                                                                    *
        ! *  3. the line from EMU+(2*NPOL*pi*i-30)*k*TK to infinity            *
        ! *              with NPNT3 integration points (Gauss-Fermi-Dirac rule)*
        ! *                                                                    *
        ! *  The total number of integration points is given by:               *
        ! *              NPNT=NPNT1+NPNT2+NPNT3+NPOL                           *
        ! *                                                                    *
        ! *  The integration points and weights on three lines are chosen      *
        ! *  according to Gauss integration rules. Only in third interval      *
        ! *  the Fermi function matters since exp(x) < 10**(-10) for x < -25.  *
        ! *                                                                    *
        ! *  There are two special cases determined by NPOL = 0 and NPOL < 0.  *
        ! *                                                                    *
        ! *  a) NPOL = 0 leads to density-of-states calculations               *
        ! *  with constant integration weights and equally distributed points  *
        ! *  between EBOT - pi*i*k*TK and EMU - pi*i*k*TK.                     *
        ! *                                                                    *
        ! *  The total number of integration points is given by:               *
        ! *              NPNT=NPNT2                                            *
        ! *                                                                    *
        ! *  b) NPOL < 0 is meant for calculations where the Fermi-Dirac       *
        ! *  function is replaced by a step function with step at EMU. When    *
        ! *  this option is used no poles of the Fermi-Dirac function are used *
        ! *  and the contour consists of the three straight lines:             *
        ! *                                                                    *
        ! *  1. the line from EBOT to EBOT-2*NPOL*pi*i*k*TK                    *
        ! *              with NPNT1 integration points (Gauss-Legendre rule)   *
        ! *                                                                    *
        ! *  2. the line from EBOT-2*NPOL*pi*i*k*TK to EMU-2*NPOL*pi*i*k*TK    *
        ! *              with NPNT2 integration points (Gauss-Legendre rule)   *
        ! *                                                                    *
        ! *  3. the line from EMU-2*NPOL*pi*i*k*TK to EMU                      *
        ! *              with NPNT3 integration points (Gauss-Legendre rule)   *
        ! *                                                                    *
        ! *  The total number of integration points is given by:               *
        ! *              NPNT=NPNT1+NPNT2+NPNT3                                *
        ! *                                                                    *
        ! **********************************************************************
        !     
            integer :: NPNT, NPNT1, NPNT2, NPNT3, NPOL, I, IEMXD
            real*8 :: EBOT, EMU, TK, ER, ETK, KB, PI, WI(128), XI(128)
            complex*16 :: DE, DF(IEMXD), EZ(IEMXD)
            !DATA RYD/13.6058D0/
        !
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        !    write (6,'(5X,A,F12.6," (Ry)",8X,A,F12.6," (Ry)")') 'E min = ', EBOT, 'Fermi energy = ', EFERMI
        !    write (6,'(5X,A,F12.6," (Ry)",8X,A,F12.6," (K )",/,5X,62(1H-))') 'E max = ', EMU, 'Temperature  = ', TK
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
        !
            ETK = PI*KB*TK ! Broadening
        ! ======================================================================
            if (NPOL == 0) then
                DE = (EMU-EBOT)
                if (NPNT2 > 1) then
                    DE = DE/(NPNT2-1)
                else
                    DE = DCMPLX(1.0D0,0.0D0)
                endif
                NPNT = 0
                do I = 1, NPNT2
                    NPNT = NPNT + 1
                    ER = EBOT + (I-1)*REAL(DE) ! Before changing it was "ER = EBOT + (I-1)*DE"
                    EZ(NPNT) = DCMPLX(ER,ETK)
                    DF(NPNT) = DE
                end do

                !write (6,FMT=9000) NPNT,ETK,ETK*RYD
            ! ------------------------------------------------------------- NPOL > 0
            else if (NPOL > 0) then
                call GAULEG(XI,WI,NPNT1)
                DE = NPOL*DCMPLX(0.0D0,ETK)
                NPNT = 0

                do I = 1, NPNT1
                    NPNT = NPNT + 1
                    EZ(NPNT) = XI(I)*DE + DE + EBOT
                    DF(NPNT) = WI(I)*DE
                end do

                call GAULEG(XI,WI,NPNT2)
                DE = (EMU-30*KB*TK-EBOT)*0.5D0
                do I = 1, NPNT2
                    NPNT = NPNT + 1
                    EZ(NPNT) = XI(I)*DE + DE + EBOT + 2*NPOL*DCMPLX(0.0D0,ETK)
                    DF(NPNT) = WI(I)*DE
                end do

                call GAUFD(XI,WI,NPNT3)
                DE = 30*KB*TK
                do I = 1, NPNT3
                    NPNT = NPNT + 1
                    EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
                    DF(NPNT) = WI(I)*DE
                end do

                do I = NPOL,1,-1
                    NPNT = NPNT + 1
                    EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
                    DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
                end do

                !write(6,9090) NPNT,NPOL,NPNT1,NPNT2,NPNT3
        ! ------------------------------------------------------------- NPOL < 0
            else
                if (NPNT1 > 0) call GAULEG(XI,WI,NPNT1)
                DE = -NPOL*DCMPLX(0.0D0,ETK)
                NPNT = 0
                do I = 1, NPNT1
                    NPNT = NPNT + 1
                    EZ(NPNT) = XI(I)*DE + DE + EBOT
                    DF(NPNT) = WI(I)*DE
                end do

                call GAULEG(XI,WI,NPNT2)
                DE = (EMU-EBOT)*0.5D0

                do I = 1, NPNT2
                    NPNT = NPNT + 1
                    EZ(NPNT) = XI(I)*DE + DE + EBOT - 2*NPOL*DCMPLX(0.0D0,ETK)
                    DF(NPNT) = WI(I)*DE
                end do

                if (NPNT3 > 0) call GAULEG(XI,WI,NPNT3)
                DE = -NPOL*DCMPLX(0.0D0,ETK)
                do I = NPNT3,1,-1
                    NPNT = NPNT + 1
                    EZ(NPNT) = XI(I)*DE + DE + EMU
                    DF(NPNT) = -WI(I)*DE
                end do
                
                !write(6,9091) NPNT,-NPOL,NPNT1,NPNT2,NPNT3

            endif
        ! ======================================================================
        !write(6,*)

        !9000 FORMAT (5X,'Density-of-States calculation',/,5X,'Number of energy points :',I4,4X,'broadening =', &
        !    &        3P,F9.3,' ( mRy )',/,48X,' =',3P,F9.3,' ( meV )')
        !9090 FORMAT (5X,'GF integration rectangular contour ( ImE > 0 )',/,5X,'Number of energy points :',I4,13X, &
        !    &       'poles =',I2,/, 23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
        !9091 FORMAT (5X,'GF integration rectangular contour ( ImE < 0 )',/, 5X,'Number of energy points :',I4,13X,'poles =',I2,/, &
        !    &       23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
     
    end subroutine EMESHT

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

end program GREEN