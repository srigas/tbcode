program BDG_IMP
    implicit none

    integer :: i, j, NUME, NUMIMP, IE, INFO, FTIMO, NUMEDOS, NPNT, NUMIT, MAXNUMIT
    integer, allocatable, dimension(:) :: IPIV
    real*8 :: PI, KB, tempval1, tempval2, epsilon
    real*8, allocatable, dimension(:) :: E0, E0IMP, nuzero, nuzeroIMP, ULCN, ULCNIMP, NEWNU, VSUPCOND, &
    &VSUPCONDIMP, NU, DIFFN, PREVNU, DIFFD
    real*8, allocatable, dimension(:,:) :: BETA, BETAIMP, densityperimpurity
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), ONE, MINUSONE, posham(2,2), EZ
    complex*16, allocatable, dimension(:) :: DELTA, energies, WEIGHTS, EZDOS, NEWDELTA, PREVDELTA
    complex*16, allocatable, dimension(:,:) :: GREEN, HAMIMP, IDENTITY
    character(len = 1) :: selfconfin, finaldosorno, printgimp

    print *, 'Initiating the impurity problem calculations.'

    call CONSTANTS(IdentityPauli,xPauli,yPauli,zPauli,CI,PI,KB) ! Sets some universal constants.

    ! Reads the configuration info from the BdG.f90 output
    open(1, file = 'impconfig.dat', action = 'read')
    read(1,*) NUMIMP
    read(1,*) NUME
    read(1,*) NUMEDOS
    read(1,*) printgimp
    read(1,*) epsilon
    close(1)

    allocate(E0(NUMIMP))
    allocate(E0IMP(NUMIMP))
    allocate(ULCN(NUMIMP))
    allocate(ULCNIMP(NUMIMP))
    allocate(BETA(3,NUMIMP))
    allocate(BETAIMP(3,NUMIMP))
    allocate(VSUPCOND(NUMIMP))
    allocate(VSUPCONDIMP(NUMIMP))
    allocate(nuzero(NUMIMP))
    allocate(nuzeroIMP(NUMIMP))

    allocate(GREEN(4*NUMIMP,4*NUMIMP))
    allocate(HAMIMP(4*NUMIMP,4*NUMIMP))
    allocate(IDENTITY(4*NUMIMP,4*NUMIMP))
    
    allocate(NU(NUMIMP))
    allocate(DELTA(NUMIMP))
    allocate(NEWNU(NUMIMP))
    allocate(NEWDELTA(NUMIMP))
    allocate(PREVNU(NUMIMP))
    allocate(PREVDELTA(NUMIMP))
    allocate(DIFFN(NUMIMP))
    allocate(DIFFD(NUMIMP))

    allocate(energies(NUME))
    allocate(EZDOS(NUMEDOS))
    allocate(WEIGHTS(NUME))
    allocate(densityperimpurity(1+NUMIMP,NUMEDOS))

    open(1, file = 'impconfig.dat', action = 'read')
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    do i = 1, NUMIMP
        read(1,*) NU(i), tempval1, tempval2
        DELTA(i) = dcmplx(tempval1,tempval2)
    end do
    close(1)

    allocate(IPIV(4*NUMIMP))

    finaldosorno = 'n' ! Starting value

    open(1, file = 'energies.dat', action = 'read')
    do i = 1, NUME
        read(1,*) tempval1, tempval2
        energies(i) = dcmplx(tempval1, tempval2)
    end do
    close(1)

    open(1, file = 'fermiweights.dat', action = 'read')
    do i = 1, NUME
        read(1,*) tempval1, tempval2
        WEIGHTS(i) = dcmplx(tempval1, tempval2)
    end do
    close(1)

    open(1, file = 'impatoms.dat', action = 'read')
    do i = 1, NUMIMP
        read (1,*) E0(i), BETA(1,i), BETA(2,i), BETA(3,i), VSUPCOND(i), ULCN(i), nuzero(i)
        read (1,*) E0IMP(i), BETAIMP(1,i), BETAIMP(2,i), BETAIMP(3,i), VSUPCONDIMP(i), ULCNIMP(i), nuzeroIMP(i)
    end do
    close(1)

    ! Initial values
    NEWNU = NU
    NEWDELTA = DELTA
    NUMIT = 1
    DIFFD = 1.0
    DIFFN = 1.0

    print *, 'Please enter the maximum number of self-consistency cyles for the impurity problem.'
    read *, MAXNUMIT

    180 HAMIMP(:,:) = (0.d0,0.d0)

    ! Constructs the ΔH matrix
    do i = 1, NUMIMP

        FTIMO = 4*(i-1)

        posham = (E0IMP(i) - E0(i) + ULCNIMP(i)*(NEWNU(i) - nuzeroIMP(i)) - ULCN(i)*(NU(i) - nuzero(i)))*IdentityPauli -&
        & (BETAIMP(1,i)-BETA(1,i))*xPauli - (BETAIMP(2,i)-BETA(2,i))*yPauli - (BETAIMP(3,i)-BETA(3,i))*zPauli

        HAMIMP(1 + FTIMO,1 + FTIMO) = posham(1,1)
        HAMIMP(1 + FTIMO,2 + FTIMO) = posham(1,2)
        !HAMIMP(1 + FTIMO,3 + FTIMO) = 0.0
        HAMIMP(1 + FTIMO,4 + FTIMO) = NEWDELTA(i) - DELTA(i)

        HAMIMP(2 + FTIMO,1 + FTIMO) = posham(2,1)
        HAMIMP(2 + FTIMO,2 + FTIMO) = posham(2,2)
        HAMIMP(2 + FTIMO,3 + FTIMO) = HAMIMP(1 + FTIMO,4 + FTIMO)
        !HAMIMP(2 + FTIMO,4 + FTIMO) = 0.0

        !HAMIMP(3 + FTIMO,1 + FTIMO) = 0.0
        HAMIMP(3 + FTIMO,2 + FTIMO) = CONJG(HAMIMP(1 + FTIMO,4 + FTIMO))
        HAMIMP(3 + FTIMO,3 + FTIMO) = -HAMIMP(1 + FTIMO,1 + FTIMO)
        HAMIMP(3 + FTIMO,4 + FTIMO) = CONJG(HAMIMP(1 + FTIMO,2 + FTIMO))

        HAMIMP(4 + FTIMO,1 + FTIMO) = CONJG(HAMIMP(2 + FTIMO,3 + FTIMO))
        !HAMIMP(4 + FTIMO,2 + FTIMO) = 0.0
        HAMIMP(4 + FTIMO,3 + FTIMO) = CONJG(HAMIMP(2 + FTIMO,1 + FTIMO))
        HAMIMP(4 + FTIMO,4 + FTIMO) = -CONJG(HAMIMP(2 + FTIMO,2 + FTIMO))

    end do

    if (finaldosorno == 'y') then
        print *, 'Commencing the DoS calculation.'
        NPNT = NUMEDOS
        open(5, file = 'greenhostfordos.txt', action = 'read')
        if (printgimp == 'y') then
            open(6, file = 'greenimpfordos.txt', action = 'write')
        endif
    else 
        print *, 'Commencing the calculation of G_imp for each energy value.'
        NPNT = NUME
        PREVDELTA = NEWDELTA
        PREVNU = NEWNU
        NEWNU = 0.d0
        NEWDELTA = (0.d0,0.d0)
        open(5, file = 'greenhost.txt', action = 'read')
        if (printgimp == 'y') then
            open(6, file = 'greenimp.txt', action = 'write')
        endif
    endif

    ! This commences the calculation of the impurity Green's function for each energy E
    do IE = 1, NPNT

        if (finaldosorno == 'y') then
            EZ = EZDOS(IE)

            densityperimpurity(1,IE) = REAL(EZ)
            do i = 1, NUMIMP
                densityperimpurity(1+i,IE) = 0.0
            end do
        else
            EZ = energies(IE)
        endif
        
        ! Constructs a 4*NUMIMP x 4*NUMIMP Identity matrix
        IDENTITY(:,:) = (0.0,0.0)
        do i = 1, 4*NUMIMP
            IDENTITY(i,i) = (1.0,0.0)
        end do

        ! This reads the G_0 matrix from the output greenimp.txt of BdG.f90 for each energy E
        do i = 1, 4*NUMIMP
            do j = 1, 4*NUMIMP
                read(5,*) tempval1, tempval2
                GREEN(i,j) = dcmplx(tempval1, tempval2)
            end do
        end do

        ! This routine returns the IDENTITY matrix as the the (I - G_0*ΔH) matrix
        ONE = (1.0,0.0)
        MINUSONE = (-1.0,0.0)
        call ZGEMM ('N', 'N', 4*NUMIMP, 4*NUMIMP, 4*NUMIMP, MINUSONE, GREEN, 4*NUMIMP, HAMIMP,&
        & 4*NUMIMP, ONE, IDENTITY, 4*NUMIMP)

        ! The following routines produce the impurity Green's function
        call ZGETRF (4*NUMIMP, 4*NUMIMP, IDENTITY, 4*NUMIMP, IPIV, INFO)
        call ZGETRS ('N', 4*NUMIMP, 4*NUMIMP, IDENTITY, 4*NUMIMP, IPIV, GREEN, 4*NUMIMP, INFO)

        ! Writes the new G elements in the corresponding file
        if (printgimp == 'y') then
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (6,'(F17.8,F17.8)') GREEN(i,j)
                end do  
            end do
        endif

        if (finaldosorno == 'y') then
            do i = 1, NUMIMP ! Calculation of density for each atom
                FTIMO = 4*(i-1)
                do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                    densityperimpurity(1+i,IE) = densityperimpurity(1+i,IE) -&
                    &(1.0/PI)*AIMAG(GREEN(j + FTIMO, j + FTIMO))
                end do
            end do
        else
            ! Calculation of the new charges and Deltas
            do i = 1, NUMIMP

                FTIMO = 4*(i-1)

                NEWNU(i) = NEWNU(i) - (1.0/PI)*AIMAG(WEIGHTS(IE)*(GREEN(1 + FTIMO, 1 + FTIMO) +&
                & GREEN(2 + FTIMO, 2 + FTIMO)))

                NEWDELTA(i) = NEWDELTA(i) + (1.0/PI)*VSUPCONDIMP(i)*AIMAG(WEIGHTS(IE)*GREEN(2 + FTIMO, 3 + FTIMO))
            end do
        endif

    end do

    close(5)
    close(6)

    if (finaldosorno /= 'y') then
        print *, 'Finished Run No.', NUMIT
        print *, '==========================='
        do i = 1, NUMIMP
            DIFFD(i) = abs(abs(NEWDELTA(i)) - abs(PREVDELTA(i)))
            DIFFN(i) = abs(NEWNU(i) - PREVNU(i))
            print *, 'Charge for atom No. ', i, '= ', NEWNU(i)
            print *, 'Delta for atom No. ', i, '= ', NEWDELTA(i)
        end do
        print *, '==========================='
    endif

    if (finaldosorno == 'y') then
        open(1, file = 'impdensities.txt', action = 'write')
            do j = 1, NUMEDOS
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

        print *, 'DoS for impurities finished.'
    endif

    ! Here an if check is inserted so that we go back and re-build the Hamiltonian using
    ! the new charges and Deltas. This is the self-consistency for the impurity problem.

    ! If the self-consistency needs to continue, we go back to 180 and reset the Hamiltonian
    ! using the NEWDELTA and NEWNU calculated above.

    ! Number of iterations
    if ((MAXVAL(DIFFN) > epsilon .or. MAXVAL(DIFFD) > epsilon) .and. NUMIT < MAXNUMIT) then  
        NUMIT = NUMIT + 1
        goto 180
    endif

    ! If the self-consistency criterion is met, then the selfconfin switch is turned to 'y'
    selfconfin = 'y'
    if (finaldosorno /= 'y') then
        print *, 'Should a DoS be performed for the converged impurity system? y/n'
        20 read *, finaldosorno
        if (finaldosorno /= 'y' .and. finaldosorno /= 'n') then
            print *, 'Invalid input, please enter y or n.'
            goto 20
        else if (finaldosorno == 'y') then

            open(1, file = 'energiesfordos.dat', action = 'read')
            do i = 1, NUMEDOS
                read(1,*) tempval1, tempval2
                EZDOS(i) = dcmplx(tempval1, tempval2)
            end do
            close(1)

            goto 180
            
        endif
    endif

    contains

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

end program BDG_IMP