program BDG_IMP
    implicit none

    integer :: i, j, NUME, NUMIMP, IE, INFO, FTIMO, NUMEDOS, NPNT, BSTEPS, BCOUNTER, num, denom
    integer, allocatable, dimension(:) :: IPIV
    real*8 :: PI, KB, tempval1, tempval2, BFINAL, magB, BINC, ROTAXIS(3), theta
    real*8, allocatable, dimension(:) :: E0, E0IMP, nuzero, nuzeroIMP, ULCN, ULCNIMP, NEWNU, VSUPCOND, &
    &VSUPCONDIMP, NU, realenergies, BMAGS
    real*8, allocatable, dimension(:,:) :: BETA, BETAIMP, EBDOS
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), ONE, MINUSONE, posham(2,2), EZ
    complex*16, allocatable, dimension(:) :: DELTA, energies, WEIGHTS, EZDOS, NEWDELTA
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
    close(1)

    open(1, file = 'greenconfig.dat', action = 'read')
        do i = 1, 9
            read(1,*)
        end do ! Skips the file's first 9 lines
        read(1,*) ROTAXIS
        read(1,*) num, denom
    close(1)

    theta = (num*PI)/denom

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

    allocate(energies(NUME))
    allocate(EZDOS(NUMEDOS))
    allocate(WEIGHTS(NUME))

    allocate(realenergies(NUMEDOS))

    allocate(IPIV(4*NUMIMP))

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

    open(1, file = 'energiesfordos.dat', action = 'read')
        do i = 1, NUMEDOS
            read(1,*) tempval1, tempval2
            EZDOS(i) = dcmplx(tempval1, tempval2)
        end do
    close(1)

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
	
	do i = 1, NUMIMP
		BETAIMP(:,i) = 0.d0
	end do

    print *, 'Enter the final magnitude for B.'
    read *, BFINAL

    print *, 'Enter the number of steps.'
    read *, BSTEPS

    BCOUNTER = 1
    magB = abs(BETAIMP(1,1))
    BINC = (BFINAL-magB)/(BSTEPS-1)
    allocate(EBDOS(BSTEPS,NUMEDOS))
    allocate(BMAGS(BSTEPS))
    BMAGS(1) = magB
	
	call MAGCHAIN(NUMIMP,ROTAXIS,theta,magB,BETAIMP)

    190 NEWNU = NU ! Initial values
    NEWDELTA = DELTA ! Initial values

    finaldosorno = 'n' ! Starting value

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
    else 
        print *, 'Commencing the calculation of G_imp for each energy value.'
        NPNT = NUME
        NEWNU = 0.d0
        NEWDELTA = (0.d0,0.d0)
        open(5, file = 'greenhost.txt', action = 'read')
    endif

    ! This commences the calculation of the impurity Green's function for each energy E
    do IE = 1, NPNT

        if (finaldosorno == 'y') then
            EZ = EZDOS(IE)

            realenergies(IE) = REAL(EZ)
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

        if (finaldosorno == 'y') then
            ! Density of first impurity
            EBDOS(BCOUNTER,IE) = -(1.0/PI)*(AIMAG(GREEN(1,1))+AIMAG(GREEN(2,2)))
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
        print *, 'G_imp calculated.'
        print *, '==========================='
        !do i = 1, NUMIMP
        !    print *, 'Charge for atom No. ', i, '= ', NEWNU(i)
        !    print *, 'Delta for atom No. ', i, '= ', NEWDELTA(i)
        !end do
    endif

    if (finaldosorno == 'y') then
        print *, 'DoS for impurities finished.'
    endif

    ! Here an if check should be inserted so that we go back and re-build the Hamiltonian using
    ! the new charges and Deltas. This is the self-consistency for the impurity problem.

    ! If the self-consistency needs to continue, we go back to 180 and reset the Hamiltonian
    ! using the NEWDELTA and NEWNU calculated above.

    ! Number of iterations
    !if (NUMIT < MAXNUMIT) then  
    !    NUMIT = NUMIT + 1
    !    goto 180
    !endif

    ! If the self-consistency criterion is met, then the selfconfin switch is turned to 'y'
    selfconfin = 'y'
    if (finaldosorno /= 'y') then
        !print *, 'Should a DoS be performed for the converged impurity system? y/n'
        !20 read *, finaldosorno
        finaldosorno = 'y'
        !if (finaldosorno /= 'y' .and. finaldosorno /= 'n') then
            !print *, 'Invalid input, please enter y or n.'
            !goto 20
        !else if (finaldosorno == 'y') then

            goto 180
            
        !endif
    endif

    if (BCOUNTER < BSTEPS) then

        magB = magB + BINC
        BCOUNTER = BCOUNTER + 1

        BMAGS(BCOUNTER) = magB

        call MAGCHAIN(NUMIMP,ROTAXIS,theta,magB,BETAIMP)

        print *, 'Finished iteration No.', BCOUNTER-1, '.'

        goto 190

    else

        !open(1, file = 'EBDOS.txt', action = 'write')
        !    do j = 1, NUMEDOS
        !        write (1,'(F17.8)',advance='no') realenergies(j), ','
        !        do i = 1, BSTEPS
        !            if (i == BSTEPS) then
        !                write (1,'(F17.8)',advance='no') EBDOS(i,j)
        !            else
        !                write (1,'(F17.8, A)',advance='no') EBDOS(i,j), ','
        !            endif
        !        end do
        !        write (1,*)
        !    end do
        !close(1)

        !open(1, file = 'BMAGS.txt', action = 'write')
        !    do i = 1, BSTEPS
        !        write (1,'(F17.8)') BMAGS(i)
        !    end do
        !close(1)

        open(1, file = 'EBDOS.txt', action = 'write')
            do i = 1, BSTEPS
                do j = 1, NUMEDOS
                    write (1,'(F17.8, A, F17.8, A, F17.8)') BMAGS(i), ',', realenergies(j), ',', EBDOS(i,j)
                end do
            end do
        close(1)

    endif


    contains

    subroutine MAGCHAIN(NUMIMP,ROTAXIS,theta,magB,NEWBETA)
        implicit none

        integer :: NUMIMP, i
        real*8 :: theta, ROTAXIS(3), NEWBETA(3,NUMIMP), magB, initheta, TOL

        TOL = 0.000001 ! The tolerance for the "if" checks.

        NEWBETA(:,:) = (0.d0,0.d0)

        if (abs(ROTAXIS(1) - 1.0) < TOL .and. abs(ROTAXIS(2)) < TOL .and. abs(ROTAXIS(3)) < TOL) then

            initheta = 0.d0

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

            initheta = 0.d0

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

            initheta = 0.d0

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
            print *, 'Invalid input for the magchain rotation axis.'

        endif

    endsubroutine MAGCHAIN

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