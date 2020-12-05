program BDG_IMP
    implicit none

    integer :: i, j, NUME, NUMIMP, IE, INFO, FTIMO
    integer, allocatable, dimension(:) :: IPIV
    real*8 :: PI, KB, tempval1, tempval2
    real*8, allocatable, dimension(:) :: E0, E0IMP
    real*8, allocatable, dimension(:,:) :: BETA, BETAIMP, densityperimpurity
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), ONE, MINUSONE, posham(2,2), EZ
    complex*16, allocatable, dimension(:) :: DELTA, DELTAIMP, energies
    complex*16, allocatable, dimension(:,:) :: GREEN, HAMIMP, IDENTITY

    print *, 'Initiating the impurity problem calculations.'

    call CONSTANTS(IdentityPauli,xPauli,yPauli,zPauli,CI,PI,KB) ! Sets some universal constants.

    ! Reads the configuration info from the BdG.f90 output
    open(1, file = 'impconfig.dat', action = 'read')
    read(1,*) NUMIMP
    read(1,*) NUME
    close(1)

    allocate(E0(NUMIMP))
    allocate(E0IMP(NUMIMP))
    allocate(BETA(3,NUMIMP))
    allocate(BETAIMP(3,NUMIMP))
    allocate(DELTA(NUMIMP))
    allocate(DELTAIMP(NUMIMP))
    allocate(GREEN(4*NUMIMP,4*NUMIMP))
    allocate(HAMIMP(4*NUMIMP,4*NUMIMP))
    allocate(IDENTITY(4*NUMIMP,4*NUMIMP))

    allocate(energies(NUME+1))
    allocate(densityperimpurity(1+NUMIMP,NUME+1))

    allocate(IPIV(4*NUMIMP))

    open(1, file = 'energies.dat', action = 'read')
    do i = 1, NUME+1
        read(1,*) tempval1, tempval2
        energies(i) = dcmplx(tempval1, tempval2)
    end do
    close(1)

    open(1, file = 'impatoms.dat', action = 'read')
    do i = 1, NUMIMP
        read (1,*) E0(i), BETA(1,i), BETA(2,i), BETA(3,i), tempval1, tempval2
        DELTA(i) = dcmplx(tempval1, tempval2)
        read (1,*) E0IMP(i), BETAIMP(1,i), BETAIMP(2,i), BETAIMP(3,i), tempval1, tempval2
        DELTAIMP(i) = dcmplx(tempval1, tempval2)
    end do
    close(1)

    HAMIMP(:,:) = (0.d0,0.d0)

    ! Constructs the ΔH matrix
    do i = 1, NUMIMP

        FTIMO = 4*(i-1)

        posham = (E0IMP(i) - E0(i))*IdentityPauli - (BETAIMP(1,i)-BETA(1,i))*xPauli - (BETAIMP(2,i)-BETA(2,i))*yPauli -&
        & (BETAIMP(3,i)-BETA(3,i))*zPauli

        HAMIMP(1 + FTIMO,1 + FTIMO) = posham(1,1)
        HAMIMP(1 + FTIMO,2 + FTIMO) = posham(1,2)
        !HAMIMP(1 + FTIMO,3 + FTIMO) = 0.0
        HAMIMP(1 + FTIMO,4 + FTIMO) = DELTAIMP(i) - DELTA(i)

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

    print *, 'Commencing the calculation of G_imp for each energy value.'

    ! This commences the calculation of the impurity Green's function for each energy E
    do IE = 1, NUME+1

        EZ = energies(IE)
        densityperimpurity(1,IE) = REAL(EZ)
        do i = 1, NUMIMP
            densityperimpurity(1+i,IE) = 0.0
        end do
        
        ! Constructs a 4*NUMIMP x 4*NUMIMP Identity matrix
        IDENTITY(:,:) = (0.0,0.0)
        do i = 1, 4*NUMIMP
            IDENTITY(i,i) = (1.0,0.0)
        end do
    
        ! This reads the G_0 matrix from the output greenimp.txt of BdG.f90 for each energy E
        open(1, file = 'greenhost.txt', action = 'read')
        if (IE /= 1) then
            do i = 1, (IE-1)*(4*NUMIMP)**2 ! Skips the proper number of lines
                read (1,*)
            end do
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    read(1,*) tempval1, tempval2
                    GREEN(i,j) = dcmplx(tempval1, tempval2)
                end do  
            end do
        else
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    read(1,*) tempval1, tempval2
                    GREEN(i,j) = dcmplx(tempval1, tempval2)
                end do  
            end do
        endif
        close(1)

        ! This routine returns the IDENTITY matrix as the the (I - G_0*ΔH) matrix
        ONE = (1.0,0.0)
        MINUSONE = (-1.0,0.0)
        call ZGEMM ('N', 'N', 4*NUMIMP, 4*NUMIMP, 4*NUMIMP, MINUSONE, GREEN, 4*NUMIMP, HAMIMP,&
        & 4*NUMIMP, ONE, IDENTITY, 4*NUMIMP)

        ! The following routines produce the impurity Green's function
        call ZGETRF (4*NUMIMP, 4*NUMIMP, IDENTITY, 4*NUMIMP, IPIV, INFO)
        call ZGETRS ('N', 4*NUMIMP, 4*NUMIMP, IDENTITY, 4*NUMIMP, IPIV, GREEN, 4*NUMIMP, INFO)

        do i = 1, NUMIMP ! Calculation of density for each atom
            FTIMO = 4*(i-1)
            do j = 1, 2 ! n = u↑u↑* + u↓u↓*
                densityperimpurity(1+i,IE) = densityperimpurity(1+i,IE) -&
                &(1.0/PI)*AIMAG(GREEN(j + FTIMO, j + FTIMO))
            end do
        end do

        open(1, file = 'greenimp.txt', action = 'readwrite')
        if (IE /= 1) then
            do i = 1, (IE-1)*(4*NUMIMP)**2 ! Skips the proper number of lines
                read (1,*)
            end do
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1,'(F17.8,F17.8)') GREEN(i,j)
                end do  
            end do
        else
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    write (1,'(F17.8,F17.8)') GREEN(i,j)
                end do  
            end do
        endif
        close(1)

    end do

    print *, 'G_imp calculated.'

    open(1, file = 'impdensities.txt', action = 'write')
        do j = 1, NUME+1 ! Energies = Intervals + 1
            do i = 1, NUMIMP+1
                write (1,'(F17.8)',advance='no') densityperimpurity(i,j)
            end do
            write (1,*)
        end do
    close(1)

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