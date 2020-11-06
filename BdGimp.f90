program BDG_IMP
    implicit none

    integer :: i, j, NUME, NUMIMP, IE, INFO
    integer, allocatable, dimension(:) :: IPIV
    real*8 :: PI, KB, tempval1, tempval2
    real*8, allocatable, dimension(:) :: E0, E0IMP
    real*8, allocatable, dimension(:,:) :: BETA, BETAIMP
    complex*16 :: CI, IdentityPauli(2,2), xPauli(2,2), yPauli(2,2), zPauli(2,2), ONE, MINUSONE, posham(2,2)
    complex*16, allocatable, dimension(:) :: DELTA, DELTAIMP
    complex*16, allocatable, dimension(:,:) :: GREEN, HAMIMP, IDENTITY

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

    allocate(IPIV(4*NUMIMP))

    open(1, file = 'impatoms.dat', action = 'read')
    do i = 1, NUMIMP
        read (1,*) E0(i), BETA(1,i), BETA(2,i), BETA(3,i), tempval1, tempval2
        DELTA(i) = cmplx(tempval1, tempval2)
        read (1,*) E0IMP(i), BETAIMP(1,i), BETAIMP(2,i), BETAIMP(3,i), tempval1, tempval2
        DELTAIMP(i) = cmplx(tempval1, tempval2)
    end do
    close(1)

    HAMIMP(:,:) = (0.0,0.0)

    ! Constructs the ΔH matrix
    do i = 1, NUMIMP

        posham = (E0IMP(i) - E0(i))*IdentityPauli - (BETAIMP(1,i)-BETA(1,i))*xPauli - (BETAIMP(2,i)-BETA(2,i))*yPauli -&
        & (BETAIMP(3,i)-BETA(3,i))*zPauli

        HAMIMP(1 + (i-1)*4,1 + (i-1)*4) = posham(1,1)
        HAMIMP(1 + (i-1)*4,2 + (i-1)*4) = posham(1,2)
        !HAMIMP(1 + (i-1)*4,3 + (i-1)*4) = 0.0
        HAMIMP(1 + (i-1)*4,4 + (i-1)*4) = DELTAIMP(i) - DELTA(i)

        HAMIMP(2 + (i-1)*4,1 + (i-1)*4) = posham(2,1)
        HAMIMP(2 + (i-1)*4,2 + (i-1)*4) = posham(2,2)
        HAMIMP(2 + (i-1)*4,3 + (i-1)*4) = HAMIMP(1 + (i-1)*4,4 + (i-1)*4)
        !HAMIMP(2 + (i-1)*4,4 + (i-1)*4) = 0.0

        !HAMIMP(3 + (i-1)*4,1 + (i-1)*4) = 0.0
        HAMIMP(3 + (i-1)*4,2 + (i-1)*4) = CONJG(HAMIMP(1 + (i-1)*4,4 + (i-1)*4))
        HAMIMP(3 + (i-1)*4,3 + (i-1)*4) = -HAMIMP(1 + (i-1)*4,1 + (i-1)*4)
        HAMIMP(3 + (i-1)*4,4 + (i-1)*4) = CONJG(HAMIMP(1 + (i-1)*4,2 + (i-1)*4))

        HAMIMP(4 + (i-1)*4,1 + (i-1)*4) = CONJG(HAMIMP(2 + (i-1)*4,3 + (i-1)*4))
        !HAMIMP(4 + (i-1)*4,2 + (i-1)*4) = 0.0
        HAMIMP(4 + (i-1)*4,3 + (i-1)*4) = CONJG(HAMIMP(2 + (i-1)*4,1 + (i-1)*4))
        HAMIMP(4 + (i-1)*4,4 + (i-1)*4) = -CONJG(HAMIMP(2 + (i-1)*4,2 + (i-1)*4))

    end do

    ! This commences the calculation of the impurity Green's function for each energy E
    do IE = 1, NUME+1
        
        ! Constructs a 4*NUMIMP x 4*NUMIMP Identity matrix
        IDENTITY(:,:) = (0.0,0.0)
        do i = 1, 4*NUMIMP
            IDENTITY(i,i) = (1.0,0.0)
        end do
    
        ! This reads the G_0 matrix from the output greenimp.txt of BdG.f90 for each energy E
        open(1, file = 'greenimp.txt', action = 'read')
        if (IE /= 1) then
            do i = 1, (IE-1)*(4*NUMIMP)**2 ! Skips the proper number of lines
                read (1,*)
            end do
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    read(1,*) tempval1, tempval2
                    GREEN(i,j) = cmplx(tempval1, tempval2)
                end do  
            end do
        else
            do i = 1, 4*NUMIMP
                do j = 1, 4*NUMIMP
                    read(1,*) tempval1, tempval2
                    GREEN(i,j) = cmplx(tempval1, tempval2)
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

        open(1, file = 'greenimpNEW.txt', action = 'readwrite')
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