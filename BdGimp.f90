program BDG_IMP
    implicit none

    integer :: i, j, NUME, NUMIMP, IE, INFO
    integer, allocatable, dimension(:) :: IPIV
    real*8 :: tempval1, tempval2, tempval3, tempval4
    real*8, allocatable, dimension(:) :: E0, BETA, E0IMP, BETAIMP
    complex*16 :: ONE, MINUSONE
    complex*16, allocatable, dimension(:) :: DELTA, DELTAIMP
    complex*16, allocatable, dimension(:,:) :: GREEN, HAMIMP, IDENTITY

    ! Reads the configuration info from the BdG.f90 output
    open(1, file = 'impconfig.dat', action = 'read')
    read(1,*) NUMIMP
    read(1,*) NUME
    close(1)

    allocate(E0(NUMIMP))
    allocate(E0IMP(NUMIMP))
    allocate(BETA(NUMIMP))
    allocate(BETAIMP(NUMIMP))
    allocate(DELTA(NUMIMP))
    allocate(DELTAIMP(NUMIMP))
    allocate(GREEN(4*NUMIMP,4*NUMIMP))
    allocate(HAMIMP(4*NUMIMP,4*NUMIMP))
    allocate(IDENTITY(4*NUMIMP,4*NUMIMP))

    allocate(IPIV(4*NUMIMP))

    open(1, file = 'impatoms.dat', action = 'read')
    do i = 1, NUMIMP
        read (1,*) E0(i), E0IMP(i), BETA(i), BETAIMP(i), tempval1, tempval2, tempval3, tempval4
        DELTA(i) = cmplx(tempval1, tempval3)
        DELTAIMP(i) = cmplx(tempval2, tempval4)
    end do
    close(1)

    ! Constructs the ΔH matrix
    do i = 1, NUMIMP

        HAMIMP(1 + (i-1)*4,1 + (i-1)*4) = E0IMP(i) - E0(i) + BETAIMP(i) - BETA(i)
        HAMIMP(1 + (i-1)*4,2 + (i-1)*4) = 0.0 ! If spin-flip effects are ignored
        HAMIMP(1 + (i-1)*4,3 + (i-1)*4) = 0.0
        HAMIMP(1 + (i-1)*4,4 + (i-1)*4) = DELTAIMP(i) - DELTA(i)

        HAMIMP(2 + (i-1)*4,1 + (i-1)*4) = 0.0 ! If spin-flip effects are ignored
        HAMIMP(2 + (i-1)*4,2 + (i-1)*4) = E0IMP(i) - E0(i) - BETAIMP(i) + BETA(i)
        HAMIMP(2 + (i-1)*4,3 + (i-1)*4) = HAMIMP(1 + (i-1)*4,4 + (i-1)*4)
        HAMIMP(2 + (i-1)*4,4 + (i-1)*4) = 0.0

        HAMIMP(3 + (i-1)*4,1 + (i-1)*4) = 0.0
        HAMIMP(3 + (i-1)*4,2 + (i-1)*4) = CONJG(HAMIMP(1 + (i-1)*4,4 + (i-1)*4))
        HAMIMP(3 + (i-1)*4,3 + (i-1)*4) = -HAMIMP(1 + (i-1)*4,1 + (i-1)*4)
        HAMIMP(3 + (i-1)*4,4 + (i-1)*4) = CONJG(HAMIMP(1 + (i-1)*4,2 + (i-1)*4))

        HAMIMP(4 + (i-1)*4,1 + (i-1)*4) = CONJG(HAMIMP(2 + (i-1)*4,3 + (i-1)*4))
        HAMIMP(4 + (i-1)*4,2 + (i-1)*4) = 0.0
        HAMIMP(4 + (i-1)*4,3 + (i-1)*4) = CONJG(HAMIMP(2 + (i-1)*4,1 + (i-1)*4))
        HAMIMP(4 + (i-1)*4,4 + (i-1)*4) = -CONJG(HAMIMP(2 + (i-1)*4,2 + (i-1)*4))

    end do

    ! This commences the calculation of the impurity Green's function for each energy E
    do IE = 1, NUME+1
        
        ! Constructs a 4*NUMIMP x 4*NUMIMP Identity matrix
        do i = 1, 4*NUMIMP
            do j = 1, 4*NUMIMP
                if (i == j) then
                    IDENTITY(i,j) = (1.0,0.0)
                else
                    IDENTITY(i,j) = (0.0,0.0)
                endif
            end do
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

end program BDG_IMP