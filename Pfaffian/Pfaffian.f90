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

    real :: ALTPFAFZERO(2), ALTPFAFPI(2)
    real, allocatable, dimension(:) :: WORK
    real, allocatable, dimension(:,:) :: ZEROALPHA, PIALPHA

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
    allocate(ZEROALPHA(4*NUMT,4*NUMT))
    allocate(PIALPHA(4*NUMT,4*NUMT))

    call RPTS(a_1,a_2,a_3,NCELLS,RLATT,IRLATTMAX) ! Here we construct the RLATT Matrix consisting of the lattice sites

    call HOPPS(RLATT,IRLATTMAX,R0,NUMCHEMTYPES,LHOPS,NUMT,CHEMTYPE,TPTS,PREFACTORS)

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
    
    kcounter = NUMK/2 ! Corresponds to k = 0
    call FINALHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)
    call MAKEALPHA(NUMT,HAMILTONIAN,ZEROALPHA)

    call ALTPFAF('U', 'P', 4*NUMT, ZEROALPHA, 4*NUMT, ALTPFAFZERO, IWORK, WORK, LWORK, INFO)
    
    kcounter = NUMK ! Corresponds to k = pi
    call FINALHAM(kcounter,NUMK,HAMILTONIANPREP,EXPONS,HOPPVALS,NEIGHBNUM,JTYPE,MAXNEIGHB,NUMT,HAMILTONIAN)
    call MAKEALPHA(NUMT,HAMILTONIAN,PIALPHA)

    call ALTPFAF('U', 'P', 4*NUMT, PIALPHA, 4*NUMT, ALTPFAFPI, IWORK, WORK, LWORK, INFO)

    TOTSIGN = ALTPFAFPI(1)*ALTPFAFZERO(1)

    if (TOTSIGN > 0.0) then
        TOTSIGN = 1.0
    else if (TOTSIGN < 0.0) then
        TOTSIGN = -1.0
    endif

    !print *, 'Pfaffian for k = 0 is', ALTPFAFZERO(1), 'times 10^', ALTPFAFZERO(2)
    !print *, 'Pfaffian for k = pi is', ALTPFAFPI(1), 'times 10^', ALTPFAFPI(2)

    if (manyruns == 'y') then

        if (BCOUNTER < BSTEPS) then

            PFSIGNS(BCOUNTER,1) = magB
            PFSIGNS(BCOUNTER,2) = TOTSIGN

            magB = magB + BINC
            BCOUNTER = BCOUNTER + 1

            print *, 'Finished iteration No.', BCOUNTER-1
            print *, 'The Majorana number was equal to', TOTSIGN

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

    subroutine PFCALC (SK,LDS,N,Pf,WORK,IH)
        implicit none
        integer,intent(in) :: N,LDS
        integer :: IS, M, Iden, IH
        Real*8,intent(inout),dimension(LDS,N) :: SK
        Real*8,intent(inout),dimension(LDS,N) :: WORK
        Real*8 :: ddot
        Real*8 :: epsln,xmod,u2,ki,al,one,zero
        Real*8,intent(inout) :: Pf
        
        one = 1.0d+00
        zero = 0.0d+00
        
        epsln = 1.0d-13   ! smallest number such that 1+epsln=1
        
        if (IH /= 1 .and. IH /= 2) then
           write(6,'("Value of IH undefined in pfaffianH",I2)') IH
        end if
        
        if (mod(N,2) == 1) then 
            ! N odd  Pf=0              
            Pf = 0.0d+00 
        else
            ! N even                               
            if (N == 2) then
                ! The pfaffian of a 2x2 matrix is Pf=Sk(1,2)
                Pf = Sk(1,2)
            else
                ! N even and > 2       
                Pf = one
                Iden = 0 ! counts the number of identity transformations
        
                do IS = 1, N-2, IH
        
                    M = N - IS
                    ! Column N-IS+1 of Sk is copied into column IS of working vector
                    ! This is vector x with dimension M of the Householder transformation
                        
                    call dcopy(M,Sk(1,N-IS+1),1,WORK(1,IS),1) 
        
                    ! Norm of x
                    xmod = dsqrt( ddot(M,WORK(1,IS),1,WORK(1,IS),1) )  ! |x|
        
                    ! if Norm of x = 0 do nothing
        
                    if (xmod > epsln) then
                        ki = sign(xmod,Sk(M,M+1))
        
                        !   u = x -+ sign(Sk(N-IS,N-IS+1))*xmod e (N-IS)
                        !   e (M) is the unit cartesian vector with 1 in position M and
                        !   zero elsewhere
                        work(M,IS) = work(M,IS) - ki
        
                        !   Norm of u
                        u2 = ddot(M,WORK(1,IS),1,WORK(1,IS),1)  ! |u|**2
        
                        !  If the norm of u is too small the negative sign in the definition
                        !  of u is taken
        
                        if (u2 < epsln**2) then 
        
                            !  ki changes sign
                            ki = -ki
        
                            work(M,IS) = work(M,IS) - 2.0d+00*ki
                            
                            u2 = ddot(M,WORK(1,IS),1,WORK(1,IS),1)  ! |u|**2
        
                        end if
        
                        ! v = A(N-1) u --> work(1,N)
                 
                        call dgemv('N',M,M,one,Sk,LDS,work(1,IS),1,zero,work(1,N),1)  
        
                        ! Update of A(N-IS)
                        al = 2.0d+00/u2 
        
                        ! Sk -> Sk + al*u v'
                        call dger(M,M, al,work(1,IS),1,work(1,N),1,Sk,LDS)
        
                        ! Sk -> Sk - al*v u'
                        call dger(M,M,-al,work(1,N),1,work(1,IS),1,Sk,LDS)
        
                    else ! xmod > epsln
        
                        ki = 0.0d+00
                        Iden = Iden + 1
               
                    end if   ! xmod > epsln
        
                    if (mod(IS,2) == 1) Pf = Pf * ki
        
                end do ! IS
        
                Pf = Pf * Sk(1,2) * dfloat(1-2*mod(Iden,2))
        
                if (IH == 2) Pf = Pf * float(1-2*mod(N/2-1,2))
        
            end if ! N==2
           
        end if ! mod(N,2) == 1

        return

    end subroutine PFCALC
    
    SUBROUTINE SMUL10(A, B)
        implicit none
        REAL A(2)
        REAL B
        REAL EXPONENT, SLAMCH
        INTEGER IEXPONENT
        EXTERNAL SLAMCH
  
        A(1) = A(1) * B
  
        IF( A(1).EQ.0.0E0 ) THEN
           A(1) = 0.0E0
           A(2) = 0.0E0
        ELSE
           EXPONENT = LOG10(ABS(A(1)))
  
           IF (EXPONENT .GE. 0E0) THEN
              IEXPONENT=INT(EXPONENT)
           ELSE
              IEXPONENT=INT(EXPONENT)-1
           END IF
  
        !     If 10**IEXPONENT is smaller than the safe minimum for inversion,
        !     the number is considered 0
           IF( SLAMCH( "S" ) > 10E0**IEXPONENT ) THEN
              A(1) = 0.0E0
              A(2) = 0.0E0
           ELSE
              A(1) = A(1)/10E0**IEXPONENT
              A(2) = A(2) + IEXPONENT
           END IF
        END IF
  
    end subroutine SMUL10
  
    subroutine DMUL10(A, B)
        implicit none
        DOUBLE PRECISION A(2)
        DOUBLE PRECISION B
        DOUBLE PRECISION EXPONENT, DLAMCH
        INTEGER IEXPONENT
        EXTERNAL DLAMCH
  
        A(1) = A(1) * B
  
        IF( A(1).EQ.0.0D0 ) THEN
           A(1) = 0.0D0
           A(2) = 0.0D0
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
              A(1) = 0.0D0
              A(2) = 0.0D0
           ELSE
              A(1) = A(1)/10D0**IEXPONENT
              A(2) = A(2) + IEXPONENT
           END IF
        END IF
  
    end subroutine DMUL10
  
    subroutine CMUL10(A, B)
        implicit none
        COMPLEX A(2)
        COMPLEX B
        REAL EXPONENT, SLAMCH
        INTEGER IEXPONENT
        EXTERNAL SLAMCH
        A(1) = A(1) * B
  
        IF( A(1).EQ.(0.0E0, 0.0E0) ) THEN
           A(1) = (0.0E0, 0.0E0)
           A(2) = (0.0E0, 0.0E0)
        ELSE
           EXPONENT = LOG10(ABS(A(1)))
  
           IF (EXPONENT .GE. 0E0) THEN
              IEXPONENT=INT(EXPONENT)
           ELSE
              IEXPONENT=INT(EXPONENT)-1
           END IF
  
        !     If 10**IEXPONENT is smaller than the safe minimum for inversion,
        !     the number is considered 0
           IF( SLAMCH( "S" ) > 10E0**IEXPONENT ) THEN
              A(1) = (0.0E0, 0.0E0)
              A(2) = (0.0E0, 0.0E0)
           ELSE
              A(1) = A(1)/10E0**IEXPONENT
              A(2) = A(2) + IEXPONENT
           END IF
        END IF
  
    end subroutine CMUL10
  
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

    subroutine SSKMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
        implicit none
        !     .. Scalar Arguments ..
        REAL ALPHA,BETA
        INTEGER INCX,INCY,LDA,N
        CHARACTER UPLO
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKMV  performs the matrix-vector  operation
        !
        !     y := alpha*A*x + beta*y,
        !
        !  where alpha and beta are scalars, x and y are n element vectors and
        !  A is an n by n skew-symmetric matrix.
        !
        !  Arguments
        !  ==========
        !
        !  UPLO   - CHARACTER*1.
        !           On entry, UPLO specifies whether the upper or lower
        !           triangular part of the array A is to be referenced as
        !           follows:
        !
        !              UPLO = 'U' or 'u'   Only the upper triangular part of A
        !                                  is to be referenced.
        !
        !              UPLO = 'L' or 'l'   Only the lower triangular part of A
        !                                  is to be referenced.
        !
        !           Unchanged on exit.
        !
        !  N      - INTEGER.
        !           On entry, N specifies the order of the matrix A.
        !           N must be at least zero.
        !           Unchanged on exit.
        !
        !  ALPHA  - REAL      .
        !           On entry, ALPHA specifies the scalar alpha.
        !           Unchanged on exit.
        !
        !  A      - REAL       array of DIMENSION ( LDA, n ).
        !           Before entry with  UPLO = 'U' or 'u', the leading n by n
        !           upper triangular part of the array A must contain the upper
        !           triangular part of the skew-symmetric matrix and the strictly
        !           lower triangular part of A is not referenced.
        !           Before entry with UPLO = 'L' or 'l', the leading n by n
        !           lower triangular part of the array A must contain the lower
        !           triangular part of the skew-symmetric matrix and the strictly
        !           upper triangular part of A is not referenced.
        !           Note that the imaginary parts of the diagonal elements need
        !           not be set and are assumed to be zero.
        !           Unchanged on exit.
        !
        !  LDA    - INTEGER.
        !           On entry, LDA specifies the first dimension of A as declared
        !           in the calling (sub) program. LDA must be at least
        !           max( 1, n ).
        !           Unchanged on exit.
        !
        !  X      - REAL       array of dimension at least
        !           ( 1 + ( n - 1 )*abs( INCX ) ).
        !           Before entry, the incremented array X must contain the n
        !           element vector x.
        !           Unchanged on exit.
        !
        !  INCX   - INTEGER.
        !           On entry, INCX specifies the increment for the elements of
        !           X. INCX must not be zero.
        !           Unchanged on exit.
        !
        !  BETA   - REAL      .
        !           On entry, BETA specifies the scalar beta. When BETA is
        !           supplied as zero then Y need not be set on input.
        !           Unchanged on exit.
        !
        !  Y      - REAL       array of dimension at least
        !           ( 1 + ( n - 1 )*abs( INCY ) ).
        !           Before entry, the incremented array Y must contain the n
        !           element vector y. On exit, Y is overwritten by the updated
        !           vector y.
        !
        !  INCY   - INTEGER.
        !           On entry, INCY specifies the increment for the elements of
        !           Y. INCY must not be zero.
        !           Unchanged on exit.
        !
        !
        !  Level 2 Blas routine.
        !
        !  -- Written on 20/22/2010
        !     Michael Wimmer, Universiteit Leiden
        !
        !     .. Parameters ..
        REAL ONE
        PARAMETER (ONE= 1.0E+0)
        REAL ZERO
        PARAMETER (ZERO= 0.0E+0)
        !     ..
        !     .. Local Scalars ..
        REAL TEMP1,TEMP2
        INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
        !     ..
        !     .. External Functions ..
        LOGICAL LSAME
        EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
        EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
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
        !
        !     Quick return if possible.
        !
        IF ((N.EQ.0) .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
        !
        !     Set up the start points in  X  and  Y.
        !
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
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through the triangular part
        !     of A.
        !
        !     First form  y := beta*y.
        !
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
        !
        !        Form  y  when A is stored in upper triangle.
        !
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
        !
        !        Form  y  when A is stored in lower triangle.
        !
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
        !
        RETURN
    end subroutine SSKMV

    subroutine SSKR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        implicit none
        !     .. Scalar Arguments ..
        REAL ALPHA
        REAL BETA
        INTEGER K,LDA,LDB,LDC,N
        CHARACTER TRANS,UPLO
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),B(LDB,*),C(LDC,*)
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKR2K  performs one of the skew-symmetric rank 2k operations
        !
        !     C := alpha*A*B^T - alpha*B*A^T + beta*C,
        !
        !  or
        !
        !     C := alpha*A^T*B - alpha*B^T*A + beta*C,
        !
        !  where  alpha and beta are scalars,  C is an  n by n
        !  skew-symmetric matrix and  A and B  are  n by k matrices in the first case
        !  and  k by n  matrices in the second case.
        !
        !  Arguments
        !  ==========
        !
        !  UPLO   - CHARACTER*1.
        !           On  entry,   UPLO  specifies  whether  the  upper  or  lower
        !           triangular  part  of the  array  C  is to be  referenced  as
        !           follows:
        !
        !              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
        !                                  is to be referenced.
        !
        !              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
        !                                  is to be referenced.
        !
        !           Unchanged on exit.
        !
        !  TRANS  - CHARACTER*1.
        !           On entry,  TRANS  specifies the operation to be performed as
        !           follows:
        !
        !              TRANS = 'N' or 'n'    C := alpha*A*B^T   -
        !                                         alpha*B*A^T   +
        !                                         beta*C.
        !
        !              TRANS = 'T' or 't'    C := alpha*A^T*B   -
        !                                         alpha*B^T*A   +
        !                                         beta*C.
        !
        !           Unchanged on exit.
        !
        !  N      - INTEGER.
        !           On entry,  N specifies the order of the matrix C.  N must be
        !           at least zero.
        !           Unchanged on exit.
        !
        !  K      - INTEGER.
        !           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
        !           of  columns  of the  matrices  A and B,  and on  entry  with
        !           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
        !           matrices  A and B.  K must be at least zero.
        !           Unchanged on exit.
        !
        !  ALPHA  - REAL         .
        !           On entry, ALPHA specifies the scalar alpha.
        !           Unchanged on exit.
        !
        !  A      - REAL       array of DIMENSION ( LDA, ka ), where ka is
        !           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
        !           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
        !           part of the array  A  must contain the matrix  A,  otherwise
        !           the leading  k by n  part of the array  A  must contain  the
        !           matrix A.
        !           Unchanged on exit.
        !
        !  LDA    - INTEGER.
        !           On entry, LDA specifies the first dimension of A as declared
        !           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
        !           then  LDA must be at least  max( 1, n ), otherwise  LDA must
        !           be at least  max( 1, k ).
        !           Unchanged on exit.
        !
        !  B      - REAL       array of DIMENSION ( LDB, kb ), where kb is
        !           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
        !           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
        !           part of the array  B  must contain the matrix  B,  otherwise
        !           the leading  k by n  part of the array  B  must contain  the
        !           matrix B.
        !           Unchanged on exit.
        !
        !  LDB    - INTEGER.
        !           On entry, LDB specifies the first dimension of B as declared
        !           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
        !           then  LDB must be at least  max( 1, n ), otherwise  LDB must
        !           be at least  max( 1, k ).
        !           Unchanged on exit.
        !
        !  BETA   - REAL          .
        !           On entry, BETA specifies the scalar beta.
        !           Unchanged on exit.
        !
        !  C      - REAL          array of DIMENSION ( LDC, n ).
        !           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
        !           upper triangular part of the array C must contain the upper
        !           triangular part  of the  skew-symmetric matrix  and the strictly
        !           lower triangular part of C is not referenced.  On exit, the
        !           upper triangular part of the array  C is overwritten by the
        !           upper triangular part of the updated matrix.
        !           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
        !           lower triangular part of the array C must contain the lower
        !           triangular part  of the  skew-symmetric matrix  and the strictly
        !           upper triangular part of C is not referenced.  On exit, the
        !           lower triangular part of the array  C is overwritten by the
        !           lower triangular part of the updated matrix.
        !           Note that the diagonal elements need
        !           not be set,  they are assumed to be zero,  and on exit they
        !           are set to zero.
        !
        !  LDC    - INTEGER.
        !           On entry, LDC specifies the first dimension of C as declared
        !           in  the  calling  (sub)  program.   LDC  must  be  at  least
        !           max( 1, n ).
        !           Unchanged on exit.
        !
        !
        !  Level 3 Blas routine.
        !
        !  -- Written on 10/22/2010
        !     Michael Wimmer, Universiteit Leiden
        !     Based on ZHER2K from BLAS (www.netlib.org)
        !
        !     .. External Functions ..
        LOGICAL LSAME
        EXTERNAL LSAME
        !     ..
        !     .. External Subroutines ..
        EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !     .. Local Scalars ..
        REAL TEMP1,TEMP2
        INTEGER I,INFO,J,L,NROWA
        LOGICAL UPPER
        !     ..
        !     .. Parameters ..
        REAL ONE
        PARAMETER (ONE= 1.0E+0)
        REAL ZERO
        PARAMETER (ZERO= 0.0E+0)
        !     ..
        !
        !     Test the input parameters.
        !
        IF (LSAME(TRANS,'N')) THEN
            NROWA = N
        ELSE
            NROWA = K
        END IF
        UPPER = LSAME(UPLO,'U')
        !
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
            CALL XERBLA('SSKR2K',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF ((N.EQ.0) .OR. (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
        !
        !     And when  alpha.eq.zero.
        !
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
        !
        !     Start the operations.
        !
        IF (LSAME(TRANS,'N')) THEN
        !
        !        Form  C := alpha*A* B^T - alpha*B*A^T + beta*
        !                   C.
        !
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
        !
        !        Form  C := alpha*A^T*B - alpha*B^T*A + beta*
        !                   C.
        !
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
        !
        RETURN
    end subroutine SSKR2K

    subroutine SSKR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
        implicit none
        !     .. Scalar Arguments ..
        REAL ALPHA
        INTEGER INCX,INCY,LDA,N
        CHARACTER UPLO
        !     ..
        !     .. Array Arguments ..
        REAL A(LDA,*),X(*),Y(*)
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKR2  performs the skew-symmetric rank 2 operation
        !
        !     A := alpha*x*y^T - alpha *y*x^T + A,
        !
        !  where alpha is a scalar, x and y are n element vectors and A is an n
        !  by n skew-symmetric matrix (A^T=-A).
        !
        !  Arguments
        !  ==========
        !
        !  UPLO   - CHARACTER*1.
        !           On entry, UPLO specifies whether the upper or lower
        !           triangular part of the array A is to be referenced as
        !           follows:
        !
        !              UPLO = 'U' or 'u'   Only the upper triangular part of A
        !                                  is to be referenced.
        !
        !              UPLO = 'L' or 'l'   Only the lower triangular part of A
        !                                  is to be referenced.
        !
        !           Unchanged on exit.
        !
        !  N      - INTEGER.
        !           On entry, N specifies the order of the matrix A.
        !           N must be at least zero.
        !           Unchanged on exit.
        !
        !  ALPHA  - REAL      .
        !           On entry, ALPHA specifies the scalar alpha.
        !           Unchanged on exit.
        !
        !  X      - REAL       array of dimension at least
        !           ( 1 + ( n - 1 )*abs( INCX ) ).
        !           Before entry, the incremented array X must contain the n
        !           element vector x.
        !           Unchanged on exit.
        !
        !  INCX   - INTEGER.
        !           On entry, INCX specifies the increment for the elements of
        !           X. INCX must not be zero.
        !           Unchanged on exit.
        !
        !  Y      - REAL       array of dimension at least
        !           ( 1 + ( n - 1 )*abs( INCY ) ).
        !           Before entry, the incremented array Y must contain the n
        !           element vector y.
        !           Unchanged on exit.
        !
        !  INCY   - INTEGER.
        !           On entry, INCY specifies the increment for the elements of
        !           Y. INCY must not be zero.
        !           Unchanged on exit.
        !
        !  A      - REAL       array of DIMENSION ( LDA, n ).
        !           Before entry with  UPLO = 'U' or 'u', the leading n by n
        !           upper triangular part of the array A must contain the upper
        !           triangular part of the skew-symmetric matrix and the strictly
        !           lower triangular part of A is not referenced. On exit, the
        !           upper triangular part of the array A is overwritten by the
        !           upper triangular part of the updated matrix.
        !           Before entry with UPLO = 'L' or 'l', the leading n by n
        !           lower triangular part of the array A must contain the lower
        !           triangular part of the skew-symmetric matrix and the strictly
        !           upper triangular part of A is not referenced. On exit, the
        !           lower triangular part of the array A is overwritten by the
        !           lower triangular part of the updated matrix.
        !           Note that the imaginary parts of the diagonal elements need
        !           not be set, they are assumed to be zero, and on exit they
        !           are set to zero.
        !
        !  LDA    - INTEGER.
        !           On entry, LDA specifies the first dimension of A as declared
        !           in the calling (sub) program. LDA must be at least
        !           max( 1, n ).
        !           Unchanged on exit.
        !
        !
        !  Level 2 Blas routine.
        !
        !  -- Written on 10/22/2010
        !     Michael Wimmer, Universiteit Leiden
        !
        !
        !     .. Parameters ..
        REAL ZERO
        PARAMETER (ZERO = 0.0E+0)
        !     ..
        !     .. Local Scalars ..
        REAL TEMP1,TEMP2
        INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
        !     ..
        !     .. External Functions ..
        LOGICAL LSAME
        EXTERNAL LSAME
        !     ..
        !     .. External routines ..
        EXTERNAL XERBLA
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC MAX
        !     ..
        !
        !     Test the input parameters.
        !
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
            CALL XERBLA('SSKR2 ',INFO)
            RETURN
        END IF
        !
        !     Quick return if possible.
        !
        IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
        !
        !     Set up the start points in X and Y if the increments are not both
        !     unity.
        !
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
        !
        !     Start the operations. In this version the elements of A are
        !     accessed sequentially with one pass through the triangular part
        !     of A.
        !
        IF (LSAME(UPLO,'U')) THEN
        !
        !        Form  A  when A is stored in the upper triangle.
        !
            IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
                DO 20 J = 1,N
                    IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                        TEMP1 = ALPHA*Y(J)
                        TEMP2 = ALPHA*X(J)
                        DO 10 I = 1,J - 1
                            A(I,J) = A(I,J) + X(I)*TEMP1 - Y(I)*TEMP2
        10              CONTINUE
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
        !
        !        Form  A  when A is stored in the lower triangle.
        !
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
        !
            RETURN
    end subroutine SSKR2

    subroutine SLASKTRF(UPLO, MODE, N, NB, A, LDA, IPIV, W, LDW, INFO)
        implicit none
        !
        !     .. Scalar Arguments ..
        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, LDW, N, NB, STEP
        !     ..
        !     .. Array Arguments ..
        INTEGER            IPIV( * )
        REAL               A( LDA, * )
        REAL               W( LDW, * )
        !
        !  Purpose
        !  =======
        !
        !  SLASKTRF computes NB steps of the factorization of a skew-symmetric
        !  matrix A using the Parlett-Reid algorithm:
        !
        !     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
        !
        !  If UPLO = 'U', SLASKTRF reduces the last NB rows and columns of a
        !  matrix, of which the upper triangle is supplied;
        !  if UPLO = 'L', SLASKTRF reduces the first NB rows and columns of a
        !  matrix, of which the lower triangle is supplied.
        !
        !  Alternatively, the routine can also be used to compute a partial
        !  tridiagonal form
        !
        !  This is an auxiliary routine called by DSKTRF.
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          Specifies whether the upper or lower triangular part of the
        !          symmetric matrix A is stored:
        !          = 'U':  Upper triangular
        !          = 'L':  Lower triangular
        !
        !  MODE    (input) CHARACTER*1
        !          = 'N':  A is fully tridiagonalized
        !          = 'P':  A is partially tridiagonalized for Pfaffian computation
        !                  (details see below)
        !
        !  N       (input) INTEGER
        !          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
        !
        !  NB      (input) INTEGER
        !          The number of rows and columns to be reduced.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the skew-symmetric matrix A.
        !            If UPLO = 'U', the leading n-by-n upper triangular part
        !               of A contains the upper triangular part of the matrix A,
        !               and the strictly lower triangular part of A is not referenced.
        !            If UPLO = 'L', the leading n-by-n lower triangular part
        !               of A contains the lower triangular part of the matrix A,
        !               and the strictly upper triangular part of A is not referenced.
        !          On exit:
        !          if UPLO = 'U', the last NB columns have been reduced to
        !            tridiagonal form, with the diagonal elements overwriting
        !            the diagonal elements of A; the elements above the diagonal
        !            represent the upper unit triangular matrix U. If MODE = 'P',
        !            only the even columns of the last NB columns contain
        !            meaningful values.
        !          if UPLO = 'L', the first NB columns have been reduced to
        !            tridiagonal form, with the diagonal elements overwriting
        !            the diagonal elements of A; the elements below the diagonal
        !            represent the lower unit triangular matrix L. If MODE = 'P',
        !            only the  odd columns of the first NB columns contain
        !            meaningful values.
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        !  IPIV    (output) INTEGER array, dimension (N)
        !          Information about the permutation matrix P: row and column
        !          i are interchanged with IPIV(i). If UPLO = 'U', those
        !          interchanges are done in the order i = N ... 1, if UPLO = 'L'
        !          in the order i = 1 ... N.
        !
        !  W       (workspace) REAL array, dimension (LDW,NB)
        !          The n-by-nb matrix W required to update the unreduced part
        !          of A. (The update is performed in this routine)
        !
        !  LDW     (input) INTEGER
        !          The leading dimension of the array W. LDW >= max(1,N).
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
   
        !     .. Local Scalars ..
        INTEGER            K, KK, KP, NPANEL, WK
        REAL               COLMAX, T
        !     ..
        !     .. External Functions ..
        LOGICAL            LSAME
        INTEGER            ISAMAX
        EXTERNAL           LSAME, ISAMAX
        !     ..
        !     .. External Subroutines ..
        EXTERNAL           SSCAL, SSWAP, SGEMV, SCOPY, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     No safety checks, it's an internal function
   
        INFO = 0
   
        IF( LSAME( MODE, 'P' ) ) THEN
            STEP = 2
        ELSE
            STEP = 1
        END IF
   
        !     double the amount of panels before the block update, if STEP == 2
        NPANEL = NB * STEP
   
        IF( LSAME( UPLO, 'U' ) ) THEN
   
        !     Factorize A as U * T * U^T using the upper triangle of A
   
        !     Consider the last NB columns of the matrix:
        !     First compute the transformations and apply them to the last NB
        !     columns (only), then apply the accumulated transformations to the rest
        !     of the matrix
   
            WK = 0
            DO 10 K=N, MAX(N-NPANEL+1, 2), -1
        !
        !     Update A(1:K,K) with all the accumulated transformations
        !     (if K<N: K=N,N-1 is updated when computing the Gauss vector,
        !     K=N-1 is not affected by the K=N transform)
        !
               KK = K-1
   
               IF( K .LT. N) THEN
   
                  IF( WK .GT. 0 ) THEN
                     A( K, K ) = ZERO
                     CALL SGEMV( 'N', K, WK, +ONE, A( 1, N-(WK-1)*STEP ), LDA*STEP, W( K, NB-WK+1 ),&
                     & LDW, ONE, A( 1, K ), 1 )
                     CALL SGEMV( 'N', K, WK, -ONE, W( 1, NB-WK+1 ), LDW, A( K, N-(WK-1)*STEP ),&
                     & LDA*STEP, ONE, A( 1, K ), 1 )
                     A( K, K ) = ZERO
                  END IF
   
        !     Store the (updated) column K in W(:,NB-WK+1)
                  IF( MOD(K, STEP).EQ.1 .OR. STEP.EQ.1 ) THEN
                     WK = WK + 1
                     CALL SCOPY(K, A(1, K), 1, W(1, NB-WK+1 ), 1)
                  END IF
               END IF
   
               IF( MOD(K, STEP) .EQ. 0) THEN
        !     For STEP == 1, process every column, but if
        !     STEP == 2, do only things for the even columns
        
        !     Find the pivot
                  KP = ISAMAX(K-1, A( 1, K ), 1)
                  COLMAX = ABS( A( KP, K ) )
   
                  IF( COLMAX.EQ.ZERO ) THEN
        !     The column is completely zero - do nothing
                     IF( INFO.EQ.0 ) THEN
                        INFO = K - 1
                     END IF
                     KP = KK
                  END IF
   
        !     swap rows and columns K+1 and KP in the
        !     full matrix A(1:N,1:N)
        !     Also, swap the first K-1 columns of W
   
                  IF( KP .NE. KK ) THEN
                     CALL SSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ),1)
                     CALL SSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
   
                     CALL SSWAP( N-K+1, A( KK, K), LDA, A( KP, K), LDA)
   
                     CALL SSCAL(KK-KP, -ONE, A(KP, KK), 1)
                     CALL SSCAL(KK-KP-1, -ONE, A(KP, KP+1), LDA)
   
        !     Swap in the last columns of W
                     IF( WK .GT. 0 ) THEN
                        CALL SSWAP( WK, W( KK, NB-WK+1 ), LDW, W( KP, NB-WK+1 ), LDW)
                     END IF
                  END IF
   
        !     (The column/row K+1 is not affected by the update)
                  IF( COLMAX .NE. ZERO ) THEN
        !     Store L(k+1) in A(k)
                     CALL SSCAL(K-2, ONE/A( K-1, K ), A(1, K), 1)
                  END IF
   
        !     Delay the update of the trailing submatrix (done at the beginning
        !     of each loop for every column of the first NB columns, and then
        !     finally in a rank-3 update for the last N-NB block)
        
        !     Store Pivot
                  IPIV( K-1 ) = KP
   
               ELSE
        !     STEP == 2 and an even column, do nothing
                  IPIV( K-1 ) = K-1
               END IF
        10      CONTINUE
   
        !     Now update the leading A(1:N-NB, 1:N-NB) submatrix in a level 3 update,
        !     if necessary
            IF( N-NPANEL+1 .GT. 2 ) THEN
   
        !     For this we have to set the N-NB,N-NB+1 entry to zero,
        !     but restore it later
               T = A( N-NPANEL, N-NPANEL+1 )
               A( N-NPANEL, N-NPANEL+1 ) = ZERO
   
               IF( WK .LT. NB) THEN
        !     Store the column N-NB in W, and update it in W
                  CALL SCOPY(N-NPANEL, A(1, N-NPANEL), 1, W(1, 1), 1)
   
                  W( N-NPANEL, 1 ) = ZERO
                  CALL SGEMV( 'N', N-NPANEL, WK, +ONE, A( 1, N-(WK-1)*STEP ), LDA*STEP,&
                  & W( N-NPANEL, NB-WK+1 ), LDW, ONE, W( 1, 1 ), 1 )
                  CALL SGEMV( 'N', N-NPANEL, WK, -ONE, W( 1, NB-WK+1 ), LDW, A( N-NPANEL, N-(WK-1)*STEP ), &
                  &LDA*STEP, ONE, W( 1, 1 ), 1 )
                  W( N-NPANEL, 1 ) = ZERO
   
                  WK = WK + 1
               END IF
   
        !     Now do the rank-3 update
               CALL SSKR2K( UPLO, "N", N-NPANEL, NB, ONE, A(1, N-(WK-1)*STEP), LDA*STEP, W(1,1), LDW,&
               & ONE, A(1, 1), LDA)
   
               A( N-NPANEL, N-NPANEL+1 ) = T
            END IF
   
         ELSE
   
        !     Factorize A as L * T * L^T using the lower triangle of A
        
        !     Consider the first NB columns of the matrix:
        !     First compute the transformations and apply them to the first NB
        !     columns (only), then apply the accumulated transformations to the rest
        !     of the matrix
   
            WK = 0
            DO 30 K=1, MIN(NPANEL, N-1)
   
        !
        !     Update A(K+1:n,K+1) with all the accumulated transformations
        !     (if K>2: K=1,2 is updated when computing the Gauss vector,
        !     K=2 is not affected by the K=1 transform)
        !
               KK = K+1
   
               IF( K .GT. 1) THEN
   
                  IF( WK .GT. 0) THEN
                     A( K, K ) = ZERO
                     CALL SGEMV( 'N', N-K+1, WK, +ONE, A( K, 1 ), LDA*STEP, W( K, 1 ), LDW, ONE, &
                     &A( K, K ), 1 )
                     CALL SGEMV( 'N', N-K+1, WK, -ONE, W( K, 1 ), LDW, A( K, 1 ), LDA*STEP, ONE, &
                     &A( K, K ), 1 )
                     A( K, K ) = ZERO
                  END IF
   
        !     Store the (updated) column K in W(:,WK)
                  IF( MOD(K, STEP) .EQ. 0 ) THEN
                     WK = WK + 1
                     CALL SCOPY(N-K+1, A(K, K), 1, W(K, WK), 1)
                  END IF
               END IF
   
               IF( MOD(K, STEP) .EQ. 1 .OR. STEP .EQ. 1) THEN
        !     For STEP == 1, process every column, but if
        !     STEP == 2, do only things for the odd columns
        
        !     Find the pivot
                  KP = K + ISAMAX(N-K, A( K+1, K ), 1)
                  COLMAX = ABS( A( KP, K ) )
   
                  IF( COLMAX.EQ.ZERO ) THEN
        !     The column is completely zero - do nothing
                     IF( INFO.EQ.0 ) THEN
                        INFO = K
                     END IF
                     KP = KK
                  END IF
   
        !     swap rows and columns K+1 and KP in the
        !     full matrix A(1:N,1:N)
        !     Also, swap the first K-1 columns of W
   
                  IF( KP .NE. KK ) THEN
                     IF( KP.LT.N ) THEN
                        CALL SSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ),1 )
                     END IF
   
                     CALL SSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )
   
                     CALL SSWAP( K, A( KK, 1), LDA, A( KP, 1), LDA)
   
        !     The matrix is anti-symmetric, hence swaps beyond the diagonal need a -1
   
                     CALL SSCAL( KP-KK, -ONE, A(KK+1, KK), 1 )
                     CALL SSCAL( KP-KK-1, -ONE, A(KP, KK+1), LDA )
   
        !     Swap in the first columns of W
                     CALL SSWAP( WK, W( KK, 1 ), LDW, W( KP, 1 ), LDW)
                  END IF
   
        !     (The column/row K+1 is not affected by the update)
                  IF( COLMAX .NE. ZERO .AND. K .LE. N-2) THEN
        !     Store L(k+1) in A(k)
                     CALL SSCAL( N-K-1, ONE/A( K+1, K ), A(K+2, K), 1 )
                  END IF
   
        !     Delay the update of the trailing submatrix (done at the beginning
        !     of each loop for every column of the first NB columns, and then
        !     finally in a rank-3 update for the last N-NB block)
        
        !     Store Pivot
                  IPIV( K+1 ) = KP
   
               ELSE
        !     STEP == 2 and an even column, do nothing
                  IPIV(K+1) = K+1
               END IF
        30      CONTINUE
   
        !     Now update the trailing A(NB+1:N, NB+1:N) submatrix in a level 3 update,
        !     if necessary
   
            IF( NPANEL .LT. N-1) THEN
        !     For this we have to set the NB+1,NB entry to zero, but restore it later
               T = A( NPANEL+1, NPANEL )
               A( NPANEL+1, NPANEL ) = ZERO
   
               IF( WK .LT. NB) THEN
        !     Store the column NB+1 in W, and update it in W
                  CALL SCOPY(N-NPANEL, A(NPANEL+1, NPANEL+1), 1, W(NPANEL+1, NB), 1)
   
                  W( NPANEL+1, NB ) = ZERO
                  CALL SGEMV( 'N', N-NPANEL, NB-1, +ONE, A( NPANEL+1, 1 ), LDA*STEP, W( NPANEL+1, 1 ), &
                  &LDW, ONE, W( NPANEL+1, NB ), 1 )
                  CALL SGEMV( 'N', N-NPANEL, NB-1, -ONE, W( NPANEL+1, 1 ), LDW, A( NPANEL+1, 1 ), &
                  &LDA*STEP, ONE, W( NPANEL+1, NB ), 1 )
                  W( NPANEL+1, NB ) = ZERO
               END IF
   
        !     Now do the rank-3 update
               CALL SSKR2K( UPLO, "N", N-NPANEL, NB, ONE, A(NPANEL+1,1), LDA*STEP, W(NPANEL+1,1), LDW, &
               & ONE, A(NPANEL+1, NPANEL+1), LDA)
   
               A(NPANEL+1, NPANEL)=T
            END IF
   
         END IF
   
    end subroutine SLASKTRF

    subroutine SSKTD2(UPLO, MODE, N, A, LDA, E, TAU, INFO )
        implicit none
        !  -- Written on 10/22/2010
        !     Michael Wimmer, Universiteit Leiden
        !
        !  Based on  the LAPACK routine ZHETD2 (www.netlib.org)
        !  -- LAPACK routine (version 3.2) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2006
        !
        !     .. Scalar Arguments ..
        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
        REAL               E(*)
        REAL               A(LDA,*), TAU(*)
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKTD2 reduces a real skew-symmetric matrix A to skew-symmetric
        !  tridiagonal form T by an orthognal similarity transformation:
        !  Q^T * A * Q = T. Alternatively, the routine can also compute
        !  a partial tridiagonal form useful for computing the Pfaffian.
        !
        !  This routine uses unblocked code.
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          Specifies whether the upper or lower triangular part of the
        !          skew-symmetric matrix A is stored:
        !          = 'U':  Upper triangular
        !          = 'L':  Lower triangular
        !
        !  MODE    (input) CHARACTER*1
        !          = 'N':  A is fully tridiagonalized
        !          = 'P':  A is partially tridiagonalized for Pfaffian computation
        !                  (details see below)
        !
        !  N       (input) INTEGER
        !          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the skew-symmetric matrix A.
        !            If UPLO = 'U', the leading N-by-N upper triangular part
        !              of A contains the upper triangular part of the matrix A,
        !              and the strictly lower triangular part of A is not referenced.
        !            If UPLO = 'L', the leading N-by-N lower triangular part
        !              of A contains the lower triangular part of the matrix A,
        !              and the strictly upper triangular part of A is not referenced.
        !          On exit, if MODE = 'N':
        !            If UPLO = 'U', the diagonal and first superdiagonal
        !              of A are overwritten by the corresponding elements of the
        !              tridiagonal matrix T, and the elements above the first
        !              superdiagonal, with the array TAU, represent the unitary
        !              matrix Q as a product of elementary reflectors;
        !            If UPLO = 'L', the diagonal and first subdiagonal of A are over-
        !              written by the corresponding elements of the tridiagonal
        !              matrix T, and the elements below the first subdiagonal, with
        !              the array TAU, represent the unitary matrix Q as a product
        !              of elementary reflectors.
        !            See Further Details, also for information about MODE = 'P'.
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        !  E       (output) REAL array, dimension (N-1)
        !          The off-diagonal elements of the tridiagonal matrix T:
        !          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
        !          If MODE = 'P', only the entries at i odd are well-defined
        !          (see Further Details)
        !
        !  TAU     (output) REAL array, dimension (N-1)
        !          The scalar factors of the elementary reflectors (see Further
        !          Details).
        !
        !  INFO    (output) INTEGER
        !          = 0:  successful exit
        !          < 0:  if INFO = -i, the i-th argument had an illegal value.
        !
        !  Further Details
        !  ===============
        !
        !  The normal use for SSKTD2 is to compute the tridiagonal form of
        !  a skew-symmetric matrix under an orthogonal similarity transformation,
        !  and chosen by setting MODE = 'N' ("normal" mode). The other
        !  use of SSKTD2 is the computation the Pfaffian of a skew-symmetric matrix,
        !  which only requires a partial tridiagonalization, this mode is chosen
        !  by setting MODE = 'P' ("Pfaffian" mode).
        !
        !  Normal mode (MODE = 'N'):
        !  ========================
        !
        !  The routine computes a tridiagonal matrix T and an orthogonal Q such
        !  that A = Q * T * Q^T .
        !
        !  If UPLO = 'U', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(n-1) . . . H(2) H(1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v'
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
        !  A(1:i-1,i+1), and tau in TAU(i).
        !
        !  If UPLO = 'L', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(1) H(2) . . . H(n-1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v'
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
        !  and tau in TAU(i).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 5:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   v2  v3  v4 )              (  0                  )
        !    (      0   e   v3  v4 )              (  e   0              )
        !    (          0   e   v4 )              (  v1  e   0          )
        !    (              0   e  )              (  v1  v2  e   0      )
        !    (                  0  )              (  v1  v2  v3  e   0  )
        !
        !  where d and e denote diagonal and off-diagonal elements of T, and vi
        !  denotes an element of the vector defining H(i).
        !
        !  The LAPACK routine SORGTR can be used to form the transformation
        !  matrix explicitely, and SORMTR can be used to multiply another
        !  matrix without forming the transformation.
        !
        !  Pfaffian mode (MODE = 'P'):
        !  ==========================
        !
        !  For computing the Pfaffian, it is enough to bring A into a partial
        !  tridiagonal form. In particular, assuming n even, it is enough to
        !  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
        !  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
        !  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
        !  only the off-diagonal entries in the odd columns (if UPLO = 'L')
        !  or in the even columns (if UPLU = 'U') are properly computed by SSKTD2.
        !
        !  A is brought into this special form pT using an orthogonal matrix Q:
        !  A = Q * pT * Q^T
        !
        !  If UPLO = 'U', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(n-1) H(n-3) . . . H(3) H(1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v^T
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
        !  A(1:i-1,i+1), and tau in TAU(i).
        !
        !  If UPLO = 'L', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(1) H(3) . . . H(n-3) H(n-1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v^T
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
        !  and tau in TAU(i).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 6:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   x   v3  x   v5 )        (  0                      )
        !    (      0   x   v3  x   v5 )        (  e   0                  )
        !    (          0   e   x   v5 )        (  v1  x   0              )
        !    (              0   x   v5 )        (  v1  x   e   0          )
        !    (                  0   e  )        (  v1  x   v3  x   0      )
        !    (                      0  )        (  v1  x   v3  x   e   0  )
        !
        !  where d and e denote diagonal and off-diagonal elements of T, vi
        !  denotes an element of the vector defining H(i), and x denotes an
        !  element not computed by SSKTD2.
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
        !     ..
        !     .. Local Scalars ..
        LOGICAL            UPPER, NORMAL
        INTEGER            I, STEP
        REAL               ALPHA, TAUI
        !     ..
        !     .. External Subroutines ..
        EXTERNAL           XERBLA, SLARFG
        !     ..
        !     .. External Functions ..
        LOGICAL            LSAME
        EXTERNAL           LSAME
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX, MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
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
            CALL XERBLA( 'SSKTD2', -INFO )
            RETURN
        END IF
        !
        !     Quick return if possible
        !
        IF( N.LE.0 ) RETURN
        !
        IF( .NOT. NORMAL ) THEN
            STEP = 2
        !     Make sure that all elements of TAU are initialized (only the
        !     odd elements are set below).
            DO 5 I = 2, N-2, 2
                TAU( I ) = ZERO
        5       CONTINUE
        ELSE
            STEP = 1
        END IF

        IF( UPPER ) THEN
        !
        !        Reduce the upper triangle of A
        !
            A( N, N ) = ZERO
            DO 10 I = N - 1, 1, -STEP
        !
        !           Generate elementary reflector H(i) = I - tau * v * v'
        !           to annihilate A(1:i-1,i+1)
        !
                ALPHA = A( I, I+1 )
                CALL SLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
                E( I ) = ALPHA
        !
                IF( TAUI.NE.ZERO ) THEN
        !
        !              Apply H(i) from both sides to A(1:i-step+1,1:i-step+1)
        !
                A( I, I+1 ) = ONE
        !
        !              Compute  x := tau * A * v  storing x in TAU(1:i)

                CALL SSKMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ),1, ZERO, TAU, 1 )
        !
        !              Apply the transformation as a rank-2 update:
        !                 A := A + v * w^T - w * v^T
        !
                CALL SSKR2( UPLO, I-STEP+1, ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA )
        !
                ELSE
                A( I, I ) = ZERO
                END IF
                A( I, I+1 ) = E( I )
                TAU( I ) = TAUI
        10    CONTINUE
        ELSE
        !
        !        Reduce the lower triangle of A
        !
            A( 1, 1 ) = ZERO
            DO 20 I = 1, N - 1, STEP
        !
        !           Generate elementary reflector H(i) = I - tau * v * v'
        !           to annihilate A(i+2:n,i)
        !
                ALPHA = A( I+1, I )
                CALL SLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
                E( I ) = ALPHA
        !
                IF( TAUI.NE.ZERO ) THEN
        !
        !              Apply H(i) from both sides to A(i+step:n,i+step:n)
        !
                A( I+1, I ) = ONE
        !
        !              Compute  x := tau^* * A * v^*  storing y in TAU(i:n-1)
        !
                CALL SSKMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 )
        !
        !              Apply the transformation as a rank-2 update:
        !                 A := A + v * x^T - x * v^T
        !
                CALL SSKR2( UPLO, N-I-STEP+1, ONE, A( I+STEP, I ), 1, TAU( I+STEP-1 ), 1, &
                &A( I+STEP, I+STEP ), LDA )
        !
                ELSE
                A( I+1, I+1 ) = ZERO
                END IF
                A( I+1, I ) = E( I )
                TAU( I ) = TAUI
        20    CONTINUE
        END IF
        !
        RETURN
    end subroutine SSKTD2

    subroutine SLASKTRD(UPLO, MODE, N, NB, A, LDA, E, TAU, W, LDW)
        implicit none
        !  -- Written on 10/22/2010
        !     Michael Wimmer, Universiteit Leiden
        !     Derived from the LAPACK routine ZLATRD (www.netlib.org)
        !  -- LAPACK auxiliary routine (version 3.2) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2006
        !
        !     .. Scalar Arguments ..
        CHARACTER          UPLO, MODE
        INTEGER            LDA, LDW, N, NB
        !     ..
        !     .. Array Arguments ..
        REAL               E( * )
        REAL               A( LDA, * ), TAU( * ), W( LDW, * )
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SLASKTRD reduces NB rows and columns of a real skew-symmetric matrix A to
        !  skew-symmetric tridiagonal form by an orthogonal similarity
        !  transformation Q^T * A * Q, and returns the matrices V and W which are
        !  needed to apply the transformation to the unreduced part of A.
        !
        !  If UPLO = 'U', SLASKTRD reduces the last NB rows and columns of a
        !  matrix, of which the upper triangle is supplied;
        !  if UPLO = 'L', SLASKTRD reduces the first NB rows and columns of a
        !  matrix, of which the lower triangle is supplied.
        !
        !  Alternatively, the routine can also be used to compute a partial
        !  tridiagonal form
        !
        !  This is an auxiliary routine called by SSKTRD.
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          Specifies whether the upper or lower triangular part of the
        !          skew-symmetric matrix A is stored:
        !          = 'U': Upper triangular
        !          = 'L': Lower triangular
        !
        !  MODE    (input) CHARACTER*1
        !          = 'N':  A is fully tridiagonalized
        !          = 'P':  A is partially tridiagonalized for Pfaffian computation
        !
        !  N       (input) INTEGER
        !          The order of the matrix A. N must be even if MODE = 'P'.
        !
        !  NB      (input) INTEGER
        !          The number of rows and columns to be reduced.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the skew-symmetric matrix A.
        !            If UPLO = 'U', the leading n-by-n upper triangular part
        !               of A contains the upper triangular part of the matrix A,
        !               and the strictly lower triangular part of A is not referenced.
        !            If UPLO = 'L', the leading n-by-n lower triangular part
        !               of A contains the lower triangular part of the matrix A,
        !               and the strictly upper triangular part of A is not referenced.
        !          On exit:
        !          if UPLO = 'U', the last NB columns have been reduced to
        !            tridiagonal form, with the diagonal elements overwriting
        !            the diagonal elements of A; the elements above the diagonal
        !            with the array TAU, represent the unitary matrix Q as a
        !            product of elementary reflectors. If MODE = 'P' only the
        !            even columns of the last NB columns contain meaningful values.
        !          if UPLO = 'L', the first NB columns have been reduced to
        !            tridiagonal form, with the diagonal elements overwriting
        !            the diagonal elements of A; the elements below the diagonal
        !            with the array TAU, represent the  unitary matrix Q as a
        !            product of elementary reflectors. If MODE = 'P' only the
        !            odd columns of the first NB columns contain meaningful values.
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        !  E       (output) REAL array, dimension (N-1)
        !          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
        !          elements of the last NB columns of the reduced matrix;
        !          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
        !          the first NB columns of the reduced matrix.
        !
        !  TAU     (output) REAL array, dimension (N-1)
        !          The scalar factors of the elementary reflectors, stored in
        !          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
        !          See Further Details.
        !
        !  W       (output) REAL array, dimension (LDW,NB)
        !          The n-by-nb matrix W required to update the unreduced part
        !          of A.
        !
        !  LDW     (input) INTEGER
        !          The leading dimension of the array W. LDW >= max(1,N).
        !
        !  =====================================================================
        !
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0)
        !     ..
        !     .. Local Scalars ..
        INTEGER            I, NW, NW2, STEP, NPANEL
        REAL               ALPHA
        !     ..
        !     .. External Subroutines ..
        EXTERNAL           SGEMV, SLARFG
        !     ..
        !     .. External Functions ..
        LOGICAL            LSAME
        EXTERNAL           LSAME
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MIN
        !     ..
        !     .. Executable Statements ..
        !
        !     Quick return if possible
        !
        IF( N.LE.0 ) RETURN
        !

        IF( LSAME( MODE, 'P' ) ) THEN
            STEP = 2
        ELSE
            STEP = 1
        END IF
        NPANEL = NB * STEP

        IF( LSAME( UPLO, 'U' ) ) THEN
        !
        !        Reduce last NPANEL columns of upper triangle
        !
            NW=0

            DO 10 I = N, MAX(N - NPANEL + 1, 2), -1

                NW2 = NW - MOD(I,STEP)
                IF( NW2 .GT. 0 ) THEN
        !
        !              Update A(1:i,i)
        !
                A( I, I ) = ZERO
                CALL SGEMV( 'No transpose', I, NW2, +ONE, A( 1, N-(NW2-1)*STEP ), LDA*STEP, &
                &W( I, NB-NW2+1 ), LDW, ONE, A( 1, I ), 1 )
                CALL SGEMV( 'No transpose', I, NW2, -ONE, W( 1, NB-NW2+1 ), LDW, A( I, N-(NW2-1)*STEP ), &
                &LDA*STEP, ONE, A( 1, I ), 1 )
                A( I, I ) = ZERO
                END IF

        !     In the Pfaffian mode, only zero all even columns
                IF( STEP.EQ.2 .AND. MOD( I, STEP ).EQ.1 ) THEN
                TAU( I-1 ) = ZERO
                GOTO 10
                END IF

                IF( I.GT.1 ) THEN
        !
        !              Generate elementary reflector H(i) to annihilate
        !              A(1:i-2,i)
        !
                ALPHA = A( I-1, I )
                CALL SLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
                E( I-1 ) = ALPHA
                A( I-1, I ) = ONE
        !
        !              Compute W(1:i-1,i)
        !

                CALL SSKMV( 'Upper', I-1, TAU(I-1), A, LDA, A( 1, I ), 1, ZERO, W( 1, NB-NW ), 1 )
                IF( NW .GT. 0 ) THEN
                    CALL SGEMV( 'Transpose', I-1, NW, ONE, W( 1, NB-NW+1 ), LDW, A( 1, I ), 1, &
                    &ZERO, W( I+1, NB-NW ), 1 )
                    CALL SGEMV( 'No transpose', I-1, NW, TAU(I-1), A( 1, N-(NW-1)*STEP ), LDA*STEP, &
                    &W( I+1, NB-NW ), 1, ONE, W( 1, NB-NW ), 1 )
                    CALL SGEMV( 'Transpose', I-1, NW, ONE, A( 1, N-(NW-1)*STEP ), LDA*STEP, A( 1, I ), &
                    &1, ZERO, W( I+1, NB-NW ), 1 )
                    CALL SGEMV( 'No transpose', I-1, NW, -TAU(I-1), W( 1, NB-NW+1 ), LDW, W( I+1, NB-NW ), &
                    & 1, ONE, W( 1, NB-NW ), 1 )
                END IF

        !     One more complete entry in W
                NW = NW + 1

                END IF

        !     Note: setting A(I-1,I) back to alpha happens in the calling routine

        10    CONTINUE
        ELSE
        !
        !        Reduce first NPANEL columns of lower triangle
        !
            NW = 0

            DO 20 I = 1, MIN(NPANEL, N-1)
        !
        !           Update A(i:n,i)
        !
                NW2 = NW - MOD(I+1,STEP)
                IF( NW2 .GT. 0 ) THEN
                A( I, I ) = ZERO
                CALL SGEMV( 'No transpose', N-I+1, NW2, +ONE, A( I, 1 ), LDA*STEP, W( I, 1 ), LDW, &
                &ONE, A( I, I ), 1 )
                CALL SGEMV( 'No transpose', N-I+1, NW2, -ONE, W( I, 1 ), LDW, A( I, 1 ), LDA*STEP, &
                &ONE, A( I, I ), 1 )
                A( I, I ) = ZERO
                END IF

        !     In the Pfaffian mode, only zero all odd columns
                IF( STEP.EQ.2 .AND. MOD( I, STEP ).EQ.0 ) THEN
                TAU( I ) = ZERO
                GOTO 20
                END IF

                IF( I.LT.N ) THEN
        !
        !              Generate elementary reflector H(i) to annihilate
        !              A(i+2:n,i)
        !
                ALPHA = A( I+1, I )
                CALL SLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAU( I ) )
                E( I ) = ALPHA
                A( I+1, I ) = ONE
        !
        !              Compute W(i+1:n,i)
        !              This is given by tau A^(i)*v^*=tau(A*v^* + VW^T v^* - WV^T v^*)


                CALL SSKMV( 'Lower', N-I, TAU( I ), A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, &
                &W( I+1, NW+1 ), 1 )
                IF( NW .GT. 0 ) THEN
                    CALL SGEMV( 'Transpose', N-I, NW, ONE, W( I+1, 1 ), LDW, A( I+1, I ), 1, &
                    &ZERO, W( 1, NW+1 ), 1 )
                    CALL SGEMV( 'No transpose', N-I, NW, TAU( I ), A( I+1, 1 ), LDA*STEP, &
                    &W( 1, NW+1 ), 1, ONE, W( I+1, NW+1 ), 1 )
                    CALL SGEMV( 'Transpose', N-I, NW, ONE, A( I+1, 1 ), LDA*STEP, A( I+1, I ), &
                    &1, ZERO, W( 1, NW+1 ), 1 )
                    CALL SGEMV( 'No transpose', N-I, NW, -TAU( I ), W( I+1, 1 ), LDW, W( 1, NW+1 ), &
                    &1, ONE, W( I+1, NW+1 ), 1 )
                END IF

        !     One more complete entry in W
                NW = NW + 1
                END IF
        !
        20    CONTINUE
        END IF
        !
        RETURN
    end subroutine SLASKTRD

    subroutine SSKTF2(UPLO, MODE, N, A, LDA, IPIV, INFO )
        implicit none
        !
        !     .. Scalar Arguments ..
        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, N
        !     ..
        !     .. Array Arguments ..
        INTEGER            IPIV( * )
        REAL               A( LDA, * )
        !
        !
        !  Purpose
        !  =======
        !
        !  SSKTF2 computes the factorization of a skew-symmetric matrix A
        !  using the Parlett-Reid algorithm:
        !
        !     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
        !
        !  where U (or L) unit upper (lower) triangular matrix (^T denotes
        !  the transpose), T is a skew-symmetric tridiagonal matrix and P
        !  is a permutation matrix. In addition to being unit triangular,
        !  U(1:n-1,n)=0 and L(2:n,1)=0.
        !  Instead of a full tridiagonalization, DFSKTF2 can also compute a
        !  partial tridiagonal form for computing the Pfaffian.
        !
        !  This is the unblocked version of the algorithm, calling Level 2 BLAS.
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          Specifies whether the upper or lower triangular part of the
        !          symmetric matrix A is stored:
        !          = 'U':  Upper triangular
        !          = 'L':  Lower triangular
        !
        !  MODE    (input) CHARACTER*1
        !          = 'N':  A is fully tridiagonalized
        !          = 'P':  A is partially tridiagonalized for Pfaffian computation
        !                  (details see below)
        !
        !  N       (input) INTEGER
        !          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the symmetric matrix A.
        !             If UPLO = 'U', the leading n-by-n upper triangular part
        !                of A contains the upper triangular part of the matrix A,
        !                and the strictly lower triangular part of A is not referenced.
        !             If UPLO = 'L', the leading n-by-n lower triangular part
        !                of A contains the lower triangular part of the matrix A,
        !                and the strictly upper triangular part of A is not referenced.
        !
        !          On exit, the tridiagonal matrix T and the multipliers used
        !          to obtain the factor U or L (see below for further details).
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        !  IPIV    (output) INTEGER array, dimension (N)
        !          Information about the permutation matrix P: row and column
        !          i are interchanged with IPIV(i). If UPLO = 'U', those
        !          interchanges are done in the order i = N ... 1, if UPLO = 'L'
        !          in the order i = 1 ... N.
        !
        !  INFO    (output) INTEGER
        !          = 0: successful exit
        !          < 0: if INFO = -k, the k-th argument had an illegal value
        !          > 0: if INFO = k, the off-diagonal entry in the k-th row
        !                            (UPLO = 'U') or k-th column (UPLO = 'L')
        !                            is exactly zero.
        !
        !  Further Details
        !  ===============
        !
        !
        !  The normal use for SSKTD2 is to compute the U T U^T or L T L^T
        !  decomposition of a skew-symmetric matrix with pivoting. This mode
        !  is chosen by setting MODE = 'N' ("normal" mode). The other
        !  use of SSKTD2 is the computation the Pfaffian of a skew-symmetric matrix,
        !  which only requires a partial computation of T, this mode is chosen
        !  by setting MODE = 'P' ("Pfaffian" mode).
        !
        !  Normal mode (MODE = 'N'):
        !  ========================
        !
        !  If UPLO = 'U', the U*T*U^T decomposition of A is computed. U is a
        !  upper triangular unit matrix with the additional constraint
        !  U(1:n-1,n) = 0, and T a tridiagonal matrix. The upper diagonal
        !  of T is stored on exit in A(i,i+1) for i = 1 .. n-1. The column
        !  U(1:i-1, i) is stored in A(1:i-1,i+1).
        !
        !  If UPLO = 'L', the L*T*L^T decomposition of A is computed. L is a
        !  lower triangular unit matrix with the additional constraint
        !  L(2:n,1) = 0, and T a tridiagonal matrix. The lower diagonal
        !  of T is stored on exit in A(i+1,i) for i = 1 .. n-1. The column
        !  L(i+1:n, i) is stored in A(i+1:n,i-1).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 5:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   u2  u3  u4 )              (  0                  )
        !    (      0   e   u3  u4 )              (  e   0              )
        !    (          0   e   u4 )              (  l2  e   0          )
        !    (              0   e  )              (  l2  l3  e   0      )
        !    (                  0  )              (  l2  l3  l4  e   0  )
        !
        !  where e denotes the off-diagonal elements of T, and ui (li)
        !  denotes an element of the i-th column of U (L).
        !
        !  Pfaffian mode (MODE = 'P'):
        !  ==========================
        !
        !  For computing the Pfaffian, it is enough to bring A into a partial
        !  tridiagonal form. In particular, assuming n even, it is enough to
        !  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
        !  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
        !  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
        !  only the off-diagonal entries in the odd columns (if UPLO = 'L')
        !  or in the even columns (if UPLU = 'U') are properly computed by SSKTF2.
        !
        !  If UPLO = 'U', the U*pT*U^T decomposition of A is computed. U is a
        !  upper triangular unit matrix with the additional constraint
        !  U(1:i-1,i) = 0 for even i, and pT a partially tridiagonal matrix.
        !  The entries in the odd rows of the upper diagonal of pT are stored
        !  on exit in A(i,i+1) for i odd. The column U(1:i-1, i) for odd i
        !  is stored in A(1:i-1,i+1).
        !
        !  If UPLO = 'L', the L*pT*L^T decomposition of A is computed. L is a
        !  lower triangular unit matrix with the additional constraint
        !  L(i+1:n,i) = 0 for odd i, and pT a partially tridiagonal matrix.
        !  The entries in odd columns in the lower diagonal of pT are stored
        !  on exit in A(i+1,i) for i odd. The column L(i+1:n, i) for i odd
        !  is stored in A(i+1:n,i-1).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 6:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   x   u3  x   u5 )              (  0                    )
        !    (      0   x   u3  x   u5 )              (  e   0                )
        !    (          0   e   x   u5 )              (  l2  x   0            )
        !    (              0   x   u5 )              (  l2  x   e   0        )
        !    (                  0   e  )              (  l2  x   l4  x   0    )
        !    (                      0  )              (  l2  x   l4  x   e  0 )
        !
        !  where e denotes the off-diagonal elements of T, ui (li)
        !  denotes an element of the i-th column of U (L), and x denotes an
        !  element not computed by SSKTF2.
        !
        !  =====================================================================
        !
    
        !     .. Parameters ..
        REAL               ZERO, ONE
        PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
        
        !     .. Local Scalars ..
        LOGICAL            UPPER, NORMAL
        INTEGER            K, KK, KP, STEP
        REAL               COLMAX
        !     ..
        !     .. External Functions ..
        LOGICAL            LSAME
        INTEGER            ISAMAX
        EXTERNAL           LSAME, ISAMAX
        !     ..
        !     .. External Subroutines ..
        EXTERNAL           SSCAL, SSWAP, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
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
        !     If STEP == 2, we need an even-dimensional matrix
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        END IF
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SSKTF2', -INFO )
            RETURN
        END IF
        
        !     Quick return if possible
        IF( N .EQ. 0 ) RETURN
        
        IF( NORMAL ) THEN
            STEP = 1
        ELSE
            STEP = 2
        END IF
        
        IF( UPPER ) THEN
        !     Factorize A as U * T * U^T using the upper triangle of A
            IPIV( N ) = N
        
            DO 10 K=N, 2, -1
        !     Either all columns or only the even ones (MODE = 'P')
                IF( MOD(K, STEP) .EQ. 0) THEN
        !     Find the pivot
                KP = ISAMAX(K-1, A( 1, K ), 1)
                COLMAX = ABS( A( KP, K ) )
        
                IF( COLMAX.EQ.ZERO ) THEN
        !     The column is completely zero - do nothing
                    IF( INFO.EQ.0 ) THEN
                        INFO = K - 1
                    END IF
                    KP = K-1
                END IF
        
        !     swap rows and columns K+1 and IMAX in the
        !     full sub-matrix A(1:N,1:N)
                       KK = K-1
        
                       IF( KP .NE. KK ) THEN
                          CALL SSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1)
                          CALL SSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )
        
                          CALL SSWAP( N-K+1, A( KK, K), LDA, A( KP, K), LDA)
        
                          CALL SSCAL(KK-KP, -ONE, A(KP, KK), 1)
                          CALL SSCAL(KK-KP-1, -ONE, A(KP, KP+1), LDA)
                       END IF
        
        !     Update the leading submatrix A(1:K-2, 1:K-2) in a rank 2 update
        !     (The column/row K-1 is not affected by the update)
                       IF( COLMAX .NE. ZERO ) THEN
                          CALL SSKR2( UPLO, K-2, ONE/A( K-1,K ), A( 1, K ), 1, A( 1, K-1 ), 1, A( 1, 1 ), LDA )
        
        !     Store L(k+1) in A(k)
                          CALL SSCAL(K-2, ONE/A( K-1, K ), A(1, K), 1)
                       END IF
        !     Store Pivot
                       IPIV( K-1 ) = KP
                    ELSE
                       IPIV( K-1 ) = K-1
                    END IF
        10      CONTINUE
        
              ELSE
        !     Factorize A as L * T * L^T using the lower triangle of A
        
                 IPIV( 1 ) = 1
        
                 DO 20 K=1, N-1
        !     Either all columns or only the odd ones (MODE = 'P')
                    IF( MOD(K, STEP).EQ.1 .OR. STEP.EQ.1 ) THEN
        !     Find the pivot
                       KP = K + ISAMAX(N-K, A( K+1, K ), 1)
                       COLMAX = ABS( A( KP, K ) )
        
                       IF( COLMAX.EQ.ZERO ) THEN
        !     The column is completely zero - do nothing
                          IF( INFO.EQ.0 ) THEN
                             INFO = K
                          END IF
                          KP = K+1
                       END IF
        
        !     swap rows and columns K+1 and IMAX in the
        !     full matrix A(1:N,1:N)
                       KK = K+1
        
                       IF( KP .NE. KK ) THEN
                          IF( KP.LT.N ) THEN
                             CALL SSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ),1 )
                          END IF
        
                          CALL SSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )
        
                          CALL SSWAP( K, A( KK, 1), LDA, A( KP, 1), LDA)
        
                          CALL SSCAL(KP-KK, -ONE, A(KK+1, KK), 1)
                          CALL SSCAL(KP-KK-1, -ONE, A(KP, KK+1), LDA)
                       END IF
        
        !     Update the trailing submatrix A(K+2:N, K+2:N) in a rank 2 update
        !     (The column/row K+1 is not affected by the update)
                       IF( COLMAX .NE. ZERO .AND. K+2 .LE. N) THEN
                          CALL SSKR2( UPLO, N-K-1, ONE/A( K+1,K ), A( K+2, K ), 1, A( K+2, K+1 ), &
                          &1, A( K+2, K+2 ), LDA )
        
        !     Store L(k+1) in A(k)
                          CALL SSCAL(N-K-1, ONE/A( K+1, K ), A(K+2, K), 1)
                       END IF
        
        !     Store Pivot
                       IPIV( K+1 ) = KP
                    ELSE
                       IPIV( K+1 ) = K+1
                    END IF
        20      CONTINUE
        
            END IF
        
    end subroutine SSKTF2

    subroutine SSKTRF(UPLO, MODE, N, A, LDA, IPIV, WORK, LWORK, INFO)
        implicit none
        !     .. Scalar Arguments ..
        CHARACTER          UPLO, MODE
        INTEGER            INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
        INTEGER            IPIV( * )
        REAL               A( LDA, * ), WORK( * )
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKTRF computes the factorization of a skew-symmetric matrix A
        !  using the Parlett-Reid algorithm:
        !
        !     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T
        !
        !  where U (or L) unit upper (lower) triangular matrix (^T denotes
        !  the transpose), T is a skew-symmetric tridiagonal matrix and P
        !  is a permutation matrix. In addition to being unit triangular,
        !  U(1:n-1,n)=0 and L(2:n,1)=0.
        !  Instead of a full tridiagonalization, SSKTRF can also compute a
        !  partial tridiagonal form for computing the Pfaffian.
        !
        !  This is the blocked version of the algorithm, calling Level 3 BLAS.
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          Specifies whether the upper or lower triangular part of the
        !          symmetric matrix A is stored:
        !          = 'U':  Upper triangular
        !          = 'L':  Lower triangular
        !
        !  MODE    (input) CHARACTER*1
        !          = 'N':  A is fully tridiagonalized
        !          = 'P':  A is partially tridiagonalized for Pfaffian computation
        !                  (details see below)
        !
        !  N       (input) INTEGER
        !          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the symmetric matrix A.
        !             If UPLO = 'U', the leading n-by-n upper triangular part
        !                of A contains the upper triangular part of the matrix A,
        !                and the strictly lower triangular part of A is not referenced.
        !             If UPLO = 'L', the leading n-by-n lower triangular part
        !                of A contains the lower triangular part of the matrix A,
        !                and the strictly upper triangular part of A is not referenced.
        !
        !          On exit, the tridiagonal matrix T and the multipliers used
        !          to obtain the factor U or L (see below for further details).
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        !  IPIV    (output) INTEGER array, dimension (N)
        !          Information about the permutation matrix P: row and column
        !          i are interchanged with IPIV(i). If UPLO = 'U', those
        !          interchanges are done in the order i = N ... 1, if UPLO = 'L'
        !          in the order i = 1 ... N.
        !
        !  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
        !          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
        !
        !  LWORK   (input) INTEGER
        !          The length of WORK.  LWORK >=1.  For best performance
        !          LWORK >= N*NB, where NB is the block size returned by ILAENV
        !          (at the moment, uses the same block size as SSYTRF from Lapack).
        !
        !          If LWORK = -1, then a workspace query is assumed; the routine
        !          only calculates the optimal size of the WORK array, returns
        !          this value as the first entry of the WORK array, and no error
        !          message related to LWORK is issued by XERBLA.
        !
        !  INFO    (output) INTEGER
        !          = 0: successful exit
        !          < 0: if INFO = -k, the k-th argument had an illegal value
        !          > 0: if INFO = k, the off-diagonal entry in the k-th row
        !                            (UPLO = 'U') or k-th column (UPLO = 'L')
        !                            is exactly zero.
        !
        !  Further Details
        !  ===============
        !
        !
        !  The normal use for SSKTRF is to compute the U T U^T or L T L^T
        !  decomposition of a skew-symmetric matrix with pivoting. This mode
        !  is chosen by setting MODE = 'N' ("normal" mode). The other
        !  use of SSKTRF is the computation the Pfaffian of a skew-symmetric matrix,
        !  which only requires a partial computation of T, this mode is chosen
        !  by setting MODE = 'P' ("Pfaffian" mode).
        !
        !  Normal mode (MODE = 'N'):
        !  ========================
        !
        !  If UPLO = 'U', the U*T*U^T decomposition of A is computed. U is a
        !  upper triangular unit matrix with the additional constraint
        !  U(1:n-1,n) = 0, and T a tridiagonal matrix. The upper diagonal
        !  of T is stored on exit in A(i,i+1) for i = 1 .. n-1. The column
        !  U(1:i-1, i) is stored in A(1:i-1,i+1).
        !
        !  If UPLO = 'L', the L*T*L^T decomposition of A is computed. L is a
        !  lower triangular unit matrix with the additional constraint
        !  L(2:n,1) = 0, and T a tridiagonal matrix. The lower diagonal
        !  of T is stored on exit in A(i+1,i) for i = 1 .. n-1. The column
        !  L(i+1:n, i) is stored in A(i+1:n,i-1).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 5:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   u2  u3  u4 )              (  0                  )
        !    (      0   e   u3  u4 )              (  e   0              )
        !    (          0   e   u4 )              (  l2  e   0          )
        !    (              0   e  )              (  l2  l3  e   0      )
        !    (                  0  )              (  l2  l3  l4  e   0  )
        !
        !  where e denotes the off-diagonal elements of T, and ui (li)
        !  denotes an element of the i-th column of U (L).
        !
        !  Pfaffian mode (MODE = 'P'):
        !  ==========================
        !
        !  For computing the Pfaffian, it is enough to bring A into a partial
        !  tridiagonal form. In particular, assuming n even, it is enough to
        !  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
        !  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
        !  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
        !  only the off-diagonal entries in the odd columns (if UPLO = 'L')
        !  or in the even columns (if UPLU = 'U') are properly computed by SSKTRF.
        !
        !  If UPLO = 'U', the U*pT*U^T decomposition of A is computed. U is a
        !  upper triangular unit matrix with the additional constraint
        !  U(1:i-1,i) = 0 for even i, and pT a partially tridiagonal matrix.
        !  The entries in the odd rows of the upper diagonal of pT are stored
        !  on exit in A(i,i+1) for i odd. The column U(1:i-1, i) for odd i
        !  is stored in A(1:i-1,i+1).
        !
        !  If UPLO = 'L', the L*pT*L^T decomposition of A is computed. L is a
        !  lower triangular unit matrix with the additional constraint
        !  L(i+1:n,i) = 0 for odd i, and pT a partially tridiagonal matrix.
        !  The entries in odd columns in the lower diagonal of pT are stored
        !  on exit in A(i+1,i) for i odd. The column L(i+1:n, i) for i odd
        !  is stored in A(i+1:n,i-1).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 6:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   x   u3  x   u5 )              (  0                    )
        !    (      0   x   u3  x   u5 )              (  e   0                )
        !    (          0   e   x   u5 )              (  l2  x   0            )
        !    (              0   x   u5 )              (  l2  x   e   0        )
        !    (                  0   e  )              (  l2  x   l4  x   0    )
        !    (                      0  )              (  l2  x   l4  x   e  0 )
        !
        !  where e denotes the off-diagonal elements of T, ui (li)
        !  denotes an element of the i-th column of U (L), and x denotes an
        !  element not computed by SSKTRF.
        !
        !  =====================================================================
        !
        !     .. Local Scalars ..
        LOGICAL            LQUERY, UPPER, NORMAL
        INTEGER            IINFO, J, K, K2, PIV, LWKOPT, NB, NBMIN, NPANEL
        !     ..
        !     .. External Functions ..
        LOGICAL            LSAME
        INTEGER            ILAENV
        EXTERNAL           LSAME, ILAENV
        !     ..
        !     .. External Subroutines ..
        EXTERNAL           SSWAP, XERBLA
        !     ..
        !     .. Intrinsic Functions ..
        INTRINSIC          MAX
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters.
        !
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
        !     If MODE = 'P', we need an even-dimensional matrix
            INFO = -3
        ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
            INFO = -5
        ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
            INFO = -8
        END IF
        !

        IF( INFO.EQ.0 ) THEN
        !
        !        Determine the block size
        !        Note: Obviously, since Lapack does not know about SSKTRF, I query the
        !        similar routine 'SSYTRF'
            NB = ILAENV( 1, 'SSYTRF', UPLO, N, -1, -1, -1 )
            LWKOPT = N*NB
            WORK( 1 ) = LWKOPT
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SSKTRF', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF
        !
        NBMIN=NB
        IF( NB.GT.1 .AND. NB.LT.N ) THEN
            IF( LWORK .LT. N*NB ) THEN
                NB = MAX( LWORK / N, 1 )
                NBMIN = MAX( 2, ILAENV( 2, 'SSYTRF', UPLO, N, -1, -1, -1 ) )
            END IF
        ELSE
            NB = N
        END IF

        IF( NB.LT.NBMIN ) NB = N
        !     Note: NB = N means do not use blocked code

        !     Quick return if possible
        IF( N .EQ. 0 ) RETURN

        IF( LSAME( MODE, 'N' ) ) THEN
            NPANEL = NB
        ELSE
            NPANEL = MIN(NB*2, N)
        END IF

        IF( UPPER ) THEN
        !
        !     Factorize A as L*T*L^T using the lower triangle of A

            IPIV( N ) = N
        !
        !     Loop throgh the system in steps of NPANEL
            DO 10 K = N, MAX(NPANEL, 1), -NPANEL
        !
                IF( K.GE.NPANEL*2 ) THEN
        !
        !     Factorize columns k-nb*step+1:k of A and use blocked code to
        !     update columns 1:k-nb*step
        !
                    CALL SLASKTRF( UPLO, MODE, K, NB, A, LDA, IPIV, WORK, N, IINFO )

                    K2 = K-NPANEL
                ELSE
        !
        !     Use unblocked code to factorize columns 1:k of A
        !
        !     IPIV( K ) is overwritten by SSKTF2, need to restore it later
                    PIV = IPIV( K )

                    CALL SSKTF2( UPLO, MODE, K, A, LDA, IPIV, IINFO )

                    IPIV( K ) = PIV

                    K2 = 1
                END IF
        !
        !        Set INFO on the first occurrence of a zero pivot
        !
                IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO

        !     Perform the missing row interchanges in the trailing N-K columns
                    IF( K .LT. N ) THEN
                        DO 20 J=K-1, K2, -1
                            CALL SSWAP( N-K, A( J, K+1 ), LDA, A( IPIV( J ), K+1 ), LDA )
            20            CONTINUE
                    END IF

            10   CONTINUE
        !
            ELSE
        !
        !        Factorize A as L*T*L^T using the lower triangle of A

                IPIV( 1 ) = 1
        !
        !        Loop throgh the system in steps of NPANEL
                DO 30 K = 1, MIN(N-NPANEL+1, N-1), NPANEL
        !
                    IF( K.LE.N-NPANEL*2+1 ) THEN
        !
        !     Factorize columns k:k+nb-1 of A and use blocked code to
        !     update columns k+nb-1:n
        !
                    CALL SLASKTRF( UPLO, MODE, N-K+1, NB, A( K, K ), LDA, IPIV( K ), WORK, N, IINFO )

                    K2 = K + NPANEL
                    ELSE
        !
        !     Use unblocked code to factorize columns k:n of A
        !
        !     IPIV( K ) is overwritten by SSKTF2, need to restore it later
                    PIV = IPIV( K )

                    CALL SSKTF2( UPLO, MODE, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )

                    IPIV( K ) = PIV

                    K2 = N
                    END IF
        !
        !     Set INFO on the first occurrence of a zero pivot
        !
                    IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + K - 1
        !
        !     Adjust IPIV
        !
                    DO 40 J = K+1, K2
                    IPIV( J ) = IPIV( J ) + K - 1
            40         CONTINUE

        !     Perform the missing row interchanges in the leading K-1 columns
                    IF( K .GT. 1 ) THEN
                    DO 50 J=K+1, K2
                        CALL SSWAP( K-1, A( J, 1 ), LDA, A( IPIV( J ), 1 ), LDA )
            50            CONTINUE
                    END IF


            30   CONTINUE
        !
            END IF
        !
            WORK( 1 ) = LWKOPT

            return

    end subroutine SSKTRF

    subroutine SSKTRD(UPLO, MODE, N, A, LDA, E, TAU, WORK, LWORK, INFO )
        implicit none
        
        !  -- Written on 10/22/2010
        !     Michael Wimmer
        !
        !  -- derived from LAPACK routine ZHETRD (version 3.2, www.netlib.org) --
        !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
        !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
        !     November 2006
        !
        !     .. Scalar Arguments ..
            CHARACTER          UPLO, MODE
            INTEGER            INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
            REAL               E( * )
            REAL               A( LDA, * ), TAU( * ), WORK( * )
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKTRD reduces a real skew-symmetric matrix A to skew-symmetric
        !  tridiagonal form T by an orthogonal similarity transformation:
        !  Q^T * A * Q = T.  Alternatively, the routine can also compute
        !  a partial tridiagonal form useful for computing the Pfaffian.
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          = 'U':  Upper triangle of A is stored;
        !          = 'L':  Lower triangle of A is stored.
        !  MODE    (input) CHARACTER*1
        !          = 'N':  A is fully tridiagonalized
        !          = 'P':  A is partially tridiagonalized for Pfaffian computation
        !                  (details see below)
        !
        !  N       (input) INTEGER
        !          The order of the matrix A.  N >= 0. N must be even if MODE = 'P'.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the skew-symmetric matrix A.
        !            If UPLO = 'U', the leading N-by-N upper triangular part
        !              of A contains the upper triangular part of the matrix A,
        !              and the strictly lower triangular part of A is not referenced.
        !            If UPLO = 'L', the leading N-by-N lower triangular part
        !              of A contains the lower triangular part of the matrix A,
        !              and the strictly upper triangular part of A is not referenced.
        !          On exit, if MODE = 'N':
        !            If UPLO = 'U', the diagonal and first superdiagonal
        !              of A are overwritten by the corresponding elements of the
        !              tridiagonal matrix T, and the elements above the first
        !              superdiagonal, with the array TAU, represent the unitary
        !              matrix Q as a product of elementary reflectors;
        !            If UPLO = 'L', the diagonal and first subdiagonal of A are over-
        !              written by the corresponding elements of the tridiagonal
        !              matrix T, and the elements below the first subdiagonal, with
        !              the array TAU, represent the unitary matrix Q as a product
        !              of elementary reflectors.
        !            See Further Details, also for information about MODE = 'P'.
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        ! E       (output) REAL array, dimension (N-1)
        !         The off-diagonal elements of the tridiagonal matrix T:
        !          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
        !         If MODE = 'P', only the entries at i odd are well-defined
        !          (see Further Details)
        !
        !  TAU     (output) REAL array, dimension (N-1)
        !          The scalar factors of the elementary reflectors (see Further
        !          Details).
        !
        !  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
        !          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
        !
        !  LWORK   (input) INTEGER
        !          The dimension of the array WORK.  LWORK >= 1.
        !          For optimum performance LWORK >= N*NB, where NB is the
        !          optimal blocksize.
        !
        !          If LWORK = -1, then a workspace query is assumed; the routine
        !          only calculates the optimal size of the WORK array, returns
        !          this value as the first entry of the WORK array, and no error
        !          message related to LWORK is issued by XERBLA.
        !
        !  INFO    (output) INTEGER
        !          = 0:  successful exit
        !          < 0:  if INFO = -i, the i-th argument had an illegal value
        !
        !  Further Details
        !  ===============
        !
        !  The normal use for SSKTRD is to compute the tridiagonal form of
        !  a skew-symmetric matrix under an orthogonal similarity transformation,
        !  and chosen by setting MODE = 'N' ("normal" mode). The other
        !  use of SSKTRD is the computation the Pfaffian of a skew-symmetric matrix,
        !  which only requires a partial tridiagonalization, this mode is chosen
        !  by setting MODE = 'P' ("Pfaffian" mode).
        !
        !  Normal mode (MODE = 'N'):
        !  ========================
        !
        !  The routine computes a tridiagonal matrix T and an orthogonal Q such
        !  that A = Q * T * Q^T .
        !
        !  If UPLO = 'U', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(n-1) . . . H(2) H(1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v'
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
        !  A(1:i-1,i+1), and tau in TAU(i).
        !
        !  If UPLO = 'L', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(1) H(2) . . . H(n-1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v'
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
        !  and tau in TAU(i).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 5:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   v2  v3  v4 )              (  0                  )
        !    (      0   e   v3  v4 )              (  e   0              )
        !    (          0   e   v4 )              (  v1  e   0          )
        !    (              0   e  )              (  v1  v2  e   0      )
        !    (                  0  )              (  v1  v2  v3  e   0  )
        !
        !  where d and e denote diagonal and off-diagonal elements of T, and vi
        !  denotes an element of the vector defining H(i).
        !
        !  The LAPACK routine DORGTR can be used to form the transformation
        !  matrix explicitely, and DORMTR can be used to multiply another
        !  matrix without forming the transformation.
        !
        !  Pfaffian mode (MODE = 'P'):
        !  ==========================
        !
        !  For computing the Pfaffian, it is enough to bring A into a partial
        !  tridiagonal form. In particular, assuming n even, it is enough to
        !  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
        !  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
        !  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
        !  only the off-diagonal entries in the odd columns (if UPLO = 'L')
        !  or in the even columns (if UPLU = 'U') are properly computed by SSKTRD.
        !
        !  A is brought into this special form pT using an orthogonal matrix Q:
        !  A = Q * pT * Q^T
        !
        !  If UPLO = 'U', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(n-1) H(n-3) . . . H(3) H(1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v^T
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
        !  A(1:i-1,i+1), and tau in TAU(i).
        !
        !  If UPLO = 'L', the matrix Q is represented as a product of elementary
        !  reflectors
        !
        !     Q = H(1) H(3) . . . H(n-3) H(n-1).
        !
        !  Each H(i) has the form
        !
        !     H(i) = I - tau * v * v^T
        !
        !  where tau is a real scalar, and v is a real vector with
        !  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
        !  and tau in TAU(i).
        !
        !  The contents of A on exit are illustrated by the following examples
        !  with n = 6:
        !
        !  if UPLO = 'U':                       if UPLO = 'L':
        !
        !    (  0   e   x   v3  x   v5 )        (  0                      )
        !    (      0   x   v3  x   v5 )        (  e   0                  )
        !    (          0   e   x   v5 )        (  v1  x   0              )
        !    (              0   x   v5 )        (  v1  x   e   0          )
        !    (                  0   e  )        (  v1  x   v3  x   0      )
        !    (                      0  )        (  v1  x   v3  x   e   0  )
        !
        !  where d and e denote diagonal and off-diagonal elements of T, vi
        !  denotes an element of the vector defining H(i), and x denotes an
        !  element not computed by SSKTRD.
        !
        !  =====================================================================
        !
        !     .. Parameters ..
            REAL               ONE
            PARAMETER          ( ONE = 1.0E+0 )
        !     ..
        !     .. Local Scalars ..
            LOGICAL            LQUERY, UPPER, NORMAL
            INTEGER            I, IINFO, IWS, J, LDWORK, LWKOPT, NB, NBMIN, NX, STEP, NPANEL, NXPANEL
        !     ..
        !     .. External Subroutines ..
            EXTERNAL           XERBLA
        !     ..
        !     .. Intrinsic Functions ..
            INTRINSIC          MAX
        !     ..
        !     .. External Functions ..
            LOGICAL            LSAME
            INTEGER            ILAENV
            EXTERNAL           LSAME, ILAENV
        !     ..
        !     .. Executable Statements ..
        !
        !     Test the input parameters
        !
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
        !
            IF( INFO.EQ.0 ) THEN
        !
        !        Determine the block size.
        !        NOTE: I keep 'SSYTRD" here, as of course ILAENV has no information
        !              about SSKTRD ...

                NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 )
                LWKOPT = N*NB
                WORK( 1 ) = LWKOPT
            END IF
        !
            IF( INFO.NE.0 ) THEN
                CALL XERBLA( 'SSKTRD', -INFO )
                RETURN
            ELSE IF( LQUERY ) THEN
                RETURN
            END IF
        !
        !     Quick return if possible
        !
            IF( N.EQ.0 ) THEN
                WORK( 1 ) = 1
                RETURN
            END IF
        !

            NX = N
            IWS = 1
            IF( NB.GT.1 .AND. NB.LT.N ) THEN
        !
        !        Determine when to cross over from blocked to unblocked code
        !        (last block is always handled by unblocked code).
        !
                NX = MAX( NB, ILAENV( 3, 'SSYTRD', UPLO, N, -1, -1, -1 ) )
                IF( NX.LT.N ) THEN
        !
        !           Determine if workspace is large enough for blocked code.
        !
                    LDWORK = N
                    IWS = LDWORK*NB
                    IF( LWORK.LT.IWS ) THEN
        !
        !              Not enough workspace to use optimal NB:  determine the
        !              minimum value of NB, and reduce NB or force use of
        !              unblocked code by setting NX = N.
        !
                    NB = MAX( LWORK / LDWORK, 1 )
                    NBMIN = ILAENV( 2, 'SSYTRD', UPLO, N, -1, -1, -1 )
                    IF( NB.LT.NBMIN .OR. NB.LE.1 ) NX = N
                    END IF
                ELSE
                    NX = N
                END IF
            ELSE
                NB = 1
            END IF
        !

            IF( .NOT.NORMAL ) THEN
                STEP = 2
            ELSE
                STEP = 1
            END IF

            NPANEL = NB * STEP
            NXPANEL = NX * STEP

            IF( UPPER ) THEN
        !
        !        Reduce the upper triangle of A.
        !        Columns 1:kk are handled by the unblocked method.
        !
                DO 20 I = N, NXPANEL + NPANEL, -NPANEL
        !
        !           Reduce columns i-npanel+1:i to tridiagonal form and form the
        !           matrix W which is needed to update the unreduced part of
        !           the matrix
        !
                    CALL SLASKTRD( UPLO, MODE, I, NB, A, LDA, E, TAU, WORK, LDWORK )
        !
        !           Update the unreduced submatrix A(1:i-npanel,1:i-npanel), using an
        !           update of the form:  A := A + V*W^T - W*V^T
        !
                    CALL SSKR2K( UPLO, 'No transpose', I-NPANEL, NB, ONE, A( 1, I-NPANEL+STEP ), &
                    &LDA*STEP, WORK, LDWORK, ONE, A, LDA )
        !
        !           Copy superdiagonal elements back into A
        !
                    DO 10 J = I-NPANEL+1+STEP-1, I, STEP
                    A( J-1, J ) = E( J-1 )
        10       CONTINUE
        20    CONTINUE
        !
        !        Use unblocked code to reduce the last or only block
        !
                CALL SSKTD2( UPLO, MODE, I, A, LDA, E, TAU, IINFO )
            ELSE
        !
        !        Reduce the lower triangle of A
        !
                DO 40 I = 1, N - NXPANEL, NPANEL
        !
        !           Reduce columns i:i+npanel-1 to tridiagonal form and form the
        !           matrix W which is needed to update the unreduced part of
        !           the matrix
        !
                    CALL SLASKTRD( UPLO, MODE, N-I+1, NB, A( I, I ), LDA, E( I ), TAU( I ), WORK, LDWORK )
        !
        !           Update the unreduced submatrix A(i+npanel:n,i+npanel:n), using
        !           an update of the form:  A := A + V*W^T - W*V^T
        !
                    CALL SSKR2K( UPLO, 'No transpose', N-I-NPANEL+1, NB, ONE, A( I+NPANEL, I ), LDA*STEP, &
                    &WORK( NPANEL+1 ), LDWORK, ONE, A( I+NPANEL, I+NPANEL ), LDA )
        !
        !           Copy subdiagonal elements back into A
        !
                    DO 30 J = I, I + NPANEL - 1, STEP
                    A( J+1, J ) = E( J )
        30       CONTINUE
        40    CONTINUE
        !
        !        Use unblocked code to reduce the last or only block
        !
                CALL SSKTD2( UPLO, MODE, N-I+1, A( I, I ), LDA, E( I ), TAU( I ), IINFO )
            END IF

            WORK( 1 ) = LWKOPT
        RETURN
    end subroutine SSKTRD

    subroutine ALTPFAF(UPLO, MTHD, N, A, LDA, PFAFF, IWORK, WORK, LWORK, INFO)
        implicit none
        !     .. Scalar Arguments ..
        CHARACTER          UPLO, MTHD
        INTEGER            INFO, LDA, LWORK, N
        !     ..
        !     .. Array Arguments ..
        INTEGER            IWORK( * )
        REAL               PFAFF( 2 )
        REAL               A( LDA, * ), WORK( * )
        !     ..
        !
        !  Purpose
        !  =======
        !
        !  SSKPF10 computes the Pfaffian of a real skew-symmetric matrix, taking
        !  special care to avoid numerical under- or overflow.
        !  (at the cost of possible additional round-off errors)
        !
        !  Arguments
        !  =========
        !
        !  UPLO    (input) CHARACTER*1
        !          = 'U':  Upper triangle of A is stored;
        !          = 'L':  Lower triangle of A is stored.
        !
        !  MTHD    (input) CHARACTER*1
        !          = 'P': Compute Pfaffian using Parlett-Reid algorithm (recommended)
        !          = 'H': Compute Pfaffian using Householder reflections
        !
        !  N       (input) INTEGER
        !          The order of the matrix A.  N >= 0.
        !
        !  A       (input/output) REAL array, dimension (LDA,N)
        !          On entry, the skew-symmetric matrix A.
        !             If UPLO = 'U', the upper triangular part of A contains
        !                the upper triangular part of the matrix A, and the
        !                strictly lower triangular part of A is not referenced.
        !             If UPLO = 'L', the lower triangular part of A contains
        !                the lower triangular part of the matrix A, and the
        !                strictly upper triangular part of A is not referenced.
        !          If the matrix size is odd, A is not referenced. If the matrix
        !          size is even, A is overwritten by values generated during
        !          the computation.
        !
        !  LDA     (input) INTEGER
        !          The leading dimension of the array A.  LDA >= max(1,N).
        !
        !  PFAFF   (output) REAL array, dimension 2
        !          The value of the Pfaffian in the form
        !          PFAFF(1)*10**PFAFF(2).
        !
        !  IWORK   (workspace) INTEGER array, dimension (N)
        !          Not referenced if MTHD = 'H'.
        !
        !  WORK    (workspace) REAL array,
        !             dimension (MAX(1, LWORK)), if MTHD = 'P';
        !             dimension (MAX(2*N-1,LWORK)), if MTHD = 'H'.
        !          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
        !
        !  LWORK   (input) INTEGER
        !          The dimension of the array WORK.
        !          If MTHD = 'P', LWORK >= 1,
        !          If MTHD = 'H', LWORK >= 2*N-1.
        !
        !          For optimum performance LWORK >= N*NB for MTHD = 'P' or
        !          LWORK >= N*NB+2*N-2 for MTHD = 'H', where NB is the
        !          optimal blocksize.
        !
        !          If LWORK = -1, then a workspace query is assumed; the routine
        !          only calculates the optimal size of the WORK array, returns
        !          this value as the first entry of the WORK array, and no error
        !          message related to LWORK is issued by XERBLA.
        !
        !  INFO    (output) INTEGER
        !          = 0:  successful exit
        !          < 0:  if INFO = -i, the i-th argument had an illegal value
        !
        !  Further Details
        !  ===============
        !
        !  The Pfaffian is computed by bringing the skew-symmetric matrix A into
        !  a partial tridiagonal form pT, either by computing a partial L pT L^T
        !  decomposition (MTHD = 'P'), or by a by a unitary congruence transformation
        !  Q^H * A * Q^* = pT (MTHD = 'H').
        !  These transformations are computed by the routines SSKTRF or SSKTRD,
        !  respectively (for further details see there).
        !
        REAL               ONE, ZERO
        PARAMETER          ( ONE = 1.0E+0 )
        PARAMETER          ( ZERO = 0.0E+0 )

        INTEGER            I

        !     .. Local Scalars ..
        LOGICAL            LQUERY, UPPER, LTL

        !     .. External Subroutines ..
        EXTERNAL           XERBLA
        !     .. External Functions ..
        LOGICAL            LSAME
        EXTERNAL           LSAME

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
        ELSE IF( MOD(N,2).NE.1 .AND. .NOT.LTL .AND. LWORK.LT.2*N-1 .AND. .NOT.LQUERY ) THEN
            INFO = -9
        END IF

        IF( INFO.EQ.0 .AND. LQUERY) THEN
            IF( MOD(N,2).EQ.1 ) THEN
                WORK(1) = 1
            ELSE IF( LTL ) THEN
        !     Defer workspace query to SSKTRF
                CALL SSKTRF( UPLO, "P", N, A, LDA, IWORK, WORK, LWORK, INFO )
            ELSE
        !    Defer workspace query to SSKTRD
                CALL SSKTRD( UPLO, "P", N, A, LDA, WORK, WORK, WORK, LWORK, INFO)
                WORK(1) = WORK(1) + 2*N - 2
            END IF
        END IF
        !
        IF( INFO.NE.0 ) THEN
            CALL XERBLA( 'SSKPF10', -INFO )
            RETURN
        ELSE IF( LQUERY ) THEN
            RETURN
        END IF

        PFAFF( 1 ) = ONE
        PFAFF( 2 ) = ZERO

        !     Quick return if possible
        IF( N.EQ.0 ) THEN
            RETURN
        ELSE IF( MOD(N,2).EQ.1 ) THEN
            PFAFF( 1 ) = ZERO
            RETURN
        END IF

        IF( LTL ) THEN
        !     Compute tridiagonal form
            CALL SSKTRF( UPLO, "P", N, A, LDA, IWORK, WORK, LWORK, INFO )

        !     In case one of the (relevant) off-diagonal elements is zero, the
        !     pfaffian is zero, too.
            IF( INFO .GT. 0 ) THEN
                PFAFF( 1 ) = ZERO
                PFAFF( 2 ) = ZERO
                INFO = 0
            ELSE
                IF( UPPER ) THEN

                DO 10 I = 1, N-1, 2
                    CALL SMUL10( PFAFF, A( I, I+1 ) )

        !     Accumulate the determinant of the permutations
                    IF( IWORK( I ) .NE. I ) PFAFF( 1 ) = -PFAFF( 1 )
        10            CONTINUE

                ELSE

                DO 20 I = 1, N-1, 2
                    CALL SMUL10( PFAFF, -A( I+1, I ) )

        !     Accumulate the determinant of the permutations
                    IF( IWORK( I+1 ) .NE. I+1 ) PFAFF( 1 ) = -PFAFF( 1 )
        20            CONTINUE

                END IF
            END IF
        ELSE

        !     Reduce to tridiagonal form
            CALL SSKTRD(UPLO, "P", N, A, LDA, WORK(1), WORK(N), WORK( 2*N-1 ), LWORK-2*N+2, INFO)

            IF( UPPER ) THEN
        !     Multiply every other entry on the superdiagonal
                DO 30 I = 1, N-1, 2
                CALL SMUL10( PFAFF, WORK( I ) )

        !     Accumulate the determinant of the Householder reflection
        !     (which in the real case can only be +1 or -1)
                IF (WORK( N-1+I ) .GT. ZERO) PFAFF( 1 ) = -PFAFF( 1 )
        30         CONTINUE

            ELSE

        !     Multiply every other entry on the superdiagonal
                DO 40 I = 1, N-1, 2
                CALL SMUL10( PFAFF, -WORK( I ) )

        !     Accumulate the determinant of the Householder reflection
        !     (which in the real case can only be +1 or -1)
                IF (WORK( N-1+I ) .GT. ZERO) PFAFF( 1 ) = -PFAFF( 1 )
        40         CONTINUE

            END IF

            !     Shift optimal workspace size to first position in the WORK array
            WORK( 1 ) = WORK( 2*N-1 ) + 2*N-2
        END IF

        RETURN
    end subroutine ALTPFAF

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

    subroutine MAKEALPHA(NUMT,HAMILTONIAN,ALPHAMATRIX)
        implicit none
        integer :: NUMT, i, j
        real :: ALPHAMATRIX(4*NUMT,4*NUMT)
        complex*16 :: HAMILTONIAN(4*NUMT,4*NUMT)

        do i = 1, NUMT
            do j = 1, NUMT

                ALPHAMATRIX(i, i) = 0.D0
                ALPHAMATRIX(i, j + NUMT) = AIMAG(HAMILTONIAN(i, j + 3*NUMT))
                ALPHAMATRIX(i, j + 2*NUMT) = AIMAG(HAMILTONIAN(i + NUMT, j)) - 0.5*REAL(HAMILTONIAN(i, j)+ &
                &HAMILTONIAN(i + NUMT, j + NUMT))
                ALPHAMATRIX(i, j + 3*NUMT) = REAL(HAMILTONIAN(i, j + NUMT) - HAMILTONIAN(i, j + 3*NUMT))

                ALPHAMATRIX(i + NUMT, j) = -1.D0*ALPHAMATRIX(i, j + NUMT)
                ALPHAMATRIX(i + NUMT, j + NUMT) = 0.D0
                ALPHAMATRIX(i + NUMT, j + 2*NUMT) = -1.D0*REAL(HAMILTONIAN(i, j + NUMT) + HAMILTONIAN(i, j + 3*NUMT))
                ALPHAMATRIX(i + NUMT, j + 3*NUMT) = AIMAG(HAMILTONIAN(i + NUMT, j)) + 0.5*REAL(HAMILTONIAN(i, j)+ &
                &HAMILTONIAN(i + NUMT, j + NUMT))

                ALPHAMATRIX(i + 2*NUMT, j) = -1.D0*ALPHAMATRIX(i, j + 2*NUMT)
                ALPHAMATRIX(i + 2*NUMT, j + NUMT) = -1.0*ALPHAMATRIX(i + NUMT, j + 2*NUMT)
                ALPHAMATRIX(i + 2*NUMT, j + 2*NUMT) = 0.D0
                ALPHAMATRIX(i + 2*NUMT, j + 3*NUMT) = AIMAG(HAMILTONIAN(i, j + 3*NUMT)) + 0.5*REAL(HAMILTONIAN(i, j)- &
                &HAMILTONIAN(i + NUMT, j + NUMT))

                ALPHAMATRIX(i + 3*NUMT, j) = -1.D0*ALPHAMATRIX(i, j + 3*NUMT)
                ALPHAMATRIX(i + 3*NUMT, j + NUMT) = -1.D0*ALPHAMATRIX(i + NUMT, j + 3*NUMT)
                ALPHAMATRIX(i + 3*NUMT, j + 2*NUMT) = -1.D0*ALPHAMATRIX(i + 2*NUMT, j + 3*NUMT)
                ALPHAMATRIX(i + 3*NUMT, j + 3*NUMT) = 0.D0

            end do
        end do

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
                    
                    KPTS(1, counter) = (b_1(1)*c_1)/flx + (b_2(1)*c_2)/fly + (b_3(1)*c_3)/flz
                    KPTS(2, counter) = (b_1(2)*c_1)/flx + (b_2(2)*c_2)/fly + (b_3(2)*c_3)/flz
                    KPTS(3, counter) = (b_1(3)*c_1)/flx + (b_2(3)*c_2)/fly + (b_3(3)*c_3)/flz
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

end program PFAFFIAN