program MISCROUTINES
    implicit none

    integer :: i, M, N, SITES
    real*8 :: nparam, mparam, surfparam

    ! Routine to print impurities on a grid
    M = 3 ! POSITIVE horizontal sites
    N = 3 ! POSITIVE vertical sites
    surfparam = 0.0 ! The layer corresponding to the surface

    SITES = (2*M+1)*(2*N+1)

    open(1, file = 'impurities.dat', action = 'write')
        mparam = -1.0*M
        nparam = -1.0*N
        do i = 1, SITES
            write (1,'(3F5.1)') mparam, nparam, surfparam
            if (nparam == 1.0*N) then
                nparam = -1.0*N
                mparam = mparam + 1.0
            else
                nparam = nparam + 1.0
            endif
        end do
        write (1,*) '-----------'
        write (1,*) ' x   y   z'
    close(1)

end program MISCROUTINES