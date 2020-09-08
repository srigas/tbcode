program fix
    implicit none

    integer :: i,j
    real*8 :: writearray(41,1000)

    open (1, file = 'numdensityperatom.txt', action = 'read')
    do j = 1, 1000
        read(1,*) (writearray(i,j), i = 1,41)
    end do
    close(1)

    open(16, file = 'fixed.txt', action = 'write')
    do j = 1, 1000
        write (16,100) (writearray(i,j), i = 1,41)
    end do
    100 format(41F17.8)
    close(16)

end program fix