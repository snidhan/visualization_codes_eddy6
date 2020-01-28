subroutine centering(file_type, var, var_centered, nr, ntheta, nx)

    implicit none

    integer (kind = 4) :: nr, ntheta, nx, k
    real    (kind = 8) :: var(nr, ntheta, nx), var_centered(nr-2, ntheta-2, nx-2)
    character (len = 160) :: file_type

    if (file_type .eq. 'w') then
        do k = 2,nx-1
            print*, 'Centering  ', k
            var_centered(:,:,k-1) = 0.5*(var(2:nr-1, 2:ntheta-1, k) + var(2:nr-1, 2:ntheta-1, k-1))
        end do
    
    endif


    return 
end subroutine centering
