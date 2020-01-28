subroutine truncation(var_centered, var_truncated, xc, zc, xc_truncated, zc_truncated, dxc_truncated, dzc_truncated, nr, ntheta, nx, nr_truncated,nx_truncated, nx_trunc_start, nx_trunc_end, nr_trunc_end)

    integer (kind = 4) :: nr_truncated, nx_truncated, nx_trunc_start, nx_trunc_end, nr_trunc_end, i
    real (kind = 8)    :: var_centered(nr-2, ntheta-2, nx-2),  var_truncated(nr_truncated, ntheta-2, nx_truncated)
    real (kind = 8)    :: xc(nr), zc(nx), xc_truncated(nr_truncated), zc_truncated(nx_truncated), dxc_truncated(nr_truncated), dzc_truncated(nx_truncated)

    xc_truncated(1:nr_truncated) = xc(2:nr_truncated+1)
    
    zc_truncated(:) = zc(nx_trunc_start+1:nx_trunc_end+1)

    do i = 1, nx_truncated
        print*, zc((i-1)+nx_trunc_start+1+1), zc((i-1)+nx_trunc_start+1-1)  
        dzc_truncated(i) = 0.50d0*(zc((i-1)+nx_trunc_start+1+1) - zc((i-1)+nx_trunc_start+1-1))
    end do

    do i = 1, nr_truncated
        print*, xc((i-1)+2+1), xc((i-1)+2-1)
        dxc_truncated(i) = 0.50d0*(xc((i-1)+2+1) - xc((i-1)+2-1))
    end do

    var_truncated(1:nr_truncated,:,1:nx_truncated) = var_centered(1:nr_truncated,:,nx_trunc_start:nx_trunc_end)

    return

end subroutine truncation
