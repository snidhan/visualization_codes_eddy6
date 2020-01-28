subroutine filtering(i_loc, j_loc, k_loc, sigma, ntheta, nr_truncated, nx_truncated, xc_truncated, zc_truncated, dxc_truncated,dzc_truncated, var_truncated, filtered_var)
    

    implicit none
    
    real (kind = 8)    :: var_truncated(nr_truncated, ntheta-2, nx_truncated), xc_truncated(nr_truncated), zc_truncated(nx_truncated), dxc_truncated(nr_truncated), dzc_truncated(nx_truncated), sigma, filtered_var
    integer (kind = 4) :: nx_truncated, nr_truncated, ntheta
    integer (kind = 4) :: i_loc, j_loc, k_loc, i, j, k

    real (kind = 8)    :: dist, theta, theta_loc, PI, cell_area, filter_val


    PI=4.D0*DATAN(1.D0)

    theta_loc = (j_loc-1)*2*PI/(ntheta-2) + (0.5)*(2*PI)/(ntheta-2)
    
    filtered_var = 0.0d0
    
    
            do i = 1, nr_truncated
                theta = (j_loc-1)*2*PI/(ntheta-2) + (0.5)*(2*PI)/(ntheta-2)
                dist = (zc_truncated(k_loc) - zc_truncated(k_loc))**2 + xc_truncated(i_loc)**2  + xc_truncated(i)**2 - 2*xc_truncated(i_loc)*xc_truncated(i)*cos(theta_loc - theta)
                filter_val = ((1/sqrt(2*PI)*sigma)**3)*exp(-dist/(2*sigma*sigma))
                cell_area = xc_truncated(i)*dxc_truncated(i)*(2*PI/(ntheta-2))*dzc_truncated(k_loc)
                filtered_var = filtered_var + cell_area*filter_val*var_truncated(i,j_loc,k_loc)
            end do

    return
end subroutine filtering
