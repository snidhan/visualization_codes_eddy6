program main

    implicit none
    
    integer (kind = 4), parameter :: nr     = 366
    integer (kind = 4), parameter :: ntheta = 258    
    integer (kind = 4), parameter :: nx     = 4610    

    integer (kind =4),  parameter :: nx_trunc_start = 300
    integer (kind =4),  parameter :: nx_trunc_end   = 1955
    integer (kind =4),  parameter :: nr_trunc_end   = 256
    
    character (len = 160), parameter :: grid_dir  = './'
    character (len = 160), parameter :: file_dir  = './'
    character (len = 160), parameter :: file_type = 'w'
    character (len = 160), parameter :: file_name = 'w_02635000.res'

    real (kind = 8), allocatable :: var(:,:,:), var_centered(:,:,:), var_truncated(:,:,:), var_filtered(:,:,:)
    integer (kind = 4) :: nx_truncated, nr_truncated
    real (kind = 8) :: ires, jres, kres, jpres, nstepres, timeres, DTM1res, gravres 
    
    real    (kind = 8)              :: xu(nr),yv(ntheta),zw(nx),zwg(nx)
    real    (kind = 8)              :: xc(nr),yc(ntheta),zc(nx),zcg(nx)
    real    (kind = 8), allocatable :: xc_truncated(:), zc_truncated(:), dxc_truncated(:), dzc_truncated(:)
    real    (kind = 8)              :: ru(nr),rp(nr)
    real    (kind = 8)              :: nxg
    integer (kind = 4)              :: tag
    integer (kind = 4)              :: i, j, k

    real    (kind = 8)              :: sigma
    real    (kind = 8)              :: filtered_var
    
    
    character (len = 160) :: fullfile



    allocate(var(nr,ntheta,nx))

    write(fullfile, '(a,a)') trim(file_dir), trim(file_name) 

    call read_restart(file_name,nr,ntheta,nx,var,ires,jres,kres,jpres,nstepres,timeres,DTM1res,gravres)
    
    !!!!!!!!!!!!!! Debug block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    print*, 'Minval of the read restart file ', minval(var(:,:,:)) 
    print*, 'Maxval of the read restart file ', maxval(var(:,:,:)) 
    print*, 'Shape of the read restart file ',  shape(var(:,:,:)) 
                             
    nxg = nx
    call readgrid(grid_dir, xu, yv, zw, zwg, xc, yc, zc, zcg, nr, ntheta, nx, nxg, ru, rp, tag)

    !!!!!!!!!!!!! Debug block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*, 'xu(1), xu(2), xu(3), xu(4) ', xu(1), xu(2), xu(3), xu(4)
    print*, 'xc(1), xc(2), xc(3), xc(4) ', xc(1), xc(2), xc(3), xc(4)
    print*, 'zw(1), zw(2), zw(3), zw(4) ', zw(1), zw(2), zw(3), zw(4)
    print*, 'zc(1), zc(2), zc(3), zc(4) ', zc(1), zc(2), zc(3), zc(4)


    allocate(var_centered(nr-2, ntheta-2, nx-2))
    call centering(file_type, var, var_centered, nr, ntheta, nx)

    !!!!!!!!!!! Debug block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Minval of the centered restart file ', minval(var_centered(:,:,:)) 
    print*, 'Maxval of the centered restart file ', maxval(var_centered(:,:,:)) 
    print*, 'Shape of the centered restart file ',  shape(var_centered(:,:,:)) 
   
    deallocate(var)

    nr_truncated = nr_trunc_end
    nx_truncated = nx_trunc_end - nx_trunc_start + 1
    
    print*, 'nr_truncated ', nr_truncated
    print*, 'nx_truncated ', nx_truncated
    
    allocate (var_truncated(nr_truncated, ntheta-2, nx_truncated))
    allocate (xc_truncated(nr_truncated), zc_truncated(nx_truncated))
    allocate (dxc_truncated(nr_truncated), dzc_truncated(nx_truncated))
    
    call truncation(var_centered, var_truncated, xc, zc, xc_truncated, zc_truncated, dxc_truncated, dzc_truncated, nr, ntheta, nx, nr_truncated, nx_truncated, nx_trunc_start, nx_trunc_end, nr_trunc_end)

    !!!!!!!!!! Debug block !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Shape of var_truncated ', shape(var_truncated(:,:,:))
    print*, 'Shape of xc_truncated ',  shape(xc_truncated(:))
    print*, 'Shape of zc_truncated ',  shape(zc_truncated(:))
    print*, 'zc_truncated(1), zc_truncated(nx_truncated) ',  zc_truncated(1), zc_truncated(nx_truncated) 
    print*, 'xc_truncated(1), xc_truncated(nr_truncated) ',  xc_truncated(1), xc_truncated(nr_truncated) 
    print*, 'dxc_truncated(1), dxc_truncated(nr_truncated) ',  dxc_truncated(1), dxc_truncated(nr_truncated) 
    print*, 'dzc_truncated(1), dzc_truncated(nr_truncated) ',  dzc_truncated(1), dzc_truncated(nx_truncated) 
    print*, 'Minval of var_truncated ', minval(var_truncated(:,:,:))
    print*, 'Maxval of var_truncated ', maxval(var_truncated(:,:,:))

    call gaussian_filter_parameters(sigma) 

    !!!!!!!! Starting the filtering operation !!!!!!!!!!!!!!!

    allocate(var_filtered(nr_truncated, ntheta-2, nx_truncated))
    do k = 1, nx_truncated
        do j = 1, ntheta-2
            print*, 'At z index ', k, 'for index j', j
            do i = 1, nr_truncated
                call filtering(i,j,k,sigma,ntheta, nr_truncated,nx_truncated,xc_truncated,zc_truncated,dxc_truncated,dzc_truncated,var_truncated,filtered_var) 
                var_filtered(i,j,k) = filtered_var         
            end do
        end do
    end do

    !!!!!!!! Debug block !!!!!!!!!!!!!!!!!!!!!!!
    print*, 'Shape  of var_filtered ', shape(var_filtered(:,:,:))
    print*, 'Maxval of var_filtered ', maxval(var_filtered(:,:,:))
    print*, 'Minval of var_filtered ', minval(var_filtered(:,:,:))


end program main
