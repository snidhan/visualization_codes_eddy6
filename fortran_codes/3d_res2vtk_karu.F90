Program Main
  implicit none
  
  
  !number of grid points with ghost cells
  integer ( kind = 4 ), parameter :: nx = 362
  integer ( kind = 4 ), parameter :: ny = 258
  integer ( kind = 4 ), parameter :: nz = 514

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,tmp,kstart,kend,s1,jend
  integer ( kind = 4 ) :: tag
  real    ( kind = 8 ) :: time0,time1,time2,DTM1,grav,Re
  character(len=128)   :: buff,filename

  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u0,v0,w0
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u1,v1,w1,p1,dens1,baro_ver
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u2,v2,w2
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: uc,vc,wc
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: uo,vo,wo
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: diffu_st,diffu_sw,diffu_vr
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: stretch_ver
  real    ( kind = 8 ) :: omg_mag(nx,ny,nz),omg_ver(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: domg_ver(nx,ny,nz),duydy(nx,ny,nz)
  real    ( kind = 8 ) :: adv_omg_ver(nx,ny,nz),lhs_omg_ver(nx,ny,nz),rhs_omg_ver(nx,ny,nz)
  real    ( kind = 8 ) :: lhs_minus_rhs(nx,ny,nz)

  real    ( kind = 8 ) :: Q(nx,ny,nz),lambda2(nx,ny,nz)
  real    ( kind = 8 ) :: int2d_hl_1(nz),int2d_hr_1(nz)
  real    ( kind = 8 ) :: int2d_hl_2(nz),int2d_hr_2(nz)
  real    ( kind = 8 ) :: int2d_hl_3(nz),int2d_hr_3(nz)
  real    ( kind = 8 ) :: int3d_hl_1,int3d_hr_1
  real    ( kind = 8 ) :: int3d_hl_2,int3d_hr_2
  real    ( kind = 8 ) :: int3d_hl_3,int3d_hr_3

  real    ( kind = 8 ) :: r(nx),theta(ny),z(nz),dtheta


  ! In order to reduce the file size, the domain can be sliced in the streamwise and the radial direction
   kstart = 2
   kend = nz-1

   jend = ny-5
 

  call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nz,ru,rp,tag)
  

  allocate(w1(nx,ny,nz),u1(nx,ny,nz),v1(nx,ny,nz))!,p1(nx,ny,nz),dens1(nx,ny,nz))
  allocate(wc(2:nx-1,2:ny-1,2:nz-1),uc(2:nx-1,2:ny-1,2:nz-1),vc(2:nx-1,2:ny-1,2:nz-1))


   Re = 150000


  !filename = '/home/sheel/Work2/projects_data/spheroid_ar6_0aoa/trip_tests/one_eighth_grid_tripping_tests/trip_x_2_intensity_50_re1.5e5_grid_eighth/visualization_files/u_00091000.res'
  !call read_restart(filename,nx,ny,nz,u1,time1)
  !filename = '/home/sheel/Work2/projects_data/spheroid_ar6_0aoa/trip_tests/one_eighth_grid_tripping_tests/trip_x_2_intensity_50_re1.5e5_grid_eighth/visualization_files/v_00091000.res'
  !call read_restart(filename,nx,ny,nz,v1,time1)
  filename = '/home/sheel/Work/projects/spheroid_ar6_0aoa/codes/interpolation_codes/version_2/w_00038500.res' 
  call read_restart(filename,nx,ny,nz,w1,time1)

  ! calculate center velocities
  do k=2,nz-1
  do j=2,ny-1
  do i=2,nx-1
     !uc(i,j,k) = 0.5d0*(u1(i,j,k)+u1(i-1,j,k))
     !vc(i,j,k) = 0.5d0*(v1(i,j,k)+v1(i,j-1,k))
     wc(i,j,k) = 0.5d0*(w1(i,j,k)+w1(i,j,k-1))
  enddo
  enddo
  enddo


 filename='/home/sheel/Work/projects/spheroid_ar6_0aoa/codes/interpolation_codes/version_2/w_00038500.vtk'
  !call write_vtk_4vars(filename,Q,omg_x,omg_y,omg_z,nx,ny,nz,xc,yc,zcg,kstart,kend)


  !call write_vtk(filename,w1,nx,ny,nz,xc,yc,zcg,kstart,kend)
  !call write_vtk_cartesian(filename,wc,nx,ny,nz,xc,yc,zcg,kstart,kend)

   
   call write_vtk(filename,wc,nx,ny,nz,xc,yc,zcg,kstart,kend,jend)


  !call write_vtk_cartesian_3vars(filename,uc,vc,wc,nx,ny,nz,xc,yc,zcg,kstart,kend)
  !call write_vtk_2vars(filename,Q,omg_mag,nx,ny,nz,xc,yc,zcg,kstart,kend)

  
 ! !call write_vtk(filename,Q,nx,ny,nz,xc,yc,zcg,kstart,kend)
  ! !call write_vtk_3vars(filename,omg_x,omg_y,omg_z,nx,ny,nz,xc,yc,zcg,kstart,kend)
  ! call write_vtk_vector(filename,omg_x,omg_y,omg_z,nx,ny,nz,xc,yc,zcg,kstart,kend)



  stop
end Program Main

subroutine read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg,ru,rp,tag)
  implicit none

  INTEGER nx,ny,nz,nzg,tag
  REAL ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nzg)
  REAL ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nzg)
  REAL ( kind = 8 ) :: ru(nx),rp(nx)

  INTEGER i,j,k

  real rdelx,rdely,rdelz, dtheta
  real, allocatable, dimension(:) :: cug,cvg

  ALLOCATE(cug(nzg),cvg(nzg))

  ! ! READ GRID

  OPEN(UNIT=1,FILE='/home/sheel/Work/projects/spheroid_ar6_0aoa/codes/interpolation_codes/version_2/x1_grid_half.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nx-1
     read(1,*) j, xu(i)
     !write(6,*) "xu(",i,") = ", xu(i)
  enddo
  close(1)
  xc(2:nx-1) = .5*(xu(1:nx-2)+xu(2:nx-1))
  xc(1 ) = 2.*xu(1  )-xc(2  )
  xc(nx) = 2.*xu(nx-1)-xc(nx-1)

  do i= 2, nx-1
     write(6,*) "xc(",i,") = ", xc(i)
  enddo


  OPEN(UNIT=1,FILE='/home/sheel/Work/projects/spheroid_ar6_0aoa/codes/interpolation_codes/version_2/x2_grid_half.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, ny-1
     read(1,*) j, yv(i)
     !write(6,*) "yv(",i,") = ", yv(i)
  enddo
  close(1)
  yc(2:ny-1) = .5*(yv(1:ny-2)+yv(2:ny-1))
  yc(1 ) = 2.*yv(1  )-yc(2  )
  yc(ny) = 2.*yv(ny-1)-yc(ny-1)

  do i= 2, ny-1
     write(6,*) "yc(",i,") = ", yc(i)
  enddo


  OPEN(UNIT=1,FILE='/home/sheel/Work/projects/spheroid_ar6_0aoa/codes/interpolation_codes/version_2/x3_grid_half.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
     !write(6,*) "zwg(",i,") = ", zwg(i)
  enddo
  close(1)
  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  )= 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

  do i= 2, nz-1
     write(6,*) "zcg(",i,") = ", zcg(i)
  enddo


  ! ru(1:nx)=xu(1:nx)
  ! rp(1:nx)=xc(1:nx)

  close(1)

  write(6,*) "READ GRID DONE"

  return
end subroutine read_grid

subroutine read_restart(filename,nx,ny,nz,var,time)
  implicit none
  
  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav

  write(6,*) "nx,ny,nz = ", nx,ny,nz

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')  
  READ(19) I,J,K,JP 
  write(6,*) "I,J,K,JP = ", I,J,K,JP 
  DO K=1,NZ
     write(6,*) " READ K = ", K
     READ(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  READ(19) nstep
  READ(19) TIME
  write(6,*) 'time=',time
  READ(19) DTM1,grav
  CLOSE(19)
  write(6,*) "READING RESTART FILE DONE"
  
  return
end subroutine read_restart

subroutine ascii_version

  integer (kind = 4) :: i,j,k,imax,jmax,kmax
  
  open(1,file="3d.vtk", form = 'formatted', status = 'unknown')

  imax = 10
  jmax = 10
  kmax = 10

  write(1,FMT='(a26)') "# vtk DataFile Version 2.0"
  write(1,FMT='(a9)') "FlowField"
  write(1,FMT='(a5)') "ASCII"
  write(1,FMT='(a23)') "DATASET STRUCTURED_GRID"
  write(1,FMT='(a10)', advance="no") "DIMENSIONS"
  write(1,"(a1, 3i7)") " ", imax, jmax, kmax
  write(1,FMT='(a6)', advance="no") "POINTS"
  write(1,"(a1, i15)", advance="no") " ", imax*jmax*kmax
  write(1,"(a6)") " float"

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1,*) real(i), real(j), real(k)
        end do
     end do
  end do

  close(1)

end subroutine ascii_version

subroutine binary_version

  integer (kind = 4) :: i,j,k,imax,jmax,kmax,s1
  character(len=128) :: buff
  real    (kind = 4), allocatable, dimension(:,:,:) :: vtkvar
  
  open(unit=1,file='3d.vtk',access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = 10
  jmax = 10
  kmax = 10

  allocate(vtkvar(imax,jmax,kmax))

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1) real(i), real(j), real(k)
           vtkvar(i,j,k) = real(i)*real(j)*real(k)
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS vtkvar float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           write(1) vtkvar(i,j,k)
        end do
     end do
  end do


  close(1)
end subroutine binary_version

subroutine center_velocity(nx,ny,nz,Ucen,Uin,dir)
  !Passed Variables
  integer,intent(in)      :: dir,nx,ny,nz
  real (kind = 8),intent(in)         :: Uin(nx,ny,nz) 
  real (kind = 8),intent(out)        :: Ucen(nx,ny,nz) 


  !Local Variables
  integer              	:: i,j,k,err
  !Zero Output Array
  Ucen=0.d0
  !*************************************************
  !********************X1***************************
  !*************************************************
  if(dir.EQ.1) then
     !U
     do k=2,nz-1
        write(6,*) " CENTER_VELOCITY DIR1, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i-1,j,k))
           enddo
     	enddo
     enddo

     !*************************************************
     !********************X2***************************
     !*************************************************
  elseif (dir.EQ.2) then 
     !V
     !U
     do k=2,nz-1
        write(6,*) " CENTER_VELOCITY DIR2, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j-1,k))
           enddo
   	enddo
     enddo
     !*************************************************
     !********************X3***************************
     !*************************************************
  elseif (dir.EQ.3) then 
     !W
     !U
     do k=2,nz-1
        write(6,*) " CENTER_VELOCITY DIR3, K = ", K
    	do j=2,ny-1
           do i=2,nx-1
              Ucen(i,j,k)=0.50d0*(uin(i,j,k)+uin(i,j,k-1))
           enddo
   	enddo
     enddo

  else
     !Invalid direction
     !write(*,'(a60,i2)') "INVALID DIRECTION IN center_velocities 
     !&                      dir must be 1,2,3.  dir= ", dir
  endif

  Ucen(:,1,:) = Ucen(:,ny-1,:)
  Ucen(:,ny,:) = Ucen(:,2,:)

  return
end subroutine center_velocity

! subroutine vorticity(u,v,w,omg_mag,Q_cri,nx,ny,nz,xc,zcg,rp)
!   implicit none
  
!   integer ( kind = 4 ) :: i,j,k,nx,ny,nz
!   real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
!   real    ( kind = 8 ) :: uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
!   real    ( kind = 8 ) :: omgt(nx,ny,nz),omgr(nx,ny,nz),omgz(nx,ny,nz)
!   real    ( kind = 8 ) :: omg_mag(nx,ny,nz),Q_Cri(nx,ny,nz)
!   real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
!   real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
!   real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

!   real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
!   real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
!   real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)


  
!   real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz),rp(nx)

!   dtheta = 2.0*3.1415926/dble(ny-2)

!   ! CALL CENTER_VELOCITY(nx,ny,nz,uo,u,1)
!   ! CALL CENTER_VELOCITY(nx,ny,nz,vo,v,2)
!   ! CALL CENTER_VELOCITY(nx,ny,nz,wo,w,3)
!   ! write(6,*) " DONE CENTER_VELOCITY"


!   do k=2,nz-1
!      write(6,*) " OMGT: K = ", K
!      do j=2,ny-1
!         do i=2,nx-1
!            dz=zcg(k)-zcg(k-1)
!            dr=xc(i)-xc(i-1)

!            dq1x3=(0.25*(uo(i,j,k)+uo(i,j,k+1)+uo(i-1,j,k)+uo(i-1,j,k+1)) &
!                  -0.25*(uo(i,j,k)+uo(i,j,k-1)+uo(i-1,j,k)+uo(i-1,j,k-1)))/dz

!            dq3x1=(0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i+1,j,k)+wo(i+1,j,k-1)) &
!                  -0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i-1,j,k)+wo(i-1,j,k-1)))/dr
!            omgt(i,j,k)=dq1x3-dq3x1
!            !write(6,*) "omgt = ", omgt(i,j,k)

!            dq1x3_g(i,j,k)=(0.25*(uo(i,j,k)+uo(i,j,k+1)+uo(i-1,j,k)+uo(i-1,j,k+1)) &
!                  -0.25*(uo(i,j,k)+uo(i,j,k-1)+uo(i-1,j,k)+uo(i-1,j,k-1)))/dz

!            dq3x1_g(i,j,k)=(0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i+1,j,k)+wo(i+1,j,k-1)) &
!                  -0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i-1,j,k)+wo(i-1,j,k-1)))/dr

!         enddo
!      enddo
!   enddo
  
!   do k=2,nz-1
!      write(6,*) " OMGR: K = ", K
!      do j=2,ny-1
!         do i=2,nx-1
!            dz=zcg(k)-zcg(k-1)
!            dq2x3=(0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i,j,k+1)+vo(i,j-1,k+1)) &
!                  -0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i,j,k-1)+vo(i,j-1,k-1)))/dz
!            dq3x2=(0.25*(wo(i,j,k)+wo(i,j+1,k)+wo(i,j,k-1)+wo(i,j+1,k-1)) &
!                  -0.25*(wo(i,j,k)+wo(i,j-1,k)+wo(i,j,k-1)+wo(i,j-1,k-1)))/dtheta
!            omgr(i,j,k)=dq3x2/rp(i)-dq2x3
!            !write(6,*) "omgr = ", omgr(i,j,k)

!            dq2x3_g(i,j,k)=(0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i,j,k+1)+vo(i,j-1,k+1)) &
!                  -0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i,j,k-1)+vo(i,j-1,k-1)))/dz
!            dq3x2_g(i,j,k)=(0.25*(wo(i,j,k)+wo(i,j+1,k)+wo(i,j,k-1)+wo(i,j+1,k-1)) &
!                  -0.25*(wo(i,j,k)+wo(i,j-1,k)+wo(i,j,k-1)+wo(i,j-1,k-1)))/dtheta

!         enddo
!      enddo
!   enddo

!   do k=2,nz-1
!      write(6,*) " OMGZ: K = ", K
!      do j=2,ny-1
!         do i=2,nx-1
!            dr=xc(i)-xc(i-1)
!            dq2x1=(0.25*((vo(i,j,k)+vo(i,j-1,k))*rp(i)+(vo(i+1,j,k)+vo(i+1,j-1,k))*rp(i+1)) &
!                  -0.25*((vo(i,j,k)+vo(i,j-1,k))*rp(i)+(vo(i-1,j,k)+vo(i-1,j-1,k))*rp(i-1)))/dz
!            dq1x2=(0.25*(uo(i,j,k)+uo(i,j+1,k)+uo(i-1,j,k)+uo(i-1,j+1,k)) &
!                  -0.25*(uo(i,j,k)+uo(i,j-1,k)+uo(i-1,j,k)+uo(i-1,j-1,k)))/dtheta
!            omgz(i,j,k)=(dq2x1-dq1x2)/rp(i)
!            !write(6,*) "omgz = ", omgz(i,j,k)

!            dq2x1_g(i,j,k)=(0.25*((vo(i,j,k)+vo(i,j-1,k))*rp(i)+(vo(i+1,j,k)+vo(i+1,j-1,k))*rp(i+1)) &
!                  -0.25*((vo(i,j,k)+vo(i,j-1,k))*rp(i)+(vo(i-1,j,k)+vo(i-1,j-1,k))*rp(i-1)))/dz
!            dq1x2_g(i,j,k)=(0.25*(uo(i,j,k)+uo(i,j+1,k)+uo(i-1,j,k)+uo(i-1,j+1,k)) &
!                  -0.25*(uo(i,j,k)+uo(i,j-1,k)+uo(i-1,j,k)+uo(i-1,j-1,k)))/dtheta
!         enddo
!      enddo
!   enddo


!   do k=2,nz-1
!      write(6,*) " DQIXI: K = ", K
!      do j=2,ny-1
!         do i=2,nx-1

!            dq1x1_g(i,j,k)=(0.25*((uo(i,j,k)+uo(i,j-1,k))*rp(i)+(uo(i+1,j,k)+uo(i+1,j-1,k))*rp(i+1)) &
!                  -0.25*((uo(i,j,k)+uo(i,j-1,k))*rp(i)+(uo(i-1,j,k)+uo(i-1,j-1,k))*rp(i-1)))/dz
!            dq2x2_g(i,j,k)=(0.25*(vo(i,j,k)+vo(i,j+1,k)+vo(i-1,j,k)+vo(i-1,j+1,k)) &
!                  -0.25*(vo(i,j,k)+vo(i,j-1,k)+vo(i-1,j,k)+vo(i-1,j-1,k)))/dtheta
!            dq3x3_g(i,j,k)=(0.25*(wo(i,j,k)+wo(i,j,k+1)+wo(i-1,j,k)+wo(i-1,j,k+1)) &
!                  -0.25*(wo(i,j,k)+wo(i,j,k-1)+wo(i-1,j,k)+wo(i-1,j,k-1)))/dz

!         enddo
!      enddo
!   enddo


!   do k=2,nz-1
!      write(6,*) " OMG_MAG: K = ", K
!      do j=2,ny-1
!         do i=2,nx-1
!            omg_mag(i,j,k) = sqrt(omgt(i,j,k)**2+omgr(i,j,k)**2+omgz(i,j,k)**2)
!            !write(6,*) "omg_mag = ", omg_mag(i,j,k)
!            Q_Cri(i,j,k)   = 0.5*( dq1x1_g(i,j,k)**2 + dq1x2_g(i,j,k)**2 + dq1x3_g(i,j,k)**2 + &
!                                   dq2x1_g(i,j,k)**2 + dq2x2_g(i,j,k)**2 + dq2x3_g(i,j,k)**2 + &
!                                   dq3x1_g(i,j,k)**2 + dq3x2_g(i,j,k)**2 + dq3x3_g(i,j,k)**2    )
!         enddo
!      enddo
!   enddo

!   return
! end subroutine vorticity

subroutine vorticity_mag(u,v,w,omg_mag,nx,ny,nz,xu,zwg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: omg_mag(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3
  real    ( kind = 8 ) :: omgt,omgr,omgz

  ! real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)
  
  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)

  !do k=2,nz-1
  do k=kstart-10,kend+10
     write(6,*) " OMG_MAG: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)

           omg_mag(i,j,k) = sqrt(omgt**2+omgr**2+omgz**2)

        enddo
     enddo
  enddo

  do i=1,10
     omg_mag(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     omg_mag(i,:,:) = 0.0
  end do

  omg_mag(:,1,:) = omg_mag(:,ny-1,:)
  omg_mag(:,ny,:) = omg_mag(:,2,:)

  return
end subroutine vorticity_mag
                        
subroutine vorticity_3comps(u,v,w,omg_x,omg_y,omg_z,omg_mag,nx,ny,nz,xu,yc,zwg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: omg_mag(nx,ny,nz)
  real    ( kind = 8 ) :: omg_x(nx,ny,nz),omg_y(nx,ny,nz),omg_z(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3
  real    ( kind = 8 ) :: omgt,omgr,omgz

  ! real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)
  
  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx),yc(ny)

  dtheta = 2.0*3.1415926/dble(ny-2)

  !do k=2,nz-1
  do k=kstart-10,kend+10
     write(6,*) " OMG_MAG: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)

           omg_mag(i,j,k) = sqrt(omgt**2+omgr**2+omgz**2)

           omg_x(i,j,k) = sin(yc(j))*omgr + cos(yc(j))*omgt
           omg_y(i,j,k) = cos(yc(j))*omgr - sin(yc(j))*omgt
           omg_z(i,j,k) = omgz

        enddo
     enddo
  enddo

  do i=1,3
     omg_x(i,:,:) = 0.0
     omg_y(i,:,:) = 0.0
     omg_z(i,:,:) = 0.0
     omg_mag(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     omg_x(i,:,:) = 0.0
     omg_y(i,:,:) = 0.0
     omg_z(i,:,:) = 0.0
     omg_mag(i,:,:) = 0.0
  end do

  omg_mag(:,1,:) = omg_mag(:,ny-1,:)
  omg_mag(:,ny,:) = omg_mag(:,2,:)

  omg_x(:,1,:) = omg_x(:,ny-1,:)
  omg_x(:,ny,:) = omg_x(:,2,:)

  omg_y(:,1,:) = omg_y(:,ny-1,:)
  omg_y(:,ny,:) = omg_y(:,2,:)

  omg_z(:,1,:) = omg_z(:,ny-1,:)
  omg_z(:,ny,:) = omg_z(:,2,:)

  return
end subroutine vorticity_3comps

subroutine vorticity_mag_sth_wrong(u,v,w,omg_mag,nx,ny,nz,xu,xc,zwg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: omg_mag(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3
  real    ( kind = 8 ) :: omgt,omgr,omgz

  real    ( kind = 8 ) :: v_tmp_i,v_tmp_im1,v_tmp_ip1

  real    ( kind = 8 ) :: weight_1,weight_2,weight_3
  real    ( kind = 8 ) :: weight_4,weight_5,weight_6
  
  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)
  real    ( kind = 8 ) :: xc(nx),zcg(nz)

  dtheta = 2.0*3.1415926/dble(ny-2)

  !do k=2,nz-1
  do k=kstart-10,kend+10
     write(6,*) " OMG_MAG: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zcg(k+1)-zcg(k)                 
           dr=xu(i)-xu(i-1)

           ! dq_x_ are at cell-center
           !---------------------------------------------------------------------!
           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1))           &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta


           weight_1 = zwg(k  )-zcg(k  )
           weight_2 = zcg(k+1)-zwg(k  )
           weight_3 = zwg(k-1)-zcg(k-1)
           weight_4 = zcg(k  )-zwg(k-1)
           dq2x3=( (  0.5*(v(i,j,k+1)+v(i,j-1,k+1))*weight_1                    &
                     +0.5*(v(i,j,k  )+v(i,j-1,k  ))*weight_2)/(zcg(k+1)-zcg(k)) &
                  -(  0.5*(v(i,j,k  )+v(i,j-1,k  ))*weight_3                    &
                     +0.5*(v(i,j,k-1)+v(i,j-1,k-1))*weight_4)/(zcg(k)-zcg(k-1)) &
                 )/dz


           weight_1 = zwg(k  )-zcg(k  )
           weight_2 = zcg(k+1)-zwg(k  )
           weight_3 = zwg(k-1)-zcg(k-1)
           weight_4 = zcg(k  )-zwg(k-1)
           dq1x3=( (  0.5*(u(i,j,k+1)+u(i-1,j,k+1))*weight_1                    &
                     +0.5*(u(i,j,k  )+u(i-1,j,k  ))*weight_2)/(zcg(k+1)-zcg(k)) &
                  -(  0.5*(u(i,j,k  )+u(i-1,j,k  ))*weight_3                    &
                     +0.5*(u(i,j,k-1)+u(i-1,j,k-1))*weight_4)/(zcg(k)-zcg(k-1)) &
                 )/dz


           weight_1 = xc(i+1)-xu(i)
           weight_2 = xu(i  )-xc(i)
           weight_3 = xu(i-1)-xc(i-1)
           weight_4 = xc(i  )-xu(i-1)
           dq3x1=( (  0.5*(w(i  ,j,k)+w(i  ,j,k-1))*weight_1                   &
                     +0.5*(w(i+1,j,k)+w(i+1,j,k-1))*weight_2)/(xc(i+1)-xc(i))  &
                  -(  0.5*(w(i  ,j,k)+w(i  ,j,k-1))*weight_3                   &
                     +0.5*(w(i-1,j,k)+w(i-1,j,k-1))*weight_4)/(xc(i)-xc(i-1))  &
                 )/dr


           weight_1 = xc(i+1)-xu(i)
           weight_2 = xu(i  )-xc(i)
           weight_3 = xu(i-1)-xc(i-1)
           weight_4 = xc(i  )-xu(i-1)
           dq2x1=( (  0.5*(v(i  ,j,k)+v(i  ,j-1,k))*rp(i)*weight_1                    &
                     +0.5*(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)*weight_2)/(xc(i+1)-xc(i)) &
                  -(  0.5*(v(i  ,j,k)+v(i  ,j-1,k))*rp(i)*weight_3                    &
                     +0.5*(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)*weight_4)/(xc(i)-xc(i-1)) &
                 )/dr


           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta
           !---------------------------------------------------------------------!

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)

           omg_mag(i,j,k) = sqrt(omgt**2+omgr**2+omgz**2)

        enddo
     enddo
  enddo

  do i=1,10
     omg_mag(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     omg_mag(i,:,:) = 0.0
  end do

  omg_mag(:,1,:) = omg_mag(:,ny-1,:)
  omg_mag(:,ny,:) = omg_mag(:,2,:)

  return
end subroutine vorticity_mag_sth_wrong

subroutine vorticity_vertical(u,v,w,omg_ver,nx,ny,nz,xu,zwg,rp,yc,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: omg_ver(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: omgt,omgr,omgz

  ! real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)
  
  real    ( kind = 8 ) :: xu(nx),yc(ny),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)

  !do k=2,nz-1
  do k=kstart-10,kend+10
     write(6,*) " OMG_VER: K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)
           
           omg_ver(i,j,k)=omgr*sin(yc(j))+omgt*cos(yc(j))    
           !write(6,*) " omg_ver(i,j,k) = ", omg_ver(i,j,k)
        enddo
     enddo
  enddo

  do i=1,10
     omg_ver(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     omg_ver(i,:,:) = 0.0
  end do

  omg_ver(:,1,:) = omg_ver(:,ny-1,:)
  omg_ver(:,ny,:) = omg_ver(:,2,:)

  return
end subroutine vorticity_vertical

subroutine QCri(u,v,w,Q,nx,ny,nz,xu,zwg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: Q(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: s_11,s_12,s_13
  real    ( kind = 8 ) :: s_21,s_22,s_23
  real    ( kind = 8 ) :: s_31,s_32,s_33

  real    ( kind = 8 ) :: o_11,o_12,o_13
  real    ( kind = 8 ) :: o_21,o_22,o_23
  real    ( kind = 8 ) :: o_31,o_32,o_33


  ! real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)
  
  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)

  do k=2,nz-1
     write(6,*) " K, Q(100,2,k-1) = ", K, Q(100,2,k-1)
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta

           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz

           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           dq1x1=(u(i,j,k)-u(i-1,j,k))/dr
           dq2x2=(v(i,j,k)-v(i,j-1,k))/dtheta
           dq3x3=(w(i,j,k)-w(i,j,k-1))/dz

           o_11 = 0.5*(dq1x1-dq1x1)
           o_12 = 0.5*(dq1x2-dq2x1)
           o_13 = 0.5*(dq1x3-dq3x1)

           o_21 = 0.5*(dq2x1-dq1x2)
           o_22 = 0.5*(dq2x2-dq2x2)
           o_23 = 0.5*(dq2x3-dq3x2)

           o_31 = 0.5*(dq3x1-dq1x3)
           o_32 = 0.5*(dq3x2-dq2x3)
           o_33 = 0.5*(dq3x3-dq3x3)

           s_11 = 0.5*(dq1x1+dq1x1)
           s_12 = 0.5*(dq1x2+dq2x1)
           s_13 = 0.5*(dq1x3+dq3x1)

           s_21 = 0.5*(dq2x1+dq1x2)
           s_22 = 0.5*(dq2x2+dq2x2)
           s_23 = 0.5*(dq2x3+dq3x2)

           s_31 = 0.5*(dq3x1+dq1x3)
           s_32 = 0.5*(dq3x2+dq2x3)
           s_33 = 0.5*(dq3x3+dq3x3)

           Q(i,j,k) = sqrt( o_11*o_11 + o_11*o_11 + o_11*o_11 + &
                            o_11*o_11 + o_11*o_11 + o_11*o_11 + &
                            o_11*o_11 + o_11*o_11 + o_11*o_11    ) &
                     -sqrt( s_11*s_11 + s_11*s_11 + s_11*s_11 + &              
                            s_11*s_11 + s_11*s_11 + s_11*s_11 + &
                            s_11*s_11 + s_11*s_11 + s_11*s_11    ) 

        enddo
     enddo
  enddo

  do i=1,10
     Q(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     Q(i,:,:) = 0.0
  end do

  Q(:,1,:) = Q(:,ny-1,:)
  Q(:,ny,:) = Q(:,2,:)

  return
end subroutine QCri

subroutine QCri_NEW(u,v,w,Q,nx,ny,nz,xu,zwg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: Q(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: s_11,s_12,s_13
  real    ( kind = 8 ) :: s_21,s_22,s_23
  real    ( kind = 8 ) :: s_31,s_32,s_33

  real    ( kind = 8 ) :: o_11,o_12,o_13
  real    ( kind = 8 ) :: o_21,o_22,o_23
  real    ( kind = 8 ) :: o_31,o_32,o_33

  real    ( kind = 8 ) :: OMG_MAG,S_MAG

  
  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)

  do k=2,nz-1
     write(6,*) " K, Q(100,2,k-1) = ", K, Q(100,2,k-1)
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/(rp(i)*dtheta)

           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz

           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

           dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
                 -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/(rp(i)*dr)

           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/(rp(i)*dtheta)

           dq1x1=(u(i,j,k)-u(i-1,j,k))/dr
           dq2x2=(v(i,j,k)-v(i,j-1,k))/(rp(i)*dtheta)
           dq3x3=(w(i,j,k)-w(i,j,k-1))/dz

           o_11 = 0.5*(dq1x1-dq1x1)
           o_12 = 0.5*(dq1x2-dq2x1)
           o_13 = 0.5*(dq1x3-dq3x1)

           o_21 = 0.5*(dq2x1-dq1x2)
           o_22 = 0.5*(dq2x2-dq2x2)
           o_23 = 0.5*(dq2x3-dq3x2)

           o_31 = 0.5*(dq3x1-dq1x3)
           o_32 = 0.5*(dq3x2-dq2x3)
           o_33 = 0.5*(dq3x3-dq3x3)

           s_11 = 0.5*(dq1x1+dq1x1)
           s_12 = 0.5*(dq1x2+dq2x1)
           s_13 = 0.5*(dq1x3+dq3x1)

           s_21 = 0.5*(dq2x1+dq1x2)
           s_22 = 0.5*(dq2x2+dq2x2)
           s_23 = 0.5*(dq2x3+dq3x2)

           s_31 = 0.5*(dq3x1+dq1x3)
           s_32 = 0.5*(dq3x2+dq2x3)
           s_33 = 0.5*(dq3x3+dq3x3)

           OMG_MAG = sqrt( o_11**2 + o_12**2 + o_13**2 + &
                           o_21**2 + o_22**2 + o_23**2 + &
                           o_31**2 + o_32**2 + o_33**2    )

             S_MAG = sqrt( s_11**2 + s_12**2 + s_13**2 + &
                           s_21**2 + s_22**2 + s_23**2 + &
                           s_31**2 + s_32**2 + s_33**2    ) 

             Q(i,j,k) = 0.5d0*(OMG_MAG**2 - S_MAG**2)
             !Q(i,j,k) = 0.5d0*((OMG_MAG**2 - S_MAG**2)/OMG_MAG**2)

        enddo
     enddo
  enddo

  do i=1,10
     Q(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     Q(i,:,:) = 0.0
  end do

  Q(:,1,:) = Q(:,ny-1,:)
  Q(:,ny,:) = Q(:,2,:)

  return
end subroutine QCri_NEW

subroutine QCri_STH_WRONG(u,v,w,Q,nx,ny,nz,xu,xc,yc,zwg,zcg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),TMP4(nx,ny,nz)

  real,dimension(:,:,:),allocatable :: UX,UY,UZ
  real,dimension(:,:,:),allocatable :: DUX_DR,DUX_DT,DUX_DZ
  real,dimension(:,:,:),allocatable :: DUY_DR,DUY_DT,DUY_DZ
  real,dimension(:,:,:),allocatable :: DUZ_DR,DUZ_DT,DUZ_DZ
  real                              :: DUX_DX,DUX_DY
  real                              :: DUY_DX,DUY_DY
  real                              :: DUZ_DX,DUZ_DY

  real,dimension(:,:,:),allocatable :: S11,S12,S13
  real,dimension(:,:,:),allocatable :: S21,S22,S23
  real,dimension(:,:,:),allocatable :: S31,S32,S33

  real,dimension(:,:,:),allocatable :: O11,O12,O13
  real,dimension(:,:,:),allocatable :: O21,O22,O23
  real,dimension(:,:,:),allocatable :: O31,O32,O33

  real    ( kind = 8 ) :: Q(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: s_11,s_12,s_13
  real    ( kind = 8 ) :: s_21,s_22,s_23
  real    ( kind = 8 ) :: s_31,s_32,s_33

  real    ( kind = 8 ) :: o_11,o_12,o_13
  real    ( kind = 8 ) :: o_21,o_22,o_23
  real    ( kind = 8 ) :: o_31,o_32,o_33

  real    ( kind = 8 ) :: v_tmp_i,v_tmp_im1,v_tmp_ip1

  real    ( kind = 8 ) :: weight_1,weight_2,weight_3
  real    ( kind = 8 ) :: weight_4,weight_5,weight_6

  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)
  real    ( kind = 8 ) :: xc(nx),zcg(nz),yc(ny)

  real    ( kind = 8 ) :: OMG_MAG,S_MAG
  integer :: AllocateStatus

  write(6,*) "HERE 1"
  write(6,*) "nz = ", nz
  allocate(UX(nx,ny,nz))
  allocate(UY(nx,ny,nz))
  allocate(UZ(nx,ny,nz))

        DO i=2,NX-1
        DO j=2,NY-1
        DO k=2,NZ-1
                 
           ! U_F(i,j,k) are at cell-center
           UX(i,j,k) =   sin(yc(j))*0.5d0*(u(i,j,k)+u(i-1,j,k)) &
                       + cos(yc(j))*0.5d0*(v(i,j,k)+v(i,j-1,k))

           UY(i,j,k) =   cos(yc(j))*0.5d0*(u(i,j,k)+u(i-1,j,k)) &
                       - sin(yc(j))*0.5d0*(v(i,j,k)+v(i,j-1,k))

           UZ(i,j,k) =   0.5d0*(w(i,j,k)+w(i,j,k-1))

        ENDDO
        ENDDO
        ENDDO
  write(6,*) "HERE 2"
  write(6,*) "nz = ", nz
        call SET_INTERP_I(UX,NX,NY,NZ)      
        call SET_PERIODIC_J(UX,NX,NY,NZ)      

        call SET_INTERP_I(UY,NX,NY,NZ)      
        call SET_PERIODIC_J(UY,NX,NY,NZ)      

        call SET_INTERP_I(UZ,NX,NY,NZ)      
        call SET_PERIODIC_J(UZ,NX,NY,NZ)      
  write(6,*) "HERE 3"       
  write(6,*) "nz = ", nz
  write(6,*) "HERE 4"   

          write(6,*) "nx,ny,nz = ", nx,ny,nz
      !----------------------------------------------------------!
        allocate(DUX_DR(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUX_DR= ", AllocateStatus
        allocate(DUX_DT(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUX_DT= ", AllocateStatus
        allocate(DUX_DZ(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUX_DZ= ", AllocateStatus
        
        allocate(DUY_DR(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUY_DR= ", AllocateStatus
        allocate(DUY_DT(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUY_DT= ", AllocateStatus
        allocate(DUY_DZ(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUY_DZ= ", AllocateStatus

        allocate(DUZ_DR(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUZ_DR= ", AllocateStatus
        allocate(DUZ_DT(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUZ_DT= ", AllocateStatus
        allocate(DUZ_DZ(nx,ny,nz), STAT = AllocateStatus)
        write(6,*) "AllocateStatus DUZ_DZ= ", AllocateStatus

        TMP4 = 0.0d0
        do i=1,nx-1
           ! TMP4 is at U(i,j,k)
           !TMP4(i,:,:)=0.5d0*(UXF(i,:,:)+UXF(i+1,:,:))

           weight_1 = xc(i+1)-xu(i  )
           weight_2 = xu(i  )-xc(i  )
           TMP4(i,:,:) = (UX(i,:,:)*weight_1 + UX(i+1,:,:)*weight_2)/(xc(i+1)-xc(i))
        enddo
        call SET_INTERP_I(TMP4,NX,NY,NZ)
        call SET_PERIODIC_J(TMP4,NX,NY,NZ)
        ! DUXF_DR is at cell-center
        call DERIVATIVE(TMP4,DUX_DR,NX,NY,NZ,XU,ZWG,1,1)
        call SET_INTERP_I(DUX_DR,NX,NY,NZ)      
        call SET_PERIODIC_J(DUX_DR,NX,NY,NZ)      
        write(6,*) "HERE 4.1"       
        TMP4 = 0.0d0
        do j=1,ny-1
           ! TMP4 is at V(i,j,k)
           TMP4(:,j,:)=0.50d0*(UX(:,j,:)+UX(:,j+1,:))
        enddo
        ! DUX_DT is at cell-center
        call DERIVATIVE(TMP4,DUX_DT,NX,NY,NZ,XU,ZWG,2,2)
        call SET_INTERP_I(DUX_DT,NX,NY,NZ)
        call SET_PERIODIC_J(DUX_DT,NX,NY,NZ)
        write(6,*) "HERE 4.2"       
        TMP4 = 0.0d0
        do k=1,nz-1
           weight_1 = zcg(k+1)-zwg(k  )
           weight_2 = zwg(k  )-zcg(k  )
           TMP4(:,:,k)=(UX(:,:,k)*weight_1 + UX(:,:,k+1)*weight_2)/(zcg(k+1)-zcg(k))
        enddo
        call SET_INTERP_I(TMP4,NX,NY,NZ)
        call SET_PERIODIC_J(TMP4,NX,NY,NZ)
        ! DUX_DZ is at cell-center
        call DERIVATIVE(TMP4,DUX_DZ,NX,NY,NZ,XU,ZWG,3,3)
        call SET_INTERP_I(DUX_DZ,NX,NY,NZ)
        call SET_PERIODIC_J(DUX_DZ,NX,NY,NZ)
        !----------------------------------------------------------!
       
        write(6,*) "HERE 4.3"       
        TMP4 = 0.0d0
        do i=1,nx-1
           weight_1 = xc(i+1)-xu(i  )
           weight_2 = xu(i  )-xc(i  )
           TMP4(i,:,:) = (UY(i,:,:)*weight_1 + UY(i+1,:,:)*weight_2)/(xc(i+1)-xc(i))
        enddo
        call SET_INTERP_I(TMP4,NX,NY,NZ)
        call SET_PERIODIC_J(TMP4,NX,NY,NZ)
        ! DUY_DR is at cell-center
        call DERIVATIVE(TMP4,DUY_DR,NX,NY,NZ,XU,ZWG,1,1)
        call SET_INTERP_I(DUY_DR,NX,NY,NZ)
        call SET_PERIODIC_J(DUY_DR,NX,NY,NZ)
        write(6,*) "HERE 4.4"       
        TMP4 = 0.0d0
        do j=1,ny-1
           ! TMP4 is at V(i,j,k)
           TMP4(:,j,:)=0.50d0*(UY(:,j,:)+UY(:,j+1,:))
        enddo
        ! DUY_DT is at cell-center
        call DERIVATIVE(TMP4,DUY_DT,NX,NY,NZ,XU,ZWG,2,2)
        call SET_INTERP_I(DUY_DT,NX,NY,NZ)
        call SET_PERIODIC_J(DUY_DT,NX,NY,NZ)
        write(6,*) "HERE 4.5"       

        TMP4 = 0.0d0
        do k=1,nz-1
           weight_1 = zcg(k+1)-zwg(k  )
           weight_2 = zwg(k  )-zcg(k  )
           TMP4(:,:,k)=(UY(:,:,k)*weight_1 + UY(:,:,k+1)*weight_2)/(zcg(k+1)-zcg(k))
        enddo
        call SET_INTERP_I(TMP4,NX,NY,NZ)
        call SET_PERIODIC_J(TMP4,NX,NY,NZ)
        ! DUY_DZ is at cell-center
        call DERIVATIVE(TMP4,DUY_DZ,NX,NY,NZ,XU,ZWG,3,3)
        call SET_INTERP_I(DUY_DZ,NX,NY,NZ)
        call SET_PERIODIC_J(DUY_DZ,NX,NY,NZ)
        !----------------------------------------------------------!
        write(6,*) "HERE 4.6"       
        TMP4 = 0.0d0
        do i=1,nx-1
           weight_1 = xc(i+1)-xu(i  )
           weight_2 = xu(i  )-xc(i  )
           TMP4(i,:,:) = (UZ(i,:,:)*weight_1 + UZ(i+1,:,:)*weight_2)/(xc(i+1)-xc(i))
        enddo
        call SET_INTERP_I(TMP4,NX,NY,NZ)
        call SET_PERIODIC_J(TMP4,NX,NY,NZ)
        ! DUZ_DR is at cell-center
        call DERIVATIVE(TMP4,DUZ_DR,NX,NY,NZ,XU,ZWG,1,1)
        call SET_INTERP_I(DUZ_DR,NX,NY,NZ)
        call SET_PERIODIC_J(DUZ_DR,NX,NY,NZ)
        write(6,*) "HERE 4.7"       
        TMP4 = 0.0d0
        do j=1,ny-1
           ! TMP4 is at V(i,j,k)
           TMP4(:,j,:)=0.50d0*(UZ(:,j,:)+UZ(:,j+1,:))
        enddo
        call DERIVATIVE(TMP4,DUZ_DT,NX,NY,NZ,XU,ZWG,2,2)
        call SET_INTERP_I(DUZ_DT,NX,NY,NZ)
        call SET_PERIODIC_J(DUZ_DT,NX,NY,NZ)
        write(6,*) "HERE 4.8"       
        TMP4 = 0.0d0
        do k=1,nz-1
           TMP4(:,:,k) = w(:,:,k)
        enddo
        write(6,*) "HERE 4.8.1"       
        !call SET_INTERP_I(TMP4,NX,NY,NZ)
        write(6,*) "HERE 4.8.2"       
        !call SET_PERIODIC_J(TMP4,NX,NY,NZ)
        write(6,*) "HERE 4.8.3"       
        ! DUZ_DZ is at cell-center
        call DERIVATIVE(TMP4,DUZ_DZ,NX,NY,NZ,XU,ZWG,3,3)
        write(6,*) "HERE 4.8.4"       
        call SET_INTERP_I(DUZ_DZ,NX,NY,NZ)
        write(6,*) "HERE 4.8.5"       
        call SET_PERIODIC_J(DUZ_DZ,NX,NY,NZ)
        !----------------------------------------------------------!
        write(6,*) "HERE 5"       

        DO i=2,NX-1
        DO j=2,NY-1
        DO k=2,NZ-1

         DUX_DX=    sin(yc(j))       *DUX_DR(i,j,k)  & 
                 + (cos(yc(j))/rp(i))*DUX_DT(i,j,k)

         DUX_DY=    cos(yc(j))       *DUX_DR(i,j,k)  & 
                 - (sin(yc(j))/rp(i))*DUX_DT(i,j,k)
         !---------------------------------------------------!
         DUY_DX=    sin(yc(j))       *DUY_DR(i,j,k)  &
                 + (cos(yc(j))/rp(i))*DUY_DT(i,j,k)

         DUY_DY=    cos(yc(j))       *DUY_DR(i,j,k)  &
                 - (sin(yc(j))/rp(i))*DUY_DT(i,j,k)
         !---------------------------------------------------!
         DUZ_DX=    sin(yc(j))       *DUZ_DR(i,j,k)  &
                 + (cos(yc(j))/rp(i))*DUZ_DT(i,j,k)

         DUZ_DY=    cos(yc(j))       *DUZ_DR(i,j,k)  &
                 - (sin(yc(j))/rp(i))*DUZ_DT(i,j,k)

        o_11 = 0.5d0*(DUX_DX        - DUX_DX       )
        o_12 = 0.5d0*(DUX_DY        - DUY_DX       )
        o_13 = 0.5d0*(DUX_DZ(i,j,k) - DUZ_DX       )

        o_21 = 0.5d0*(DUY_DX        - DUX_DY       )
        o_22 = 0.5d0*(DUY_DY        - DUY_DY       )
        o_23 = 0.5d0*(DUY_DZ(i,j,k) - DUZ_DY       )

        o_31 = 0.5d0*(DUZ_DX        - DUX_DZ(i,j,k))
        o_32 = 0.5d0*(DUZ_DY        - DUY_DZ(i,j,k))
        o_33 = 0.5d0*(DUZ_DZ(i,j,k) - DUZ_DZ(i,j,k))

        s_11 = 0.5d0*(DUX_DX        + DUX_DX       )
        s_12 = 0.5d0*(DUX_DY        + DUY_DX       )
        s_13 = 0.5d0*(DUX_DZ(i,j,k) + DUZ_DX       )

        s_21 = 0.5d0*(DUY_DX        + DUX_DY       )
        s_22 = 0.5d0*(DUY_DY        + DUY_DY       )
        s_23 = 0.5d0*(DUY_DZ(i,j,k) + DUZ_DY       )

        s_31 = 0.5d0*(DUZ_DX        + DUX_DZ(i,j,k))
        s_32 = 0.5d0*(DUZ_DY        + DUY_DZ(i,j,k))
        s_33 = 0.5d0*(DUZ_DZ(i,j,k) + DUZ_DZ(i,j,k))

           OMG_MAG = sqrt(  o_11**2 + o_12**2 + o_13**2 &
                          + o_21**2 + o_22**2 + o_23**2 &
                          + o_31**2 + o_32**2 + o_33**2 &
                         )

             S_MAG = sqrt(  s_11**2 + s_12**2 + s_13**2 &
                          + s_21**2 + s_22**2 + s_23**2 &
                          + s_31**2 + s_32**2 + s_33**2 &
                         )

           Q(i,j,k) = OMG_MAG - S_MAG

           write(6,*) "OMG_MAG(",i,",",j,",",k,") = ", OMG_MAG
           write(6,*) "S_MAG(",i,",",j,",",k,") = ", S_MAG
           write(6,*) "Q(",i,",",j,",",k,") = ", Q(i,j,k)
         
        ENDDO
        ENDDO
        ENDDO
        write(6,*) "HERE 6 DONE"       
  do i=1,10
     Q(i,:,:) = 0.0
  end do

  do i=nx-50,nx
     Q(i,:,:) = 0.0
  end do

  Q(:,1,:) = Q(:,ny-1,:)
  Q(:,ny,:) = Q(:,2,:)

  return
end subroutine QCri_STH_WRONG

subroutine Lamb2(u,v,w,Lambda2,nx,ny,nz,xu,zwg,rp,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: Lambda2(nx,ny,nz), MAT_L(1:3,1:3)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)

  real    ( kind = 8 ) :: s_11,s_12,s_13
  real    ( kind = 8 ) :: s_21,s_22,s_23
  real    ( kind = 8 ) :: s_31,s_32,s_33

  real    ( kind = 8 ) :: o_11,o_12,o_13
  real    ( kind = 8 ) :: o_21,o_22,o_23
  real    ( kind = 8 ) :: o_31,o_32,o_33

  ! FOR LAPACK
  CHARACTER               JOBZ, UPLO
  INTEGER ( kind = 4)  :: INFO, LDA
  INTEGER, parameter   :: LWORK = 2*(3*3-1)
  REAL    ( kind = 8)  :: WORK(LWORK), EGV(3)

  dtheta = 2.0*3.1415926/dble(ny-2)

  do k=2,nz-1
     write(6,*) " K, Lambda2(100,2,k-1) = ", K, Lambda2(100,2,k-1)
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           dq1x1=(u(i,j,k)-u(i-1,j,k))/dr
           dq2x2=(v(i,j,k)-v(i,j-1,k))/dtheta
           dq3x3=(w(i,j,k)-w(i,j,k-1))/dz

           o_11 = 0.5*(dq1x1-dq1x1)
           o_12 = 0.5*(dq1x2-dq2x1)
           o_13 = 0.5*(dq1x3-dq3x1)

           o_21 = 0.5*(dq2x1-dq1x2)
           o_22 = 0.5*(dq2x2-dq2x2)
           o_23 = 0.5*(dq2x3-dq3x2)

           o_31 = 0.5*(dq3x1-dq1x3)
           o_32 = 0.5*(dq3x2-dq2x3)
           o_33 = 0.5*(dq3x3-dq3x3)

           s_11 = 0.5*(dq1x1+dq1x1)
           s_12 = 0.5*(dq1x2+dq2x1)
           s_13 = 0.5*(dq1x3+dq3x1)

           s_21 = 0.5*(dq2x1+dq1x2)
           s_22 = 0.5*(dq2x2+dq2x2)
           s_23 = 0.5*(dq2x3+dq3x2)

           s_31 = 0.5*(dq3x1+dq1x3)
           s_32 = 0.5*(dq3x2+dq2x3)
           s_33 = 0.5*(dq3x3+dq3x3)

           MAT_L(1,1) = s_11**2 + o_11**2
           MAT_L(1,2) = s_12**2 + o_12**2
           MAT_L(1,3) = s_13**2 + o_13**2

           MAT_L(2,1) = s_21**2 + o_21**2
           MAT_L(2,2) = s_22**2 + o_22**2
           MAT_L(2,3) = s_23**2 + o_23**2

           MAT_L(3,1) = s_31**2 + o_31**2
           MAT_L(3,2) = s_32**2 + o_32**2
           MAT_L(3,3) = s_33**2 + o_33**2


           !Check symmetric
           !write(6,*) "MAT_L(1,2), MAT_L(2,1) = ", MAT_L(1,2), MAT_L(2,1)
           !!!!!!!!!!!!!!!!!

           ! STEP 2: SOLVE for Eigenvalue and Eigenvector
           !
           !      2.1) Call simple driver in LAPACK, SSYEV

           !JOBZ = 'V' ! Compute eigenvalue and eigenvector
           JOBZ = 'N' ! Compute eigenvalue only

           UPLO = 'U' ! Upper triangle of C is stored
           LDA  =  3  ! The leading dimension of array MAT_L, LDA >= max(1,N)

           ! On exit, if JOBZ='V', then if INFO = 0, C contains the orthonormal eigenvectors 
           ! of the matrix C, W contains corresponding eigenvalue
           !call SSYEV(JOBZ,UPLO,N,C,LDA,W,WORK,LWORK,INFO)
           call DSYEV(JOBZ,UPLO,3,MAT_L,LDA,EGV,WORK,LWORK,INFO)
           
           !write(6,*) "EGV(1),EGV(2),EGV(3) = ",EGV(1),EGV(2),EGV(3)

           Lambda2(i,j,k) = EGV(2)
           !write(6,*) "Lambda2(i,j,k) = ", Lambda2(i,j,k) 
        enddo
     enddo
  enddo

  ! do i=1,10
  !    Q(i,:,:) = 0.0
  ! end do

  ! do i=nx-50,nx
  !    Q(i,:,:) = 0.0
  ! end do

  Lambda2(:,1,:) = Lambda2(:,ny-1,:)
  Lambda2(:,ny,:) = Lambda2(:,2,:)

  return
end subroutine Lamb2

subroutine write_vtk(filename,var,nx,ny,nz,xc,yc,zcg,kstart,kend,jend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend,jend
  character (len = 128)   :: buff,filename
  real    ( kind = 8 ) :: var(2:nx-1,2:ny-1,2:nz-1)
!  real    ( kind = 8 ) :: xc(2:nx-1),yc(2:ny-1),zcg(2:nz-1)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  !erase previous output file
  open(unit=1, iostat=s1, file=filename, status='old')
  if (s1 == 0) close(1, status='delete')


  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx-2
  jmax = ny-2
!  kmax = kend-kstart+1
  kmax = nz-2

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,(jmax+1),kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*(jmax+1)*kmax, " float"
  write(1) buff//char(10)

!  do k = kstart, kend
   do k = 2, nz-1  
    write(6,*) "GRID: WRITE K = ", kstart, k, kend
     do j = 2, ny
        do i = 2, nx-1 
           write(1) real(xc(i)*cos(yc(j))), real(xc(i)*sin(yc(j))), real(zcg(k))   
        end do
     end do
   end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*(jmax+1)*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS vtkvar float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 2, ny
        do i = 2, nx-1
           write(1) real(var(i,j,k))
        end do
     end do
  end do


  close(1)

end subroutine write_vtk

subroutine write_vtk_cartesian(filename,var,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx-2
  jmax = ny-2
  kmax = nz-2

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, "float"
  write(1) buff//char(10)

  do k = 2,nz-1
     write(6,*) "GRID: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(xc(i)), real(yc(j)), real(zcg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS vtkvar float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var(i,j,k))
        end do
     end do
  end do


  close(1)


end subroutine write_vtk_cartesian

subroutine write_vtk_cartesian_3vars(filename,var1,var2,var3,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: var2(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: var3(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx-2
  jmax = ny-2
  kmax = nz-2

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, "float"
  write(1) buff//char(10)

  do k = 2,nz-1
     write(6,*) "GRID: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(xc(i)), real(yc(j)), real(zcg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS U float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var1(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS V float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var2(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS W float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = 2,nz-1
     write(6,*) "DATA: WRITE K = ", 2, k, nz-1
     do j = 2,ny-1
        do i = 2,nx-1
           write(1) real(var3(i,j,k))
        end do
     end do
  end do


  close(1)


end subroutine write_vtk_cartesian_3vars

subroutine write_vtk_2vars(filename,var1,var2,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(nx,ny,nz),var2(nx,ny,nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx
  jmax = ny-1
  kmax = kend-kstart+1

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = kstart, kend
     write(6,*) "GRID: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
!           write(1) real(r(i)*cos(theta(j))), real(r(i)*sin(theta(j))), real(z(k))
           write(1) real(xc(i)*cos(yc(j))), real(xc(i)*sin(yc(j))), real(zcg(k))
           !write(6,*) real(xu(i)*cos(yv(j))), real(xu(i)*sin(yv(j))), real(zwg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS var1 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var1(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS var2 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var2(i,j,k))
        end do
     end do
  end do


  close(1)


end subroutine write_vtk_2vars

subroutine write_vtk_3vars(filename,var1,var2,var3,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(nx,ny,nz),var2(nx,ny,nz),var3(nx,ny,nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx
  jmax = ny-1
  kmax = kend-kstart+1

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = kstart, kend
     write(6,*) "GRID: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
!           write(1) real(r(i)*cos(theta(j))), real(r(i)*sin(theta(j))), real(z(k))
           write(1) real(xc(i)*cos(yc(j))), real(xc(i)*sin(yc(j))), real(zcg(k))
           !write(6,*) real(xu(i)*cos(yv(j))), real(xu(i)*sin(yv(j))), real(zwg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS var1 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var1(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS var2 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var2(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS var3 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var3(i,j,k))
        end do
     end do
  end do


  close(1)

end subroutine write_vtk_3vars

subroutine write_vtk_4vars(filename,var1,var2,var3,var4,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(nx,ny,nz),var2(nx,ny,nz)
  real    ( kind = 8 ) :: var3(nx,ny,nz),var4(nx,ny,nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx
  jmax = ny-1
  kmax = kend-kstart+1

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = kstart, kend
     write(6,*) "GRID: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
!           write(1) real(r(i)*cos(theta(j))), real(r(i)*sin(theta(j))), real(z(k))
           write(1) real(xc(i)*cos(yc(j))), real(xc(i)*sin(yc(j))), real(zcg(k))
           !write(6,*) real(xu(i)*cos(yv(j))), real(xu(i)*sin(yv(j))), real(zwg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "SCALARS var1 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA_VAR1: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var1(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS var2 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA_VAR2: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var2(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS var3 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA_VAR3: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var3(i,j,k))
        end do
     end do
  end do

  write(1) "SCALARS var4 float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA_VAR4: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var4(i,j,k))
        end do
     end do
  end do

  close(1)

end subroutine write_vtk_4vars

subroutine write_vtk_vector(filename,var1,var2,var3,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var1(nx,ny,nz),var2(nx,ny,nz),var3(nx,ny,nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx
  jmax = ny-1
  kmax = kend-kstart+1

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,jmax,kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*jmax*kmax, " float"
  write(1) buff//char(10)

  do k = kstart, kend
     write(6,*) "GRID: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
!           write(1) real(r(i)*cos(theta(j))), real(r(i)*sin(theta(j))), real(z(k))
           write(1) real(xc(i)*cos(yc(j))), real(xc(i)*sin(yc(j))), real(zcg(k))
           !write(6,*) real(xu(i)*cos(yv(j))), real(xu(i)*sin(yv(j))), real(zwg(k))
        end do
     end do
  end do

  write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
  write(1) char(10)//buff//char(10)

  write(1) "VECTORS var float 1"//char(10)
  write(1) "LOOKUP_TABLE default"//char(10)
  do k = kstart, kend
     write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 1, jmax
        do i = 1, imax
           write(1) real(var2(i,j,k)),real(var1(i,j,k)),real(var3(i,j,k))
        end do
     end do
  end do
  close(1)

end subroutine write_vtk_vector

subroutine write_vtk_2d(filename,var,nx,ny,nz,xc,yc,zcg)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend,np1,np2
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: var(nx,ny,nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)
  real(4)              :: g1vtk(nx),g2vtk(nz),vtkvar(nx,nz)
  character(len=25)    :: ss
  character(len=80)    :: basename

  np1 = nx
  np2 = nz
  g1vtk = real(xc)
  g2vtk = real(zcg)
  vtkvar(:,:) = real(var(:,3,:))

  basename = "VAR"
  
  open(unit=13,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  !HEADER: note termination with char(10)
  write(13) "# vtk DataFile Version 3.0"//char(10)
  write(13) trim(basename)//char(10)
  write(13) "BINARY"//char(10)
  write(13) "DATASET RECTILINEAR_GRID"//char(10)

  ! if (dir.EQ.1) then 
  !  write(ss,fmt='(A10,3I5)') "DIMENSIONS",1,np1,np2
  !  write(13) ss//char(10)
  !  !X-grid
  !  write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",1," float"
  !  write(13) char(10)//ss//char(10)
  !  if (iu.EQ.1) then
  !   write(13) eL
  !  else
  !   write(13) cL
  !  endif
  !  !Y-grid
  !  write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np1," float"
  !  write(13) char(10)//ss//char(10)
  !  do j = 1, np1
  !   write(13) g1vtk(j)
  !  enddo
  !  !Z-grid
  !  write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
  !  write(13) char(10)//ss//char(10)
  !  do k = 1, np2
  !   write(13) g2vtk(k)
  !  enddo
  ! elseif (dir.EQ.2) then 
   write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,1,np2
   write(13) ss//char(10)
   !X-grid
   write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
   write(13) char(10)//ss//char(10)
   do i = 1, np1
    write(13) g1vtk(i)
   enddo
   !Y-grid
   write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",1," float"
   write(13) char(10)//ss//char(10)
   !if (iv.EQ.1) then
   !write(13) eL
   !else
   write(13) real(2)
   ! write(13) cL
   !endif
   !Z-grid
   write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",np2," float"
   write(13) char(10)//ss//char(10)
   do k = 1, np2
    write(13) g2vtk(k)
   enddo
  ! elseif (dir.EQ.3) then 
  !  write(ss,fmt='(A10,3I5)') "DIMENSIONS",np1,np2,1
  !  write(13) ss//char(10)
  !  !X-grid
  !  write(ss,fmt='(A13,I6,A6)') "X_COORDINATES",np1," float"
  !  write(13) char(10)//ss//char(10)
  !  do i = 1, np1
  !   write(13) g1vtk(i)
  !  enddo
  !  !Y-grid
  !  write(ss,fmt='(A13,I6,A6)') "Y_COORDINATES",np2," float"
  !  write(13) char(10)//ss//char(10)
  !  do j = 1, np2
  !   write(13) g2vtk(j)
  !  enddo
  !  !Z-grid
  !  write(ss,fmt='(A13,I6,A6)') "Z_COORDINATES",1," float"
  !  write(13) char(10)//ss//char(10)
  !  if (iw.EQ.1) then
  !   write(13) eL
  !  else
  !   write(13) cL
  !  endif
  ! endif

   !Field
   write(ss,fmt='(A10,I15)') "POINT_DATA",np1*np2
   write(13) char(10)//ss//char(10)
   write(13) "SCALARS Fraction float 1"//char(10)
   write(13) "LOOKUP_TABLE default"//char(10)
   write(13) vtkvar
   !Close VTK File
  close(13)

!endif


end subroutine write_vtk_2d

SUBROUTINE DERIVATIVE(VAR,DERI,NX,NY,NZ,XU,ZWG,DIR,IVAR)
  implicit none

  INTEGER ( kind = 4 ) :: NX,NY,NZ
  REAL    ( kind = 8 ) :: XU(NX),ZWG(NZ)
  INTEGER ( kind = 4 ) :: I,J,K,DIR,IVAR
  REAL    ( kind = 8 ) :: VAR(NX,NY,NZ),DERI(NX,NY,NZ)
  REAL    ( kind = 8 ) :: dtheta,dr,dz

  dtheta = 2.0*3.1415926/dble(ny-2)  

  IF(DIR.EQ.1) THEN

     do k=2,NZ-1
        !write(6,*) " DIR 1: K = ", K
        do j=2,NY-1
           do i=2,NX-1
              dr=xu(i)-xu(i-1)
              
              !write(6,*) "xu(",i,"), xu(",i-1,"), dr = ", xu(i),xu(i-1), dr


              IF(IVAR.EQ.1) THEN
                 DERI(I,J,K) = (VAR(I,J,K)-VAR(I-1,J,K))/dr

                 !if(DERI(I,J,K) .gt. 1000000) write(6,*) "IVAR1,DIR1: I,J,K,dr = ", I,J,K,dr

              ELSEIF(IVAR.EQ.2) THEN
                 DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i+1,j,k)+VAR(i+1,j-1,k)) &
                              -0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i-1,j,k)+VAR(i-1,j-1,k)))/dr

                 !if(DERI(I,J,K) .gt. 1000000) write(6,*) "IVAR2,DIR1: I,J,K = ", I,J,K

              ELSEIF(IVAR.EQ.3) THEN
                 DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j,k-1)+VAR(i+1,j,k)+VAR(i+1,j,k-1)) &
                              -0.25*(VAR(i,j,k)+VAR(i,j,k-1)+VAR(i-1,j,k)+VAR(i-1,j,k-1)))/dr
              ENDIF
           enddo
        enddo
     enddo

  ELSEIF(DIR.EQ.2) THEN

     do k=2,NZ
        !write(6,*) " DIR 2: K = ", K
        do j=2,NY
           do i=2,NX-1
              IF(IVAR.EQ.1) THEN
                 DERI(I,J,K) = (0.25*(VAR(i,j,k)+VAR(i,j+1,k)+VAR(i-1,j,k)+VAR(i-1,j+1,k)) &
                               -0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i-1,j,k)+VAR(i-1,j-1,k)))/dtheta

                 !if(DERI(I,J,K) .gt. 1000) write(6,*) "IVAR1,DIR2: I,J,K = ", I,J,K

              ELSEIF(IVAR.EQ.2) THEN
                 DERI(I,J,K) = (VAR(I,J,K)-VAR(I,J-1,K))/dtheta

                 !if(DERI(I,J,K) .gt. 1000) write(6,*) "IVAR2,DIR2: I,J,K = ", I,J,K
              ELSEIF(IVAR.EQ.3) THEN
                 DERI(I,J,K) = (0.25*(VAR(i,j,k)+VAR(i,j+1,k)+VAR(i,j,k-1)+VAR(i,j+1,k-1)) &
                               -0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i,j,k-1)+VAR(i,j-1,k-1)))/dtheta
              ENDIF
           enddo
        enddo
     enddo

  ELSEIF(DIR.EQ.3) THEN

     do k=2,NZ
        write(6,*) " DIR 3: K = ", K
        do j=2,NY
           do i=2,NX-1
              dz=zwg(k)-zwg(k-1)
              IF(IVAR.EQ.1) THEN
                 DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j,k+1)+VAR(i-1,j,k)+VAR(i-1,j,k+1)) &
                              -0.25*(VAR(i,j,k)+VAR(i,j,k-1)+VAR(i-1,j,k)+VAR(i-1,j,k-1)))/dz
              ELSEIF(IVAR.EQ.2) THEN
                 DERI(I,J,K) =(0.25*(VAR(i,j,k)+VAR(i,j-1,k)+VAR(i,j,k+1)+VAR(i,j-1,k+1)) &
                              -0.25*(VAR(i,j,k)+ VAR(i,j-1,k)+VAR(i,j,k-1)+VAR(i,j-1,k-1)))/dz
              ELSEIF(IVAR.EQ.3) THEN
                 if(k .eq. 2306) write(6,*) " DIR 3: I,J,K = ", I,J,K
                 DERI(I,J,K) = (VAR(I,J,K)-VAR(I,J,K-1))/dz
              ENDIF
           enddo
        enddo
     enddo


  ENDIF

  RETURN
END SUBROUTINE DERIVATIVE

! subroutine DUYDY_cal(u,v,DUYDY,nx,ny,nz,xc,yc,xu,zwg)
!   implicit none

!   integer ( kind = 4 ) :: i,j,k,nx,ny,nz

!   real    ( kind = 8 ) :: xu(nx),zwg(nz)
!   real    ( kind = 8 ) :: xc(nx),yc(ny)

!   real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
!   real    ( kind = 8 ) :: DUYDY(nx,ny,nz)
!   real    ( kind = 8 ) :: DUDR(nx,ny,nz),DVDR(nx,ny,nz)
!   real    ( kind = 8 ) :: DUDTHETA(nx,ny,nz),DVDTHETA(nx,ny,nz)


!   call DERIVATIVE(u,DUDR    ,NX,NY,NZ,XU,ZWG,1,1)
!   call DERIVATIVE(v,DVDR    ,NX,NY,NZ,XU,ZWG,1,2)
!   call DERIVATIVE(u,DUDTHETA,NX,NY,NZ,XU,ZWG,2,1)
!   call DERIVATIVE(v,DVDTHETA,NX,NY,NZ,XU,ZWG,2,2)


!   DO k=2,NZ
!      write(6,*) " K, DUYDY(100,2,k-1) = ", K, DUYDY(100,2,k-1)
!      DO j=2,NY
!         DO i=2,NX-1

!         DUYDY(i,j,k) = sin(yc(j))**2.0d0*DUDR(i,j,k)+ &
!                        sin(yc(j))*cos(yc(j))*DVDR(i,j,k)+ &
!                        sin(yc(j))*cos(yc(j))/xc(i)*DUDTHETA(i,j,k) + &
!                        u(i,j,k)/xc(i)*cos(yc(j))**2.0d0 + &
!                        cos(yc(j))**2.0d0/xc(i)*DVDTHETA(i,j,k) - &
!                       -v(i,j,k)/xc(i)*sin(yc(j))*cos(yc(j))

        
!         !if(DUYDY(I,J,K) .gt. 100) write(6,*) "DUYDY: I,J,K = ", I,J,K

!         ENDDO
!      ENDDO
!   ENDDO


!   return
! end subroutine DUYDY_cal

subroutine gradient(var,gv_r,gv_t,gv_z,nx,ny,nz,xu,zwg,rp,kstart,kend)
  ! var is at cell-center
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: var(nx,ny,nz)
  real    ( kind = 8 ) :: gv_r(nx,ny,nz),gv_t(nx,ny,nz),gv_z(nx,ny,nz)
  real    ( kind = 8 ) :: dz,dr,dtheta
  real    ( kind = 8 ) :: xu(nx),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)
  
  do k=2,nz-1
     write(6,*) " GRADIENT: K = ", K
     do j=2,ny-1
        do i=2,nx-1
  
           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           gv_r(i,j,k)=(0.5*(var(i,j,k)+var(i+1,j,k)) &
                       -0.5*(var(i,j,k)+var(i-1,j,k))  )/dr

           gv_t(i,j,k)=((0.5*(var(i,j,k)+var(i,j+1,k)) &
                        -0.5*(var(i,j,k)+var(i,j-1,k)))/dtheta)/rp(i)
       
           gv_z(i,j,k)=(0.5*(var(i,j,k)+var(i,j,k+1)) &
                       -0.5*(var(i,j,k)+var(i,j,k-1))  )/dz

           !write(6,*) "gv_r, gv_t, gv_z(",i,",",j,",",k,") = ", &
           ! write(6,*) "gv_r, gv_t, gv_z = ", &
           !             gv_r(i,j,k), gv_t(i,j,k), gv_z(i,j,k)
                
           ! if(gv_r(i,j,k) > 100.) write(6,*) i,j,k,"gv_r = ", gv_r(i,j,k)
           ! if(gv_t(i,j,k) > 100.) write(6,*) i,j,k,"gv_t = ", gv_t(i,j,k)
           ! if(gv_z(i,j,k) > 100.) write(6,*) i,j,k,"gv_z = ", gv_z(i,j,k)
           ! if(gv_r(i,j,k) < -100.) write(6,*) i,j,k,"gv_r = ", gv_r(i,j,k)
           ! if(gv_t(i,j,k) < -100.) write(6,*) i,j,k,"gv_t = ", gv_t(i,j,k)
           ! if(gv_z(i,j,k) < -100.) write(6,*) i,j,k,"gv_z = ", gv_z(i,j,k)

        enddo
     enddo
  enddo
  return
end subroutine gradient

subroutine domg_ver_dt(u0,v0,w0,u2,v2,w2,time0,time2,domg_ver,nx,ny,nz,xu,zwg,rp,yc)
  ! calculate (partial vorticity_ver / partial t)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u0(nx,ny,nz),v0(nx,ny,nz),w0(nx,ny,nz)
  real    ( kind = 8 ) :: u2(nx,ny,nz),v2(nx,ny,nz),w2(nx,ny,nz)
  real    ( kind = 8 ) :: domg_ver(nx,ny,nz),omg_ver0(nx,ny,nz),omg_ver2(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3
  real    ( kind = 8 ) :: time0,time2,t2mt0

  real    ( kind = 8 ) :: omgt,omgr,omgz

  ! real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)
  
  real    ( kind = 8 ) :: xu(nx),yc(ny),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)
  t2mt0  = time2 - time0

  do k=2,nz-1
     write(6,*) " DVORT_VER_DT(OMG_VER0): K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w0(i,j,k)+w0(i,j+1,k)+w0(i,j,k-1)+w0(i,j+1,k-1)) &
                 -0.25*(w0(i,j,k)+w0(i,j-1,k)+w0(i,j,k-1)+w0(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v0(i,j,k)+v0(i,j-1,k)+v0(i,j,k+1)+v0(i,j-1,k+1)) &
                 -0.25*(v0(i,j,k)+v0(i,j-1,k)+v0(i,j,k-1)+v0(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u0(i,j,k)+u0(i,j,k+1)+u0(i-1,j,k)+u0(i-1,j,k+1)) &
                 -0.25*(u0(i,j,k)+u0(i,j,k-1)+u0(i-1,j,k)+u0(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w0(i,j,k)+w0(i,j,k-1)+w0(i+1,j,k)+w0(i+1,j,k-1)) &
                 -0.25*(w0(i,j,k)+w0(i,j,k-1)+w0(i-1,j,k)+w0(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v0(i,j,k)+v0(i,j-1,k))*rp(i)+(v0(i+1,j,k)+v0(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v0(i,j,k)+v0(i,j-1,k))*rp(i)+(v0(i-1,j,k)+v0(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u0(i,j,k)+u0(i,j+1,k)+u0(i-1,j,k)+u0(i-1,j+1,k)) &
                 -0.25*(u0(i,j,k)+u0(i,j-1,k)+u0(i-1,j,k)+u0(i-1,j-1,k)))/dtheta

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)
           
           omg_ver0(i,j,k)=omgr*sin(yc(j))+omgt*cos(yc(j))    

        enddo
     enddo
  enddo

  omg_ver0(:,1,:) = omg_ver0(:,ny-1,:)
  omg_ver0(:,ny,:) = omg_ver0(:,2,:)

  do k=2,nz-1
     write(6,*) " DVORT_VER_DT(OMG_VER2): K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w2(i,j,k)+w2(i,j+1,k)+w2(i,j,k-1)+w2(i,j+1,k-1)) &
                 -0.25*(w2(i,j,k)+w2(i,j-1,k)+w2(i,j,k-1)+w2(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v2(i,j,k)+v2(i,j-1,k)+v2(i,j,k+1)+v2(i,j-1,k+1)) &
                 -0.25*(v2(i,j,k)+v2(i,j-1,k)+v2(i,j,k-1)+v2(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u2(i,j,k)+u2(i,j,k+1)+u2(i-1,j,k)+u2(i-1,j,k+1)) &
                 -0.25*(u2(i,j,k)+u2(i,j,k-1)+u2(i-1,j,k)+u2(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w2(i,j,k)+w2(i,j,k-1)+w2(i+1,j,k)+w2(i+1,j,k-1)) &
                 -0.25*(w2(i,j,k)+w2(i,j,k-1)+w2(i-1,j,k)+w2(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v2(i,j,k)+v2(i,j-1,k))*rp(i)+(v2(i+1,j,k)+v2(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v2(i,j,k)+v2(i,j-1,k))*rp(i)+(v2(i-1,j,k)+v2(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u2(i,j,k)+u2(i,j+1,k)+u2(i-1,j,k)+u2(i-1,j+1,k)) &
                 -0.25*(u2(i,j,k)+u2(i,j-1,k)+u2(i-1,j,k)+u2(i-1,j-1,k)))/dtheta

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)
           
           omg_ver2(i,j,k)=omgr*sin(yc(j))+omgt*cos(yc(j))    

        enddo
     enddo
  enddo

  omg_ver2(:,1,:) = omg_ver2(:,ny-1,:)
  omg_ver2(:,ny,:) = omg_ver2(:,2,:)

  do k=2,nz-1
     write(6,*) " DVORT_VER_DT(DOMG_VER): K = ", K
     do j=2,ny-1
        do i=2,nx-1
           domg_ver(i,j,k) = (omg_ver2(i,j,k)-omg_ver0(i,j,k))/t2mt0
           !write(6,*) "omg_ver2,omg_ver1,domg_ver = ", omg_ver2(i,j,k),omg_ver0(i,j,k),domg_ver(i,j,k)
        enddo
     enddo
  enddo

  domg_ver(:,1,:) = domg_ver(:,ny-1,:)
  domg_ver(:,ny,:) = domg_ver(:,2,:)

  ! do i=1,10
  !    omg_ver0(i,:,:) = 0.0
  ! end do

  ! do i=nx-50,nx
  !    omg_ver0(i,:,:) = 0.0
  ! end do

  return
end subroutine domg_ver_dt

subroutine advection_omg_ver(u,v,w,adv_omg_ver,nx,ny,nz,xu,zwg,rp,yc,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: adv_omg_ver(nx,ny,nz),omg_ver(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: omgt,omgr,omgz
  real    ( kind = 8 ) :: gov_r(nx,ny,nz),gov_t(nx,ny,nz),gov_z(nx,ny,nz)
                          ! gov = gradient omega vertical


  ! real    ( kind = 8 ) :: dq1x1_g(nx,ny,nz),dq1x2_g(nx,ny,nz),dq1x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq2x1_g(nx,ny,nz),dq2x2_g(nx,ny,nz),dq2x3_g(nx,ny,nz)
  ! real    ( kind = 8 ) :: dq3x1_g(nx,ny,nz),dq3x2_g(nx,ny,nz),dq3x3_g(nx,ny,nz)
  
  real    ( kind = 8 ) :: xu(nx),yc(ny),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)

  do k=2,nz-1
     write(6,*) " ADV_OMG_VER(OMG_VER): K = ", K
     do j=2,ny-1
        do i=2,nx-1

           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)
           
           omg_ver(i,j,k)=omgr*sin(yc(j))+omgt*cos(yc(j))    

        enddo
     enddo
  enddo

  omg_ver(:,1,:) = omg_ver(:,ny-1,:)
  omg_ver(:,ny,:) = omg_ver(:,2,:)

  call gradient(omg_ver,gov_r,gov_t,gov_z,nx,ny,nz,xu,zwg,rp,kstart,kend)

  do k=2,nz-1
     write(6,*) " ADV_OMG_VER(ADV_OMG_VER): K = ", K
     do j=2,ny-1
        do i=2,nx-1
           adv_omg_ver(i,j,k) = 0.5*(u(i,j,k)+u(i-1,j,k))*gov_r(i,j,k) + &
                                0.5*(v(i,j,k)+v(i,j-1,k))*gov_t(i,j,k) + &
                                0.5*(w(i,j,k)+w(i,j,k-1))*gov_z(i,j,k) 
           !write(6,*) "adv_omg_ver(i,j,k) = ", adv_omg_ver(i,j,k)
        enddo
     enddo
  enddo

  adv_omg_ver(:,1,:) = adv_omg_ver(:,ny-1,:)
  adv_omg_ver(:,ny,:) = adv_omg_ver(:,2,:)

  return
end subroutine advection_omg_ver

subroutine baroclinic_ver(ov_y,pres,dens,nx,ny,nz,xu,zwg,rp,yc,kstart,kend)
  implicit none
  
  ! ov is (d x p)/dens^2 term 
  ! p  is pressure gradient
  ! d  is density gradient

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: ov_r,ov_t,ov_z(nx,ny,nz)
  real    ( kind = 8 ) :: ov_y(nx,ny,nz)
  real    ( kind = 8 ) :: pres(nx,ny,nz),dens(nx,ny,nz)

  real    ( kind = 8 ) :: p_r(nx,ny,nz),p_t(nx,ny,nz),p_z(nx,ny,nz)
  real    ( kind = 8 ) :: d_r(nx,ny,nz),d_t(nx,ny,nz),d_z(nx,ny,nz)

  real    ( kind = 8 ) :: xu(nx),yc(ny),zwg(nz),rp(nx)

  call gradient(pres,p_r,p_t,p_z,nx,ny,nz,xu,zwg,rp,kstart,kend)
  call gradient(dens,d_r,d_t,d_z,nx,ny,nz,xu,zwg,rp,kstart,kend)

  do k=2,nz-1
     write(6,*) " BAROCLINIC_VER: K = ", K
     do j=2,ny-1
        do i=2,nx-1
  
           ov_r        = d_t(i,j,k)*p_z(i,j,k) - p_t(i,j,k)*d_z(i,j,k)
           ov_t        = p_r(i,j,k)*d_z(i,j,k) - d_r(i,j,k)*p_z(i,j,k)
           !ov_z(i,j,k) = d_r(i,j,k)*p_t(i,j,k) - p_r(i,j,k)*d_t(i,j,k)

           ov_y(i,j,k) = (ov_r*sin(yc(j))+ov_t*cos(yc(j)))!/dens(i,j,k)**2
           
        enddo
     enddo
  enddo

  ov_y(:,1 ,:) = ov_y(:,ny-1,:)
  ov_y(:,ny,:) = ov_y(:,2   ,:)

  return
end subroutine baroclinic_ver

subroutine stretch_vertical(u,v,w,stretch_ver,nx,ny,nz,xu,zwg,rp,yc,kstart,kend)
  implicit none
  
  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  real    ( kind = 8 ) :: uo(nx,ny,nz),vo(nx,ny,nz)
  real    ( kind = 8 ) :: stretch_ver(nx,ny,nz)
  real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta,sine,cosine
  real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
  real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3

  real    ( kind = 8 ) :: omgt,omgr,omgz
  
  real    ( kind = 8 ) :: xu(nx),yc(ny),zwg(nz),rp(nx)

  dtheta = 2.0*3.1415926/dble(ny-2)

  CALL CENTER_VELOCITY(nx,ny,nz,uo,u,1)
  CALL CENTER_VELOCITY(nx,ny,nz,vo,v,2)

  do k=2,nz-1
     write(6,*) " STRETCH_VER: K = ", K
     do j=2,ny-1
        do i=2,nx-1
           
           dz=zwg(k)-zwg(k-1)
           dr=xu(i)-xu(i-1)

           dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
           dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
                 -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

           dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
                 -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
           dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
                 -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

        dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
              -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
           dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
                 -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

           dq1x1=(u(i,j,k)-u(i-1,j,k))/dr
           dq2x2=(v(i,j,k)-v(i,j-1,k))/dtheta
           dq3x3=(w(i,j,k)-w(i,j,k-1))/dz

           omgr = dq3x2/rp(i)-dq2x3
           omgt = dq1x3-dq3x1
           omgz = (dq2x1-dq1x2)/rp(i)
           
           sine   = sin(yc(j))
           cosine = cos(yc(j))

           stretch_ver(i,j,k) = omgr*sine*dq1x1 + omgr*cosine*dq2x1 + &
                       (omgt/rp(i))*sine*dq1x2 + (omgt/rp(i))*cosine*uo(i,j,k) - &
                       (omgt/rp(i))*vo(i,j,k)*sine + (omgt/rp(i))*cosine*dq2x2 + &
                       omgz*sine*dq1x3 + omgz*cosine*dq2x3

        enddo
     enddo
  enddo

  ! do i=1,10
  !    stretch_ver(i,:,:) = 0.0
  ! end do

  do i=nx-50,nx
     stretch_ver(i,:,:) = 0.0
  end do

  stretch_ver(:,1,:) = stretch_ver(:,ny-1,:)
  stretch_ver(:,ny,:) = stretch_ver(:,2,:)

  return
end subroutine stretch_vertical

! subroutine omg_diffusion(u,v,w,diffu_st,diffu_sw,diffu_vr,nx,ny,nz,xu,zwg,rp,yc,xc,zcg,kstart,kend,Re)
!   ! diffu_vr = diffusion term vertical
!   ! diffu_sw = diffusion term horizontal cross stream (spanwise)
!   ! diffu_st = diffusion term horizontal streamwise
!   implicit none
  
!   integer ( kind = 4 ) :: i,j,k,nx,ny,nz,kstart,kend
!   real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
!   real    ( kind = 8 ) :: uo(nx,ny,nz),vo(nx,ny,nz)
!   real    ( kind = 8 ) :: stretch_ver(nx,ny,nz)
!   real    ( kind = 8 ) :: dq1x3,dq3x1,dz,dr,dtheta,sine,cosine,Re
!   real    ( kind = 8 ) :: dq2x3,dq3x2,dq2x1,dq1x2
!   real    ( kind = 8 ) :: dq1x1,dq2x2,dq3x3
!   real    ( kind = 8 ) :: dz_l,dz_r,dr_d,dr_u

!   real    ( kind = 8 ) :: dor_drdr,dot_drdr,dor_dt,dor_dtdt
!   real    ( kind = 8 ) :: dot_dt,dot_dtdt,dor_dzdz,dot_dzdz
!   real    ( kind = 8 ) :: domg_ver_drdr, domg_ver_dtdt, domg_ver_dzdz

!   real    ( kind = 8 ) :: omgt(nx,ny,nz),omgr(nx,ny,nz),omgz(nx,ny,nz)

!   real    ( kind = 8 ) :: diffu_omgr,diffu_omgt,diffu_omgz
!   real    ( kind = 8 ) :: diffu_sw(nx,ny,nz),diffu_vr(nx,ny,nz),diffu_st(nx,ny,nz)
  
!   real    ( kind = 8 ) :: xu(nx),yc(ny),zwg(nz),rp(nx),xc(nx),zcg(nz)

!   dtheta = 2.0*3.1415926/dble(ny-2)

!   do k=2,nz-1
!      write(6,*) " diffusion: K = ", K
!      do j=2,ny-1
!         do i=2,nx-1
           
!            dz=zwg(k)-zwg(k-1)
!            dr=xu(i)-xu(i-1)

!            dz_l=zcg(k  )-zcg(k-1)
!            dz_r=zcg(k+1)-zcg(k  )

!            dr_d=xc(i  )-xc(i-1)
!            dr_u=xc(i+1)-xc(i  )

!            dq3x2=(0.25*(w(i,j,k)+w(i,j+1,k)+w(i,j,k-1)+w(i,j+1,k-1)) &
!                  -0.25*(w(i,j,k)+w(i,j-1,k)+w(i,j,k-1)+w(i,j-1,k-1)))/dtheta
!            dq2x3=(0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1)) &
!                  -0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1)))/dz

!            dq1x3=(0.25*(u(i,j,k)+u(i,j,k+1)+u(i-1,j,k)+u(i-1,j,k+1)) &
!                  -0.25*(u(i,j,k)+u(i,j,k-1)+u(i-1,j,k)+u(i-1,j,k-1)))/dz
!            dq3x1=(0.25*(w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1)) &
!                  -0.25*(w(i,j,k)+w(i,j,k-1)+w(i-1,j,k)+w(i-1,j,k-1)))/dr

!         dq2x1=(0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i+1,j,k)+v(i+1,j-1,k))*rp(i+1)) &
!               -0.25*((v(i,j,k)+v(i,j-1,k))*rp(i)+(v(i-1,j,k)+v(i-1,j-1,k))*rp(i-1)))/dr
!            dq1x2=(0.25*(u(i,j,k)+u(i,j+1,k)+u(i-1,j,k)+u(i-1,j+1,k)) &
!                  -0.25*(u(i,j,k)+u(i,j-1,k)+u(i-1,j,k)+u(i-1,j-1,k)))/dtheta

!            dq1x1=(u(i,j,k)-u(i-1,j,k))/dr
!            dq2x2=(v(i,j,k)-v(i,j-1,k))/dtheta
!            dq3x3=(w(i,j,k)-w(i,j,k-1))/dz

!            omgr(i,j,k) = dq3x2/rp(i)-dq2x3
!            omgt(i,j,k) = dq1x3-dq3x1
!            omgz(i,j,k) = (dq2x1-dq1x2)/rp(i)

!         enddo
!      enddo
!   enddo

!   diffu_vr = 0.0
!   ! diffu_sw = 0.0
!   ! diffu_st = 0.0

!   do k=3,nz-2
!      write(6,*) " diffusion: K = ", K
!      do j=3,ny-2
!         do i=3,nx-2
           
!            dz=zwg(k)-zwg(k-1)
!            dr=xu(i)-xu(i-1)

!            dz_l=zcg(k  )-zcg(k-1)
!            dz_r=zcg(k+1)-zcg(k  )

!            dr_d=xc(i  )-xc(i-1)
!            dr_u=xc(i+1)-xc(i  )

!            sine   = sin(yc(j))
!            cosine = cos(yc(j))

!            dor_dt   =(  0.5*(omgr(i,j,k)+omgr(i,j+1,k)) &
!                       - 0.5*(omgr(i,j,k)+omgr(i,j-1,k))  )/dtheta

!            dot_dt   =(  0.5*(omgt(i,j,k)+omgt(i,j+1,k)) &
!                       - 0.5*(omgt(i,j,k)+omgt(i,j-1,k))  )/dtheta

!            dor_drdr =(  (omgr(i+1,j,k)-omgr(i  ,j,k))/dr_u &
!                       - (omgr(i  ,j,k)-omgr(i-1,j,k))/dr_d  )/dr

!            dot_drdr =(  (omgt(i+1,j,k)-omgt(i  ,j,k))/dr_u &
!                       - (omgt(i  ,j,k)-omgt(i-1,j,k))/dr_d  )/dr

!            dor_dtdt =(  (omgr(i,j+1,k)-omgr(i,j  ,k))/dtheta &
!                       - (omgr(i,j  ,k)-omgr(i,j-1,k))/dtheta  )/dtheta

!            dot_dtdt =(  (omgt(i,j+1,k)-omgt(i,j  ,k))/dtheta &
!                       - (omgt(i,j  ,k)-omgt(i,j-1,k))/dtheta  )/dtheta

!            dor_dzdz =(  (omgr(i,j,k+1)-omgr(i,j,k  ))/dz_r &
!                       - (omgr(i,j,k  )-omgr(i,j,k-1))/dz_l &  )/dz

!            dot_dzdz =(  (omgt(i,j,k+1)-omgt(i,j,k  ))/dz_r &
!                       - (omgt(i,j,k  )-omgt(i,j,k-1))/dz_l &  )/dz

!            domg_ver_drdr = sine*dor_drdr + cosine*dot_drdr

!            domg_ver_dtdt = - omgr(i,j,k)*sine + cosine*dor_dt &
!                            + sine*dor_dtdt + dor_dt*cosine    &
!                            - omgt(i,j,k)*cosine - sine*dot_dt &
!                            + cosine*dot_dtdt - dot_dt*sine

!            domg_ver_dzdz = sine*dor_dzdz + cosine*dot_dzdz

!            diffu_vr(i,j,k) = domg_ver_drdr + domg_ver_dtdt + domg_ver_dzdz

!            !write(6,*) domg_ver_drdr, domg_ver_dtdt, domg_ver_dzdz, diffu_vr(i,j,k)

!            ! diffu_omgr = (((omgr(i+1,j,k)-omgr(i  ,j,k))/dr_u) - &
!            !               ((omgr(i  ,j,k)-omgr(i-1,j,k))/dr_d)    )/dr
!            ! diffu_omgr = diffu_omgr/Re


!            ! diffu_omgt = ((((omgt(i,j+1,k)-omgt(i,j  ,k))/dtheta) - &
!            !                ((omgt(i,j  ,k)-omgt(i,j-1,k))/dtheta)    )/dtheta)/(rp(i)**2)
!            ! diffu_omgt = diffu_omgt/Re


!            ! diffu_omgz = (((omgz(i,j,k+1)-omgz(i,j,k  ))/dz_r) - &
!            !               ((omgz(i,j,k  )-omgz(i,j,k-1))/dz_l)    )/dz
!            ! diffu_omgz = diffu_omgz/Re


!            ! diffu_vr(i,j,k) = diffu_omgr*sine   + diffu_omgt*cosine
!            ! diffu_sw(i,j,k) = diffu_omgr*cosine - diffu_omgt*sine
!            ! diffu_st(i,j,k) = diffu_omgz

!         enddo
!      enddo
!   enddo

!   ! do i=nx-50,nx
!   !    diffu_vr(i,:,:) = 0.0
!   ! end do

!   ! diffu_vr(:,1,:) = diffu_vr(:,ny-1,:)
!   ! diffu_vr(:,ny,:) = diffu_vr(:,2,:)

!   ! do i=nx-50,nx
!   !    diffu_sw(i,:,:) = 0.0
!   ! end do

!   ! diffu_sw(:,1,:) = diffu_sw(:,ny-1,:)
!   ! diffu_sw(:,ny,:) = diffu_sw(:,2,:)

!   ! do i=nx-50,nx
!   !    diffu_st(i,:,:) = 0.0
!   ! end do

!   ! diffu_st(:,1,:) = diffu_st(:,ny-1,:)
!   ! diffu_st(:,ny,:) = diffu_st(:,2,:)

!   return
! end subroutine omg_diffusion


SUBROUTINE INT2DFIELD(INPUT,OUTPUT_HL,OUTPUT_HR,NX,NY,NZ,XU,XC)
  ! OUTPUT_HL is output half-left
  ! OUTPUT_HL is output half-right
  implicit none
  INTEGER ( kind = 4) :: NX,NY,NZ,nymod
  INTEGER ( kind = 4) :: I,J,K
  REAL    ( kind = 8) :: XU(NX),XC(NX)
  REAL    ( kind = 8) :: INPUT(NX,NY,NZ),OUTPUT_HL(NZ),OUTPUT_HR(NZ)
  REAL    ( kind = 8) :: dtheta,dr

  nymod=ny-2
  dtheta = 2.0*3.1415926/dble(ny-2)
  OUTPUT_HL = 0.0
  OUTPUT_HR = 0.0

  DO K=2,NZ-1
     write(6,*) " INT2DFIELD: K = ", K
     DO I=2,NX-1
        dr=xu(i)-xu(i-1)

        DO J=2,NY
           if(j >= 2 .and. j <= (nymod/4 + 1)) then
              OUTPUT_HL(K)=OUTPUT_HL(K)+INPUT(I,J,K)*dr*XC(I)*dtheta
              !write(6,*) "OUTPUT_HL HERE1"
           elseif(j >= (nymod/4 + 2) .and. j <= (nymod/2)+(nymod/4)+1) then
              OUTPUT_HR(K)=OUTPUT_HR(K)+INPUT(I,J,K)*dr*XC(I)*dtheta
              !write(6,*) "OUTPUT_HR"
           elseif(j >= (nymod/2)+(nymod/4)+2) then
              OUTPUT_HL(K)=OUTPUT_HL(K)+INPUT(I,J,K)*dr*XC(I)*dtheta
              !write(6,*) "OUTPUT_HL HERE2"
           end if
        ENDDO

     ENDDO
  ENDDO

  RETURN
END SUBROUTINE INT2DFIELD

SUBROUTINE INT3DFIELD(INPUT,OUTPUT_HL,OUTPUT_HR,NX,NY,NZ,XU,XC,ZWG,kstart,kend)
  ! OUTPUT_HL is output half-left
  ! OUTPUT_HL is output half-right
  implicit none
  INTEGER ( kind = 4) :: NX,NY,NZ,nymod,kstart,kend
  INTEGER ( kind = 4) :: I,J,K
  REAL    ( kind = 8) :: XU(NX),XC(NX),ZWG(NZ)
  REAL    ( kind = 8) :: INPUT(NX,NY,NZ),OUTPUT_HL,OUTPUT_HR
  REAL    ( kind = 8) :: dtheta,dr,dz

  nymod=ny-2
  dtheta = 2.0*3.1415926/dble(ny-2)
  OUTPUT_HL = 0.0
  OUTPUT_HR = 0.0

  !DO K=2,NZ-1
  do k = kstart,kend
     write(6,*) " INT3DFIELD: K = ", K
     DO J=2,NY-1
     DO I=2,NX-1

        dz=zwg(k)-zwg(k-1)
        dr=xu(i)-xu(i-1)

           if(j >= 2 .and. j <= (nymod/4 + 1)) then
              OUTPUT_HL=OUTPUT_HL+INPUT(I,J,K)*dr*XC(I)*dtheta*dz
              !write(6,*) "OUTPUT_HL HERE1"
           elseif(j >= (nymod/4 + 2) .and. j <= (nymod/2)+(nymod/4)+1) then
              OUTPUT_HR=OUTPUT_HR+INPUT(I,J,K)*dr*XC(I)*dtheta*dz
              !write(6,*) "OUTPUT_HR"
           elseif(j >= (nymod/2)+(nymod/4)+2) then
              OUTPUT_HL=OUTPUT_HL+INPUT(I,J,K)*dr*XC(I)*dtheta*dz
              !write(6,*) "OUTPUT_HL HERE2"
           end if
        ENDDO

     ENDDO
  ENDDO

  RETURN
END SUBROUTINE INT3DFIELD

subroutine write_int2dout(filename,var,nz,zcg)
  implicit none

  integer ( kind = 4 ) :: k,nz
  character(len=128)   :: filename
  real    ( kind = 8 ) :: var(nz),zcg(nz)

  open(1,file=filename,form='formatted',status='new')

  !do k = 1, nz
  do k = 913, nz-2
     write(1,*) zcg(k), var(k)
  end do
  close(1)
  return
end subroutine write_int2dout

subroutine write_int3dout(filename,var)
  implicit none

  character(len=128)   :: filename
  real    ( kind = 8 ) :: var

  open(1,file=filename,form='formatted',status='new')

  write(1,*) var

  close(1)
  return
end subroutine write_int3dout

SUBROUTINE SET_PERIODIC_J(VAR,NX,NY,NZ)
  implicit none
  INTEGER :: NX,NY,NZ
  REAL    :: VAR(NX,NY,NZ)

  VAR(:,1,:) = VAR(:,NY-1,:)
  VAR(:,NY,:) = VAR(:,2,:)


  RETURN
END SUBROUTINE SET_PERIODIC_J

SUBROUTINE SET_INTERP_I(VAR,NX,NY,NZ)
  implicit none
  INTEGER :: NX,NY,NZ
  REAL    :: VAR(NX,NY,NZ)

  VAR(1,:,:) = 2*VAR(2,:,:)-VAR(3,:,:)

  RETURN
END SUBROUTINE SET_INTERP_I


subroutine write_plt_2d_3vars(filename,u,v,w,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:),ALLOCATABLE :: xc_2d,yc_2d,zc_2d,u_car,v_car,w_car

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx
  jmax_plt = 1
  kmax_plt = kend-kstart+1

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_2d(1:imax_plt,1:kmax_plt))
  ALLOCATE(yc_2d(1:imax_plt,1:kmax_plt))
  ALLOCATE(zc_2d(1:imax_plt,1:kmax_plt))
  ALLOCATE(u_car(1:imax_plt,1:kmax_plt))           
  ALLOCATE(v_car(1:imax_plt,1:kmax_plt))           
  ALLOCATE(w_car(1:imax_plt,1:kmax_plt))           
  
  ! VERTICAL PLANE
  j = (NY-2)/4

  kmod = 0
  DO k = kstart, kend
     kmod = kmod + 1

     write(6,*) "grid setup + uv cal: k = ", k
     DO i = 1, imax_plt

        xc_2d(i,kmod) = xc(i)
        zc_2d(i,kmod) = zcg(k)

        u_car(i,kmod) = u(i,j,k)
        w_car(i,kmod) = w(i,j,k)

     ENDDO
  ENDDO
  write(6,*) "done1"
  ier = tecini(trim(title)//nulchar,'x,y,u,v'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,1,kmax_plt, &
              'BLOCK'//nulchar,nulchar)
  write(6,*) "done2"

  ! Write out the field data.
  itot = imax_plt*kmax_plt
  ier = tecdat(itot,zc_2d,disdouble)
  ier = tecdat(itot,xc_2d,disdouble)
  ier = tecdat(itot,w_car,disdouble)
  ier = tecdat(itot,u_car,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_2d_3vars"
  return 
end subroutine write_plt_2d_3vars

subroutine write_plt_3d_3vars(filename,u,v,w,nx,ny,nz,xc,yc,zcg,kstart,kend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend

  integer ( kind = 4 ) :: kmod
  character(len=128)   :: buff,filename
  real    ( kind = 8 ) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
  
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  REAL ( kind = 8 ), DIMENSION(:,:,:),ALLOCATABLE :: xc_3d,yc_3d,zc_3d,u_car,v_car,w_car

  CHARACTER(len=33) :: title
  INTEGER  imax_plt,jmax_plt,kmax_plt,debug,ier,itot
  INTEGER  tecini,tecdat,teczne,tecend
  INTEGER  visdouble,disdouble,K_PLT
  CHARACTER*1 nulchar      

  write(title,'(a,i1)') "UVW"
  nulchar = char(0)
  debug   = 0
  visdouble = 0
  disdouble = 1

  imax_plt = nx-1
  jmax_plt = ny-1
  kmax_plt = kend-kstart+1

  write(6,*) "imax_plt, jmax_plt, kmax_plt = ", imax_plt, jmax_plt, kmax_plt

  ALLOCATE(xc_3d(1:imax_plt,1:jmax_plt,1:kmax_plt))
  ALLOCATE(yc_3d(1:imax_plt,1:jmax_plt,1:kmax_plt))
  ALLOCATE(zc_3d(1:imax_plt,1:jmax_plt,1:kmax_plt))
  ALLOCATE(u_car(1:imax_plt,1:jmax_plt,1:kmax_plt))           
  ALLOCATE(v_car(1:imax_plt,1:jmax_plt,1:kmax_plt))           
  ALLOCATE(w_car(1:imax_plt,1:jmax_plt,1:kmax_plt))           
  
  kmod = 0
  DO k = kstart, kend
     kmod = kmod + 1

     write(6,*) "grid setup + uv cal: k = ", k
     DO i = 1, imax_plt
        DO j = 1, jmax_plt

           xc_3d(i,j,kmod) = xc(i)*cos(yc(j))
           yc_3d(i,j,kmod) = xc(i)*sin(yc(j))
           zc_3d(i,j,kmod) = zcg(k)

           u_car(i,j,kmod) = u(i,j,k)*sin(yc(j)) + v(i,j,k)*cos(yc(j)) 
           v_car(i,j,kmod) = u(i,j,k)*cos(yc(j)) - v(i,j,k)*sin(yc(j)) 
           w_car(i,j,kmod) = w(i,j,k)

        ENDDO
     ENDDO
  ENDDO
  write(6,*) "done1"
  ier = tecini(trim(title)//nulchar,'x,y,z,u,v,w'//nulchar, &
               trim(filename)//nulchar,'.'//nulchar,debug,visdouble)

  ier = teczne(trim(title)//nulchar,imax_plt,jmax_plt,kmax_plt, &
              'BLOCK'//nulchar,nulchar)
  write(6,*) "done2"

  ! Write out the field data.
  itot = imax_plt*jmax_plt*kmax_plt
  ier = tecdat(itot,xc_3d,disdouble)
  ier = tecdat(itot,yc_3d,disdouble)
  ier = tecdat(itot,zc_3d,disdouble)
  ier = tecdat(itot,u_car,disdouble)
  ier = tecdat(itot,v_car,disdouble)
  ier = tecdat(itot,w_car,disdouble)
  ! Close the file
  ier = tecend()

  write(6,*) "done write_plt_3d_3vars"
  return 
end subroutine write_plt_3d_3vars
