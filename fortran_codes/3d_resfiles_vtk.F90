Program Main
  implicit none


  !number of grid points with ghost cells
  integer ( kind = 4 ), parameter :: nx = 368
  integer ( kind = 4 ), parameter :: ny = 514
  integer ( kind = 4 ), parameter :: nz = 514

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imax,jmax,kmax,tmp,kstart,kend,s1,jend
  integer ( kind = 4 ) :: tag
  real    ( kind = 8 ) :: time0,time1,time2,DTM1,grav,Re
  character(len=128)   :: buff,filename_in,filename_out,basename

  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u0,v0,w0
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u1,v1,w1,p1,dens1,baro_ver
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u2,v2,w2,p2,desn2
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: uc,vc,wc
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: uo,vo,wo
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: diffu_st,diffu_sw,diffu_vr
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: stretch_ver

  real    ( kind = 8 ) :: int2d_hl_1(nz),int2d_hr_1(nz)
  real    ( kind = 8 ) :: int2d_hl_2(nz),int2d_hr_2(nz)
  real    ( kind = 8 ) :: int2d_hl_3(nz),int2d_hr_3(nz)
  real    ( kind = 8 ) :: int3d_hl_1,int3d_hr_1
  real    ( kind = 8 ) :: int3d_hl_2,int3d_hr_2
  real    ( kind = 8 ) :: int3d_hl_3,int3d_hr_3

  real    ( kind = 8 ) :: r(nx),theta(ny),z(nz),dtheta


  ! In order to reduce the file size, the domain can be sliced in the streamwise and the radial direction

   !!! kstart and kend are the axial start and end points

   basename = 'w_00024000_in'
   
   filename_in = trim(basename) // '.res'
   filename_out= trim(basename) // '.vtk'
   
   kstart = 2
   kend = nz-1

   jend = ny-1


   call read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nz,ru,rp,tag)


  allocate(w1(nx,ny,nz),u1(nx,ny,nz),v1(nx,ny,nz))!,p1(nx,ny,nz),dens1(nx,ny,nz))
  allocate(wc(2:nx-1,2:ny-1,2:nz-1),uc(2:nx-1,2:ny-1,2:nz-1),vc(2:nx-1,2:ny-1,2:nz-1))

  call read_restart(filename_in,nx,ny,nz,w1,time1)

  ! Centering the face variables

  do k=2,nz-1
    do j=2,ny-1
        do i=2,nx-1
            !wc(i,j,k) = 0.5d0*(w1(i,j,k)+w1(i,j,k-1))
            wc(i,j,k) = w1(i,j,k)
        enddo
    enddo
  enddo

  call write_vtk(filename_out,wc,nx,ny,nz,xc,yc,zcg,kstart,kend,jend)

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

  OPEN(UNIT=1,FILE='x1_grid.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nx-1
     read(1,*) j, xu(i)
  end do
  close(1)
 
  xc(2:nx-1) = .5*(xu(1:nx-2)+xu(2:nx-1))
  xc(1 ) = 2.*xu(1  )-xc(2  )
  xc(nx) = 2.*xu(nx-1)-xc(nx-1)

  do i= 2, nx-1
     write(6,*) "xc(",i,") = ", xc(i)
  enddo

  dtheta = 2.0*3.1415926/dble(ny-2)
  do i= 1, ny
     yv(i)=real(i-1)*dtheta
     !write(6,*) "yv(",i,") = ", yv(i)
  enddo
         
  do j=2,ny-1
     yc(j)=0.5*(yv(j)+yv(j-1))
     write(6,*) "yc(",j,") = ", yc(j)
  enddo
  yc(1)=2.*yv(1)-yc(2)
  yc(ny)=2.*yv(ny)-yc(ny)

  OPEN(UNIT=1,FILE='x3_grid.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
  enddo
  close(1)
  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  )= 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

  do i= 2, nz-1
     write(6,*) "zcg(",i,") = ", zcg(i)
  enddo

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

subroutine write_vtk(filename,var,nx,ny,nz,xc,yc,zcg,kstart,kend,jend)
  implicit none

  integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend,jend
  character (len = 128)   :: buff,filename
  real    ( kind = 8 ) :: var(2:nx-1,2:ny-1,2:nz-1)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zcg(nz)

  !erase previous output file
  open(unit=1, iostat=s1, file=filename, status='old')
  if (s1 == 0) close(1, status='delete')


  open(unit=1,file=filename,access='stream',form='unformatted',status='new',&
         convert='big_endian',iostat=s1)

  imax = nx-2
  jmax = ny-2
  kmax = kend-kstart+1

  write(1) "# vtk DataFile Version 3.0"//char(10)
  write(1) "FlowField"//char(10)
  write(1) "BINARY"//char(10)
  write(1) "DATASET STRUCTURED_GRID"//char(10)
  write(buff,FMT='(A10,3I5)') "DIMENSIONS",imax,(jmax+1),kmax
  write(1) buff//char(10)
  write(buff,FMT='(A6,I15,A6)') "POINTS",imax*(jmax+1)*kmax, " float"
  write(1) buff//char(10)

  do k = kstart, kend
    !write(6,*) "GRID: WRITE K = ", kstart, k, kend
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
  do k = kstart,kend
    !write(6,*) "DATA: WRITE K = ", kstart, k, kend
     do j = 2, ny
        do i = 2, nx-1
           write(1) real(var(i,j,k))
        end do
     end do
  end do

  close(1)

end subroutine write_vtk
