subroutine readgrid(grid_dir,xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg,ru,rp,tag)
  implicit none

  INTEGER nx,ny,nz,nzg,tag
  REAL ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nzg)
  REAL ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nzg)
  REAL ( kind = 8 ) :: ru(nx),rp(nx)
  CHARACTER (len = 160) :: grid_dir, filename
  INTEGER i,j,k
  real rdelx,rdely,rdelz, dtheta
  real, allocatable, dimension(:) :: cug,cvg

  ALLOCATE(cug(nzg),cvg(nzg))

  write(filename,'(a,a)') trim(grid_dir), 'x1_grid.in'
  OPEN(UNIT=1,FILE=filename,STATUS='OLD',FORM='FORMATTED')

  read(1,*) j
  do i= 1, nx-1
     read(1,*) j, xu(i)
     write(6,*) "xu(",i,") = ", xu(i)
  enddo
  close(1)
  xc(2:nx-1) = .5*(xu(1:nx-2)+xu(2:nx-1))
  xc(1 ) = 2.*xu(1  )-xc(2  )
  xc(nx) = 2.*xu(nx-1)-xc(nx-1)

  dtheta = 2.0*3.1415926/dble(ny-2)
  do i= 1, ny
     yv(i)=real(i-1)*dtheta
  enddo

  do j=2,ny-1
     yc(j)=0.5*(yv(j)+yv(j-1))
     write(6,*) "yc(",j,") = ", yc(j)
  enddo
  yc(1)=2.*yv(1)-yc(2)
  yc(ny)=2.*yv(ny)-yc(ny)

  close(1)

  write(filename,'(a,a)') trim(grid_dir), 'x3_grid.in'
  OPEN(UNIT=1,FILE=filename,STATUS='OLD',FORM='FORMATTED')

  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
     write(6,*) "zwg(",i,") = ", zwg(i)
  enddo
  close(1)

  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  ) = 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

  zw = zwg
  zc = zcg
  
  ru(1:nx)=xu(1:nx)
  rp(1:nx)=xc(1:nx)

  close(1)

  write(6,*) "READ GRID DONE"

  return

end subroutine readgrid
