 subroutine write_vtk(filename,var,nx,ny,nz,xc,yc,zcg,kstart,kend)
   implicit none
   integer ( kind = 4 ) :: i,j,k,nx,ny,nz,imax,jmax,kmax,s1,kstart,kend
   character(len=128)   :: buff,filename
   real    ( kind = 8 ) :: var(nx,ny,nz)
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
            write(1) real(xc(i)*cos(yc(j))), real(xc(i)*sin(yc(j))), real(zcg(k))
         end do
      end do
   end do
 
   write(buff,fmt='(A10,I15)') "POINT_DATA",imax*jmax*kmax
   write(1) char(10)//buff//char(10)
 
   write(1) "SCALARS vtkvar float 1"//char(10)
   write(1) "LOOKUP_TABLE default"//char(10)
   do k = kstart, kend
      write(6,*) "DATA: WRITE K = ", kstart, k, kend
      do j = 1, jmax
         do i = 1, imax
            write(1) real(var(i,j,k))
         end do
      end do
   end do
   close(1)
 
 end subroutine write_vtk

