 subroutine read_restart(filename,nx,ny,nz,var,i,j,k,jp,nstep,time,DTM1,grav)
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
