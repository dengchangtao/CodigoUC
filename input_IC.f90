!Rutina que lee la batimetria y la guarda en x_global, y_global, z_global

SUBROUTINE input_IC

USE global_variables
USE geometries
USE time0
USE senales

implicit none

integer::i,j,k

!allocate(qold_global(3,Nbx,Nby))

allocate (qnew_global(3,Nbx,Nby), qold_global(3,Nbx,Nby), &
	  qreal_global(3,Nbx,Nby),q0_global(3,Nbx,Nby), V_global(Nbx,Nby), C_global(Nbx,Nby), & VC(Nbx,Nby),S1_global(Nbx,Nby),S2_global(Nbx,Nby))
allocate(qA1(3,Nby),qA2(3,Nby),qA3(3,Nbx),qA4(3,Nbx),zA1(Nby),zA2(Nby),zA3(Nbx),zA4(Nbx))
	  

print*, Nbx, Nby
 open	(unit=2, file ='IC.dat')!, form='unformatted')
! read	(unit=2) ((qold_global1(1,i,j),i=1,Nbx),j=1,Nby), ((qold_global1(2,i,j),i=1,Nbx),j=1,Nby), ((qold_global1(3,i,j),i=1,Nbx),j=1,Nby)
! close(unit=2)
! print*, qold_global1

DO i=1,Nbx
  DO j=1,Nby
  read(2,*) qold_global(1,i,j), qold_global(2,i,j), qold_global(3,i,j)  
  !print*, qold_global(1,i,j), qold_global(2,i,j), qold_global(3,i,j) 
  END DO
END DO

q0_global=qold_global

END SUBROUTINE input_IC