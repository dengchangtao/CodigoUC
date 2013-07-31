SUBROUTINE outputmat(time)

USE global_variables
USE geometries

implicit none

integer	:: i,j, ent, itent
real (kind=8)	:: ind, dec, time, dit1
logical	:: lexist, lexistT
character(len=1100):: filename, filenameT
character(len=10)::number, numbercaso
character(len=1000)::intchar, ncaso

!-------------------------
logical	:: lexistP, lexistT_eta
character(len=1100):: filenameP,filenameT_eta, path

! Blend Interpolation for times series at a point
real (kind=8), dimension(3)::x0, y0, i0, j0, H0
real (kind=8), dimension(1)::m1,m2
real (kind=8), dimension(3)::r, s, x, x00, x01, x10, x11 ! Interpolation by blending (see blend.f90, function blend_102)
! Coordenadas Puntos

x0=(/334513.54, 381251.0, 402755.0 /)
y0=(/292620.1, 282966.0, 292633.0 /)


do i=1,3
		!m=x0(i)
		m1=minloc(x_global(:,1)-x0(i),1,mask=(x_global(:,1)-x0(i)).ge.0)
		i0(i)=m1(1)-1
		r(i)=(x0(i)-x_global(i0(i),1))/(x_global(2,1)-x_global(1,1))
		m2=minloc(y_global(1,:)-y0(i),1,mask=(y_global(1,:)-y0(i)).ge.0)
		j0(i)=m2(1)-1
		s(i)=(y0(i)-y_global(1,j0(i)))/(y_global(1,2)-y_global(1,1))
		x00(i)=qreal_global(1,i0(i),j0(i))    +z_global(i0(i),j0(i))
		x10(i)=qreal_global(1,i0(i)+1,j0(i))  +z_global(i0(i)+1,j0(i))
		x01(i)=qreal_global(1,i0(i),j0(i)+1)  +z_global(i0(i),j0(i)+1)
		x11(i)=qreal_global(1,i0(i)+1,j0(i)+1)+z_global(i0(i)+1,j0(i)+1)
end do




!------------------------
if (time==0.0D0) then
qnew_global=qold_global
end if

!Dimensionalization of the variables
do i=1,Nbx; do j=1,Nby
		qreal_global(1,i,j)=qnew_global(1,i,j)*H
		qreal_global(2,i,j)=qnew_global(2,i,j)*U
		qreal_global(3,i,j)=qnew_global(3,i,j)*U
end do; end do

do i=1,3
		!H0(i)=qreal_global(1,i0(i),j0(i))+z_global(i0(i),j0(i))
		! function blend_102 by John Burkardt (see blend.f90)
		H0(i)=            	      + x00(i) &
				+ r(i) *        ( - x00(i) + x10(i) ) & 
				+ s(i) *        ( - x00(i)       + x01(i) ) &
				+ r(i) * s(i) * ( + x00(i) - x10(i) - x01(i) + x11(i) )
			
end do

!Writing files
!Write the initial conditions in the first results file
!Escribir tb las condiciones iniciales en un archivo
ind=it/dit
ent=ind
dec=ind-ent


IF (dec==0.0D0.OR.it==0.0D0) THEN
  
  itent=it
  
  write(intchar,*) itent
  intchar=adjustl(intchar)
  number=intchar
  write(ncaso,*) caso
  ncaso=adjustl(ncaso)
  numbercaso=ncaso

  path='Results/run1/'!/media/JGM_/Resultados simulacion/Malla_fija_5mm/' !data2/maru/EQ89/' !Cambiar Path segun caso
  !filename='SOL2D.'//trim(number)//'.dat'
  !filename='/Resultados/Case'//trim(numbercaso)//'/SOL2D.'//trim(number)//'.dat!Invunche:
  
  
  filename=trim(path)//'SOL2D.'//trim(number)//'.dat'  
  
  !filenameT='Time'//trim(numbercaso)//'.dat'
  !filename='\Resultados\Case'//trim(numbercaso)//'\Time'//trim(numbercaso)//'.dat'
  !Invunche:
  filenameT=trim(path)//'Time'//trim(numbercaso)//'.dat'
  
  inquire(FILE=filename, EXIST=lexist)
  
  inquire(FILE=filenameT, EXIST=lexistT)


 

  IF (.NOT. lexist) THEN
  open(10,file=filename,status='new',action='write')
  ELSE
  open(10,file=filename,status='replace',action='write')
  END IF
  
  IF (.NOT. lexistT) THEN
  open(20,file=filenameT,status='new',action='write')
  ELSE
    
    IF (it==0.0D0) THEN
    open(20,file=filenameT,status='replace',action='write')
    ELSE
    open(20,file=filenameT,status='old',action='write', position='append')
    END IF
    
 END IF

!Escribo

!Tiempo
  write(20,180) treal
  180 format(F10.3) !Formato para el vector tiempo
!Resultados

  DO j=1,Nby; DO i=1,Nbx
    
    if (abs(qreal_global(1,i,j))<=kappa) then
    qreal_global(1,i,j)=0.0D0 
    end if
    if (abs(qreal_global(2,i,j))<=kappa) then
    qreal_global(2,i,j)=0.0D0 
    end if
    if (abs(qreal_global(3,i,j))<=kappa) then
    qreal_global(3,i,j)=0.0D0 
    end if
    
    write(10,170) x_global(i,j),y_global(i,j),z_global(i,j),qreal_global(1,i,j), qreal_global(2,i,j),qreal_global(3,i,j)
    170 format(F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5, TR1, F15.5 / )
    
  END DO; END DO;

  close(10)
  close(20)
  
END IF

! -----------------------------------------------------------
! Escribir a cada paso de tiempo (o cada it1 pasos de tiempo)

dit1=10
!~ ind1=it/dit1
!~ ent1=ind1
!~ dec1=ind1-ent1
!~ print*,'NOT SAVING...', mod(it,dit1)
!~ IF (dec1==0.0D0.OR.it==0.0D0) THEN
IF (mod(it,dit1)==0.0D0.OR.it==0.0D0) THEN

!~ print*,'SAVING...', mod(it,dit1), it, dit1
filenameP=trim(path)//'Eta.dat'
filenameT_eta=trim(path)//'Time_eta.dat'


  inquire(FILE=filenameP, EXIST=lexistP)
  inquire(FILE=filenameT_eta, EXIST=lexistT_eta)

  
  IF (.NOT. lexistP) THEN  
	open(30,file=filenameP,status='new',action='write')  
  ELSE
  
  IF (it==0.0D0) THEN
		open(30,file=filenameP,status='replace',action='write')
    ELSE
		open(30,file=filenameP,status='old',action='write', position='append')
    END IF

  END IF
  
  IF (.NOT. lexistT_eta) THEN
	open(40,file=filenameT_eta,status='new',action='write')
  ELSE
    
    IF (it==0.0D0) THEN
		open(40,file=filenameT_eta,status='replace',action='write')
    ELSE
		open(40,file=filenameT_eta,status='old',action='write', position='append')
    END IF
    
  END IF
  
  
  
  !Tiempo
  write(40,190) treal
  190 format(F10.3) !Formato para el vector tiempo
  

  write(30,200) H0
    200 format(4(F15.5, TR1) / )


close(30)
close(40)
close(50)

END IF


END SUBROUTINE outputmat



