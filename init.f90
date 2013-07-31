!Rutina que inicia las variables,lee datos de entrada, lee batimetria, condiciones iniciales, condiciones de borde y asigna estos datos a variables
!Llama a tranformacion de coordenadas
SUBROUTINE init

USE global_variables
USE geometries
USE senales
USE coords

implicit none

integer :: i,j 
real (kind=8) :: U1, U2

  
!Leer informacion de la simulacion
call input_control

FR2=U**2.0D0/(g*H)
print*, 'Fr2= ', FR2

! Leer batimetria o Crear
SELECT CASE (batiopt)
  CASE(1)
  !Lee Batimetria
  call input_geom
  
  CASE(2)
  !Crea Batimetria usando bati.f90
  call batimetria
  
END SELECT


SELECT CASE (inputopt)
  CASE(1)!
  !Lee condicion inicial
  call input_IC
  print*,'input_IC_OK'
  
  CASE(2)
  !Condiciones Inciales
  call init_flowfield
  print*,'InitFlowfieldOK'
  
END SELECT


!Coeficiente de friccion: debe venir adimensionalizado
allocate(MCoef(Nbx,Nby))

IF (fopt==1) THEN
  
  IF (fM==1) THEN !Completa matriz de friccion con el valor ingresado en input.dat
    DO i=1,Nbx 
    DO j=1,Nby
    MCoef(i,j)=Coef
    END DO
    END DO
  print*, 'friccionOK'
  ELSE !Lee matriz de friccion de un archivo binario
  open	(unit=99, file ='friccion.dat', form='unformatted')
  read	(unit=99) ((MCoef(i,j),i=1,Nbx),j=1,Nby)
  close(unit=99)
  print*, 'friccionOK'
  END IF

ELSE
  Cf=0
  DO i=1,Nbx 
  DO j=1,Nby
  MCoef(i,j)=0.0D0
  END DO
  END DO 
END IF
      

!Adimensionalizar CI y Bathy
call adimension

!Calcular las metricas
call metrics
!Crea Coordenadas Xi,Eta

allocate(coordxi(Nbx),coordeta(Nby))
 coordxi(1)=dxi/2.0D0
 coordeta(1)=deta/2.0D0
do i=2,Nbx
 coordxi(i)=coordxi(i-1)+dxi
end do
do i=2,Nby
 coordeta(i)=coordeta(i-1)+deta
end do

!Angulo normales a los bordes con respecto a un eje

call angulo

!Velocidades Contravariantes para CFL
do i=1,Nbx; do j=1,Nby

	U1=qold_global(2,i,j)*xi_global(1,i,j)+qold_global(3,i,j)*xi_global(2,i,j)
	U2=qold_global(2,i,j)*eta_global(1,i,j)+qold_global(3,i,j)*eta_global(2,i,j)
	S1_global(i,j)=U1+C_global(i,j)*sqrt(xi_global(1,i,j)**2+xi_global(2,i,j)**2)
	S2_global(i,j)=U2+C_global(i,j)*sqrt(eta_global(1,i,j)**2+eta_global(2,i,j)**2)

end do; end do


print*, 'Simluacion Incializada'


END SUBROUTINE init


SUBROUTINE init_flowfield
!Condiciones Iniciales: Diferentes para cada caso
!Se podrian bajar de archivos, pero por mientras mejor hacerlas aqui

USE global_variables
USE geometries
USE senales
USE time0
implicit none

integer :: mitad, i, j, error
real (kind=8)::	zaux, hz, xmed, xo, yo, hc, hs,r,d, Haux, sigma, hmean, hestanque, &
		p3, c3, yaux, Am, a, ro, ho, eta0, omega,tau,p,S, B, uo, Dsyn, gama, x1, m, z85, u0

allocate (qnew_global(3,Nbx,Nby), qold_global(3,Nbx,Nby), &
	  qreal_global(3,Nbx,Nby),q0_global(3,Nbx,Nby), V_global(Nbx,Nby), C_global(Nbx,Nby), &
	  VC(Nbx,Nby),S1_global(Nbx,Nby),S2_global(Nbx,Nby), STAT = error)
!GA
allocate(qA1(3,Nby),qA2(3,Nby),qA3(3,Nbx),qA4(3,Nbx),zA1(Nby),zA2(Nby),zA3(Nbx),zA4(Nbx))


SELECT CASE (caso)

	CASE(1)
	print *, 'Steady State, Flat Bottom, Uniform Grid'
	!Input dimensionalized data
	!Ingresar Datos Dimensionales, mas adelante se adimensionalizan!
	do i=1,Nbx; do j=1,Nby
			
			qold_global(1,i,j)=0.5D0
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	end do; end do
	
!--------------------------------------------------------------------------------------------------------------

	CASE(2)

	print *, 'Steady State, Over a Bump, Uniform Grid '
	
	hz=0.15D0
	do i=1,Nbx; do j=1,Nby
			  zaux=z_global(i,j)
 
			  if (zaux<hz) then
			    qold_global(1,i,j)=hz-zaux
			  else
			    qold_global(1,i,j)=0.0D0
			  end if
			  qold_global(2,i,j)=0.0D0
			  qold_global(3,i,j)=0.0D0
	end do; end do
	
!--------------------------------------------------------------------------------------------------------------	
	
	CASE(3)

	print *, 'Steady State, Over a Bump, Uniform Grid, Water over the bump'
	
	do i=1,Nbx; do j=1,Nby
			zaux=z_global(i,j)
			qold_global(1,i,j)=max(0.25D0-zaux,0.0D0)
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	end do; end do
	
!--------------------------------------------------------------------------------------------------------------

	CASE(4)

	print *, '1D Dam-Break'
	xmed=24.5D0
	do i=1,Nbx; do j=1,Nby
			
			if (x_global(i,j).le.xmed) then
			
			    qold_global(1,i,j)=1.0D0
			else
			    qold_global(1,i,j)=0.0D0
			end if
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	end do; end do
	
!--------------------------------------------------------------------------------------------------------------	

	CASE(5)
	print *, '2D Dam-Break'
	xmed=100.0D0
	
	do i=1,Nbx; do j=1,Nby
			if (x_global(i,j).le.100.0D0) then
			    hz=10.0D0;
			else
			    hz=5.0D0;
			end if
        
			if (z_global(i,j).le.hz) then
			    qold_global(1,i,j)=hz-z_global(i,j)
			else
			    qold_global(1,i,j)=0.0D0
			end if
			
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	end do; end do

!--------------------------------------------------------------------------------------------------------------

	CASE(6)
	print *, '2D Dam-Break with an obstacle'
	xmed=100.0D0
	
	do i=1,Nbx; do j=1,Nby
			if (x_global(i,j).le.50.0D0) then
			    hz=10.0D0;
			else
			    hz=5.0D0;
			end if
        
			if (z_global(i,j).le.hz) then
			    qold_global(1,i,j)=hz-z_global(i,j)
			else
			    qold_global(1,i,j)=0.0D0
			end if
			
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	end do; end do	
	
!--------------------------------------------------------------------------------------------------------------
	
	CASE(7)
	print *, 'Cilindrical Dam-Break, Cartesian'
	
	xo=25.0D0
	yo=25.0D0
	hc=10.0D0
	hs=1.0D0
	r=11.0D0
	do i=1,Nbx; do j=1,Nby
			
			d=sqrt((x_global(i,j)-xo)**2+(y_global(i,j)-yo)**2);
			
			if (d<r) then
        
			qold_global(1,i,j)=hc
			else
			qold_global(1,i,j)=hs
			end if
			
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	end do; end do		

!--------------------------------------------------------------------------------------------------------------

	CASE(8)
	print *, 'Cilindrical Dam-Break, Curvilinear'
	
	r=11.000D0
	!"print*,r
	!pause
	hc=10.0D0
	hs=1.0D0
	do i=1,Nbx
	  do j=1,Nby
			
			d=sqrt((x_global(i,j))**2+(y_global(i,j))**2)
			
			
			if (d<=r) then
			qold_global(1,i,j)=hc
			else
			qold_global(1,i,j)=hs
			end if
			
			!print*, d, qold_global(1,i,j), i, j !, x_global(i,j), y_global(i,j)
			
			qold_global(2,i,j)=0.0D0
			qold_global(3,i,j)=0.0D0
	  end do
	!pause	
	end do	

!--------------------------------------------------------------------------------------------------------------

	CASE(9)
	print *, '2D Flow at rest'
	
	xo=5.0D0
	yo=5.0D0
	sigma=0.5D0
	
	hz=0.3001D0
	
	do i=1,Nbx
	  do j=1,Nby
			
			r=sqrt((x_global(i,j)-xo)**2.0D0+(y_global(i,j)-yo)**2.0D0)
			d=0.5D0*exp(-(r/sigma)**2.0D0)
			
			qold_global(1,i,j)=max(hz-d,0.0D0)
			
! 			  if (zaux<=hz) then
! 			    qold_global(1,i,j)=hz-zaux
! 			  else
! 			    qold_global(1,i,j)=0.0D0
! 			  end if
			  qold_global(2,i,j)=0.0D0
			  qold_global(3,i,j)=0.0D0
			
			
	  end do
	
	end do	
	
!--------------------------------------------------------------------------------------------------------------

	CASE(10)
	print *, 'Vaciamiento Estanque Pascua'
		
	
	hmean=0.56D0
	hestanque=0.75D0
	
	do i=1,Nbx
	  do j=1,Nby
	      
	p3=16.0D0
	c3=1.6D0-16.0D0*1.475D0
	      
	      !z_global(i,j)=z_global(i,j)*0.9D0
	      
             IF( x_global(i,j)>=0.0D0 .AND. x_global(i,j) <= 1.5D0 .AND. y_global(i,j)>=-0.1D0 .AND. y_global(i,j)<= 4.6D0) THEN
                qold_global(1,i,j) = MAX(hestanque - z_global(i,j), 0.0D0)
                qold_global(2,i,j) = 0.0D0
                qold_global(3,i,j) = 0.0D0
             
	     ELSE IF( x_global(i,j)>1.375D0 .AND. x_global(i,j)<1.475D0 .AND. y_global(i,j)>=0.0D0 .AND. y_global(i,j)<1.6D0) THEN
		yaux=p3*x_global(i,j)+c3
		IF( y_global(i,j)>=yaux) THEN
		qold_global(1,i,j) = MAX(hestanque - z_global(i,j), 0.0D0)
                qold_global(2,i,j) = 0.0D0
                qold_global(3,i,j) = 0.0D0

		ELSE
		qold_global(1,i,j) = MAX(hmean - z_global(i,j), 0.0D0)
                qold_global(2,i,j) = 0.0D0
                qold_global(3,i,j) = 0.0D0
		END IF

	     ELSE IF(x_global(i,j)>1.3750D0 .AND. x_global(i,j)<1.4750D0 .AND. y_global(i,j)>=1.60D0 .AND. y_global(i,j)<=4.6D0) THEN

                qold_global(1,i,j) = MAX(hestanque - z_global(i,j), 0.0D0) 
                qold_global(2,i,j) = 0.0D0
                qold_global(3,i,j) = 0.0D0
             ELSE 
                qold_global(1,i,j) = MAX(hmean - z_global(i,j), 0.0D0) 
                qold_global(2,i,j) = 0.0D0
                qold_global(3,i,j) = 0.0D0
             
	     end IF
	     
			
	  end do
	
	end do	
	
!--------------------------------------------------------------------------------------------------------------	

	CASE(11)
	print*, 'Thacker`s Curved Solution'
	
	xo=2.0D0
	yo=2.0D0
	a=1.0D0
	ro=0.80D0
	ho=0.10D0
	Am=(a**2.0D0-ro**2.0D0)/(a**2.0D0+ro**2.0D0)
	
	do i=1,Nbx; do j=1,Nby
	
	r=sqrt((x_global(i,j)-xo)**2.0D0+(y_global(i,j)-yo)**2.0D0)
	hz=ho*((1.0D0-Am**2)**0.5D0/(1.0D0-Am)-1.0D0-(r**2.0D0/a**2.0D0)*((1.0D0-Am**2)/(1.0D0-Am)**2.0D0-1.0D0))
	
	if (hz<z_global(i,j)) then
	
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	end if
	
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end do; end do
	
!--------------------------------------------------------------------------------------------------------------	
	CASE(12)
	print*, 'Thacker`s Planar Solution'
	
	xo=2.0D0
	yo=2.0D0
	a=1.0D0

	ho=0.1D0
	eta0=0.5D0
	omega=sqrt(2.0D0*g*ho)/a
	
	do i=1,Nbx; do j=1,Nby
	
	
	hz=(eta0*ho/a**2.0D0)*(2.0D0*(x_global(i,j)-xo)-eta0)
	
	if (hz<z_global(i,j)) then
	
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	end if
	
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=eta0*omega
	end do; end do	
	
!--------------------------------------------------------------------------------------------------------------	

	CASE(13)
	print*, 'Sampson, parabolic bathymetry'
	
	ho=10.0D0
	a=3000.0D0
	B=5.0D0
	tau=Coef
	p=sqrt(8.0D0*g*ho/a**2.0D0)
	S=sqrt(p**2.0D0-tau**2.0D0)/2.0D0
		
	do i=1,Nbx; do j=1,Nby
		
	hz=ho+(a**2.0D0)*(B**2.0D0)/(8.0D0*g**2.0D0*ho)*(tau**2.0D0/4.0D0-S**2.0D0)-B**2.0D0/(4.0D0*g)-(1.0D0/g)*B*S*(x_global(i,j))
	
	if (hz<=z_global(i,j)) then
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	end if
	
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end do; end do	
! 	print*, 'Zb=', z_global(:,1)
! 	print*, 'h=', qold_global(1,:,1)
! 	pause
!--------------------------------------------------------------------------------------------------------------	
	CASE(14)
	print*,'Onda Solitaria Synolakis'
	
	hz=h01
	
	do i=1,Nbx; do j=1,Nby
	
	if (hz<=z_global(i,j)) then
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	
	end if
	
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end do; end do	
	
	
!--------------------------------------------------------------------------------------------------------------	
	CASE(15)
	print*,'Sudden gate closure'
	
	hz=0.6D0
	
	do i=1,Nbx; do j=1,Nby
	
	if (hz<=z_global(i,j)) then
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	end if
	
	qold_global(2,i,j)=0.5D0
	qold_global(3,i,j)=0.0D0
	end do; end do	
	
	
!--------------------------------------------------------------------------------------------------------------	
	CASE(17)
	print*,'Curvilinear Converging Channel'
	
	hz=0.03D0
	
	do i=1,Nbx; do j=1,Nby
	
	if (hz<=z_global(i,j)) then
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	end if
	
	qold_global(2,i,j)=2.17D0
	qold_global(3,i,j)=0.0D0
	end do; end do	

!--------------------------------------------------------------------------------------------------------------	

	CASE(18)
	print*,'Tsunami Mataquito'
	
	hz=0.0D0+126.6737D0
	
	do i=1,Nbx; do j=1,Nby
	
	if (hz<=z_global(i,j)) then
	qold_global(1,i,j)=0.0D0
	else
	qold_global(1,i,j)=hz-z_global(i,j)
	end if
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end do; end do	

!--------------------------------------------------------------------------------------------------------------	
	CASE(19)
	print*,'Synolakis 2'
	
	Haux=0.019D0
	Dsyn=1.0D0
	gama=sqrt(3.0D0*Haux/(4.0D0*D))
	!x1=sqrt(4.0D0*D/(3.0D0*Haux))*cosh
	x1=18.2476D0
	do i=1,Nbx; do j=1,Nby
	
	ho=Haux/Dsyn*(1.0D0/(cosh(gama*(80.0D0-x_global(i,j)-x1)))**2.0D0)+Dsyn
	uo=sqrt(g/Dsyn)*(ho-D)
	
	if (z_global(i,j)>ho) then
	qold_global(1,i,j)=0.0D0
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	else
	qold_global(1,i,j)=ho-z_global(i,j)
	qold_global(2,i,j)=uo
	qold_global(3,i,j)=0.0D0
	end if
	end do; end do	
	
	!print*, qold_global(1,:,1)+z_global(:,1)
	!pause

!--------------------------------------------------------------------------------------------------------------	

!--------------------------------------------------------------------------------------------------------------	
	CASE(20)
	print*,'Pichilemu 1'
	
	ho=-0.69D0+3.6789D0
	
	do i=1,Nbx; do j=1,Nby
	if (z_global(i,j)>ho) then
	qold_global(1,i,j)=0.0D0
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	else
	qold_global(1,i,j)=ho-z_global(i,j)
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end if
	end do; end do	
	
	!print*, qold_global(1,:,1)+z_global(:,1)
	!pause

!--------------------------------------------------------------------------------------------------------------	
!--------------------------------------------------------------------------------------------------------------	
	CASE(21)
	print*,'Hydraulic Jump in a Diverging Channel'
	
	ho=0.0976D0
	
	open	(unit=100, file ='hinicial.dat', form='unformatted')

	read	(unit=100) ((qold_global(1,i,j),i=1,Nbx),j=1,Nby)
	close(unit=100)

	open	(unit=200, file ='uinicial.dat', form='unformatted')
	read	(unit=200) ((qold_global(2,i,j),i=1,Nbx),j=1,Nby)
	close(unit=200)

	do i=1,Nbx; do j=1,Nby
! 	qold_global(1,i,j)=0.0D0
! 	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end do; end do	
	
	!print*, qold_global(1,:,1)+z_global(:,1)
	!pause

!--------------------------------------------------------------------------------------------------------------	
	CASE(22)
	print*,'SandersBC1'
	
	ho=2.0D0
	
	do i=1,Nbx; do j=1,Nby
	qold_global(1,i,j)=ho
	qold_global(2,i,j)=0.7937D0
	qold_global(3,i,j)=0.0D0
	end do; end do	
	
	!print*, qold_global(1,:,1)+z_global(:,1)
	!pause

!--------------------------------------------------------------------------------------------------------------	
!--------------------------------------------------------------------------------------------------------------	
	CASE(23)
	print*,'Hydraulic Jump over a Bumb, Transcritical with Shock'
	
	ho=0.33D0
	
	do i=1,Nbx; do j=1,Nby
	if (z_global(i,j)>ho) then
	qold_global(1,i,j)=0.0D0
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	else
	qold_global(1,i,j)=ho-z_global(i,j)
	qold_global(2,i,j)=0.0D0
	qold_global(3,i,j)=0.0D0
	end if
	end do; end do	
!--------------------------------------------------------------------------------------------------------------	
!--------------------------------------------------------------------------------------------------------------	
	CASE(24)
	print*,'Dam break in a Converging-Diverging Flume'
	
	m=0.002D0
	z85=z_global(1,1)-m*8.5D0
	
	ho=0.3D0+z85
	
	do i=1,Nbx; do j=1,Nby	
	    if (x_global(i,j)<=8.5D0) then
	
		if (z_global(i,j)>ho) then
		qold_global(1,i,j)=0.0D0
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		else
		qold_global(1,i,j)=ho-z_global(i,j)
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		end if
	    else
		qold_global(1,i,j)=0.00001D0
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
	    end if
	    
	end do; end do	
!--------------------------------------------------------------------------------------------------------------	
!--------------------------------------------------------------------------------------------------------------	
	CASE(25)
	print*,'Dam break in a closed channel with topo and friction'
	
	
	ho=1.875D0
	
	do i=1,Nbx; do j=1,Nby	
	    if (x_global(i,j)<=16.0D0) then
	
		if (z_global(i,j)>ho) then
		qold_global(1,i,j)=0.0D0
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		else
		qold_global(1,i,j)=ho-z_global(i,j)
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		end if
	    else
		qold_global(1,i,j)=1.0e-6
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
	    end if
	    
	end do; end do	
!--------------------------------------------------------------------------------------------------------------	
	CASE(99) ! Ensayo
	print *, 'Ensayo perfil de playa'
		!Input dimensionalized data
		!Ingresar Datos Dimensionales, mas adelante se adimensionalizan!
		ho=0.765D0
	do i=1,Nbx; do j=1,Nby
			
		qold_global(1,i,j)=ho-z_global(i,j)
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
	end do; end do
		
!	do j=1,Nby
!		qold_global(1,100,j)=ho+0.5-z_global(100,j)
!		qold_global(2,1,j)=3.0D0
!	end do

!--------------------------------------------------------------------------------------------------------------
	CASE(100) ! Tsunami 27F, San Antonio
	print *, 'Tsunami 27F, San Antonio'
	!Ingresar maxima marea con respecto a WGS84
!  	ho=0.0D0 !nivel 27F
! 	!ho=1.4D0 !nivel 8.6 y 8.8
! 	
! 	ho=ho+207.9032D0!205.83248D0!255.6058D0 !traslado igual como trasladé la batimetría
! 	
! 	do i=1,Nbx; do j=1,Nby
! 		
! 		if (z_global(i,j)>ho) then
! 		qold_global(1,i,j)=0.0D0
! 		qold_global(2,i,j)=0.0D0
! 		qold_global(3,i,j)=0.0D0
! 		else
! 		qold_global(1,i,j)=ho-z_global(i,j)
! 		qold_global(2,i,j)=0.0D0
! 		qold_global(3,i,j)=0.0D0
! 		end if
! 	end do; end do
	  open	(unit=2, file ='h0.dat', form='unformatted')
	  read	(unit=2) ((qold_global(1,i,j),i=1,Nbx),j=1,Nby)
	  close(unit=2)

	  open	(unit=3, file ='U0.dat', form='unformatted')
	  read	(unit=3) ((qold_global(2,i,j),i=1,Nbx),j=1,Nby)
	  close(unit=3)

	  open	(unit=4,file='V0.dat', form='unformatted')
	  read	(unit=4) ((qold_global(3,i,j),i=1,Nbx),j=1,Nby)
	  close(unit=4)
	  
	  print*, qold_global(1,50,40)
	  print*, qold_global(2,50,40)
	  print*, qold_global(3,50,40)

!--------------------------------------------------------------------------------------------------------------
	CASE(200) ! Tsunami in a Large River
	print *, 'Tsunami Columbia River'
	!Ingresar maxima marea con respecto al 0 
  	ho=1.5D0 
! 
! 	ho=ho+39.309606D0!39.837508D0!39.899902D0!!!!SanAntonio207.9032D0!205.83248D0!255.6058D0 !traslado igual como trasladé la batimetría
! 	
	do i=1,Nbx; do j=1,Nby
		
		if (z_global(i,j)>ho) then
		qold_global(1,i,j)=0.0D0
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		else
		qold_global(1,i,j)=ho-z_global(i,j)
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		end if
	end do; end do
	
	! Leo Condicion Inicial
	
	  open	(unit=2, file ='h0.dat', form='unformatted')
	  read	(unit=2) ((qold_global(1,i,j),i=1,Nbx),j=1,Nby)
	  close(unit=2)

! 	  open	(unit=3, file ='U0.dat', form='unformatted')
! 	  read	(unit=3) ((qold_global(2,i,j),i=1,Nbx),j=1,Nby)
! 	  close(unit=3)
! 
! 	  open	(unit=4,file='V0.dat', form='unformatted')
! 	  read	(unit=4) ((qold_global(3,i,j),i=1,Nbx),j=1,Nby)
! 	  close(unit=4)

!--------------------------------------------------------------------------------------------------------------
	CASE(300) ! Tsunami in a Large River
	print *, 'Tsunami El Quisco'
	
  	ho=1.4D0  !Marea maxima
! 
 	ho=ho+103.416330D0 !traslado igual como trasladé la batimetría, sumar abs(zmin)
! 	
	do i=1,Nbx; do j=1,Nby
		
		if (z_global(i,j)>ho) then
		qold_global(1,i,j)=0.0D0
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		else
		qold_global(1,i,j)=ho-z_global(i,j)
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		end if
	end do; end do
! 	
! 	! Leo Condicion Inicial
! 	
! 	  open	(unit=2, file ='h0.dat', form='unformatted')
! 	  read	(unit=2) ((qold_global(1,i,j),i=1,Nbx),j=1,Nby)
! 	  close(unit=2)
! 
! ! 	  open	(unit=3, file ='U0.dat', form='unformatted')
! ! 	  read	(unit=3) ((qold_global(2,i,j),i=1,Nbx),j=1,Nby)
! ! 	  close(unit=3)
! ! 
! ! 	  open	(unit=4,file='V0.dat', form='unformatted')
! ! 	  read	(unit=4) ((qold_global(3,i,j),i=1,Nbx),j=1,Nby)
! ! 	  close(unit=4)
!--------------------------------------------------------------------------------------------------------------
	CASE(400) ! Dead Zone con Paredes Laterales
	print *, 'Canal con Dead Zonex, 15x15cm'
	ho=0.01D0/1.0D0  !Profundidad Inicial: 1cm, dividido por el factor de escala L
	!u0=25.50D0/ho/10000.0D0 !Dividido por 100*100 para que M sea en m2/s y nocm2/s y u0 quede en m/s
	u0=(255.0D0/(100*100*100))*(1.0D0/0.1D0)*(1/0.01D0)!Q*1/b*1/h=Q/(b*h) div
	u0=u0/1.0D0!Factor de escala U
	!Q en cm3/s, dividido 100*100*100 para que quede en m3/s
	!b y h en m
	!ho=0.0202D0
	!u0=0.3698D0
! 
 	do i=1,Nbx; do j=1,Nby		
	
		if (z_global(i,j)>0.002D0) then
		
		qold_global(1,i,j)=0.0D0!h
		qold_global(2,i,j)=0.0D0!u
		qold_global(3,i,j)=0.0D0!v
		else		
		qold_global(1,i,j)=ho	
		
		if (y_global(i,j)<0.1D0) then		
		qold_global(2,i,j)=u0
		end if
		
		end if
	end do; end do
	
 	q0_global=qold_global
	
! 	! Leo Condicion Inicial
! 	
! 	  open	(unit=2, file ='h0.dat', form='unformatted')
! 	  read	(unit=2) ((qold_global(1,i,j),i=1,Nbx),j=1,Nby)
! 	  close(unit=2)
! 
! ! 	  open	(unit=3, file ='U0.dat', form='unformatted')
! ! 	  read	(unit=3) ((qold_global(2,i,j),i=1,Nbx),j=1,Nby)
! ! 	  close(unit=3)
! ! 
! ! 	  open	(unit=4,file='V0.dat', form='unformatted')
! ! 	  read	(unit=4) ((qold_global(3,i,j),i=1,Nbx),j=1,Nby)
! ! 	  close(unit=4)
	!  CASE(500)
	  
!--------------------------------------------------------------------------------------------------------------
	CASE(401) ! Dead Zone con Paredes Laterales
	print *, 'Canal con Pendiente 10x50cm'
	ho=0.01D0  !Profundidad Inicial: 1cm
	!u0=25.50D0/ho/10000.0D0 !Dividido por 100*100 para que M sea en m2/s y nocm2/s y u0 quede en m/s
	u0=0.255D0 !(255.0D0/(100*100*100))*(1.0D0/0.1D0)*(1/0.01D0)!Q*1/b*1/h=Q/(b*h) div
	
	!Q en cm3/s, dividido 100*100*100 para que quede en m3/s
	!b y h en m
	
! 
 	do i=1,Nbx; do j=1,Nby		
	
		if (z_global(i,j)<0.5D0) then
		
		qold_global(1,i,j)=ho!h
		qold_global(2,i,j)=0.0D0!u
		qold_global(3,i,j)=0.0D0!v
		else
		
		qold_global(1,i,j)=0.0D0
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		
	
		end if
	end do; end do
	
 	q0_global=qold_global
!-----------------------------------------------------
	CASE(402) ! Caso de E.Mignot
	!cavidad cuadrada 30 cm en canal de 4.3 metros,con pequeña 'perturbaciòn' en la altura de la cavidad 
	!Sin pendiente
	
	print *, 'Canal E.Mignot'
	ho=0.07D0!0.01D0/1.0D0  !Profundidad Inicial: 1cm, dividido por el factor de escala L
	!u0=25.50D0/ho/10000.0D0 !Dividido por 100*100 para que M sea en m2/s y nocm2/s y u0 quede en m/s
	u0=0.1667D0!(255.0D0/(100*100*100))*(1.0D0/0.1D0)*(1/0.01D0)!Q*1/b*1/h=Q/(b*h) div
	u0=u0/1.0D0!Factor de escala U
	!Q en cm3/s, dividido 100*100*100 para que quede en m3/s
	!b y h en m
	!ho=0.0202D0
	!u0=0.3698D0
! 
 	do i=1,Nbx; do j=1,Nby		
	
		if (z_global(i,j)>0.002D0) then
		
		qold_global(1,i,j)=0.0D0!h
		qold_global(2,i,j)=0.0D0!u
		qold_global(3,i,j)=0.0D0!v
		else
		
		qold_global(1,i,j)=0.99*ho
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		
		if (y_global(i,j)<0.3D0) then
		qold_global(1,i,j)=ho
		qold_global(2,i,j)=u0
		end if
		
		end if
	end do; end do
	
 	q0_global=qold_global
	
	CASE(403) 

	print *, 'Canal con dead zone y --perturbacion-- inicial (Kimura & Hosoda)'
	ho=0.01D0/1.0D0  !Profundidad Inicial: 1cm, dividido por el factor de escala L
	!u0=25.50D0/ho/10000.0D0 !Dividido por 100*100 para que M sea en m2/s y nocm2/s y u0 quede en m/s
	u0=0.255 !Q*1/b*1/h=Q/(b*h) div
	u0=u0/1.0D0!Factor de escala U
	!Q en cm3/s, dividido 100*100*100 para que quede en m3/s
	!b y h en m
	!ho=0.0202D0
	!u0=0.3698D0
! 
 	do i=1,Nbx; do j=1,Nby		
	
		if (z_global(i,j)>0.002D0) then
		
		qold_global(1,i,j)=0.0D0!h
		qold_global(2,i,j)=0.0D0!u
		qold_global(3,i,j)=0.0D0!v
		else
		
		qold_global(1,i,j)=0.99*ho
		qold_global(2,i,j)=0.0D0
		qold_global(3,i,j)=0.0D0
		
		if (y_global(i,j)<0.1D0) then
		qold_global(1,i,j)=ho
		qold_global(2,i,j)=u0
		end if		
		end if
	end do; end do	
 	q0_global=qold_global
	
	CASE(500)
	print *, 'Playa oblicua'
	ho=15.0D0  !Profundidad Inicial Màxima 
 	do i=1,Nbx; do j=1,Nby			
		if (ho-z_global(i,j)>0.0D0) then		
		  qold_global(1,i,j)=ho-z_global(i,j)!h
		  qold_global(2,i,j)=0.0D0!u
		  qold_global(3,i,j)=0.0D0!v
		else		
		  qold_global(1,i,j)=0.0D0
		  qold_global(2,i,j)=0.0D0
		  qold_global(3,i,j)=0.0D0	
		end if
	end do; end do	
 	q0_global=qold_global
	CASE(600)
	open	(unit=2, file ='initq1.dat')!, form='unformatted')
	read	(2,*) ((qold_global(1,i,j),i=1,Nbx),j=1,Nby)
	close(unit=2)

	open	(unit=3, file ='initq2.dat')!, form='unformatted')
	read	(3,*) ((qold_global(2,i,j),i=1,Nbx),j=1,Nby)
	close(unit=3)

	open	(unit=4,file='initq3.dat')!, form='unformatted')
	read	(4,*) ((qold_global(3,i,j),i=1,Nbx),j=1,Nby)
	close(unit=4)
	
END SELECT



END SUBROUTINE init_flowfield


SUBROUTINE ADIMENSION
!Function that aplies the adimensionalization to the initial conditions
!Funcion que aplica la adimensionalizacion a las condiciones inciales y a todo

USE global_variables
USE geometries
implicit none
integer :: i,j
real (kind=8)::U1,U2
	do i=1,Nbx; do j=1,Nby
			x_global(i,j)=x_global(i,j)/L
			y_global(i,j)=y_global(i,j)/L
			z_global(i,j)=z_global(i,j)/H
			
			
			
			qold_global(1,i,j)=qold_global(1,i,j)/H
			qold_global(2,i,j)=qold_global(2,i,j)/U
			qold_global(3,i,j)=qold_global(3,i,j)/U
			
			V_global(i,j)=sqrt((qold_global(2,i,j))**2.0D0+(qold_global(3,i,j))**2.0D0)
			
			C_global(i,j)=sqrt(qold_global(1,i,j)/FR2)
			VC(i,j)=V_global(i,j)+C_global(i,j)			
	end do; end do
! 	
! 			print*, 'hij',qold_global(1,i,j),'H',H
! 			print*, 'uij',qold_global(2,i,j),'U',U
! 			print*, 'vij',qold_global(3,i,j),'V',U
			

END SUBROUTINE ADIMENSION
