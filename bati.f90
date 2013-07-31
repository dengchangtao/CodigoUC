SUBROUTINE batimetria

!Programa que crea batimetria para ir probando el codigo, 
!la idea es usar una malla regular para las pruebas
USE global_variables
USE geometries
USE senales
implicit none

real (kind=8), dimension(:), allocatable :: xaux, yaux, zaux

real (kind=8):: dx,dy,zf,Lx, Ly, dx1, dx2, ri, rext, ds, xo, yo, sigma,r, ho, a, beta

integer	:: i,j,Ntot

SELECT CASE (caso)
  
  CASE(1)
	  !Nbx=20
	  !Nby=10
	  zf=0.0D0
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby))

	  dx=0.1D0
	  dy=0.2D0
	  !print *,dx

	  xaux(1)=dx
	  yaux(1)=dy


	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do

	  !print *, xaux

	  do j=2,Nby,1
	  yaux(j)=dy+yaux(j-1)
	  end do

	  !print*,yaux

	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zf
	  end do; end do
  
  CASE(2)	!Steady State over a bump
	  !Nbx=300
	  !Nby=60
	  Lx=25.0D0
	  Ly=5.0D0

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  !dx=0.0836D0
	  !dy=0.0847D0
	  dy=Ly/(Nby-1)

	  xaux(1)=0.0D0
	  yaux(1)=0.0D0
	  zaux(1)=0.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	      
	      if (xaux(i)>8.0D0.AND.xaux(i)<12.0D0) then
	      zaux(i)=0.2D0-0.05D0*(xaux(i)-10.0D0)**2.0D0
	      else
	      zaux(i)=0.0D0
	      end if
	      
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zaux(i)
	  end do; end do
	  !print*,'x=',x_global(:,1)
	  !print*,'y=',y_global(1,:)
	  !print*,'z=',z_global(:,1)
	  !stop
  CASE(3)	!Steady State over a bump
	  !Nbx=300
	  !Nby=60
	  Lx=25.0D0
	  Ly=10.0D0

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  !dx=0.0836D0
	  !dy=0.0847D0
	  dy=Ly/(Nby-1)

	  xaux(1)=0.0D0
	  yaux(1)=0.0D0
	  zaux(1)=0.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	      
	      if (xaux(i)>8.0D0.AND.xaux(i)<12.0D0) then
	      zaux(i)=0.2D0-0.05D0*(xaux(i)-10.0D0)**2.0D0
	      else
	      zaux(i)=0.0D0
	      end if
	      
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zaux(i)
	  end do; end do
	  !print*,'x=',x_global(:,1)
	  !print*,'y=',y_global(1,:)
	  !print*,'z=',z_global(:,1)
	  !stop 
  
  CASE(4)	!1D Dam-break
	 
	  !Nbx=201
	  !Nby=41
	  zf=0.0D0
	  
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby))
	  
	  Lx=50.0D0
	  Ly=5.0D0
	  
! 	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	  dx1=0.05D0
	  !dx2=0.05D0
	  xaux(1)=0
	  yaux(1)=0


	  do i=2,Nbx,1
! 	      if (xaux(i-1)<20) then
! 	      xaux(i)=dx1+xaux(i-1)
! 	      else if (xaux(i-1)>=20.AND.xaux(i-1)<30) then
! 	      xaux(i)=dx2+xaux(i-1)
! 	      else
	      xaux(i)=dx1+xaux(i-1)
! 	      end if
	  end do

	  

	  do j=2,Nby,1
	  yaux(j)=dy+yaux(j-1)
	  end do

	  

	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zf
	  end do; end do
	  
  CASE(5)	!2D Dam-Break
  	 

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby))
	  
	  Lx=200.0D0
	  Ly=200.0D0
	  
	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	  
	  xaux(1)=0.0D0
	  yaux(1)=0.0D0


	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do
	  do j=2,Nby,1
	  yaux(j)=dy+yaux(j-1)
	  end do
	  do i=1,Nbx; do j=1,Nby
        
            !if (xaux(i)>=98.0D0.AND.xaux(i)<=102.0D0.AND.yaux(j)>=0.0D0.AND.yaux(j)<=95.0D0) then
            !if (xaux(i)==100.0D0.AND.yaux(j)>=0.0D0.AND.yaux(j)<=95.0D0) then
	    if (i>=19.AND.i<=21.AND.j>=35.AND.j<=41) then
	    z_global(i,j)=10.0D0
            !else if (xaux(i)>=98.0D0.AND.xaux(i)<=102.0D0.AND.yaux(j)>=170.0D0.AND.yaux(j)<=200.0D0) then
	    !else if (xaux(i)==100.0D0.AND.yaux(j)>=170.0D0.AND.yaux(j)<=200.0D0) then
	    else if (i>=19.AND.i<=21.AND.j>=1.AND.j<=20) then
	    z_global(i,j)=10.0D0
            else 
            z_global(i,j)=0.0D0;
            end if
	  
	  end do; end do
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  end do; end do
  
  CASE(6)	!2D Dam-Break with an obstacle
  	 

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby))
	  
	  Lx=200.0D0
	  Ly=50.0D0
	  
	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	  
	  xaux(1)=0.0D0
	  yaux(1)=0.0D0


	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do
	  do j=2,Nby,1
	  yaux(j)=dy+yaux(j-1)
	  end do
	  do i=1,Nbx; do j=1,Nby
        
            if (xaux(i)>=50.0D0.AND.xaux(i)<=150.0D0.AND.yaux(j)>=20.0D0.AND.yaux(j)<=30.0D0) then
            z_global(i,j)=10.0D0
 
            else 
            z_global(i,j)=0.0D0;
            end if
	  
	  end do; end do
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  end do; end do  

  CASE(7)	!2D Dam-Break Cilindro Cartesiano
  	 

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby))
	  
	  Lx=50.0D0
	  Ly=50.0D0
	  
	  dx=1.0D0
	  dy=1.0D0
	  
	  xaux(1)=dx/2.0D0
	  yaux(1)=dx/2.0D0


	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do
	  do j=2,Nby,1
	  yaux(j)=dy+yaux(j-1)
	  end do
	  
	  
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=0.0D0
	  end do; end do
	  
  CASE(9)	!Flow at rest 2D
  
	xo=5.0D0
	yo=5.0D0
	sigma=0.5D0
	Lx=9.0D0
	Ly=9.0D0

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  !dx=0.0836D0
	  !dy=0.0847D0
	  dy=Ly/(Nby-1)

	  xaux(1)=0.0D0
	  yaux(1)=0.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1

	  end do; end do
	
	do i=1,Nbx; do j=1,Nby
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)	
	  
	  r=sqrt((x_global(i,j)-xo)**2+(y_global(i,j)-yo)**2)
	  z_global(i,j)=0.5D0*exp(-(r/sigma)**2) 
	
	end do; end do
  
  CASE(10)
  
  ! call pascua(Nbx,Nby)
  
  
  CASE(11) !Thacker's Curved Solution
  
	  xo=2.0D0
	  yo=2.0D0
	  Lx=4.0D0
	  Ly=4.0D0
	  ho=0.1D0
	  a=1.0D0
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)

	  xaux(1)=0.0D0
	  yaux(1)=0.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1

	  end do; end do
	
	  do i=1,Nbx; do j=1,Nby
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)	
	  
	  r=sqrt((x_global(i,j)-xo)**2+(y_global(i,j)-yo)**2)
	  
	  z_global(i,j)=-ho*(1.0D0-r**2/a**2)
	
	end do; end do 

	
   CASE(12) !Thacker's Planar Solution
  
	  xo=2.0D0
	  yo=2.0D0
	  Lx=4.0D0
	  Ly=4.0D0
	  ho=0.1D0
	  a=1.0D0
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)

	  xaux(1)=0.0D0
	  yaux(1)=0.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1

	  end do; end do
	
	  do i=1,Nbx; do j=1,Nby
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)	
	  
	  r=sqrt((x_global(i,j)-xo)**2+(y_global(i,j)-yo)**2)
	  
	  z_global(i,j)=-ho*(1.0D0-r**2/a**2)
	
	end do; end do

    CASE(13) !Friction parabolic bathy
    
	  Lx=10000.0D0
	  Ly=4.0D0
	  a=3000.0D0
	  ho=10.0D0

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  !dx=0.0836D0
	  !dy=0.0847D0
	  dy=Ly/(Nby-1)
	    
	  xaux(1)=-5000.0D0
	  yaux(1)=0.0D0
	  zaux(1)=ho*((xaux(1))/a)**2.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	      
	      zaux(i)=ho*((xaux(i))/a)**2.0D0
      
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zaux(i)
	  end do; end do
	  !print*,'x=',x_global(:,1)
	  !print*,'y=',y_global(1,:)
	  !print*,'z=',z_global(:,1)
	  !stop 
	  
    CASE(14) !Synolakis Wave
	  Ly=4.0D0
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))
	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	  open(unit=99,file='bathy_synolakis3.dat')
	  Do i=1,Nbx
	  read(99,*) xaux(i), zaux(i)
	  End Do
	  close(unit=99)
	  
	  yaux(1)=0.0D0
	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)+3.87
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zaux(i)+h01
	  end do; end do
! 	  print*,'x=',x_global(:,1)
! 	  print*,'y=',y_global(1,:)
!  	  print*,'z=',z_global(:,1)
! 	  pause 

     CASE(15) !Sudden gate closure
	  
	  Lx=25.0D0
	  Ly=1.0D0
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))
	  
	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	 

	  xaux(1)=0.0D0
	  
	  do i=2,Nbx
	      xaux(i)=dx+xaux(i-1)
	  end do
	  yaux(1)=0.0D0
	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=0.0D0
	  end do; end do
	  
      CASE(19) !Synolakis 2 (Marche, 2007)
      
	  Lx=100.0D0
	  Ly=1.0D0
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))
	  beta=0.0503D0
	  xo=60.15D0
	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	  
	  xaux(1)=0.0D0
	  
	  do i=2,Nbx
	      xaux(i)=dx+xaux(i-1)
	  end do
	  yaux(1)=0.0D0
	  
	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do
	 
	  do i=1,Nbx; do j=1,Nby
	 	 if (xaux(i)>xo) then
		 zaux(i)=(xaux(i)-xo)*tan(beta)
		 else    
		 zaux(i)=0.0D0
		 end if
		 end do; end do
    	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zaux(i)
	  end do; end do
	  
!        CASE(20)
! 	  
! 	  call pichbat(Nbx,Nby)
   
   CASE(23)	!Hydraulic Jump over a bump

	  Lx=25.0D0
	  Ly=1.0D0

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx),yaux(Nby), zaux(Nbx))

	  dx=Lx/(Nbx-1)
	  !dx=0.0836D0
	  !dy=0.0847D0
	  dy=Ly/(Nby-1)

	  xaux(1)=0.0D0
	  yaux(1)=0.0D0
	  zaux(1)=0.0D0

	  do i=2,Nbx,1
	      xaux(i)=dx+xaux(i-1)
	      
	      if (xaux(i)>8.0D0.AND.xaux(i)<12.0D0) then
	      zaux(i)=0.2D0-0.05D0*(xaux(i)-10.0D0)**2.0D0
	      else
	      zaux(i)=0.0D0
	      end if
	      
	  end do

	  do j=2,Nby
	      yaux(j)=dy+yaux(j-1)
	  end do

	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(i)
	  y_global(i,j)=yaux(j)
	  z_global(i,j)=zaux(i)
	  end do; end do
	  
	  CASE(99)	!Ensayo para un perfil de playa

	  Lx=19.0D0
	  Ly=24.8D0

	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  allocate(xaux(Nbx*Nby),yaux(Nbx*Nby), zaux(Nbx*Nby))

	  dx=Lx/(Nbx-1)
	  dy=Ly/(Nby-1)
	  Ntot=Nbx*Nby
	  print *, 'dx dy', dx, dy, Ntot
	  !Lee batimetria Leandro
	  open(unit=99,file='bathy_run34.dat')
	  Do i=1,Ntot
	  read(99,*) xaux(i), yaux(i), zaux(i)
	  End Do
	  close(unit=99)
!	  print *, 'bathy', xaux
!	  yaux(1)=0.0D0
!	  do j=2,Nby
!	      yaux(j)=dy+yaux(j./t-1)
!	  end do
	  do i=1,Nbx,1; do j=1,Nby,1
	  x_global(i,j)=xaux(j+(i-1)*Nby)+0.00
	  y_global(i,j)=yaux(j+(i-1)*Nby)
	  z_global(i,j)=zaux(j+(i-1)*Nby)+0.00
	  end do; end do
	  print *, 'bathy', x_global(1,1)
	  ! Ensayo para caso 2D
	  
	  
	  !Sismo 27F
	  CASE(100)
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  open	(unit=2, file ='gridX.dat', form='unformatted')
	  read	(unit=2) ((x_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=2)

	  open	(unit=3, file ='gridY.dat', form='unformatted')
	  read	(unit=3) ((y_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=3)

	  open	(unit=4,file='gridZ.dat', form='unformatted')
	  read	(unit=4) ((z_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=4)
	  
	  ! Traslado batimetria en 207.9032 m.s.n.m
	  print*,'zmin=',minval(z_global)
	  z_global=z_global+abs(minval(z_global))
	  

!-------------------------------------------------------------------------------
	  !Columbia River
	  CASE(200)
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  open	(unit=2, file ='gridX.dat', form='unformatted')
	  read	(unit=2) ((x_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=2)

	  open	(unit=3, file ='gridY.dat', form='unformatted')
	  read	(unit=3) ((y_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=3)

	  open	(unit=4,file='gridZ.dat', form='unformatted')
	  read	(unit=4) ((z_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=4)
	  
	  ! Traslado batimetria en 207.9032 m.s.n.m
	  print*,'zmin=',minval(z_global)
	  z_global=z_global+abs(minval(z_global))

!-------------------------------------------------------------------------------
	  !El Quisco
	  CASE(300)
	  allocate (x_global(Nbx,Nby), y_global(Nbx,Nby),z_global(Nbx,Nby))
	  open	(unit=2, file ='gridX.dat', form='unformatted')
	  read	(unit=2) ((x_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=2)

	  open	(unit=3, file ='gridY.dat', form='unformatted')
	  read	(unit=3) ((y_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=3)

	  open	(unit=4,file='gridZ.dat', form='unformatted')
	  read	(unit=4) ((z_global(i,j),i=1,Nbx),j=1,Nby)
	  close(unit=4)
	  
	  ! Traslado batimetria en  m.s.n.m
	  print*,'zmin=',minval(z_global)
	  z_global=z_global+abs(minval(z_global))
!
	  
 END SELECT
 


END SUBROUTINE batimetria
