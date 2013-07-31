! Programa para leer batimetrias y señales de entrada

Program Mallas


!ifort -o Mallas LeeBatiSenales.f90 /usr/local/tecplot/lib/tecio64.a -lstdc++

implicit none

!----------------------------------------------------------------
!No modificar
real (kind=8), dimension(:,:),allocatable:: xtec,ytec,ztec,h0tec, x_global, y_global, z_global

character (len = 256) :: filename
! variables to enable writing of TecPlot binary (*.plt) files

  integer (kind = 4)           :: TecIni, TecDat, TecZne
  integer (kind = 4)           :: TecEnd
  integer (kind = 4)           :: VIsDouble = 0
  integer (kind = 4)           :: Debug = 1
  
  integer (kind = 4)           :: III

  character (len = 1)          :: nullchr = char(0)
!--------------------------------------------------------------

real(kind=8), dimension (:),allocatable::ts, ts2, xt, yt, zt, h0,u0,v0,SL
real(kind=8), dimension (:),allocatable::qinx, qiny, qinxE,qinyE, hin1, hin2, etain1,etain2
integer::Nx,Ny,Ns,i,j,k, Ns2
real (kind=8):: dt, Tf

Nx=301 !Numero de nodos en X--
Ny=102 !Numero de nodos en Y--

Ns=301 !Tamano de la senal--
Ns2=361 !Tamano de la senal 2

allocate(xt(Nx*Ny),yt(Nx*Ny),zt(Nx*Ny),ts(Ns),qinx(Ny*Ns),qiny(Ny*Ns2),qinxE(Nx*Ns2),qinyE(Nx*Ns2),h0(Nx*Ny),u0(Nx*Ny),v0(Nx*Ny),SL(Nx*Ny),hin1(Ny*Ns),hin2(Ny*Ns),etain1(Ny*Ns),etain2(Ny*Ns2),ts2(Ns2))
allocate(x_global(Nx,Ny),y_global(Nx,Ny),z_global(Nx,Ny))

!Batimetria
open	(unit=1, file ='XYZconparedes.dat') !Abrir archivo de batimetria creado en matlab, lista de x, y ,z

Do i=1,Nx*Ny
read	(1,*) xt(i), yt(i), zt(i)
end Do
close(unit=1)

open	(unit=5, file ='gridX.dat', form='unformatted')
write	(unit=5) (xt(i),i=1,Nx*Ny)
close(unit=5)

open	(unit=6, file ='gridY.dat', form='unformatted')
write	(unit=6) (yt(i),i=1,Nx*Ny)
close(unit=6)

open	(unit=7,file='gridZ.dat', form='unformatted')
write	(unit=7) (zt(i),i=1,Nx*Ny)
close(unit=7)
print*, Nx, Ny, Ny*Nx

pause

k=1
Do j=1,Ny
Do i=1,Nx
x_global(i,j)=xt(k)
y_global(i,j)=yt(k)
z_global(i,j)=zt(k)
k=k+1
end Do
end Do



! !Señal canal con deadzone, caso 15x15
Tf=300.D0 !2mins
dt=1.0D0 !datos cada 1 sec
ts(1)=0.0D0
Do i=2,Ns
ts(i)=ts(i-1)+dt
end Do

print*, ts
!Borde xi=0
open	(unit=7, file ='Suh_xi0.dat')
Do i=1,Ny*Ns
read	(7,*) qinx(i)
end Do
close(unit=7)

open	(unit=8, file ='Svh_xi0.dat')
Do i=1,Ny*Ns
read	(8,*) qiny(i)
end Do
close(unit=8)

open	(unit=81, file ='Seta_xi0.dat')
Do i=1,Ny*Ns
read	(81,*) etain1(i)
end Do
close(unit=81)

open	(unit=112,file='timeS1.dat', form='unformatted')
 write	(unit=112) (ts(i),i=1,Ns)
 close(unit=11)
 
 open	(unit=9, file ='qsx1.dat', form='unformatted')

write	(unit=9) (qinx(i),i=1,Ny*Ns)
close(unit=9)

open	(unit=10, file ='qsy1.dat', form='unformatted')
write	(unit=10) (qiny(i),i=1,Ny*Ns)
close(unit=10)

open	(unit=101, file ='etas1.dat', form='unformatted')
write	(unit=101) (etain1(i),i=1,Ny*Ns)
close(unit=101)





! !Condicion Inicial
! open	(unit=800, file ='CI_H.dat')
! Do i=1,Nx*Ny
! read	(800,*) SL(i), h0(i), u0(i)
! end Do
! close(unit=800)
! 
! open	(unit=900, file ='h0.dat', form='unformatted')
! 
! write	(unit=900) (h0(i),i=1,Nx*Ny)
! close(unit=900)
! 
! open	(unit=1000, file ='U0.dat', form='unformatted')
! write	(unit=1000) (u0(i),i=1,Nx*Ny)
! close(unit=1000)
! 
! open	(unit=1100,file='V0.dat', form='unformatted')
! write	(unit=1100) (0.0D0,i=1,Nx*Ny)
! close(unit=1100)
! 
! 
! ! Señal
! Tf= 2.880e3 !~8 horas
! dt=10.0D0 !datos cada 10 sec
! ts(1)=0.0D0
! 
! Do i=2,Ns
! ts(i)=ts(i-1)+dt
! end Do
! 
! !Borde xi=0
! 
! ! Estado Estacionario
! 
! open	(unit=81, file ='SetaSS_xi0.dat')
! Do i=1,Ny*Ns
! read	(81,*) etain1(i)
! end Do
! close(unit=81)
! 
! open	(unit=101, file ='etas1H.dat', form='unformatted')
! write	(unit=101) (etain1(i),i=1,Ny*Ns)
! close(unit=101)
! 
! ! TSUNAMI
! 
! ! open	(unit=7, file ='Su_xi0.dat')
! ! Do i=1,Ny*Ns
! ! read	(7,*) qinx(i)
! ! end Do
! ! close(unit=7)
! ! 
! ! open	(unit=8, file ='Sv_xi0.dat')
! ! Do i=1,Ny*Ns
! ! read	(8,*) qiny(i)
! ! end Do
! ! close(unit=8)
! ! 
! ! open	(unit=81, file ='Seta_xi0.dat')
! ! Do i=1,Ny*Ns
! ! read	(81,*) etain1(i)
! ! end Do
! ! close(unit=81)
! ! 
! ! open	(unit=9, file ='qsx1.dat', form='unformatted')
! ! 
! ! write	(unit=9) (qinx(i),i=1,Ny*Ns)
! ! close(unit=9)
! ! 
! ! open	(unit=10, file ='qsy1.dat', form='unformatted')
! ! write	(unit=10) (qiny(i),i=1,Ny*Ns)
! ! close(unit=10)
! ! 
! ! open	(unit=101, file ='etas1.dat', form='unformatted')
! ! write	(unit=101) (etain1(i),i=1,Ny*Ns)
! ! close(unit=101)
! ! 
! ! 
! open	(unit=11,file='timeS1.dat', form='unformatted')
! write	(unit=11) (ts(i),i=1,Ns)
! close(unit=11)
! 
! 
! !Borde Xi=N
! 
! 
! 
! Tf= 3.6e3 !~8 horas
! dt=100.0D0 !datos cada 100 sec
! ts2(1)=0.0D0
! 
! Do i=2,Ns2
! ts2(i)=ts2(i-1)+dt
! end Do
! 
! 
! open	(unit=72, file ='Su_xiN.dat')
! Do i=1,Ny*Ns2
! read	(72,*) qinxE(i)
! end Do
! close(unit=72)
! 
! open	(unit=82, file ='Sv_xiN.dat')
! Do i=1,Ny*Ns2
! read	(82,*) qinyE(i)
! end Do
! close(unit=82)
! 
! open	(unit=812, file ='Seta_xiN.dat')
! Do i=1,Ny*Ns2
! read	(812,*) etain2(i)
! end Do
! close(unit=812)
! 
! open	(unit=92, file ='qsx2.dat', form='unformatted')
! 
! write	(unit=92) (qinxE(i),i=1,Ny*Ns2)
! close(unit=92)
! 
! open	(unit=102, file ='qsy2.dat', form='unformatted')
! write	(unit=102) (qinyE(i),i=1,Ny*Ns2)
! close(unit=102)
! 
! open	(unit=1012, file ='etas2.dat', form='unformatted')
! write	(unit=1012) (etain2(i),i=1,Ny*Ns2)
! close(unit=1012)
! 
! 
! open	(unit=112,file='timeS2.dat', form='unformatted')
! write	(unit=112) (ts2(i),i=1,Ns2)
! close(unit=11)


!Writing files
! H0 para ver en Tecplot


!Batimetria en TECPLOT
  
      !filename='SOL2D.'//trim(number)
      
      write(filename, fmt = '(a,i6.6)')'Bati'
      
      I = TecIni('Bati'//NULLCHR,			&   ! title of file
		'X, Y, Z, H'//NULLCHR,	&   ! list of variables
		'BatiCR.plt'//NULLCHR,	&   ! output file name
		'.'//NULLCHR,						&
		Debug,								&
		VIsDouble)
  

allocate(xtec(Nx,Ny),ytec(Nx,Ny),ztec(Nx,Ny),h0tec(Nx,Ny))

  DO j=1,Ny; DO i=1,Nx
    
    xtec(i,j)=x_global(i,j)
    ytec(i,j)=y_global(i,j)
    ztec(i,j)=z_global(i,j)
    h0tec(i,j)=1.5D0

  END DO; END DO;

	I = TecZne('Zone'//NULLCHR,    &     
			Nx,			           &
			Ny,			           &
			1, &
			'BLOCK'//NULLCHR,		   &
			NULLCHR//NULLCHR)
			! total number of points
		III = Nx * Ny
			! write each variable
			! the last argument in the following calls indicates
			! the precision of the the variable,
			! 0 = single, 1 = double
		I   = TecDat(III,xtec,1)
		I   = TecDat(III,ytec,1)
		I   = TecDat(III,ztec,1)
		I   = TecDat(III,h0tec,1)

deallocate(xtec,ytec,ztec,h0tec)
    I   = TecEnd()
close(1)




END Program Mallas
