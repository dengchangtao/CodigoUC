SUBROUTINE input_control

USE global_variables
USE senales
implicit none

!integer 	::n,nz,i,j,ii
! 1	!Caso  
! 1	!(1)Lee Bati o (2)Crea bati
! 10	!Tiempo Final (s)
! 0.8	!CFL
! 10	!Nbx
! 5	!Nby
! 1	!dXi
! 1	!dEta
! 1	!L
! 1	!H
! 1	!U
! 1	!CB 1
! 1	!CB 2
! 1	!CB 3
! 1	!CB 4
! 1	dit
! 1e-6	kappa
! Runge Kutta Method: 1=RK4, 2=RK2
! mmopt: 1=MinMod, 2=Superbee, 3=MC
!Outopt: 1=Matlab, 2=Tecplot

open(1,file='input.dat')
read(1,*) caso  
read(1,*) batiopt !1= Read, 2=Create
read(1,*) inputopt !1=Read, 2=Create
read(1,*) tfinal !Tiempo final de simulacion en segundos
read(1,*) CFL
read(1,*) Nbx
read(1,*) Nby
read(1,*) dxi
read(1,*) deta
read(1,*) L
read(1,*) H
read(1,*) U
read(1,*) CB(1)	!Boundary Condition for xi_0,j
!Diferentes tipos de CB


IF (CB(1)==4) THEN !GENABS 1, 2, 3 o 9
read(1,*) GA1
read(1,*) Nsenal1
call readGA(1,GA1,Nsenal1)
END IF

IF (CB(1)==5) THEN !Outflow en 1, se fija una altura o Inflow en 1
read(1,*) IO1
read(1,*) Nsenal1
call readIO(1,IO1,Nsenal1)
END IF

! Esto es equivalente a CB(1)==4, ya que ahi entra a los genabs, 
! IF (CB(1)==9) THEN
! read(1,*) GA1
! read(1,*) Nsenal1
! call readGA(1,GA1,Nsenal1)
! END IF

read(1,*) CB(2)	!Boundary Condition for xi_Nbx,j
IF (CB(2)==4) THEN
read(1,*) GA2
read(1,*) Nsenal2
call readGA(2,GA2,Nsenal2)
END IF

IF (CB(2)==5) THEN !Outflow or Inflow en 2, se fija una altura o Inflow en 1
read(1,*) IO2
read(1,*) Nsenal2

call readIO(2,IO2,Nsenal2)

END IF


read(1,*) CB(3)	!Boundary Condition for eta_i,0
IF (CB(3)==4) THEN
read(1,*) GA3
read(1,*) Nsenal3
call readGA(3,GA3,Nsenal3)
END IF
IF (CB(3)==5) THEN !Outflow or Inflow en 3, se fija una altura o Inflow en 1
read(1,*) IO3
read(1,*) Nsenal3
call readIO(3,IO3,Nsenal3)
END IF
read(1,*) CB(4)	!Boundary Condition for eta_i,Nby
IF (CB(4)==4) THEN
read(1,*) GA4
read(1,*) Nsenal4
call readGA(4,GA4,Nsenal4)
END IF
IF (CB(4)==5) THEN !Outflow or Inflow en 4, se fija una altura o Inflow en 1
read(1,*) IO4
read(1,*) Nsenal4
call readIO(4,IO4,Nsenal4)
END IF

read(1,*) dit	!dit to print results, write files every pdt iterations
read(1,*) kappa !To consider a 0.0 value
read(1,*) rk !Runge Kutta method 1=Rk4, 2=Rk2
read(1,*) mmopt !1=Minmod, 2=Superbee Limiters
read(1,*) fopt !Con o sin friccion, No == 0, Si ==1
IF (fopt==0) THEN
 Cf=0
 Coef=0.0D0
END IF
IF (fopt==1) THEN
read(1,*) fM !Si fM=1, un coeficiente de friccion, si fM=2, se usará la matrzi de coeficientes
read(1,*) Cf !Tipo Fricción: Manning ==1, Chezy ==2, Sampson==3, 0=ninguno
  if (fM==1) then !Un coeficiente para todo el dominio
  read(1,*) Coef !Coeficiente friccion segun sea el caso ya adimensionalizado!!!!
  end if
END IF
read(1,*) outopt !1=Matlab, 2=Tecplot files

close(1)

write(*,100) caso
100 FORMAT ('Caso: ', T25, I4)

write(*,110) tfinal
110 FORMAT ('T_final (s): ', T25, F8.2)

write(*,120) CFL
120 FORMAT ('CFL: ', T25, F5.2)

write(*,130) Nbx
130 FORMAT ('Nbx: ', T25,  I6)

write(*,140) Nby
140 FORMAT ('Nby: ', T25, I6)


write(*,150)  L, H, U
150 FORMAT ('Escalas (L, H, U): ', T25, 3F5.2)

write(*,160) CB(1), CB(2), CB(3), CB(4)
160 FORMAT ('Condiciones de Borde: ', T25, 4I2)

if (mmopt==1) then
write(*,161) mmopt
161 FORMAT ('Usando limitador MINMOD, mmopt= ', T30, 4I2)
else if (mmopt==2) then
write(*,162) mmopt
162 FORMAT ('Usando limitador SUPERBEE, mmopt= ', T30, 4I2)
else
write(*,163) mmopt
163 FORMAT ('Usando limitador MC, mmopt= ', T30, 4I2)
end if
pause

END SUBROUTINE input_control