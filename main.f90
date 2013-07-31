!NSWE Finite Volume Model
!Curvilinear coordinate system
!Well-Balanced, Shock-Capturing
!Approximate Riemann Solver VFRoe-ncv plus 2nd order Hydrostatic Reconstruction
!Semi-implicit friction step

!Developed by Maricarmen Guerra P.

PROGRAM MAIN !Ponerle un nombre decente

USE global_variables !Use the global variables module
USE geometries	!Use geometries module
USE senales	!Use senales module for inflow or outflow BC
implicit none 

!Declare variables that are used here and only here if it is necessary

integer::i,j,k
real (kind=8) :: mindxdy,minxieta,maxUC, dtreal, maxV, maxC, zmax, maxS1, maxS2, pVol1
real (kind=8), dimension(2) :: xieta, loc
real (kind=8), dimension(:), allocatable:: dxdy


g=9.812D0
hmin=1.0E-7
it=0.0D0
dt=0.0D0
t=0.0D0!0.0D0
treal=0.0D0!960.0D0 !0.0D0

! Initialize calculation

call init !init function: reads input data and bathymethry, initialize flowfield, calculates metrics


if (outopt==1) then
call outputmat(treal)	!write initial conditions into a file
pause
else
call outputtec(treal)
end if
 !print*, 'impreso'


call massbalance	!Calculates the initial volume and mass of water

write(*,100)
100 FORMAT ('Start Solver')
print*,'------------------------------'

zmax=maxVal(z_global(:,:))

allocate(dxdy((Nbx-1)*Nby+Nbx*(Nby-1)))
k=1
do i=1,Nbx-1; do j=1,Nby
dxdy(k)=x_global(i+1,j)-x_global(i,j)
k=k+1
end do; end do

do i=1,Nbx; do j=1,Nby-1
dxdy(k)=y_global(i,j+1)-y_global(i,j)
k=k+1
end do; end do

!MAIN CYCLE!

!Do while (it<=1) !to find errors

DO while(t<=tfinal)
      
!1. Calculate time step with the CFL condition

xieta=(/dxi,deta/)
minxieta=minval(xieta)
maxS1=maxval(S1_global)
maxS2=maxval(S2_global)

maxUC=maxval((/maxS1,maxS2/))
!maxUC=maxval(VC)

dt=CFL*minxieta/maxUC		!This dt is adimensional
dtreal=dt*L/U

!If the new time is bigger than final time, exit the cycle

IF ((treal+dtreal).GE.tfinal) THEN
!Termino del programa
print *, 'Simulation Ended'
print *, 'Final Time of Computation= ', t
print *, 'Iteraciones', it
stop

ELSE

!2.Calling main_solver, which solves the 4 stages of RK method, and calculates qnew_global(h,u,v)

if (fopt==0) then

  if (rk==1) then
  call solver1 !Solver RK4, sin fricción
  else
  call solver2 !Solver RK2, sin fricción
  end if

else

  if (rk==1) then
  call solverf4 !Solver RK4, con fricción
  else
  call solverf2 !Solver RK2, con fricción
  end if

end if



t=t+dt
treal=t*L/U
it=it+1.0D0


!Print iteration information on the screen

print*, 'dt= ', dt
write(*,199) it
199 FORMAT ('Iteracion= ',F20.0)
print*, 'Time= ', t

!Adimensional Mass Balance Verification
call massbalance
write(*,170) Pvol
170 FORMAT ('%Volumen= ',F7.3 ,'%')
print*,'------------------------------'
!Pvol1=Pvol



!3.Dimensionalize results and write results into a file
if (outopt==1) then
call outputmat(treal)
else
call outputtec(treal)
end if


!Actualization of the global variables (adimensionalized)
qold_global=qnew_global


END IF

END DO

!FIN!

END PROGRAM MAIN
