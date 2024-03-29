SUBROUTINE solverf4

! Solver includes the friction step before de SWE RK steps
! 4th order RK scheme
!4th order time integration
USE global_variables
USE geometries
USE senales
implicit none

!Local variables definition
integer	:: i,j
real (kind=8) :: dtf
real (kind=8),dimension(2):: locu1, locv1, DFloc
real (kind=8),dimension(:,:,:),allocatable :: q1,q2,q3,q4, &
				     q1T,q2T,q3T,q4T, &
				     xi_T, eta_T,zT, &
				     F1mas,F2mas,F3mas,F4mas, &
				     G1mas,G2mas,G3mas,G4mas,  &
				     F1menos,F2menos,F3menos,F4menos, &
				     G1menos,G2menos,G3menos,G4menos,  &
				     qaux, qaux1, qaux2, qaux3, qaux4, &
				     SC, k1, k2, k3, k4, suma, &				!T Includes ghost cells
				     q1F,q2F,q3F,q4F
				     
!Allocate intermediate state variables
allocate( q1(3,Nbx,Nby),q2(3,Nbx,Nby),q3(3,Nbx,Nby),q4(3,Nbx,Nby), &
	  F1mas(3,Nbx,Nby),F2mas(3,Nbx,Nby),F3mas(3,Nbx,Nby),F4mas(3,Nbx,Nby), &
	  G1mas(3,Nbx,Nby),G2mas(3,Nbx,Nby),G3mas(3,Nbx,Nby),G4mas(3,Nbx,Nby), &
	  F1menos(3,Nbx,Nby),F2menos(3,Nbx,Nby),F3menos(3,Nbx,Nby),F4menos(3,Nbx,Nby), &
	  G1menos(3,Nbx,Nby),G2menos(3,Nbx,Nby),G3menos(3,Nbx,Nby),G4menos(3,Nbx,Nby), &
	  qaux(3,Nbx,Nby),qaux1(3,Nbx,Nby),qaux2(3,Nbx,Nby),qaux3(3,Nbx,Nby),qaux4(3,Nbx,Nby), &
	  q1F(3,Nbx,Nby),q2F(3,Nbx,Nby),q3F(3,Nbx,Nby),q4F(3,Nbx,Nby), &
	  SC(3,Nbx,Nby), k1(3,Nbx,Nby), k2(3,Nbx,Nby), k3(3,Nbx,Nby), k4(3,Nbx,Nby), suma(3,Nbx,Nby))
	  
	  !F1(3,Nbx,Nby),F2(3,Nbx,Nby),F3(3,Nbx,Nby),F4(3,Nbx,Nby), &
	  !G1(3,Nbx,Nby),G2(3,Nbx,Nby),G3(3,Nbx,Nby),G4(3,Nbx,Nby), &

!Allocate intermediate state variables with ghost cells
allocate(q1T(3,Nbx+4,Nby+4),q2T(3,Nbx+4,Nby+4),q3T(3,Nbx+4,Nby+4),q4T(3,Nbx+4,Nby+4),zT(1,Nbx+4,Nby+4))
allocate(xi_T(2,Nbx+4,Nby+4),eta_T(2,Nbx+4,Nby+4))

!-----------------------------
!qold,q1,q2,q3,q4,qnew=(h,u,v)
!-----------------------------




! First Friction Step, qold_global-->q1F
dtf=0.5D0*dt

call friccion(qold_global,Nbx,Nby,dtf,FR2,hmin,kappa,q1F,Cf,MCoef)


!I need qaux=(h,hu,hv) for the time integration, and q1 must have (h,u,v) for the bcs and fluxes functions:

DO i=1,Nbx; DO j=1,Nby
  qaux(1,i,j)=q1F(1,i,j)
  qaux(2,i,j)=q1F(1,i,j)*q1F(2,i,j)
  qaux(3,i,j)=q1F(1,i,j)*q1F(3,i,j)
END DO; END DO

!Cleaning qaux

Do i=1,Nbx; Do j=1,Nby
  if (qaux(1,i,j).le.hmin) then
  qaux(1,i,j)=0.0D0
  qaux(2,i,j)=0.0D0
  qaux(3,i,j)=0.0D0
  end if
  
  if (abs(qaux(2,i,j)).le.kappa) then
  qaux(2,i,j)=0.0D0
  end if
  
  if (abs(qaux(3,i,j)).le.kappa) then
  qaux(3,i,j)=0.0D0
  end if
  
End do; End do


!First RK Stage

!Calculates q(n+1/2*)

!Boundary Conditions

call bcs(fopt,Cf,MCoef,1,t,0.5D0*dt,FR2,caso,q1F,z_global,Nbx,Nby,CB,xi_global,eta_global,aj_global,dxi,deta,q1T,xi_T,eta_T,zT)

call fluxes(CB,mmopt,hmin,q1T,zT,xi_T,eta_T,dxi,deta,Nbx,Nby,FR2,F1mas,F1menos,G1mas,G1menos,SC)

DO i=1,Nbx; DO j=1,Nby
k1(:,i,j)=(SC(:,i,j)-(aj_global(i,j)/dxi)*(F1menos(:,i,j)-F1mas(:,i,j))-(aj_global(i,j)/deta)*(G1menos(:,i,j)-G1mas(:,i,j)))

qaux1(:,i,j)=qaux(:,i,j)+0.5D0*dt*k1(:,i,j)

END DO; END DO


!q1
DO i=1,Nbx; DO j=1,Nby

! !Prueba
! if (i==2.and.j==5) then
! print*, 'solver', qaux1(:,i,j)
! end if

  IF (qaux1(1,i,j).le.hmin) THEN	!Dry Cell
    qaux1(1,i,j)=0.0D0
    qaux1(2,i,j)=0.0D0
    qaux1(3,i,j)=0.0D0
  END IF
  
  q1(1,i,j)=qaux1(1,i,j)
  
  IF (q1(1,i,j)==0.0D0) THEN
  q1(2,i,j)=qaux1(2,i,j)
  q1(3,i,j)=qaux1(3,i,j)
  
  ELSE
  q1(2,i,j)=qaux1(2,i,j)/qaux1(1,i,j)
  q1(3,i,j)=qaux1(3,i,j)/qaux1(1,i,j)
  END IF
  
  IF (abs(q1(2,i,j)).le.kappa) THEN
  q1(2,i,j)=0.0D0
  END IF
  
  IF (abs(q1(3,i,j)).le.kappa) THEN
  q1(3,i,j)=0.0D0
  END IF
  
END DO; END DO

! -----END of First Step-----



!Second RK Stage
!Calculates q(n+1/2**)

call bcs(fopt,Cf,MCoef,2,t,0.5D0*dt,FR2,caso,q1,z_global,Nbx,Nby,CB,xi_global,eta_global,aj_global,dxi,deta,q2T,xi_T,eta_T,zT)

call fluxes(CB,mmopt,hmin,q2T,zT,xi_T,eta_T,dxi,deta,Nbx,Nby,FR2,F2mas,F2menos,G2mas,G2menos,SC)

DO i=1,Nbx; DO j=1,Nby
k2(:,i,j)=(SC(:,i,j)-(aj_global(i,j)/dxi)*(F2menos(:,i,j)-F2mas(:,i,j))-(aj_global(i,j)/deta)*(G2menos(:,i,j)-G2mas(:,i,j)))
qaux2(:,i,j)=qaux(:,i,j)+0.5D0*dt*k2(:,i,j)

END DO; END DO

!q2
DO i=1,Nbx; DO j=1,Nby

  IF (qaux2(1,i,j)<=hmin) THEN	!Dry Cell
    qaux2(1,i,j)=0.0D0
    qaux2(2,i,j)=0.0D0
    qaux2(3,i,j)=0.0D0
   
  END IF
  
  q2(1,i,j)=qaux2(1,i,j)
  
  IF (q2(1,i,j)==0.0D0) THEN
  q2(2,i,j)=0.0D0
  q2(3,i,j)=0.0D0
  
  ELSE
  q2(2,i,j)=qaux2(2,i,j)/qaux2(1,i,j)
  q2(3,i,j)=qaux2(3,i,j)/qaux2(1,i,j)
  END IF
  
  IF (abs(q2(2,i,j)).le.kappa) THEN
  q2(2,i,j)=0.0D0
  END IF
  
  IF (abs(q2(3,i,j)).le.kappa) THEN
  q2(3,i,j)=0.0D0
  END IF
  
END DO; END DO


!Third RK Stage
!Calculates q(n+1*)

call bcs(fopt,Cf,MCoef,3,t+0.5D0*dt,0.5D0*dt,FR2,caso,q2,z_global,Nbx,Nby,CB,xi_global,eta_global,aj_global,dxi,deta,q3T,xi_T,eta_T,zT)

call fluxes(CB,mmopt,hmin,q3T,zT,xi_T,eta_T,dxi,deta,Nbx,Nby,FR2,F3mas,F3menos,G3mas,G3menos,SC)

DO i=1,Nbx; DO j=1,Nby
k3(:,i,j)=(SC(:,i,j)-(aj_global(i,j)/dxi)*(F3menos(:,i,j)-F3mas(:,i,j))-(aj_global(i,j)/deta)*(G3menos(:,i,j)-G3mas(:,i,j)))

qaux3(:,i,j)=qaux(:,i,j)+dt*k3(:,i,j)

END DO; END DO

!q3
DO i=1,Nbx; DO j=1,Nby

  IF (qaux3(1,i,j)<=hmin) THEN	!Dry Cell
    qaux3(1,i,j)=0.0D0
    qaux3(2,i,j)=0.0D0
    qaux3(3,i,j)=0.0D0
  END IF
  
  q3(1,i,j)=qaux3(1,i,j)
  
  IF (q3(1,i,j)==0.0D0) THEN
  q3(2,i,j)=0.0D0
  q3(3,i,j)=0.0D0
  ELSE
   
  q3(2,i,j)=qaux3(2,i,j)/qaux3(1,i,j)
  q3(3,i,j)=qaux3(3,i,j)/qaux3(1,i,j)
  END IF
  
  IF (abs(q3(2,i,j)).le.kappa) THEN
  q3(2,i,j)=0.0D0
  END IF
  
  IF (abs(q3(3,i,j)).le.kappa) THEN
  q3(3,i,j)=0.0D0
  END IF
  
  
END DO; END DO


!4th RK Stage, 
!Calculates q(new)

call bcs(fopt,Cf,MCoef,4,t+0.5D0*dt,0.5D0*dt,FR2,caso,q3,z_global,Nbx,Nby,CB,xi_global,eta_global,aj_global,dxi,deta,q4T,xi_T,eta_T,zT)
call fluxes(CB,mmopt,hmin,q4T,zT,xi_T,eta_T,dxi,deta,Nbx,Nby,FR2,F4mas,F4menos,G4mas,G4menos,SC)

DO i=1,Nbx; DO j=1,Nby

k4(:,i,j)=(SC(:,i,j)-(aj_global(i,j)/dxi)*(F4menos(:,i,j)-F4mas(:,i,j))-(aj_global(i,j)/deta)*(G4menos(:,i,j)-G4mas(:,i,j)))
suma(:,i,j)=dt/6.0D0*(k1(:,i,j)+2.0D0*k2(:,i,j)+2.0D0*k3(:,i,j)+k4(:,i,j))
qaux4(:,i,j)=qaux(:,i,j)+suma(:,i,j)

END DO; END DO


!q_new
DO i=1,Nbx; DO j=1,Nby

  IF (qaux4(1,i,j)<=hmin) THEN	!Dry Cell
    qaux4(1,i,j)=0.0D0
    qaux4(2,i,j)=0.0D0
    qaux4(3,i,j)=0.0D0
  END IF
  
  q4(1,i,j)=qaux4(1,i,j)
  
  IF (q4(1,i,j)==0.0D0) THEN
  q4(2,i,j)=0.0D0
  q4(3,i,j)=0.0D0
  
  ELSE
  
  q4(2,i,j)=qaux4(2,i,j)/qaux4(1,i,j)
  q4(3,i,j)=qaux4(3,i,j)/qaux4(1,i,j)
  END IF
  
  IF (abs(q4(2,i,j)).le.kappa) THEN
  q4(2,i,j)=0.0D0
  END IF
  
  IF (abs(q4(3,i,j)).le.kappa) THEN
  q4(3,i,j)=0.0D0
  END IF
  

END DO; END DO

! Second Friction Step, qold_global-->q1F
dtf=0.5D0*dt
call friccion(q4,Nbx,Nby,dtf,FR2,hmin,kappa,q4F,Cf,MCoef)

!New Result (n+1)
qnew_global=q4F


!V y C para calcular dt
call VyC

!-----------------------------------------

!-----------------------------------------
END SUBROUTINE solverf4




