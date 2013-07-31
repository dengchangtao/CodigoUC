MODULE global_variables

!Global Variables MGP

!Variables de Estado
real (kind=8), dimension(:,:),save,allocatable	:: V_global, C_global, VC, S1_global, S2_global

!Metrics
real (kind=8), dimension(:,:,:),save,allocatable	::xi_global, eta_global, &
					  qnew_global, qold_global, &
					  qreal_global


!q real es la variable dimensionalizada					  
real (kind=8),dimension(:,:),save,allocatable	::xc,yc,xe,ye, aj_global,MCoef


!Input control
real (kind=8)	::CFL,tfinal,L,H,U, t, iteration, dt, FR2, g,dxi,deta, hmin, kappa, it, Pvol, vol0, treal, Coef
integer,dimension(4) :: CB
integer	::caso,batiopt, inputopt, Nbx, Nby, dit, mmopt, rk, outopt, fopt, fM, Cf

END MODULE global_variables

MODULE coords
!Coordenadas Curvilíneas
real (kind=8),dimension(:), save, allocatable :: coordxi, coordeta
real (kind=8),dimension(:),save,allocatable :: angulo1,angulo2,angulo3,angulo4

END MODULE coords


MODULE geometries
  real (kind=8),dimension(:,:),save,allocatable	::x_global,y_global,z_global
END MODULE geometries

MODULE senales
  integer:: GA1,GA2,GA3,GA4,Nsenal1,Nsenal2,Nsenal3,Nsenal4, IO1,IO2,IO3,IO4
  real (kind=8):: h01,h02,h03,h04
  real (kind=8),dimension(:,:),save,allocatable:: etaL1,qs1,us1,hs1, &
						  etaR2,qs2,us2,hs2,&
						  etaL3,qs3,us3,hs3, &
						  etaR4,qs4,us4,hs4, &
						  qA1,qA2,qA3,qA4, &
						  qsx1,qsy1, qsx2, qsy2, &
						  qsx3,qsy3, qsx4, qsy4, &
						  etaL9, hs9, &
						  etas1,etas2,etas3,etas4
						  
  real (kind=8),dimension(:),save,allocatable:: zA1,zA2,zA3,zA4,timeS1,timeS2,timeS3,timeS4,timeS9
  
 
END MODULE senales


MODULE time0
!Modulo que guarda las condiciones iniciales
real (kind=8), dimension(:,:,:),save,allocatable:: q0_global
real (kind=8):: epxR0, epyR0
END MODULE time0

!-----------------------
MODULE Jacobianos
!Jacobianos en interfaces de la celda
real (kind=8), dimension(:,:),save,allocatable:: Jac_global_xi, Jac_global_eta, Jac_global

END MODULE Jacobianos
!---------------------------
! MODULE Pich
!  real (kind=8),dimension(:),save,allocatable:: AreaPich
!  real (kind=8),dimension(:),save,allocatable:: VelPich
! END MODULE Pich

