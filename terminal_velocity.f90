!##############################################################
!###
!###   	Tool to calculate the free fall terminal velocity
!###    of an isolated particle    
!### 	using Schiller Naumann
!###                                     F. Audard 19/07/2016
!##############################################################
Program terminal_velocity
!
!## Variables 
  integer :: i,itemax,compo_liq,compo_gaz,Test_autre_x
  real*8 :: diam,rayon,Vol,A,gz, ro_liq, ro_gaz 	
  real*8 :: mu_gaz,pi,nu_gaz,increment,Web
  real*8 :: Cd,x_1,x_2,x_0,xmin,Rep,Taux_p
  logical :: arret_boucle

!## Parameters 
  pi=4.D0*atan(1.D0)
  diam  =  0.0005*2.
  rayon = diam/2.
  gz =  9.81
!
  compo_gaz = 9
  compo_liq = 9 

! test in report: 
! Evolution of the free-fall water drop in air 
! from Belkelfa (2011).
!  diam  =  454.E-6  ! or !  364.E-06 ! or ! 590E-6 ! or 615E-6
!  compo_gaz = 8
!  compo_liq = 1 !9
!--------
! If we want to modify the gravity: 
! write(*,*)'gz'
! read(*,*)gz
!
  !______________________________________________________________
  ! Liquid composition:
  ! 1 -> water
  ! 2 -> heptane 300 k 
  ! 3 -> heptane 350 k
  ! 5 -> décane 320 k
  ! 99 -> test
  ! 9 -> sand particle 
  !______________________________________________________________
  if (compo_liq.eq.1) then 
     ro_liq = 1000.
  elseif (compo_liq.eq.2) then 
     ro_liq = 680.06
  elseif (compo_liq.eq.3) then 
     ro_liq = 629.11
  elseif (compo_liq.eq.4) then 
     ro_liq =  677.7  
  elseif (compo_liq.eq.5) then 
     ro_liq =  709.10
  elseif (compo_liq.eq.99) then 
     ro_liq =  677.70
  elseif (compo_liq.eq.6) then 
     ro_liq =  648
  elseif (compo_liq.eq.7) then 
     ro_liq =  1000.
  elseif (compo_liq.eq.8) then 
     ro_liq =  4600.
  elseif (compo_liq.eq.9) then 
     ro_liq =  1500.
  endif

  !______________________________________________________________
  ! Gaz composition: 
  ! 1 -> air
  ! 2 -> diazote (466k ) 
  ! 3 -> diazote (766k ) 
  ! 5 -> air (766k ) 
  !______________________________________________________________
  if (compo_gaz.eq.1) then 
     ro_gaz = 1.226	
     mu_gaz = 0.0000178  
  elseif (compo_gaz.eq.2) then 
     ro_gaz = 0.7326
     mu_gaz = 2.47E-5
  elseif (compo_gaz.eq.3) then 
     ro_gaz = 0.4457
     mu_gaz = 3.313E-5
    elseif (compo_gaz.eq.4) then 
     ro_gaz =  0.460404
     mu_gaz = 3.31E-005  
    elseif (compo_gaz.eq.5) then 
     ro_gaz =  0.351516
     mu_gaz = 4.21E-005 
    elseif (compo_gaz.eq.6) then 
     ro_gaz = 0.775
     mu_gaz = 2.5E-005
    elseif (compo_gaz.eq.7) then 
     ro_gaz = 0.77
     mu_gaz = 2.5E-005
    elseif (compo_gaz.eq.8) then 
     ro_gaz = 1.2
     mu_gaz = 1.85e-05
    elseif (compo_gaz.eq.9) then 
     ro_gaz = 1000.
     mu_gaz = 5.e-03

  endif
  nu_gaz = mu_gaz / ro_gaz
  !______________________________________________________________
  ! Terminal velocity from the trajectory equation using  
  ! Newton-Raphson method
  !______________________________________________________________
  !---- Correlation used: Schiller-Neumann
  !
  !---- Galilée number: 
  Ga = (diam**3.)*gz*ro_gaz* (ro_liq-ro_gaz)/(mu_gaz**2.)
  A = 4./3. * Ga * nu_gaz /diam 
  B  = 24.*0.15*(diam**0.687)/(nu_gaz**0.687) 
Print*,'Galile number, B coeff', Ga , B
  !  Newton-Raphson
  x_0 = 0.1!1.0
  write(*,*) 'x_0'
  read(*,*) x_0  
!read(*,*) x
  itemax = 100000.
  xmin = 1.E-0010
  x_2 = x_0
  x_1 = 1. 
  i = 1  
  arret_boucle = .true.
  increment = 0.
  Test_autre_x = 0
  x_0 = 0.001
  f = 2.
 !   val1 = sign (abs(x_0)**0.313, x_0)
! A = -4./3. *gz * nu_gaz
! B = 3.6*(diam**0.687)/(nu_gaz**0.687) 

14  do while (i.lt.itemax.and.f.gt.xmin)
     !

    f =  A/x_0 - B*(x_0**0.687)-24.
    g = -A/(x_0*x_0) - B*0.687/(x_0**0.313)
!   f = (4./3.)*Ga/(x_0**2.) - (3.6**(x_0**0.687)  - 24)/(x_0)
!g = -(8./3.)*Ga/(x_0**3.)  - 24*(1+0.15*(x_0**0.687))/(x_0**2.)+&
!&        2.47320/(x_0**1.313)

   !  
     x_2 = x_0 - f/g    
     x_1 = abs(x_0-x_2)
     x_0 = x_2
     i = i+1
     Rep = (x_2*diam)/nu_gaz
     Cd  = (24./Rep)*(1.+0.15*(Rep)**(0.687))

         Print *,'i,x,x_1',i,f,g, x_0,x_1
         Print *,'Cote 1',Cd*(pi*diam**2. / 8. )*ro_gaz*x_2**2. 
         Print *,'Cote 2',gz*(4./3.)*(diam/2.)**3.*(ro_liq-ro_gaz)
         if (x_1.le.xmin) then
         arret_boucle = .false.
         endif
  end do 	
  !
  Rep = (x_2*diam)/nu_gaz
  if (Rep.lt.400) then
     Cd  = (24./Rep)*(1.+0.15*(Rep)**(0.687))
  else
     Cd = 0.044
  endif
  Taux_p = ((ro_liq/ro_gaz)+0.5)*4*diam/(3.*Cd*x_2)
  !
  Web = (x_2**2.)*ro_gaz*diam/0.021
  if (x_1.le.xmin) then
     Print *,B,Ga
     Print *, 'This values are an estimation without taking account mass transfer'
     Print *, 'Terminal velocity is:',x_2
     Print *, 'Particulate Reynolds number:',Rep
     Print *, 'Drag coefficient  with Schiller & Naumann :',Cd
     Print *, 'Time of the relaxation Temps:',Taux_p
     Print *, 'Weber number:',Web
  else 
  Print*, Test_autre_x
  Test_autre_x = Test_autre_x + 1
    if (Test_autre_x.eq.1) then
	x_0 = 0.1
         goto 14
    elseif (Test_autre_x.eq.2) then 
       x_0 = 1.
       goto 14
    elseif (Test_autre_x.eq.3) then 
       x_0 = 10.
       goto 14
    else
     Print *, 'not converge'
     endif
  endif
	
	do i =1,1000 

    f =  A/(i/1000.) - B*((i/1000.)**0.687)-24.
 !   Print*,i,f
    write(99,*) (i/1000.),f 


	end do 


end Program terminal_velocity
