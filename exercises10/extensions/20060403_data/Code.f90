! This file was created on December 09, 2004.
! Modified on March 28, 2005.
! This file is part of BCUV project.
! It computes equilibrium of the Mortensen-Pissarides model with agg. uncertainty.
!*************************************************************************************
MODULE PARAMETERS
Implicit NONE
Integer, Parameter :: dd = Selected_real_kind(15,307)

!*********
!Switches:
Integer, Parameter :: one_run=0 ! Set indicator one_run to 1 in order to compute just one point
								! Set indicator one_run to 0 in order to calibrate the model

!***************************
!Fundamentals:
Real(Kind=dd), Parameter :: p = 1
Real(Kind=dd), Parameter :: r = 0.99916
Real(Kind=dd), Parameter :: s = 0.0081
Real(Kind=dd), Parameter :: costK = 0.474
Real(Kind=dd), Parameter :: costW = 0.110
Real(Kind=dd), Parameter :: Smoothing=1600 
Real(Kind=dd), Parameter :: sd = 0.0034
Real(Kind=dd), Parameter :: pers = 0.9895
Integer, Parameter :: ZZ=35

!********
!Targets:
Real(Kind=dd), Parameter :: target_elasticity=0.4525
Real(Kind=dd), Parameter :: target_mean_theta=0.634
Real(Kind=dd), Parameter :: target_arrival_rate=0.139

!*********************
!Technical Parameters:
Integer, Parameter :: Sim_Periods=756000
Integer, Parameter :: Sim_Periods_Dropped=12000
Integer, Parameter :: IOPT=0 
Real(Kind=dd), Parameter :: eps_theta=0.0000004

!*************************
!Defining Other Variables:
Real(Kind=dd) :: chi,barg_power,out_option,mean_theta,arrival_rate,alpha,profit,elasticity,lowest_prod 
Real(Kind=dd), Dimension(ZZ) :: shock,steady_dist,theta, theta_prime, temp
Real(Kind=dd), Dimension(ZZ, ZZ) :: q
Real(Kind=dd), Dimension(:,:), Allocatable :: Cumq
Real(Kind=dd), Dimension (:), Allocatable :: Cum_steady_dist
Integer, Dimension (:), Allocatable :: seed
Integer size, count, economy_call_number, calib_attempt, counter
Integer :: convergence  
Integer, Dimension(1) :: coordinate, store

END MODULE PARAMETERS


MODULE Ziggurat
! Downloaded from http://users.bigpond.net.au/amiller/random.html on
! January 23, 2005 
 
! Marsaglia & Tsang generator for random normals.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001
   IMPLICIT NONE

   PRIVATE
   INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
   REAL(DP), PARAMETER  ::  m1=2147483648.0_DP,   m2=2147483648.0_DP,      &
                            half=0.5_DP
   REAL(DP)             ::  dn=3.442619855899_DP, tn=3.442619855899_DP,    &
                            vn=0.00991256303526217_DP,                     &
                            q,                    de=7.697117470131487_DP, &
                            te=7.697117470131487_DP,                       &
                            ve=0.003949659822581572_DP
   INTEGER,  SAVE       ::  iz, jz, jsr=123456789, kn(0:127),              &
                            ke(0:255), hz
   REAL(DP), SAVE       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
   LOGICAL,  SAVE       ::  initialized=.FALSE.

   PUBLIC  :: zigset, shr3, uni, rnor

CONTAINS

SUBROUTINE zigset( jsrseed )
   INTEGER, INTENT(IN)  :: jsrseed
   INTEGER  :: i

   !  Set the seed
   jsr = jsrseed

   !  Tables for RNOR
   q = vn*EXP(half*dn*dn)
   kn(0) = (dn/q)*m1
   kn(1) = 0
   wn(0) = q/m1
   wn(127) = dn/m1
   fn(0) = 1.0_DP
   fn(127) = EXP( -half*dn*dn )
   DO  i = 126, 1, -1
      dn = SQRT( -2.0_DP * LOG( vn/dn + EXP( -half*dn*dn ) ) )
      kn(i+1) = (dn/tn)*m1
      tn = dn
      fn(i) = EXP(-half*dn*dn)
      wn(i) = dn/m1
   END DO
   initialized = .TRUE.
   RETURN
END SUBROUTINE zigset

!  Generate random 32-bit integers
FUNCTION shr3( ) RESULT( ival )
   INTEGER  ::  ival

   jz = jsr
   jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
   jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
   jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
   ival = jz + jsr
   RETURN
END FUNCTION shr3

!  Generate uniformly distributed random numbers
FUNCTION uni( ) RESULT( fn_val )
   REAL(DP)  ::  fn_val

   fn_val = half + 0.2328306e-9_DP * shr3( )
   RETURN
END FUNCTION uni

!  Generate random normals
FUNCTION rnor( ) RESULT( fn_val )
   REAL(DP)             ::  fn_val

   REAL(DP), PARAMETER  ::  r = 3.442620_DP
   REAL(DP)             ::  x, y

   IF( .NOT. initialized ) CALL zigset( jsr )
   hz = shr3( )
   iz = IAND( hz, 127 )
   IF( ABS( hz ) < kn(iz) ) THEN
      fn_val = hz * wn(iz)
   ELSE
      DO
         IF( iz == 0 ) THEN
            DO
               x = -0.2904764_DP* LOG( uni( ) )
               y = -LOG( uni( ) )
               IF( y+y >= x*x ) EXIT
            END DO
            fn_val = r+x
            IF( hz <= 0 ) fn_val = -fn_val
            RETURN
         END IF
         x = hz * wn(iz)
         IF( fn(iz) + uni( )*(fn(iz-1)-fn(iz)) < EXP(-half*x*x) ) THEN
            fn_val = x
            RETURN
         END IF
         hz = shr3( )
         iz = IAND( hz, 127 )
         IF( ABS( hz ) < kn(iz) ) THEN
            fn_val = hz * wn(iz)
            RETURN
         END IF
      END DO
   END IF
   RETURN
END FUNCTION rnor

END MODULE ziggurat
!**********************************************************

Program full_simplex_calibration
Use PARAMETERS
Use Ziggurat
Implicit NONE
External Economy
 
Integer iter,i,j
Character(len=10) date, time
Integer, Parameter :: ndim=3, mp=4, np=3
Real(Kind=dd), Dimension(mp,np) :: pp, kk
Real(Kind=dd), Dimension(np) :: x
Real(Kind=dd), Dimension(np) :: best_point
Real(Kind=dd), Dimension(mp) :: y

Real(Kind=dd) funk
Real(Kind=dd), Parameter :: ftol=1.e-4

economy_call_number=0
Open(15,file='Output.txt',STATUS='REPLACE',IOSTAT=i)
Call DATE_AND_TIME(Date, Time)
Write(15,'(2X,12A)') 'This file was created on ', Date(1:4),'-',Date(5:6),'-',Date(7:10), &
 'at ', Time(1:2),':',Time(3:4),':',Time(5:7)
Write(15,'(2X,A)') 'Computing equilibrium of the MP model with agg. uncertainty.'

Call SYSTEM_CLOCK(Count)
CALL zigset(count)
Call Random_Seed(size)
Allocate(seed(size))
Call SYSTEM_CLOCK(Count)
seed=count

Call TI(sd,pers,ZZ,shock,q,steady_dist)
Do i=1,ZZ
	Do j=1,ZZ
		If (q(i,j) < 0) q(i,j)=0 
	Enddo
Enddo

Do i=1, ZZ
	If (Sum(q(i,:)) < 1) q(i,ZZ)=1-sum(q(i,1:ZZ-1))
Enddo

WHERE (steady_dist<0) steady_dist=0
If (Sum(steady_dist) < 1) steady_dist(ZZ)=1-sum(steady_dist(1:ZZ-1))

Allocate(Cumq(ZZ,ZZ))
Do i=1,ZZ
	Cumq(i,1)=q(i,1)
	Do j=2,ZZ
		Cumq(i,j)=Cumq(i,j-1)+q(i,j)
	Enddo
Enddo

Do i=1,ZZ
	If (Cumq(i,ZZ)<1) Cumq(i,ZZ)=1
Enddo
WHERE (Cumq>1) Cumq=1

Allocate(Cum_steady_dist(ZZ))
Cum_steady_dist(1)=steady_dist(1)
Do i=2,ZZ
	Cum_steady_dist(i)=Cum_steady_dist(i-1)+steady_dist(i)
Enddo
If (Cum_steady_dist(ZZ)<1) Cum_steady_dist(ZZ)=1
WHERE (Cum_steady_dist>1) Cum_steady_dist=1

theta_prime=shock
lowest_prod=Minval(Shock)
counter=0
calib_attempt=1

If (one_run==1) then
	pp(1,1)=log(0.40715)
	pp(1,2)=log(((0.97012*lowest_prod)/lowest_prod)/(1-(0.97012*lowest_prod)/lowest_prod))
	pp(1,3)=log(0.0524/(1-0.0524))
Else
	pp(1,1)=log(0.40715)
	pp(1,2)=log(((0.97012*lowest_prod)/lowest_prod)/(1-(0.97012*lowest_prod)/lowest_prod))
	pp(1,3)=log(0.0524/(1-0.0524))
Endif

pp(2,1)=log(0.405)
pp(2,2)=log(((0.97*lowest_prod)/lowest_prod)/(1-(0.97*lowest_prod)/lowest_prod))
pp(2,3)=log(0.055/(1-0.055))

pp(3,1)=log(0.39)
pp(3,2)=log(((0.968*lowest_prod)/lowest_prod)/(1-(0.968*lowest_prod)/lowest_prod))
pp(3,3)=log(0.05/(1-0.05))

pp(4,1)=log(0.395)
pp(4,2)=log(((0.969*lowest_prod)/lowest_prod)/(1-(0.969*lowest_prod)/lowest_prod))
pp(4,3)=log(0.048/(1-0.048))

kk=pp

1212 convergence=0

If (calib_attempt>1 .and. counter<calib_attempt) then
	kk(store(1),:)=best_point*.995
	pp=kk
Endif

56 Do j=1,mp
	Do i=1,np
		x(i)=pp(j,i)
	Enddo
	If (one_run==0) then 
		Call func(x,ndim,funk)
		y(j)=funk
	Else 
		convergence=1
		Call func(x,ndim,funk)
		Stop
	Endif
Enddo

If (counter<calib_attempt) then
	counter=calib_attempt
	store=maxloc(y)
Endif

call simplex(pp,y,mp,np,ndim,ftol,iter)

If (iter>=100) then
	Write(15,'(A)') 'Best points after the last 100 iterations'
	Write(15,'(1X,A,11(F20.11, 1X))')  'First point=', exp(pp(1,1)),(exp(pp(1,2))/(1+exp(pp(1,2))))*lowest_prod, exp(pp(1,3))/(1+exp(pp(1,3)))
	Write(15,'(1X,A,11(F20.11, 1X))')  'Second point=', exp(pp(2,1)),(exp(pp(2,2))/(1+exp(pp(2,2))))*lowest_prod, exp(pp(2,3))/(1+exp(pp(2,3)))
	Write(15,'(1X,A,11(F20.11, 1X))')  'Third point=',exp(pp(3,1)),(exp(pp(3,2))/(1+exp(pp(3,2))))*lowest_prod, exp(pp(3,3))/(1+exp(pp(3,3)))
	Write(15,'(1X,A,11(F20.11, 1X))')  'Fourth point=',exp(pp(4,1)),(exp(pp(4,2))/(1+exp(pp(4,2))))*lowest_prod, exp(pp(4,3))/(1+exp(pp(4,3)))
	goto 56
Endif

If (minval(y)>0.05 .and. calib_attempt<mp) then
	Write(15,'(A)') ' '
	Write(15,'(A)') 'Converged, but not too well. Restarting.'
	coordinate=minloc(y)
	Do i=1,np
		best_point(i)=pp(coordinate(1),i)
	Enddo
	calib_attempt=calib_attempt+1
	goto 1212 
Endif

Write(15,'(A)') 'CONVERGED.'
	Write(15,'(1X,A,11(F20.11, 1X))')  'First point=', exp(pp(1,1)),(exp(pp(1,2))/(1+exp(pp(1,2))))*lowest_prod, exp(pp(1,3))/(1+exp(pp(1,3)))
	Write(15,'(1X,A,11(F20.11, 1X))')  'Second point=', exp(pp(2,1)),(exp(pp(2,2))/(1+exp(pp(2,2))))*lowest_prod, exp(pp(2,3))/(1+exp(pp(2,3)))
	Write(15,'(1X,A,11(F20.11, 1X))')  'Third point=',exp(pp(3,1)),(exp(pp(3,2))/(1+exp(pp(3,2))))*lowest_prod, exp(pp(3,3))/(1+exp(pp(3,3)))
	Write(15,'(1X,A,11(F20.11, 1X))')  'Fourth point=',exp(pp(4,1)),(exp(pp(4,2))/(1+exp(pp(4,2))))*lowest_prod, exp(pp(4,3))/(1+exp(pp(4,3)))

convergence=1
If(convergence==1) then
	coordinate=minloc(y)
	Do i=1,np
		x(i)=pp(coordinate(1),i)
	Enddo
	Call func(x,ndim,funk)
Endif

Close(15)

contains       

Subroutine simplex(pp,y,mp,np,ndim,ftol,iter)
Integer, Intent(In) :: ndim,mp,np
Integer, Intent(Out) :: iter
Real(Kind=dd), Intent(In) :: ftol
Real(Kind=dd), Dimension(mp,np), Intent(InOut) :: pp
Real(Kind=dd), Dimension(mp), Intent(InOut) :: y
Real(Kind=dd), Parameter :: TINY=1.e-10
Integer, Parameter :: itmax=50000
Real(Kind=dd) funk, fac
Integer i,ihi,ilo,inhi,j,m,n
Real(Kind=dd) rtol,sum,swap,ysave,ytry,psum(ndim), try_value
iter=0       
1 Do n=1,ndim
      sum=0.0
      Do m=1,ndim+1
            sum=sum+pp(m,n)
      Enddo
         psum(n)=sum
  Enddo
2 ilo=1
If (iter>=100) return 
if (y(1).gt.y(2)) then
   ihi=1
   inhi=2
else
   ihi=2
   inhi=1
end if
Do i=1,ndim+1
   if (y(i).le.y(ilo)) ilo=i
   if (y(i).gt.y(ihi)) then
      inhi=ihi
      ihi=i
   else if (y(i).gt.y(inhi)) then
        if (i.ne.ihi) inhi=i
   end if
Enddo
rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi)) + abs(y(ilo))+TINY)
If (rtol.lt.ftol) then
   swap=y(1)
   y(1)=y(ilo)
   y(ilo)=swap
      Do n=1,ndim
         swap=pp(1,n)
         pp(1,n)=pp(ilo,n)
         pp(ilo,n)=swap
      Enddo
   return
Endif
If ((iter.ge.itmax).and.(iter.gt.0)) then
	Write(15,'(A)') 'ITMAX exceeded'
	Pause 
Endif
iter=iter+2
fac=-1.0
Call try(pp,y,psum,mp,np,ndim,ihi,fac,try_value)
ytry=try_value
If (ytry.le.y(ilo)) then
	fac=2.0
	Call try(pp,y,psum,mp,np,ndim,ihi,fac,try_value)
	ytry=try_value
Elseif (ytry.ge.y(inhi)) then
    ysave=y(ihi)
	fac=0.5
    Call try(pp,y,psum,mp,np,ndim,ihi,fac,try_value)
	ytry=try_value
		If (ytry.ge.ysave) then
		   Do i=1,ndim+1
             If (i.ne.ilo) then
               Do j=1,ndim
                 psum(j)=0.5*(pp(i,j)+pp(ilo,j))
                 pp(i,j)=psum(j)
               Enddo
			   	Call func(psum, ndim, funk)
			    y(j)=funk
             Endif
           Enddo
          iter=iter+ndim
          goto 1
        Endif
Else
   iter=iter-1
Endif
goto 2
End Subroutine simplex  
      
Subroutine try(pp,y,psum,mp,np,ndim,ihi,fac,yy)
      
Integer, Intent(In) :: ihi,mp,ndim,np
Real(Kind=dd) :: funk
Real(Kind=dd), Intent(In) :: fac
Real(Kind=dd), Intent(Out) :: yy
Real(Kind=dd), Dimension(mp,np), Intent(InOut) :: pp
Real(Kind=dd), Dimension(np), Intent(InOut) :: psum
Real(Kind=dd), Dimension(mp), Intent(InOut) :: y
Integer j
Real(Kind=dd) fac1,fac2,ytry
Real(Kind=dd), Dimension(ndim) :: ptry
fac1=(1.-fac)/ndim
fac2=fac1-fac
Do j=1,ndim
  ptry(j)=psum(j)*fac1-pp(ihi,j)*fac2
Enddo
Call func(ptry, ndim, funk)
ytry=funk
If (ytry.lt.y(ihi)) then
    y(ihi)=ytry
		 Do j=1,ndim
           psum(j)=psum(j)-pp(ihi,j)+ptry(j)
           pp(ihi,j)=ptry(j)
         Enddo
Endif
yy=ytry
Return
End Subroutine try

!*****************************************
Subroutine func(x, ndim, xx)
      
Integer, Intent(In) :: ndim
Real(Kind=dd), Dimension(ndim), Intent(In) :: x
Real(Kind=dd), Intent(Out) :: xx
Real(Kind=dd) :: x1, x2, x3
chi=x(1)
chi=exp(chi)
alpha =chi
chi=1
out_option=x(2)
out_option=(exp(out_option)/(1+exp(out_option)))*lowest_prod+0.03
barg_power=x(3) 
barg_power=exp(barg_power)/(1+exp(barg_power))
Call Economy
x1=((elasticity-target_elasticity)/target_elasticity)**2.0 
x2=((mean_theta-target_mean_theta)/target_mean_theta)**2.0 
x3=((arrival_rate-target_arrival_rate)/target_arrival_rate)**2.0 
xx=0
xx=x1+x2+x3
Write(15,'(1X,A,F15.10)')  'The calibration function value is ', xx
Return
End Subroutine func

Subroutine TI(STDE,RHO,ZZ,SZ,PZEE,pi)
Implicit NONE
Real(Kind=dd), Intent(In) :: STDE,RHO
Integer, Intent(In) :: ZZ
Real(Kind=dd), Intent(Out) :: SZ(ZZ),PZEE(ZZ,ZZ),pi(ZZ)
INTEGER i,j,k
Integer, parameter :: draws_number=10000000
Real(Kind=dd)  STDW, step, M, temp1
Real(Kind=dd), Dimension(draws_number) :: normal_draws
M=2.0
STDW=STDE/sqrt(1-RHO*RHO)
SZ(ZZ)=M*STDW
SZ(1)=-SZ(ZZ)
step=2*STDW*M/FLOAT(ZZ-1)
DO i=2,ZZ-1
   SZ(i)=SZ(i-1)+step
Enddo
DO  i = 1, draws_number
	normal_draws(i) = rnor( )
END DO
normal_draws=normal_draws*stde
PZEE=0
Do j=1,ZZ
	temp1=rho*SZ(j)
	Do i=1,draws_number
		If (normal_draws(i)<=SZ(1)+step/2-temp1) PZEE(j,1)=PZEE(j,1)+1
		If (normal_draws(i)> SZ(ZZ)-step/2-temp1) PZEE(j,ZZ)=PZEE(j,ZZ)+1
	Enddo
Enddo
Do j=1,ZZ
	temp1=rho*SZ(j)
	Do i=1,draws_number
		Do k=2, ZZ-1
			If ((normal_draws(i)> SZ(k)-step/2-temp1) .and. (normal_draws(i)<=SZ(k)+step/2-temp1)) &
				 PZEE(j,k)=PZEE(j,k)+1
		Enddo
	Enddo
Enddo
PZEE=PZEE/draws_number
DO i=1,ZZ
   SZ(i)=EXP(SZ(i))
Enddo
pi=1.0/ZZ
Do i=1,400
	pi=Matmul(pi,PZEE)
Enddo
RETURN
END Subroutine TI
End Program full_simplex_calibration

Subroutine ECONOMY
Use PARAMETERS
Implicit NONE
Real(Kind=dd) :: mean_log_wage,mean_log_shock,numerator,denominator,sdev,var,corr
Real(Kind=dd), Dimension (:), Allocatable :: rn,Trend,Deviations,Stats,Stats1,Stats2,Unemp,Productivity,Tightness
Real(Kind=dd), Dimension (:,:), Allocatable :: HP_Temp
Integer, Dimension (:), Allocatable :: Shocks, Shocks_match
Integer :: b,i,j,k,conv_theta

economy_call_number=economy_call_number+1

Write(15,'(A)') ' '
Write(15,'(A)') ' '
Write(15,'(A)') '***************************************'

	If (economy_call_number==1) then
		Write(15,'(A)') '*Economy is called for the first time*'
	Else If (economy_call_number==2) then
		Write(15,'(A)') '*Economy is called for the second time*'
	Else If (economy_call_number==3) then
		Write(15,'(A)') '*Economy is called for the third time*'
	Else
Write(15,'(1X,A,I5,A)') '*Economy is called for the',economy_call_number,'th time'
	Endif
Write(15,'(1X,A,F16.12)')  '* alpha is ', alpha
Write(15,'(1X,A,F16.12)')  '* out_option is ', out_option
Write(15,'(1X,A,F16.12)')  '* barg_power is ', barg_power

Write(15,'(A)') '****************************************'
Write(15,'(A)') ' '

conv_theta=0
b=1
Do While (conv_theta==0 .and. b<=50000)
	temp=0
	Do i=1, ZZ
   		temp(i)=(r*chi)*( (1-barg_power)*(shock(i)-out_option)- &
			(costK*shock(i)+costW*shock(i)**(target_elasticity))*barg_power*theta_prime(i)+ &
			(costK*shock(i)+costW*shock(i)**(target_elasticity))*(1-s)*((1+theta_prime(i)**(alpha))**(1/alpha))/chi )
	Enddo
	theta=0
	Do i=1,ZZ
		Do k=1,ZZ
			theta(i)=theta(i)+temp(k)*q(i,k)/(costK*shock(i)+costW*shock(i)**(target_elasticity))
		Enddo
			theta(i)=(theta(i)**(alpha)-1)**(1/alpha)
	Enddo
		If (Maxval(abs(theta-theta_prime)/theta_prime)<eps_theta) conv_theta=1
		theta_prime=0.95*theta_prime+0.05*theta
	b=b+1
Enddo
Write(15,'(A,I6,1X,A)') 'Theta Function Converged After', b-1, 'Iterations.'

If (convergence==1) then
Allocate(Unemp(Sim_periods))
Allocate(Productivity(Sim_periods))
Allocate(Tightness(Sim_periods))
Allocate(Shocks(Sim_periods))
Allocate(Shocks_match(Sim_periods))
Allocate(rn(Sim_periods))

Call random(seed,size,Sim_periods,rn)
Shocks=1
Do i=2,ZZ
	If (rn(1)>=Cum_steady_dist(i-1).and. rn(1)<Cum_steady_dist(i)) Shocks(1)=i
Enddo
Do i=2,Sim_periods
	Do j=2,ZZ
		If (rn(i)>=Cumq(Shocks(i-1),j-1) .and. rn(i)<Cumq(Shocks(i-1),j)) Shocks(i)=j
	Enddo 
Enddo

Do i=1,Sim_periods
	Tightness(i)=theta(Shocks(i))
Enddo

Unemp(1)=0.05
Do i=2,Sim_periods
	Unemp(i)=s+Unemp(i-1)*(1-s-chi*((Tightness(i-1)**(-alpha)+1)**(-1/alpha)))
Enddo

Do i=1,Sim_periods
	Productivity(i)=shock(Shocks(i))
Enddo


!********************
!Compute statistics *
!********************
Allocate(Stats(Sim_periods))

k=Sim_periods-Sim_Periods_Dropped

Stats=Unemp
Deallocate(Unemp)
Allocate(Unemp(k))
Unemp=Stats(Sim_Periods_Dropped+1:Sim_Periods)

Stats=Tightness
Deallocate(Tightness)
Allocate(Tightness(k))
Tightness=Stats(Sim_Periods_Dropped+1:Sim_Periods)

Stats=Productivity
Deallocate(Productivity)
Allocate(Productivity(k))
Productivity=Stats(Sim_Periods_Dropped+1:Sim_Periods)

Deallocate(Stats)
Allocate(Stats(Sim_periods-Sim_Periods_Dropped))

k=(Sim_periods-Sim_Periods_Dropped)/12

Stats=Unemp
Deallocate(Unemp)
Allocate(Unemp(k))
Do i=1,k
	Unemp(i)=Sum(Stats(i*12-11:i*12))/12
Enddo

Stats=Tightness
Deallocate(Tightness)
Allocate(Tightness(k))
Do i=1,k
	Tightness(i)=Sum(Stats(i*12-11:i*12))/12
Enddo

Stats=Productivity
Deallocate(Productivity)
Allocate(Productivity(k))
Do i=1,k
	Productivity(i)=Sum(Stats(i*12-11:i*12))/12
Enddo


Deallocate(Stats)

k=k-1
Allocate(Stats1(k))
Allocate(Stats2(k))
Allocate(Trend(k))
Allocate(Deviations(k))
Allocate(HP_Temp(k,3))

Write(15,'(A)') ' '
Write(15,'(2X,A)') '*********************************************************'
Write(15,'(2X,A)') '*             HP Filtered LOG Unemployment              *'
Write(15,'(2X,A)') '*********************************************************'
Stats1=Unemp(2:k+1)
Stats2=Unemp(1:k)
Stats1=log(Stats1)
Stats2=log(Stats2)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Call Moments(Stats1,k,sdev)
Write(15,'(2X,A,F14.9)') 'Standard Deviation=', sdev
Write(15,'(2X,A,F14.9)') 'Autocorrelation=', corr

Write(15,'(A)') ' '
Write(15,'(2X,A)') '*********************************************************'
Write(15,'(2X,A)') '*              HP Filtered LOG Vacancies                *'
Write(15,'(2X,A)') '*********************************************************'
Do i=1,k
	Stats1(i)=Tightness(i+1)*Unemp(i+1)
	Stats2(i)=Tightness(i)*Unemp(i)
Enddo
Stats1=log(Stats1)
Stats2=log(Stats2)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Call Moments(Stats1,k,sdev)
Write(15,'(2X,A,F14.9)') 'Standard Deviation=', sdev
Write(15,'(2X,A,F14.9)') 'Autocorrelation=', corr


Write(15,'(A)') ' '
Write(15,'(2X,A)') '*********************************************************'
Write(15,'(2X,A)') '*           HP Filtered LOG  Market  Tightness          *'
Write(15,'(2X,A)') '*********************************************************'
Stats1=Tightness(2:k+1)
Stats2=Tightness(1:k)
Stats1=log(Stats1)
Stats2=log(Stats2)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Call Moments(Stats1,k,sdev)
Write(15,'(2X,A,F14.9)') 'Standard Deviation=', sdev
Write(15,'(2X,A,F14.9)') 'Autocorrelation=', corr

Write(15,'(A)') ' '
Write(15,'(2X,A)') '*********************************************************'
Write(15,'(2X,A)') '*             HP Filtered LOG Productivity              *'
Write(15,'(2X,A)') '*********************************************************'
Stats1=Productivity(2:k+1)
Stats2=Productivity(1:k)
Stats1=log(Stats1)
Stats2=log(Stats2)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Call Moments(Stats1,k,sdev)
Write(15,'(2X,A,F14.9)') 'Standard Deviation=', sdev
Write(15,'(2X,A,F14.9)') 'Autocorrelation=', corr
Write(15,'(A)') ' '



Write(15,'(A)') ' '
Write(15,'(2X,A)') '*********************************************************'
Write(15,'(2X,A)') '*    Correlations (of HP Filtered Loged Variables)      *'
Write(15,'(2X,A)') '*********************************************************'

! Unemployment
Stats1=Unemp(1:k)
Stats1=log(Stats1)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations

Do i=1,k
	Stats2(i)=Tightness(i)*Unemp(i)
Enddo
Stats2=log(Stats2)
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Write(15,'(2X,A,F14.9)') 'Unemployment and Vacancies =', corr

Stats2=Tightness(1:k)
Stats2=log(Stats2)
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Write(15,'(2X,A,F14.9)') 'Unemployment and Mrkt Tightness =', corr

Stats2=Productivity(1:k)
Stats2=log(Stats2)
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Write(15,'(2X,A,F14.9)') 'Unemployment and Productivity =', corr


! Vacancies
Do i=1,k
	Stats1(i)=Tightness(i)*Unemp(i)
Enddo
Stats1=log(Stats1)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations

Stats2=Tightness(1:k)
Stats2=log(Stats2)
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Write(15,'(2X,A,F14.9)') 'Vacancies and Mrkt Tightness =', corr

Stats2=Productivity(1:k)
Stats2=log(Stats2)
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Write(15,'(2X,A,F14.9)') 'Vacancies and Productivity =', corr


! Tightness
Stats1=Tightness(1:k)
Stats1=log(Stats1)
Call HPFILT(Stats1,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats1=Deviations

Stats2=Productivity(1:k)
Stats2=log(Stats2)
Call HPFILT(Stats2,Trend,Deviations,HP_Temp,k,Smoothing,IOPT)
Stats2=Deviations
Call Correlation(Stats1,Stats2,k,corr)
Write(15,'(2X,A,F14.9)') 'Mrkt Tightness and Productivity =', corr


Deallocate(HP_Temp); Deallocate(Trend); Deallocate(Deviations)
Deallocate(Stats1); Deallocate(Stats2)
Deallocate(Unemp);
Deallocate(Tightness); Deallocate(Productivity)
Deallocate(Shocks); Deallocate(rn)

Endif 



!**************************************************
! *******    Compute Calibration Targets   ********
!**************************************************

Write(15,'(A)') ' '
Write(15,'(2X,A)') 'Computing Calibration Targets.'

!**************************************************
!1. Compute Elasticity of Wages wrt Productivity
temp=0	
elasticity=0
mean_log_wage=0
mean_log_shock=0
numerator=0
denominator=0

Do i=1,ZZ
	temp(i)=barg_power*shock(i)+(1-barg_power)*out_option+(costK*shock(i)+costW*shock(i)**(target_elasticity))*barg_power*theta_prime(i)
Enddo
Do i=1,ZZ
	mean_log_wage=mean_log_wage+log(temp(i))*steady_dist(i)
	mean_log_shock=mean_log_shock+log(shock(i))*steady_dist(i)
Enddo
Do i=1,ZZ
	numerator=numerator+steady_dist(i)*(log(temp(i))-mean_log_wage)*(log(shock(i))-mean_log_shock)
	denominator=denominator+steady_dist(i)*(log(shock(i))-mean_log_shock)**2
Enddo
elasticity=numerator/denominator
Write(15,'(1X,A,F16.12)')  '* Elasticity wage/prod is ', elasticity


!**************************************************
!2. Compute Average Vacancy-Unemployment Ratio
mean_theta=0
Do i=1,ZZ
	mean_theta=mean_theta+theta_prime(i)*steady_dist(i)
Enddo
Write(15,'(1X,A,F16.12)')  '* Mean theta is ', mean_theta


!**************************************************
!3. Compute Average Arrival Rate of Offers
arrival_rate=0
Do i=1,ZZ
	arrival_rate=arrival_rate+steady_dist(i)*chi*(1+theta_prime(i)**(-alpha))**(-1/alpha)
Enddo
Write(15,'(1X,A,F16.12)')  '* Arrival rate is ', arrival_rate


return

contains       

Subroutine random(seed,size,Size1,rn)
	Integer, Intent(In) :: Size
	Integer, Intent(In) :: Size1
	Integer, Dimension(size), Intent(InOut) :: Seed
	Real(Kind=dd), Dimension(Size1), Intent(Out) :: rn
	Call Random_SEED(Put=Seed)
	Call Random_Number(rn)
	Call Random_SEED(Get=Seed)
End Subroutine random

Subroutine Moments(Array,n,sdev)
Integer, Intent(In) :: n
Real(Kind=dd), Intent(In) :: Array(n)
Real(Kind=dd), Intent(Out) :: sdev
Integer j
Real(Kind=dd) p,s,ep,ave,var

s=SUM(Array)
ave=s/n
var=0.; ep=0.
do j=1,n
	s=Array(j)-ave
	ep=ep+s
	p=s*s
	var=var+p
Enddo
var=(var-ep**2/n)/(n-1)
sdev=sqrt(var)
End Subroutine Moments

Subroutine Correlation(Array1,Array2,n,corr)
Integer, Intent(In) :: n
Real(Kind=dd), Intent(In) :: Array1(n),Array2(n)
Real(Kind=dd), Intent(Out) :: corr
Integer j
Real(Kind=dd) ave1,ave2,sdev1,sdev2,s,ep

s=SUM(Array1)
ave1=s/n
var=0.; ep=0.
do j=1,n
	s=Array1(j)-ave1
	ep=ep+s
	var=var+s*s
Enddo
var=(var-ep**2/n)/(n-1)
sdev1=sqrt(var)

s=SUM(Array2)
ave2=s/n
var=0.; ep=0.
do j=1,n
	s=Array2(j)-ave2
	ep=ep+s
	var=var+s*s
Enddo
var=(var-ep**2/n)/(n-1)
sdev2=sqrt(var)

corr=0.
do j=1,n
	corr=corr+(Array1(j)-ave1)*(Array2(j)-ave2)
Enddo
corr=corr/(sdev1*sdev2*(n-1))

End Subroutine Correlation

SUBROUTINE HPFILT(Y,T,D,V,N,S,IOPT)
! ----------------------------------------------------------------------
!  SR: hpfilt
!  Kalman smoothing routine for HP filter written by E Prescott.
!   y=data series, d=deviations from trend, t=trend, n=no. obs,
!   s=smoothing parameter (eg, 1600 for std HP).
!   Array v is scratch area and must have dimension at least 3n.
!   If IOPT=1 and n and s are the same as for the previous call,
!   the numbers in v are not recomputed.  This reduces execution
!   time by about 30 percent.  Note that if this option is exercised,
!   v cannot be used for other purposes between calls.
! ----------------------------------------------------------------------

Integer, Intent(In) :: N,IOPT
Real(Kind=dd), Intent(In) :: Y(N),S
Real(Kind=dd), Intent(Out) :: T(N), D(N)
Real(Kind=dd), Intent(InOut) :: V(N,3)

Real(Kind=dd) :: SS,M1,M2,V11,V12,V22,X,Z,B11,B12,B22,DET,E1,E2
Integer :: NN,I,I1,IB

DATA SS,NN/0.D0,0/

!     Compute sequences of covariance matrix for f[x(t),x(t-1) | y(<t)]
      IF(IOPT.NE.1.OR.NN.NE.N.OR.S.NE.SS)  THEN
        SS=S
        NN=N
        V11=1.D0
        V22=1.D0
        V12=0.D0
       DO 5 I=3,N
        X=V11
        Z=V12
        V11=1.D0/S + 4.D0*(X-Z) + V22
        V12=2.D0*X - Z
        V22=X
        DET=V11*V22-V12*V12
        V(I,1)=V22/DET
        V(I,3)=V11/DET
        V(I,2)=-V12/DET
        X=V11+1.D0
        Z=V11
        V11=V11-V11*V11/X
        V22=V22-V12*V12/X
        V12=V12-Z*V12/X
  5    CONTINUE
                                       ENDIF
!
!     this is the forward pass
!
      M1=Y(2)
      M2=Y(1)
      DO 10 I=3,N
        X=M1
        M1=2.0*M1-M2
        M2=X
        T(I-1)= V(I,1)*M1+V(I,2)*M2
        D(I-1)= V(I,2)*M1+V(I,3)*M2
        DET=V(I,1)*V(I,3)-V(I,2)*V(I,2)
        V11=V(I,3)/DET
        V12=-V(I,2)/DET
        Z=(Y(I)-M1)/(V11+1.D0)
        M1=M1+V11*Z
        M2=M2+V12*Z
 10     CONTINUE
      T(N)=M1
      T(N-1)=M2
!
!       this is the backward pass
!
         M1=Y(N-1)
         M2=Y(N)
        DO 15 I=N-2,1,-1
           I1=I+1
           IB=N-I+1
         X=M1
         M1=2.D0*M1 - M2
         M2=X
!
!           combine info for y(.lt.i) with info for y(.ge.i)
!
         IF(I.GT.2)                 THEN
           E1=V(IB,3)*M2 + V(IB,2)*M1 + T(I)
           E2=V(IB,2)*M2 + V(IB,1)*M1 + D(I)
           B11=V(IB,3)+V(I1,1)
           B12=V(IB,2)+V(I1,2)
           B22=V(IB,1)+V(I1,3)
           DET=B11*B22-B12*B12
           T(I)=(-B12*E1+B11*E2)/DET
                                    ENDIF
!
!           end of combining
!
        DET=V(IB,1)*V(IB,3)-V(IB,2)*V(IB,2)
        V11=V(IB,3)/DET
        V12=-V(IB,2)/DET
        Z=(Y(I)-M1)/(V11+1.D0)
         M1=M1+V11*Z
         M2=M2+V12*Z
 15     CONTINUE
       T(1)=M1
       T(2)=M2
        DO 20 I=1,N
 20      D(I)=Y(I)-T(I)
        RETURN
END SUBROUTINE HPFILT


!************************************
! The End of SIMULATION SUBROUTINES *
!************************************

End Subroutine Economy


