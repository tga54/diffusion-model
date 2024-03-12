PROGRAM main
!$use omp_lib
implicit none
integer*4, parameter :: NumT=501, NumE=21, NumR=5000, NumNu=7, NumX=51, NumTheta=99, NumZeta=91 
integer*4 i, j, k, w, Ndiv
real*8 c, pi, m_e, tauc, tau0, We, s, tage, braking, p0, period, pdot, L0, Ls, Q0, count1, MA, eta
real*8 gme(NumE), gme_temp(NumE), gmet(NumE, NumT), Neinj(NumE), Ne0, t(NumT), St(NumT), Qinj(NumT), r(NumR),  DE(NumE), be(NumE), gmetemp, temp, tempt(NumE), tempn(NumE), ndis(NumE,NumR), lambda(NumE, NumT), lambda_temp(NumE), gmax, gmc, dne(NumE), ncr_temp(NumE), cosmu, Ztemp(NumE), D0, delta, y(NumE), Zy(NumE)
real*8 B0, ucmb, ufir, unir, uvis, Tcmb, Tfir, Tnir, Tvis, cfir, cnir, cvis, dE0todE(NumE), rdis, rc, zc, time1, time2, theta(NumTheta)
real*8 dpul, rmax, xx(NumX), xmax, xmin, xfrac(NumX), eg1(NumNu), nu(NumNu), flux_syn(NumNu), flux_ic(NumNu), flux_ic_cmb(NumNu), flux_ic_fir(NumNu), flux_ic_nir(NumNu), flux_ic_vis(NumNu), single_spec(NumNu), intensity_syn(NumNu, NumZeta, NumTheta), intensity_ic(NumNu, NumZeta, NumTheta), intensity_ic_psf(NumNu, NumZeta, NumTheta), psf(NumNu), sigma(NumNu), da, cosda, flux_syn_tot(NumNu), flux_ic_tot(NumNu)
real*8 zeta(NumZeta), dzeta, cosalpha, alpha0, phi, tantheta_cmin, tantheta_cmax, opacity(NumNu), km2a_psf_data(2,5)
character(100) filename
!$ call omp_set_num_threads(40)
c=2.998d10
pi=3.1416d0
m_e=9.1d-28
braking=3d0
period=0.246  !!!
pdot=1.54d-14  !!!
p0=0.05d0       !!!
!Ls=1d45*(2*3.14/period**2*pdot)*(2*3.14/period)*0.5     !!!
Ls=4.1d34
tauc=period/pdot/2
tage=2*tauc/(braking-1)*(1-(p0/period)**(braking-1))
tau0=tage/((p0/period)**(1-braking)-1)
do i=1, NumE
    gme(i)=10**((i-1d0)/10+13)/0.511d6
enddo
filename='J0543'
eta=10**-1.67
MA=0.356
s=1.02d0      !!!
phi=18.92d0/180*pi
D0=1d28  
delta=0.5
B0=5d-6    !!!
gmax=6d8  !0.5muG-1e9, s=2.3, D=2e27

gmc=1d6
Neinj=(gme/gmc)**(-s)*dexp(-gme/gmax)
!WHERE(gme .GT. gmc) Neinj=(gme/gmc)**(-2.4)!*dexp(-gme/gmax)
WHERE(gme .LT. 2D3) Neinj=0


temp=0
do i=1, NumE-1
    temp=temp+(gme(i)*Neinj(i)+gme(i+1)*Neinj(i+1))*0.511d6*(gme(i+1)-gme(i))/2
enddo

!do i=1, NumR  !!!
!    r(i)=((i-1)+0.1)*3.086d18*0.2
!enddo
do i=1, NumR  !!!
    r(i)=((i-1d0)*1d0+1d0)*3.086d18
enddo

do i=1, NumT
    t(i)=10**((i-1d0)/100-5)*tage
enddo


Q0=Ls/temp/1.6d-12
Qinj=eta*Q0*(1+tage/tau0)**((braking+1d0)/(braking-1d0))/(1+(tage-t)/tau0)**((braking+1d0)/(braking-1d0))/MA**4


!print*, tage/3.1536d7, tau0/3.1536d7, tauc/3.1536d7

Tcmb=2.73d0/11600
Tfir=27d0/11600
Tnir=4d2/11600
Tvis=5034d0/11600

ucmb=4.2d-13
ufir=3.4d-13
unir=1.2d-13
uvis=3.4d-13

cfir=ufir/(7.56d-15*(Tfir*11600)**4)
cnir=unir/(7.56d-15*(Tnir*11600)**4)
cvis=uvis/(7.56d-15*(Tvis*11600)**4)



!DO i=1, NumNu
!    nu(i)=10**((i-1d0)/10+13.4)
!ENDDO
!eg1=6.63d-27*nu/1.6d-12

do i=1, NumNu
    eg1(i)=10**((i-1d0)/10+13.4) !!! 25 TeV -- 100 TeV
enddo
nu=eg1*1.6d-12/6.63d-27 !Hz

DO i=1, NumTheta  !!! 10 deg
    theta(i)=((i-1)*1d0+1d0)/10*pi/180
ENDDO
do i=1, NumZeta ! 180 deg (1 pi)
    zeta(i)=(i-1d0)/(NumZeta-1d0)*pi
enddo
dzeta=zeta(2)-zeta(1)

!=========================
!=========================

!实际上是d\gamma/dt
Ndiv=10
lambda_temp=0
gmet(:,1)=gme
do i=2, NumT
    if (t(i) .gt. tage) cycle
    gme_temp=gmet(:,i-1)
    do j=1, Ndiv
        DE=D0*(gme_temp*0.511d6/1d9)**delta 
        lambda_temp=lambda_temp+DE*(t(i)-t(i-1))/Ndiv
        be=4d0/3*6.65d-25*3d10*(gme_temp**2-1)/m_e/c**2*(B0**2/8/pi+ucmb/(1+(gme_temp*Tcmb*2.82/0.511d6)**0.6)**(1.9/0.6)+ufir/(1+(gme_temp*Tfir*2.82/0.511d6)**0.6)**(1.9/0.6)+unir/(1+(gme_temp*Tnir*2.82/0.511d6)**0.6)**(1.9/0.6)+uvis/(1+(gme_temp*Tvis*2.82/0.511d6)**0.6)**(1.9/0.6))
        gme_temp=gme_temp-be*(t(i)-t(i-1))/Ndiv
    enddo
    lambda(:,i)=lambda_temp
    gmet(:,i)=gme_temp
enddo


!$omp parallel private(dE0todE, tempt), REDUCTION(+:tempn)
!$omp do schedule(guided)  
do k=1, NumR
    tempn=0
    do i=2, NumT
        if ( (t(i) .gt. tage)  .or. (r(k) .ge. c*t(i)) ) cycle !
!        if (t(i) .gt. tage) cycle !

        do j=2, NumE-1
            dE0todE(j)=(gme(j+1)-gme(j))/(gmet(j+1,i)-gmet(j,i))/2+(gme(j)-gme(j-1))/(gmet(j,i)-gmet(j-1,i))/2
        enddo
        dE0todE=(gme**2-1)*(B0**2/8/pi+ucmb/(1+(gme*Tcmb*2.82/0.511d6)**0.6)**(1.9/0.6)+ufir/(1+(gme*Tfir*2.82/0.511d6)**0.6)**(1.9/0.6))/( (gmet(:,i)**2-1)*(B0**2/8/pi+ucmb/(1+(gmet(:,i)*Tcmb*2.82/0.511d6)**0.6)**(1.9/0.6)+ufir/(1+(gmet(:,i)*Tfir*2.82/0.511d6)**0.6)**(1.9/0.6)) )
 
        tempt=Qinj(i)/(4*pi*lambda(:,i))**(3d0/2)*DEXP( -r(k)**2/4/lambda(:,i) )*Neinj*dE0todE*(t(i)-t(i-1))

        tempt=dexp(interpol(dlog(gmet(:,i)),dlog(tempt+1d-300), dlog(gme), NumE, NumE))

!  tempt=interpol(gmet(:,i),tempt, gme, NumE, NumE)
        where(gme .gt. gmet(NumE,i)) tempt=0

        tempn=tempn+tempt

    enddo

    ndis(:,k)=tempn
!    print*, k
enddo
!$omp end do
!$omp end parallel

ndis=ndis+1d-300
!open(unit=10, file='ne_'//adjustl(trim(filename))//'.dat')
!write(10,*) ndis
!close(10)
!stop
!====================================================================================
DE=D0*(gme*0.511d6/1d9)**delta


DO i=1, NumX
    xfrac(i)=(i-1d0)/(NumX-1)
ENDDO


dpul=1565d0    !!!
rmax=125d0

r=r/3.086d18

!intensity_syn=0
intensity_ic=0
CALL CPU_TIME(time1)
!$omp parallel private(xmin, xmax, xx, zc, rc, flux_ic_tot, rdis, count1, ncr_temp, dne, flux_ic, single_spec)
!$omp do schedule(guided)
DO i=1, NumTheta     !theta
!DO i=1, 3     !theta
    if (dpul*sin(theta(i)) .gt. rmax) cycle

    xmin=dpul*dcos(theta(i))-sqrt(rmax**2-dpul**2*dsin(theta(i))**2)
    xmax=dpul*dcos(theta(i))+sqrt(rmax**2-dpul**2*dsin(theta(i))**2)

    xx=(xmax-xmin)*xfrac+xmin

    do w=1, NumZeta !zeta
!        flux_syn_tot=0
        flux_ic_tot=0

        do j=1, NumX-1 !xfrac
            zc=dabs( (dpul-xx(j)*dcos(theta(i)))*dcos(phi) + xx(j)*dsin(theta(i))*dcos(zeta(w))*dsin(phi) )
            rc=sqrt( ( (dpul-xx(j)*dcos(theta(i)))*dsin(phi) -  xx(j)*dsin(theta(i))*dcos(zeta(w))*dcos(phi) )**2 + ( xx(j)*dsin(theta(i))*dsin(zeta(w)) )**2 )
            rdis = (rc**2 / MA**4 + zc**2)**0.5

            count1=COUNT(rdis .gt. r)
            if (count1 .eq. 0) count1=1
            if (count1 .eq. NumR) count1=NumR-1

            ncr_temp=exp(dlog( ndis(:,count1+1)/ndis(:,count1))/(r(count1+1)-r(count1))*(rdis-r(count1))+dlog(ndis(:,count1)) )/4/pi

!            flux_syn=syn_spec(nu, gme, m_e, ncr_temp, B0)*(xx(j+1)-xx(j))*3.086d18

            flux_ic=0   
            do k=1, NumE-1 !gme
                single_spec=IC_anal(Tcmb/0.511d6, gme(k), Eg1/0.511d6)+IC_anal(Tfir/0.511d6, gme(k), Eg1/0.511d6)*cfir+IC_anal(Tnir/0.511d6, gme(k), Eg1/0.511d6)*cnir+IC_anal(Tvis/0.511d6, gme(k), Eg1/0.511d6)*cvis
                flux_ic=flux_ic+single_spec*ncr_temp(k)*(gme(k+1)-gme(k))*(xx(j+1)-xx(j))*3.086d18
            enddo
!            flux_syn_tot=flux_syn_tot+flux_syn*nu
            flux_ic_tot=flux_ic_tot+flux_ic*eg1/0.511d6*Eg1*1.6d-12 
        enddo
!        intensity_syn(:,w, i)=flux_syn_tot
        intensity_ic(:,w, i)=flux_ic_tot
    enddo
    CALL CPU_TIME(time2)
!    print*, i, (time2-time1)/60
enddo
!$omp enddo
!$omp end parallel




!open(unit=10, file='syn_'//adjustl(trim(filename))//'.dat')
!write(10, *) intensity_syn
!close(10)

!open(unit=10, file='ic_'//adjustl(trim(filename))//'.dat')
open(unit=10, file=adjustl(trim(filename))//'.dat')
write(10, *) intensity_ic
close(10)


CONTAINS

FUNCTION interpol(x_dat, y_dat, xx, num_dat, num_x)
integer*4 :: num_dat, num_xx, i, num, num_x
real*8 x_dat(num_dat), y_dat(num_dat), xx(num_x), yy(num_x)
real*8 interpol(num_x)

DO i=1, num_x
num=COUNT(xx(i) .GT. x_dat)
IF (num .EQ. 0) yy(i)= (y_dat(2)-y_dat(1))/(x_dat(2)-x_dat(1))*(xx(i)-x_dat(1))+y_dat(1)
IF (num .GT. 0 .AND. num .LT. num_dat)  yy(i)=(y_dat(num+1)-y_dat(num))/(x_dat(num+1)-x_dat(num))*(xx(i)-x_dat(num))+y_dat(num)
IF (num .EQ. num_dat) yy(i)=(y_dat(num_dat)-y_dat(num_dat-1))/(x_dat(num_dat)-x_dat(num_dat-1))*(xx(i)-x_dat(num_dat))+y_dat(num_dat)

ENDDO

interpol=yy

RETURN

END FUNCTION


FUNCTION IC_anal(Tem, Ee, Eg1)
real*8 pi, r0, h, hbar, m_e, c, alpha3, beta3, a3, b3, c3, alpha4, beta4, a4, b4, c4
real*8 Tem, Ee, Eg1(NumNu), z(NumNu), t, x0(NumNu), gg3(NumNu), gg4(NumNu), G3(NumNu), G4(NumNu), IC_anal(NumNu)

pi=3.1416d0
r0=2.81d-13
h=6.63d-27
hbar=h/2/pi
m_e=9.1d-28
c=2.998d10

alpha3=0.606
beta3=1.481
a3=0.443
b3=0.54
c3=0.319

alpha4=0.461
beta4=1.457
a4=0.726
b4=0.382
c4=6.62


z=Eg1/Ee
t=4*Ee*Tem
x0=z/(1-z)/t


gg3=1/( 1+a3*x0**alpha3/(1+b3*x0**beta3) )
gg4=1/( 1+a4*x0**alpha4/(1+b4*x0**beta4) )


G3=pi**2/6d0*(1+c3*x0)/(1+pi**2*c3/6*x0)*dexp(-x0)*gg3
G4=pi**2/6d0*(1+c4*x0)/(1+pi**2*c4/6*x0)*dexp(-x0)*gg4
IC_anal=2*r0**2*m_e**3*c**4*Tem**2/pi/hbar**3/Ee**2*( z**2/2d0/(1-z)*G3+G4 )

WHERE(z .ge. 1) IC_anal=0

RETURN
END FUNCTION


END PROGRAM
