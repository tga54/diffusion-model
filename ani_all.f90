PROGRAM main
!$use_omp_lib
implicit none

integer*4, parameter :: NumE=61, NumT=40000, NumR=181, NumZ=51
integer*4 i,j,k,l,w, tt,count1, count2, istart, jstart
real*8 c, pi, kb, m_e, m_p,e, sigmaT, braking, spin_index,p0,period,pdot, time1, time2,h, sigmaTh,Di
real*8 r(NumR),z(NumZ),Qe(NumE), NinjE(NumE), gam(NumE), pcr(NumE), dr(NumR-1), dz(NumZ-1), lnp(NumE), dlnp, dt, r0, MA, delta
real*8 pcr1(151),gam1(151),Ninj1(151)
real*8 DEpara(NumE), DEperp(NumE), ncrE(NumR, NumZ, NumE), ncrr(NumR, NumZ, NumE), ncrz(NumR, NumZ, NumE), ncr1(NumR, NumZ, NumE), ncr2(NumR, NumZ, NumE), bp(NumE)
real*8 aaz(NumZ,NumE), bbz(NumZ,NumE), ccz(NumZ,NumE), ddz(NumZ,NumE), eez(NumZ,NumE), ffz(NumZ,NumE)
real*8 aar(NumR,NumE), bbr(NumR,NumE), ccr(NumR,NumE), ddr(NumR,NumE), eer(NumR,NumE), ffr(NumR,NumE), p, temp, Le, dpcr(NumE-1)
real*8 Tstar, Tir, Tcmb, Tuv, B0, B1, ustar, uir, ucmb, uuv, Qt(NumT), tauc, t0, t, tage,  timei, timef
real*8 ndis0(NumR, NumZ,NumE)   !NumR, NumZ, NumE
real*8 Ek(NumE), rdis(NumR), zdis(43), rt, zt, ncr(NumR, 43, NumE), ncr0(NumR, 43, NumE), ncr_temp1(43, NumE), ncr_temp2(NumE)
real*8 dpul, rmax, xx(201),xmax, xmin, xfrac(201), slope1, slope2
real*8 TTstar, TTir, TTcmb, TTuv, zeta(120), phi, theta(101), sinalpha, cosalpha, dzeta, Eg1(18), dne(NumE), cir, cstar, cuv, falpha_data(100,3), falpha(100), cosdalpha2(100), cosdelta
real*8 flux_ic(18), intensity_ic(18, 120, 101), flux_ic_tot(18), single_spec(18)
character(60) filename, filename1, ma_name, phi_name, d_name
COMMON /constant/ pi, c, e, h, sigmaTh
integer :: omp_get_thread_num
!$ call omp_set_num_threads(1)

c=2.998d10
pi=3.1416d0
kb=1.38d-16
m_e=9.1e-28
sigmaT=6.65d-25
e=4.8d-10

B0=3d-6
MA=0.25
Di=1e28 !1GeV
p=1.6
delta=0.5

Tcmb=2.73d0
Tir=30d0
Tstar=5d3
Tuv=2d4

ucmb=4.2d-13
uir=4.8d-13
ustar=4.8d-13
uuv=1.6d-13*0



do i=1, 40
r(i)=(i-1)*3.086d18*0.5
enddo

do i=41,120
!r(i)=20*10**( (i-41d0)/100d0)*3.086d18  !最大277 pc
r(i)=20*3.086d18+(i-41)*3.086d18  !20-99 pc
enddo

!do i=81,120
!r(i)=20*10**( (i-41d0)/100d0)*3.086d18  !最大277 pc
!r(i)=60*3.086d18+(i-81)*3.086d18*2  !20-99 pc
!enddo

do i=121, NumR
!r(i)=90*10**( (i-101d0)/110d0)*3.086d18  !最大300 pc
r(i)=100*3.086d18+(i-121)*3*3.086d18
enddo


do i=1, NumR-1
dr(i)=r(i+1)-r(i)
enddo

do i=1, 20
z(i)=(i-1)*3.086d18
enddo
do i=21, NumZ
z(i)=20*10**((i-21d0)/18d0)*3.086d18  !1000pc
enddo

do i=1, NumZ-1
dz(i)=z(i+1)-z(i)
enddo
   !original:dz=2pc,dr=0.3pc,dt=40yr

do i=1, 151
pcr1(i)=10**((i-1)/30d0+5)*1d6*1.6d-12/c           ! momentum in unit of erg*s/cm, eneryg from 1e6eV - 1e18eV,     1/1000能量精度+1yr时间精度->不稳定, 1/100+10yr -> 稳定
enddo
gam1=sqrt((pcr1)**2*c**2+m_e**2*c**4)/m_e/c**2


do i=1, NumE
pcr(i)=10**((i-1)/20d0+7)*1d6*1.6d-12/c           ! momentum in unit of erg*s/cm, eneryg from 1e6eV - 1e18eV,     1/1000能量精度+1yr时间精度->不稳定, 1/100+10yr -> 稳定
enddo
do i=1, NumE-1
dpcr(i)=pcr(i+1)-pcr(i)
enddo
gam=sqrt((pcr)**2*c**2+m_e**2*c**4)/m_e/c**2


!r0=50
!B1=3d-6
!DO i=1, NumR
! IF (r(i) .LT. r0*3.086D18) THEN
!   DE(i,:)=4.5d27*(pcr*c/1.6d-12/1d14)**(1d0/3)      ! diffusion coefficient, normalized at 1GeV
!   bp(i,:)=4d0/3*sigmaT*(gam**2-1)*(B1**2/8/pi+ucmb/(1+4*gam*3*kb*Tcmb/m_e/c**2)**1.5+uir/(1+4*gam*3*kb*Tir/m_e/c**2)**1.5+ustar/(1+4*gam*3*kb*Tstar/m_e/c**2)**1.5+uuv/(1+4*gam*3*kb*Tuv/m_e/c**2)**1.5)
! ELSE
!   DE(i,:)=3.86d28*(pcr*c/1.6d-12/1d9)**(1d0/3) 
!   bp(i,:)=4d0/3*sigmaT*(gam**2-1)*(B0**2/8/pi+ucmb/(1+4*gam*3*kb*Tcmb/m_e/c**2)**1.5+uir/(1+4*gam*3*kb*Tir/m_e/c**2)**1.5+ustar/(1+4*gam*3*kb*Tstar/m_e/c**2)**1.5+uuv/(1+4*gam*3*kb*Tuv/m_e/c**2)**1.5)
!  ENDIF
!ENDDO

DEpara=Di*(pcr*c/1.6d-12/1d9)**(delta)
DEperp=DEpara*MA**4
bp=4d0/3*sigmaT*(gam**2-1)*(B0**2/8/pi+ucmb/(1+(gam*2.82*kb*Tcmb/m_e/c**2)**0.6)**(1.9/0.6)+uir/(1+(gam*2.82*kb*Tir/m_e/c**2)**0.6)**(1.9/0.6)+ustar/(1+(gam*2.82*kb*Tstar/m_e/c**2)**0.6)**(1.9/0.6)+uuv/(1+(gam*2.82*kb*Tuv/m_e/c**2)**0.6)**(1.9/0.6))

lnp=DLOG(pcr)
dlnp=lnp(2)-lnp(1)

dt=3.1536d7*4
t=dt


braking=3d0
spin_index=(braking+1)/(braking-1)
p0=50d-3
period=0.237
pdot=1.1d-14                  !0.1, 1d-14; 0.3, 3d-15
tauc=period/pdot/2
tage=2*tauc/(braking-1)*(1-(p0/period)**(braking-1))               !puslar年龄
t0=2*tauc/(braking-1)*(p0/period)**(braking-1)             !pulsar spindown timescale

Le=3.25d34   !I=1d45 * Omega * Omegadot * 50%  current spin down luminosity

Ninj1=gam1**(-p)*exp(-gam1/4d8) 
where(gam1*0.511d6 .lt. 1d11) Ninj1=0

Ninje=gam**(-p)*exp(-gam/4d8)    ! 1.4 4d8  , 1.6 4d8
where(gam*0.511d6 .lt. 1d11) Ninje=0

temp=0
DO i=1, 151-1
temp=temp+( gam1(i)*m_e*c**2*Ninj1(i)+gam1(i+1)*m_e*c**2*Ninj1(i+1) )*( gam1(i+1)-gam1(i) )/2d0
ENDDO


Qe=Le/temp*Ninje

print*, tage/3.1536d7

ncr1=0
ncr1(1,1,:)=Qe/6d0/pi/dr(1)**2/dz(1)
ncr1(1,2,:)=Qe/6d0/pi/dr(1)**2/dz(1)
ncr1(2,1,:)=Qe/6d0/pi/dr(1)**2/dz(1)
ncr1(2,2,:)=Qe/6d0/pi/dr(1)**2/dz(1)
ncr1(:,:,NumE)=0



tt=1     !时间，用于计数，把感兴趣的时间点的数据记录下来
ndis0=0
ncr2=0
ncrr=0
ncrz=0
ncrE=0

DO l=1, NumT   !时间
if (t .GT. 2d4*3.1536d7) cycle
!call cpu_time(timei)
!$omp parallel 
!$omp do schedule(guided)
   DO j=1, NumZ-1
    DO k=NumE-1, 1, -1

     IF (k .EQ. NumE-1) THEN
      ncrE(:,j,k)=( (4*bp(k+1)*ncr1(:,j,k+1)-3*bp(k)*ncr1(:,j,k)+4*bp(k+1)*ncrE(:,j,k+1))/(4*dpcr(k))*(dt/2)+ncr1(:,j,k) )/( 1+3*bp(k)*(dt/2)/(4*dpcr(k)) )
     ELSE
      ncrE(:,j,k)=( (-bp(k+2)*ncr1(:,j,k+2)+4*bp(k+1)*ncr1(:,j,k+1)-3*bp(k)*ncr1(:,j,k)-bp(k+2)*ncrE(:,j,k+2)+4*bp(k+1)*ncrE(:,j,k+1))/(4*dpcr(k))*(dt/2)+ncr1(:,j,k) )/( 1+3*bp(k)*(dt/2)/(4*dpcr(k)) )
     ENDIF
    ENDDO

   ENDDO
!$omp enddo
!$omp end parallel
!call cpu_time(timef)
!print*, 'finish E', timef-timei
  ncrE(:,:,NumE)=0
  ncrE(NumR,:,:)=0
  ncrE(:,NumZ,:)=0


!$omp parallel private(aar,bbr,ccr,ddr,eer,ffr)
!$omp do schedule(guided)
    
   DO j=1, NumZ-1   !l+1/2

    DO i=1, NumR-1

    IF (i .EQ. 1) THEN
     aar(i,:)=0
     bbr(i,:)=1+DEperp*dt/(dr(1)**2)
     ccr(i,:)=-DEperp*dt/(dr(1)**2)
     IF (j .EQ. 1) ddr(i,:)=DEpara*dt/(dz(j)**2)*( ncrE(i,j+1,:) - ncrE(i,j,:) ) + ncrE(i,j,:)
     IF (j .NE. 1) ddr(i,:)=DEpara*dt*( ncrE(i,j+1,:)/dz(j)/(dz(j)+dz(j-1))- ncrE(i,j,:)/dz(j-1)/dz(j)+ ncrE(i,j-1,:)/dz(j-1)/(dz(j)+dz(j-1))) + ncrE(i,j,:)
     eer(i,:)=ccr(i,:)/bbr(i,:)
     ffr(i,:)=ddr(i,:)/bbr(i,:)
    ELSE
     aar(i,:)=DEperp*dt/(2*(dr(i)+dr(i-1)))*(1d0/r(i)-2/dr(i-1))         ! r(i+1) -r(i-1)= dr(i)+dr(i-1), dr(i)=r(i+1)-r(i) 
     bbr(i,:)=1+DEperp*dt/(2*(dr(i)+dr(i-1)))*(2/dr(i-1)+2/dr(i))   !
     ccr(i,:)=-DEperp*dt/(2*(dr(i)+dr(i-1)))*(1d0/r(i)+2/dr(i))  !
     IF (j .EQ. 1) ddr(i,:)=DEpara*dt/(dz(j)**2)*( ncrE(i,j+1,:) - ncrE(i,j,:) ) + ncrE(i,j,:)  ! N_-1=N_1
     IF (j .NE. 1) ddr(i,:)=DEpara*dt*( ncrE(i,j+1,:)/dz(j)/(dz(j)+dz(j-1))- ncrE(i,j,:)/dz(j-1)/dz(j)+ ncrE(i,j-1,:)/dz(j-1)/(dz(j)+dz(j-1))) + ncrE(i,j,:)  !
     eer(i,:)=ccr(i,:)/(bbr(i,:)-aar(i,:)*eer(i-1,:))
     ffr(i,:)=(ddr(i,:)-aar(i,:)*ffr(i-1,:))/(bbr(i,:)-aar(i,:)*eer(i-1,:))
    ENDIF
    ENDDO
 

    DO i=NumR-1, 1, -1
    IF (i .EQ. NumR-1) THEN
     ncrr(i,j,:)=ffr(i,:)
    ELSE
     ncrr(i,j,:)=ffr(i,:)-eer(i,:)*ncrr(i+1,j,:)
    ENDIF
    ENDDO

   ENDDO
!$omp enddo
!$omp end parallel  
 

 ncrr(:,:,NumE)=0
 ncrr(NumR,:,:)=0
 ncrr(:,NumZ,:)=0


!$omp parallel private(aaz,bbz,ccz,ddz,eez,ffz)
!$omp do schedule(guided)
    
   DO i=1, NumR-1  !l+1

    DO j=1, NumZ-1   !空间

    IF (j .EQ. 1) THEN
     aaz(j,:)=0
     bbz(j,:)=( 1+DEpara*dt/dz(1)**2 ) 
     ccz(j,:)=-DEpara*dt/(dz(1)**2)  
     IF (i .EQ. 1) ddz(j,:)=DEperp*dt/(dr(1)**2)*( ncrr(i+1,j,:) - ncrr(i,j,:) ) + ncrr(i,j,:) 
     IF (i .NE. 1) ddz(j,:)=DEperp*dt/(2*r(i))*( ncrr(i+1,j,:) - ncrr(i-1,j,:))/(dr(i)+dr(i-1)) + DEperp*dt*( ncrr(i+1,j,:)/(dr(i)+dr(i-1))/dr(i)-ncrr(i,j,:)/dr(i)/dr(i-1)+ncrr(i-1,j,:)/(dr(i-1)+dr(i))/dr(i-1))+ ncrr(i,j,:)
     eez(j,:)=ccz(j,:)/bbz(j,:)
     ffz(j,:)=ddz(j,:)/bbz(j,:)
    ELSE
     aaz(j,:)=-DEpara*dt/(dz(j-1)+dz(j))/dz(j-1) !
     bbz(j,:)=1+DEpara*dt/dz(j-1)/dz(j)   !
     ccz(j,:)=-DEpara*dt/( dz(j-1)+dz(j))/dz(j)   !
     IF (i .EQ. 1) ddz(j,:)=DEperp*dt/(dr(1)**2)*( ncrr(i+1,j,:) - ncrr(i,j,:) ) + ncrr(i,j,:)  !
     IF (i .NE. 1) ddz(j,:)=DEperp*dt/(2*r(i))*( ncrr(i+1,j,:) - ncrr(i-1,j,:))/(dr(i)+dr(i-1)) + DEperp*dt*( ncrr(i+1,j,:)/(dr(i)+dr(i-1))/dr(i)-ncrr(i,j,:)/dr(i)/dr(i-1)+ncrr(i-1,j,:)/(dr(i-1)+dr(i))/dr(i-1))+ ncrr(i,j,:)  !
     eez(j,:)=ccz(j,:)/(bbz(j,:)-aaz(j,:)*eez(j-1,:))
     ffz(j,:)=(ddz(j,:)-aaz(j,:)*ffz(j-1,:))/(bbz(j,:)-aaz(j,:)*eez(j-1,:))
    ENDIF
    ENDDO


    DO j=NumZ-1, 1, -1
    IF (j .EQ. NumZ-1) THEN
     ncrz(i,j,:)=ffz(j,:)
    ELSE
     ncrz(i,j,:)=ffz(j,:)-eez(j,:)*ncrz(i,j+1,:)
    ENDIF
    ENDDO

   ENDDO
!$omp enddo
!$omp end parallel   

 ncrz(:,:,NumE)=0
 ncrz(NumR,:,:)=0
 ncrz(:,NumZ,:)=0



!$omp parallel 
!$omp do schedule(guided)
   DO j=1, NumZ-1
    DO k=NumE-1, 1, -1  
     IF (k .EQ. NumE-1) THEN
      ncr2(:,j,k)=( (4*bp(k+1)*ncrz(:,j,k+1)-3*bp(k)*ncrz(:,j,k)+4*bp(k+1)*ncr2(:,j,k+1))/(4*dpcr(k))*(dt/2)+ncrz(:,j,k) )/( 1+3*bp(k)*(dt/2)/(4*dpcr(k)) )
     ELSE
      ncr2(:,j,k)=( (-bp(k+2)*ncrz(:,j,k+2)+4*bp(k+1)*ncrz(:,j,k+1)-3*bp(k)*ncrz(:,j,k)-bp(k+2)*ncr2(:,j,k+2)+4*bp(k+1)*ncr2(:,j,k+1))/(4*dpcr(k))*(dt/2)+ncrz(:,j,k))/( 1+3*bp(k)*(dt/2)/(4*dpcr(k)) )
     ENDIF
    ENDDO

  ENDDO
!$omp enddo
!$omp end parallel   


ncr2(:, :, NumE)=0
ncr2(NumR,:,:)=0
ncr2(:,NumZ,:)=0

ncr1=ncr2
  
!IF (t .LE. tage) THEN
! DO k=1, 151
   ndis0=ndis0+(1+tage/t0)**spin_index/(1+(tage-t)/t0)**spin_index*ncr2(:,:, :)*dt
! ENDDO
!ENDIF
ncr2=0
!print*, ncr2(100,100,30)
IF (MOD(l, NumT/1000) .EQ. 0) THEN
print*, tt, t/3.1536d7,'yr' , ndis0(1,1,1)
tt=tt+1
ENDIF
t=t+dt
ENDDO

!OPEN(UNIT=10, FILE='ncr0_ani_3muG_1.6_1d28_MA0.25_p050_n3_200TeV_dt4yr_NumE61_NumR181_NumZ51log_1wyr_all.dat')
!WRITE(10,*) ndis0(:,1:43,:)
!CLOSE(10)

!辐射






TTcmb=2.73d0/11600
TTir=30d0/11600
TTstar=5d3/11600
TTuv=2d4/11600

cir=uir/(7.56d-15*(TTir*11600)**4)
cstar=ustar/(7.56d-15*(TTstar*11600)**4)
cuv=0

Ek=pcr*c
DO i=1, 201
xfrac(i)=(i-1d0)/200
ENDDO

zdis=z(1:43)/3.086d18
rdis=r/3.086d18


DO i=1, 120
zeta(i)=(i-1d0)*3d0/180*pi  !zeta就是天球上的phi
ENDDO
dzeta=zeta(2)-zeta(1)

DO i=1, 18
eg1(i)=10**( (i-1)/10d0+13d0)  !10 TeV-500TeV !出射能量
ENDDO

DO i=1, 100
cosdalpha2(i)=(i-1)*0.01d0+0.005d0
ENDDO

phi=10d0
dpul=250
cosdelta=dsqrt(1/(1+MA**2))   !用Lazarian Pogosyan 2012的近似

write(ma_name, '(f15.2)') MA
write(phi_name, '(I2)') NINT(phi)
write(d_name, '(I2)') NINT(dpul)

phi=phi/180*pi

ncr=ndis0(:,1:43,:)
where(ncr .le. 0) ncr=1d-300

rmax=250d0

DO i=1, 101
theta(i)=(i-1d0)/(10d0)/180*pi  !0-10°
ENDDO

print*, 'start to calculate the intensity map'
intensity_ic=0
CALL CPU_TIME(time1)
!$omp parallel private(xmin, xmax, xx,  flux_ic_tot,zt, rt, count1, count2, ncr_temp1, ncr_temp2, dne, sinalpha, cosalpha, flux_ic, single_spec)
!$omp do schedule(guided)

DO i=1, 100                     !theta
   if (dpul*dsin(theta(i)) .gt. rmax) cycle
   xmin=dpul*dcos(theta(i))-sqrt(rmax**2-dpul**2*dsin(theta(i))**2)
   xmax=dpul*dcos(theta(i))+sqrt(rmax**2-dpul**2*dsin(theta(i))**2)
   xx=(xmax-xmin)*xfrac+xmin

   do w=1, 120 !zeta
     !flux_syn_tot=0
     flux_ic_tot=0    
     !print*, w
     do j=1, 200   !xfrac
    
      zt=dabs( (dpul-xx(j)*dcos(theta(i)))*dcos(phi) + xx(j)*dsin(theta(i))*dcos(zeta(w))*dsin(phi) )
      rt=sqrt( ( (dpul-xx(j)*dcos(theta(i)))*dsin(phi) -  xx(j)*dsin(theta(i))*dcos(zeta(w))*dcos(phi) )**2 + ( xx(j)*dsin(theta(i))*dsin(zeta(w)) )**2 )
      count1=COUNT(rt .gt. rdis)
      count2=COUNT(zt .gt. zdis) 
      if (count1 .eq. 0) count1=1
      if (count2 .eq. 0) count2=1
      if (count1 .eq. 181) count1=180
      if (count2 .eq. 43) count2=42
      ncr_temp1=exp(dlog( ncr(count1+1, :, :)/ncr(count1, :,:))/(rdis(count1+1)-rdis(count1))*(rt-rdis(count1))+dlog(ncr(count1,:,:)) )  
      where (ncr_temp1 .le. 0) ncr_temp1=1d-300
      ncr_temp2=exp(dlog( ncr_temp1(count2+1,:)/ncr_temp1(count2,:))/(zdis(count2+1)-zdis(count2))*(zt-zdis(count2))+dlog(ncr_temp1(count2,:)) )
      where (ncr_temp2 .le. 0) ncr_temp2=1d-300  
!      dne=ncr_temp2*dzeta*xx(j)*dsin((theta(i)+theta(i+1))/2)*(xx(j)*(theta(i+1)-theta(i)))*(xx(j+1)-xx(j))*3.086d18**3
      dne=ncr_temp2*(xx(j+1)-xx(j))*3.086d18
      cosalpha=dcos(theta(i))*dcos(phi)-dsin(theta(i))*dcos(zeta(w))*dsin(phi)

     !flux_syn=0
!     do k=1, 100
!     if (falpha(k) .eq. 0) cycle
!     sinalpha=sqrt(1-(cosalpha*sqrt(cosdalpha2(k)))**2)
      sinalpha=sqrt(1-(cosalpha*cosdelta)**2)
      !flux_syn=flux_syn*0+syn_spec(nu, gme, m_e, dne, B0*sinalpha)/4/3.14!/(xx(j)*3.086d18)**2!*falpha(k)
!     enddo
      flux_ic=0
     do k=1, 60 !gme
      single_spec=IC_anal(TTcmb/0.511d6, gam(k), Eg1/0.511d6)+IC_anal(TTir/0.511d6, gam(k), Eg1/0.511d6)*cir+IC_anal(TTstar/0.511d6, gam(k), Eg1/0.511d6)*cstar+IC_anal(TTuv/0.511d6, gam(k), Eg1/0.511d6)*cuv
      flux_ic=flux_ic+single_spec*dne(k)*(gam(k+1)-gam(k))/4/3.14!/(xx(j)*3.086d18)**2
     enddo
     !flux_syn_tot=flux_syn_tot+flux_syn*nu
     flux_ic_tot=flux_ic_tot+flux_ic*eg1/0.511d6*Eg1*1.6d-12   
 
  enddo

  !intensity_syn(:,w,i)=flux_syn_tot!/(dzeta*(dcos(theta(i))-dcos(theta(i+1))))
  intensity_ic(:,w,i)=flux_ic_tot!/(dzeta*(dcos(theta(i))-dcos(theta(i+1))))
  enddo
  CALL CPU_TIME(time2)
  print*, i, (time2-time1)/60
enddo
!$omp enddo
!$omp end parallel
open(unit=10, file='intensity_ic_1e28_3muG_p1.6_phi10_Simp_dt4yr_NumE61_NumR181_2wyr_1col.dat')
write(10, *) intensity_ic
close(10)

CONTAINS


FUNCTION IC_anal(Tem, Ee, Eg1)
real*8 pi, r0, h, hbar, m_e, c, alpha3, beta3, a3, b3, c3, alpha4, beta4, a4, b4, c4
real*8 Tem, Ee, Eg1(18), zz(18), tz, x0(18), gg3(18), gg4(18), G3(18), G4(18), IC_anal(18)

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


zz=Eg1/Ee
tz=4*Ee*Tem
x0=zz/(1-zz)/tz


gg3=1/( 1+a3*x0**alpha3/(1+b3*x0**beta3) )
gg4=1/( 1+a4*x0**alpha4/(1+b4*x0**beta4) )


G3=pi**2/6d0*(1+c3*x0)/(1+pi**2*c3/6*x0)*dexp(-x0)*gg3
G4=pi**2/6d0*(1+c4*x0)/(1+pi**2*c4/6*x0)*dexp(-x0)*gg4
IC_anal=2*r0**2*m_e**3*c**4*Tem**2/pi/hbar**3/Ee**2*( zz**2/2d0/(1-zz)*G3+G4 )

WHERE(zz .ge. 1) IC_anal=0

RETURN
END FUNCTION







END PROGRAM


