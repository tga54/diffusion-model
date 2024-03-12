PROGRAM main
!$use_omp_lib
implicit none

integer*4, parameter :: NumE=151, NumT=150000, NumR=601, NumZ=1001
integer*4 i,j,k,l,w, tt
real*8 c, pi, kb, m_e, m_p, sigmaT, braking, spin_index,p0,period,pdot
real*8 r(NumR), Qe(NumE), NinjE(NumE), gam(NumE), pcr(NumE), dr, dz, lnp(NumE), dlnp, dt, r0, MA
real*8 DEpara(NumE), DEperp(NumE), ncrE(NumR, NumZ, NumE), ncrr(NumR, NumZ, NumE), ncrz(NumR, NumZ, NumE), ncr1(NumR, NumZ, NumE), ncr2(NumR, NumZ, NumE), bp(NumE)
real*8 aaz(NumZ,NumE), bbz(NumZ,NumE), ccz(NumZ,NumE), ddz(NumZ,NumE), eez(NumZ,NumE), ffz(NumZ,NumE)
real*8 aar(NumR,NumE), bbr(NumR,NumE), ccr(NumR,NumE), ddr(NumR,NumE), eer(NumR,NumE), ffr(NumR,NumE), p, temp, Le, dpcr(NumE-1)
real*8 Tstar, Tir, Tcmb, Tuv, B0, B1, ustar, uir, ucmb, uuv, Qt(NumT), tauc, t0, t, tage,  timei, timef
real*8 ndis0(501, 301, 151)  ,ndis(600,151) !NumR, NumZ, NumE
integer :: omp_get_thread_num
!$ call omp_set_num_threads(251)

c=2.998d10
pi=3.1416d0
kb=1.38d-16
m_e=9.1e-28
sigmaT=6.65d-25

B0=3d-6

Tcmb=2.73d0
Tir=30d0
Tstar=5d3
Tuv=2d4

ucmb=4.2d-13
uir=4.8d-13
ustar=4.8d-13
uuv=1.6d-13*0


do i=1, NumR
r(i)=(i-1)*3.086d18*0.5
enddo

dr=r(2)-r(1)
dz=3.086d18*1    !original:dz=2pc,dr=0.3pc,dt=40yr

do i=1, NumE
pcr(i)=10**((i-1)/30d0+5)*1d6*1.6d-12/c           ! momentum in unit of erg*s/cm, eneryg from 1e6eV - 1e18eV,     1/1000能量精度+1yr时间精度->不稳定, 1/100+10yr -> 稳定
enddo

do i=1, NumE-1
dpcr(i)=pcr(i+1)-pcr(i)
enddo

gam=sqrt((pcr)**2*c**2+m_e**2*c**4)/m_e/c**2

r0=50
B1=3d-6
!DO i=1, NumR
! IF (r(i) .LT. r0*3.086D18) THEN
!   DE(i,:)=4.5d27*(pcr*c/1.6d-12/1d14)**(1d0/3)      ! diffusion coefficient, normalized at 1GeV
!   bp(i,:)=4d0/3*sigmaT*(gam**2-1)*(B1**2/8/pi+ucmb/(1+4*gam*3*kb*Tcmb/m_e/c**2)**1.5+uir/(1+4*gam*3*kb*Tir/m_e/c**2)**1.5+ustar/(1+4*gam*3*kb*Tstar/m_e/c**2)**1.5+uuv/(1+4*gam*3*kb*Tuv/m_e/c**2)**1.5)
! ELSE
!   DE(i,:)=3.86d28*(pcr*c/1.6d-12/1d9)**(1d0/3) 
!   bp(i,:)=4d0/3*sigmaT*(gam**2-1)*(B0**2/8/pi+ucmb/(1+4*gam*3*kb*Tcmb/m_e/c**2)**1.5+uir/(1+4*gam*3*kb*Tir/m_e/c**2)**1.5+ustar/(1+4*gam*3*kb*Tstar/m_e/c**2)**1.5+uuv/(1+4*gam*3*kb*Tuv/m_e/c**2)**1.5)
!  ENDIF
!ENDDO

MA=0.25
DEpara=1d28*(pcr*c/1.6d-12/1d9)**(1d0/2)
DEperp=DEpara*MA**4
bp=4d0/3*sigmaT*(gam**2-1)*(B0**2/8/pi+ucmb/(1+(gam*2.82*kb*Tcmb/m_e/c**2)**0.6)**(1.9/0.6)+uir/(1+(gam*2.82*kb*Tir/m_e/c**2)**0.6)**(1.9/0.6)+ustar/(1+(gam*2.82*kb*Tstar/m_e/c**2)**0.6)**(1.9/0.6)+uuv/(1+(gam*2.82*kb*Tuv/m_e/c**2)**0.6)**(1.9/0.6))

lnp=DLOG(pcr)
dlnp=lnp(2)-lnp(1)

dt=3.1536d7*3
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
p=1.6
Ninje=gam**(-p)*exp(-gam/4d8)    ! 1.4 4d8  , 1.6 4d8
where(gam*0.511d6 .lt. 1d11) Ninje=0
temp=0
DO i=1, NumE-1
temp=temp+( gam(i)*m_e*c**2*Ninje(i)+gam(i+1)*m_e*c**2*Ninje(i+1) )*( gam(i+1)-gam(i) )/2
ENDDO

Qe=Le/temp*Ninje

print*, tage/3.1536d7

ncr1=0
ncr1(1,1,:)=Qe/6/pi/dr**2/dz
ncr1(1,2,:)=Qe/6/pi/dr**2/dz
ncr1(2,1,:)=Qe/6/pi/dr**2/dz
ncr1(2,2,:)=Qe/6/pi/dr**2/dz
ncr1(:,:,NumE)=0



tt=1     !时间，用于计数，把感兴趣的时间点的数据记录下来
ndis0=0
ncr2=0
ncrr=0
ncrz=0
ncrE=0

DO l=1, NumT   !时间
if (t .GT. tage) cycle
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
     bbr(i,:)=( 1+DEperp*dt/dr**2 )
     ccr(i,:)=-DEperp*dt/(dr**2)
     IF (j .EQ. 1) ddr(i,:)=DEpara*dt/(dz**2)*(ncrE(i,j+1,:)-ncrE(i,j,:))+ncrE(i,j,:)
     IF (j .NE. 1) ddr(i,:)=DEpara*dt/(2*dz**2)*(ncrE(i,j+1,:)-2*ncrE(i,j,:)+ncrE(i,j-1,:))+ncrE(i,j,:)
     eer(i,:)=ccr(i,:)/bbr(i,:)
     ffr(i,:)=ddr(i,:)/bbr(i,:)
    ELSE
     aar(i,:)=-DEperp*dt/(2*dr**2)*( 1-dr/r(i)/2 )       
     bbr(i,:)=1+DEperp*dt/dr**2
     ccr(i,:)=-DEperp*dt/(2*dr**2)*( 1+dr/r(i)/2 )
     IF (j .EQ. 1) ddr(i,:)=DEpara*dt/(dz**2)*( ncrE(i,j+1,:) - ncrE(i,j,:) ) + ncrE(i,j,:)
     IF (j .NE. 1) ddr(i,:)=DEpara*dt/(2*dz**2)*( ncrE(i,j+1,:) - 2*ncrE(i,j,:)+ ncrE(i,j-1,:) ) + ncrE(i,j,:)
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
     bbz(j,:)=( 1+DEpara*dt/dz**2 )
     ccz(j,:)=-DEpara*dt/(dz**2)
     IF (i .EQ. 1) ddz(j,:)=DEperp*dt/(dr**2)*(ncrr(i+1,j,:)-ncrr(i,j,:))+ncrr(i,j,:)
     IF (i .NE. 1) ddz(j,:)=DEperp*dt/(4*r(i)*dr)*( ncrr(i+1,j,:)-ncrr(i-1,j,:) )+ DEperp*dt/(2*dr**2)*( ncrr(i+1,j,:)-2*ncrr(i,j,:)+ncrr(i-1,j,:) ) + ncrr(i,j,:)
     eez(j,:)=ccz(j,:)/bbz(j,:)
     ffz(j,:)=ddz(j,:)/bbz(j,:)
    ELSE
     aaz(j,:)=-DEpara*dt/(2*dz**2)     
     bbz(j,:)=1+DEpara*dt/dz**2
     ccz(j,:)=-DEpara*dt/(2*dz**2)
     IF (i .EQ. 1) ddz(j,:)=DEperp*dt/(dr**2)*( ncrr(i+1,j,:) - ncrr(i,j,:) ) + ncrr(i,j,:)
     IF (i .NE. 1) ddz(j,:)=DEperp*dt/(4*r(i)*dr)*( ncrr(i+1,j,:) - ncrr(i-1,j,:) ) + DEperp*dt/(2*dr**2)*( ncrr(i+1,j,:)-2*ncrr(i,j,:)+ncrr(i-1,j,:) ) + ncrr(i,j,:)
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
   ndis0=ndis0+(1+tage/t0)**spin_index/(1+(tage-t)/t0)**spin_index*ncr2(1:501, 1:301, :)*dt
! ENDDO
!ENDIF
ncr2=0


!print*, ncr2(100,100,30)
IF (MOD(l, 200) .EQ. 0) THEN
ndis(tt,:)=ndis0(1,250,:)
print*, tt, t/3.1536d7,'yr'
tt=tt+1
ENDIF
t=t+dt
ENDDO



OPEN(UNIT=10, FILE='ncr0_ani_3muG_1.6_1d28_MA0.25_p050_n3_200TeV_dt3yr.dat')
WRITE(10,*) ndis0
CLOSE(10)

OPEN(UNIT=10, FILE='ani_positron_t.dat')
WRITE(10,*) ndis
CLOSE(10)





END
