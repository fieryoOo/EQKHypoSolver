      subroutine rad_pattern_r(tm,nper,dper,t,eigH,deigH,eigV,deigV,ampr,wvr,
     +                         perlist,nperlst,azimuth,groupT,phaseT,amplitude)
c To calculate group_delay as a function of azimuth and period
      integer*4 nper, nperlst
C      real*4 strike, dip, rake
C      parameter (ntmax=500)
      integer*4 pos1, pos2
      real*4 eigH(nper), deigH(nper), eigV(nper), deigV(nper)
      real*4 t(nper),fr(nper),wvr(nper),ampr(nper),ampl_max(nper)
C      real*4 v(3,nper),dvdz(3,nper),
C      v(eigen), dvdz(derivative_eigen): 1=R_horizontal, 2=R_vertical, 3=L
      real*4 tm(6),du(3),vu(3),wvn(2)
      complex*8 br(6),bl(6),sumr,step
      character*40 bred
C      character*2 symbik
      character*1 sigR,sigL
      real*4 pq(181,nper),ph(181,nper),gr_time(181,nper),aml(181,nper)
      real*4 perlist(nperlst),azimuth(181),groupT(181,nperlst),phaseT(181,nperlst),amplitude(181,nperlst)
C      real*4 coef_amp(nperlst), wavenum(nperlst)
      real*4 temp_ph(nper),unph(nper),grt(nper),stepr,dper
      data marg/6/,pi/3.1415927/,oo2pi/0.1591549431/,r/2./,eps/0.0001/
      data const/1.E+20/,const2/5013256549262000.0/

C----------- Initiation------------
      drad=180./pi
      sigR='+'
      sigL='-'
C-----------read m and dper from .phv file ------
c!$OMP CRITICAL
c      open(1,file=phvfile(1:phvlen),status='OLD')

c      m=0
c      pos1=1
c      pos2=pos1+40
c      do i=1,1000
c         read(fphv_buff(pos1:pos2),'(a)', end=9988) bred
c         pos2 = INDEX(fphv_buff(pos1:),NEW_LINE('a'))
c         pos1 = pos1+pos2
c         pos2 = pos1+40
c         if(i.eq.1)read(bred,*)per1
c         if(i.eq.2)read(bred,*)per2
cc         if(bred(6:10).eq.'     ')go to 9988
c         if(pos1.ge.phvlen)goto 9988
c         m=m+1
c      enddo
c9988  close(1)
c!$OMP END CRITICAL
c9988     dper=per2-per1
c-----------reading OLD_SURF_DEEP output------------c
C      symbik='1 '
C      nd=1000 
      nt=nper
c      if(nt.ge.ntmax) then
c         STOP"(rad_pattern_r): num of pers in .phv exceeds the limit!"
c      endif
C      call surfreadRad(feig_buff(1:eiglen),eiglen,sigR,sigL,symbik,nt,nd,
C     +              depth,t,cr,ur,wvr,cl,ul,wvl,v,dvdz,ampr,ampl)
c----------Source term calculations-----------------c
C      call angles2tensorRad(strike,dip,rake,tm)
C      write(*,*) tm(1), tm(2), tm(3), tm(4), tm(5), tm(6)

c    period loop
      DO j=1,nt
         ampl_max(j)=0.0
         fr(j)=1./t(j)
         vu(1)=eigH(j)
         vu(2)=eigV(j)
C         vu(3)=v(3,j)
         du(1)=deigH(j)
         du(2)=deigV(j)
C         du(3)=dvdz(3,j)
         w=pi*2.0*fr(j)
         wvn(1)=wvr(j)
c         wvn(2)=wvl(j)
c         step=cmplx(0.,-1./w)
         stepr=1./w
c        azimuthal loop
         Do  jkl=1,181
            AZI=2.*float(jkl-1)
            AZ_rad=AZI/drad
            cs=cos(AZ_rad)
            sc=sin(AZ_rad)
c           convolution with moment tensor
            call sourceRad(sigR,sigL,cs,sc,wvn,vu,du,br,bl)
            sumr=(0.0,0.0)
            do m=1,6
               sumr= sumr+tm(m)*br(m)*stepr
            end do
            aml(jkl,j)=cabs(sumr)
            sume=real(sumr)
            sumi=aimag(sumr)
            pq(jkl,j)=phaRad(sumi,sume)
c      if(t(j).eq.10.) write(*,*) "rad_pattern_Rayl new:",t(j),pq(jkl,j),sumi,sume,cabs(sumr)
cYT            aml(jkl,j)=aq*ampr(j)*const
cYT            aml(jkl,j)=aq*sqrt(ampr(j)*const2)
cYT         source amp norm term : ampr = 1./(2.*c*ugr*sumi0)/sqrt(6.28318)*1.e-15
cYT            aml(jkl,j)=aq* sqrt( ampr(j)*1.0e15*sqrt(8.0*pi) )
         EndDo
1     ENDDO

c-------------unwrap phase and get group time---------
      do jkl=1,181
         do j=1,nt
            temp_ph(j)=pq(jkl,j)
            if(aml(jkl,j).gt.ampl_max(j)) ampl_max(j)=aml(jkl,j)
         enddo
         call unwrap(dper,t,nt,temp_ph,unph,grt,r)
         do j=1,nt
            pq(jkl,j)=unph(j)
            gr_time(jkl,j)=grt(j)
            ph(jkl,j)=unph(j)
         enddo
      enddo
c----azimuth-dependent output for a set of requested periods------
      do jkl=1,181
         azimuth(jkl) = 2. * (jkl-1)
      enddo
      DO jop=1,nperlst
C         write(*,*) "nperlst = ",nperlst
C---------------selection and output -----------------------S
         Do jpa=1,nt 
C        write(*,*) "nper = ",jpa,perlist(jop),t(jpa)
         if(abs(perlist(jop)-t(jpa)).lt.eps)then
C            coef_amp(jop)=ampr(jpa); wavenum(jop)=wvr(jpa)
            do jkl=1,181
               groupT(jkl,jop)=gr_time(jkl,jpa)
               phaseT(jkl,jop)=pq(jkl,jpa)*oo2pi*t(jpa)
               amplitude(jkl,jop)=aml(jkl,jpa)
            enddo
         endif
         endDo
      endDO
      end
