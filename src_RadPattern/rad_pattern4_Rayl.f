      subroutine rad_pattern_r(feig_buff,eiglen, phvnper,dperin,
     +                         strike,dip,rake,depth, period,nper,
     +                         azimuth,groupT,phaseT,amplitude)
c To calculate group_delay as a function of azimuth and period
      integer*4 ntmax, phvnper, dperin
      parameter (ntmax=200)
      integer*4 pos1, pos2
      integer*4 eiglen
      real*4 v(3,ntmax),dvdz(3,ntmax),ampr(ntmax),ampl(ntmax)
      real*4 tm(6),du(3),vu(3),wvn(2),period(20),ampl_max(ntmax)
      complex*8 br(6),bl(6),sumr
      character*20000000 feig_buff
c(phvlen)
      character*40 bred
      character*2 symbik
      character*1 sigR,sigL
      real*4 pq(181,ntmax),ph(181,ntmax),gr_time(181,ntmax),aml(181,ntmax)
      real*4 azimuth(181),groupT(181,nper),phaseT(181,nper),amplitude(181,nper)
      real*4 cr(ntmax),ur(ntmax),wvr(ntmax),t(ntmax),fr(ntmax)
      real*4 cl(ntmax),ul(ntmax),wvl(ntmax)
      real*4 temp_ph(ntmax),unph(ntmax),grt(ntmax)
      data marg/6/,pi/3.1415927/,oo2pi/0.1591549431/,r/2./,eps/0.0001/
      data const/1.E+20/

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
        m=phvnper
        dper=dperin
c-----------reading OLD_SURF_DEEP output------------c
      symbik='1 '
      nd=1000 
      nt=m
      if(nt.ge.ntmax) then
         STOP"(rad_pattern_r): num of pers in .phv exceeds the limit!"
      endif
      call surfreadRad(feig_buff(1:eiglen),eiglen,sigR,sigL,symbik,nt,nd,
     +              depth,t,cr,ur,wvr,cl,ul,wvl,v,dvdz,ampr,ampl)
c----------Source term calculations-----------------c
      call angles2tensorRad(strike,dip,rake,tm)

c    period loop
      DO j=1,nt
         ampl_max(j)=0.0
         fr(j)=1./t(j)
         vu(1)=v(1,j)
         vu(2)=v(2,j)
         vu(3)=v(3,j)
         du(1)=dvdz(1,j)
         du(2)=dvdz(2,j)
         du(3)=dvdz(3,j)
         w=pi*2.0*fr(j)
         wvn(1)=wvr(j)
         wvn(2)=wvl(j)
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
               sumr= sumr+tm(m)*br(m)
            end do
            aq=cabs(sumr)
            sume=real(sumr)
            sumi=aimag(sumr)
            pq(jkl,j)=pha(sumi,sume)
            aml(jkl,j)=aq*ampr(j)*const
         EndDo
1     ENDDO

c-------------unwrap phase and get group time---------
      do jkl=1,181
         do j=1,nt
            temp_ph(j)=pq(jkl,j)
            if(aml(jkl,j).gt.ampl_max(j)) ampl_max(j)=aml(jkl,j)
         enddo
         call unwrapR(dper,t,nt,temp_ph,unph,grt,r)
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
      DO jop=1,nper 
C---------------selection and output -----------------------S
         Do jpa=1,nt 
         if(abs(period(jop)-t(jpa)).lt.eps)then
            do jkl=1,181
               groupT(jkl,jop)=gr_time(jkl,jpa)
               phaseT(jkl,jop)=pq(jkl,jpa)*oo2pi*t(jpa)
               amplitude(jkl,jop)=aml(jkl,jpa)
            enddo
         endif
         endDo
      endDO
      end
