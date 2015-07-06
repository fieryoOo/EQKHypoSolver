       subroutine surfread(feig_buff,eiglen,sigR,sigL,symb,nt,nd,depth,
     +            fr,cr,ur,wvr,cl,ul,wvl,v,dvdz,ampr,ampl,ratio,qR,qL,I0)
C------to read SURFLEV output for producing synthetic seismograms--
C----------------INPUT ARGUMENTS-----------------------------------
C------infile -file with spectral eigenfunctions, wvnumbers ,etc---
C------sign is R or L----------------------------------------------
C------nt is number of periods-------------------------------------
C------nd is number of depths for which these functions exist------
C------depth is a wanted depth of a point source-------------------
C---------------OUTPUT ARGUMENTS-----------------------------------
C------c is a phase velocity;u -group velocity;wvn -wavenumber;
C------v(3,2000)---eigenfunction; dvdz(3,2000) - depth derivative; 
C------amp -ampitude factor; ratio is ellipticity------
      integer*4 pos1, pos2, eiglen
      character*80 bred,vzdor,mura
      character*300 linetmp
      character*20000000 feig_buff
      character*1 sigR,sigL
      character*2 symb
      real*4 depold, vaold, vdold, dfactor
      real*4 v(3,2000),dvdz(3,2000),ratio(2000)
      real*4 t(2000),fr(2000),ur(2000),cr(2000),wvr(2000),ampr(2000),qR(2000)
      real*4 ul(2000),cl(2000),wvl(2000),ampl(2000),qL(2000),I0(2000)
      data pi2/6.28318/,tlim/10000.0/,eps/1.0/
C-----------------------------------------------------------------
      lsy=lnblnk(symb)
      wvl(1)=0.0
      wvr(1)=0.0
      t(1)=0.0
      t(nt+2)=10000.
      ires=0
      nd_real=nd
C----------Reading Love stuff----------------------S
      if (sigL.eq.'+') THEN
c         open(1,file=infile,STATUS='OLD')
         pos1=1
         do m=1,2000000
             call readline80(feig_buff, pos1, pos2, bred)
             if(pos1.ge.eiglen) goto 11
c            read(feig_buff(pos1:pos2),'(a)',end =11),bred
c            pos2 = INDEX(feig_buff(pos1:),NEW_LINE('a'))
c            pos1 = pos1+pos2
c            pos2 = pos1+80
            do j=1,80
               if(bred(j:j+3).eq.'Love')then
                  if(lsy.eq.1.and.bred(j+11:j+11).eq.symb(1:1))go to 1
                  if(lsy.eq.2.and.bred(j+11:j+12).eq.symb(1:2))go to 1
               endif
            end do
         end do
         STOP 'NO LOVE'
C----------Reading Love stuff----------------------S
1       ires=1
        Do k=2,nt+1
            if(nd_real.eq.nd)  then
               call readline80(feig_buff, pos1, pos2, vzdor)
               if(pos1.ge.eiglen) goto 11
c               read(feig_buff(pos1:pos2),'(a)',end=11)vzdor
c               pos2 = INDEX(feig_buff(pos1:),NEW_LINE('a'))
c               pos1 = pos1+pos2
c               pos2 = pos1+80
            endif
c            read(1, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
            call readline300(feig_buff, pos1, pos2, linetmp)
            read(linetmp, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
c            read(1,'(a)') ,vzdor
            call readline80(feig_buff, pos1, pos2, vzdor)
c           Love component---S
            depold=0.
            vaold=0.
            vdold=0.
            do l=1,nd
c               read(1,'(a)')mura
               call readline80(feig_buff, pos1, pos2, mura)
               if(pos1.ge.eiglen) goto 11
               if(mura(2:4).eq.'@@@')then
                  nd_real=l-1 
                  go to 5454
               endif
               read(mura,'(3(E14.7,2X))',end=11),dep,va,vd
               if(depold.le.depth.and.depth.lt.dep) then
c               if(abs(dep-depth).lt.eps) then
                  dfactor=(depth-depold)/(dep-depold)
                  v(3,k)=vaold+(va-vaold)*dfactor
                  dvdz(3,k)=vdold+(vd-vdold)*dfactor
               end if
               depold=dep
               vaold=va
               vdold=vd
            end do
6666        nd_real=nd
c           Love component---E
5454        continue
         enddo
C         ires=1
c         close(1)
         wvl(1)=pi2/cl(2)/tlim
         wvl(nt+2)=0.0
         ampl(nt+2)=0.0

      END if
C----------Reading Love stuff----------------------E
c      if(ires.eq.0) print *,'NO LOVE'
11    continue
C----------Reading Rayleigh stuff----------------------S
      if (sigR.eq.'+') Then
c         open(1,file=infile,STATUS='OLD')
C----------------------Search for Rayleigh--------S
         pos1=1
         do m=1,200000
            call readline80(feig_buff, pos1, pos2, bred)
            if(pos1.ge.eiglen) goto 99
c            read(1,'(a)',end =99),bred
            do j=1,80
            if(bred(j:j+3).eq.'Rayl') then
               if(lsy.eq.1.and.bred(j+15:j+15).eq.symb(1:1))go to 2
               if(lsy.eq.2.and.bred(j+15:j+16).eq.symb(1:2))go to 2
            endif
            end do
         end do
         STOP 'NO RAYLEIGH'
C----------------------------------PERIOD LOOP------S
2        DO k=2,nt+1
            nd_real=nd
c            read(1,'(a)',end=9797)vzdor
            call readline80(feig_buff, pos1, pos2, vzdor)
            if(pos1.ge.eiglen) goto 9797
c            read(1, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
            call readline300(feig_buff, pos1, pos2, linetmp)
c        write(*,*) "a new line with k=",k," : ",linetmp
            read(linetmp, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
            fr(k) = 1./t(k)
C           PRint*,k,t(k),ampr(k)
c            read(1,'(a)') ,vzdor
            call readline80(feig_buff, pos1, pos2, vzdor)
            read(vzdor, '(E14.7)') I0(k)
C----------Rayl. Horizontal component------S
            depold=0.
            vaold=0.
            vdold=0.
            do l=1,nd
c               read(1,'(a)')mura
               call readline80(feig_buff, pos1, pos2, mura)
               if(mura(2:4).eq.'$$$')then 
                  nd_real=l-1
                  go to 5554
               endif
               read(mura,'(3(E14.7,2X))'),dep,va,vd
               if(depold.le.depth.and.depth.lt.dep) then
c               if(abs(dep-depth).lt.eps) then
                  dfactor=(depth-depold)/(dep-depold)
                  v(1,k)=vaold+(va-vaold)*dfactor
                  dvdz(1,k)=vdold+(vd-vdold)*dfactor
               end if
               depold=dep
               vaold=va
               vdold=vd
            end do
C----------Rayl. Horizontal component------E
C5554        continue   
5554        if(nd_real.eq.nd) call readline80(feig_buff, pos1, pos2, vzdor)
C----------Rayl. Vertical   component------S
            depold=0.
            vaold=0.
            vdold=0.
            do l=1,nd_real
c               read(1,*),dep,va,vd
               call readline300(feig_buff, pos1, pos2, linetmp)
               read(linetmp,*),dep,va,vd
               if(depold.le.depth.and.depth.lt.dep) then
c               if(abs(dep-depth).lt.eps) then
                  dfactor=(depth-depold)/(dep-depold)
                  v(2,k)=vaold+(va-vaold)*dfactor
                  dvdz(2,k)=vdold+(vd-vdold)*dfactor
               end if
               depold=dep
               vaold=va
               vdold=vd
            end do
C----------Rayl. Vertical   component------E
            wvr(1)=pi2/cr(2)/tlim
            wvr(nt+2)=0.0
            ampr(nt+2)=0.0
         end DO
C----------------------------------PERIOD LOOP------E
C----------Reading Rayleigh stuff----------------------E
c9797     close(1)
C         nt=k-1
9797     ires=ires+1
      End if
99    if(ires.gt.0)  return
      STOP 'NO RAYLEIGH OR LOVE'
      END


C      subroutine readline80( fbuff, pos1, pos2, linetmp )
C      integer*4 pos1, pos2
C      character*20000000 fbuff
C      character*80 linetmp
C      pos2 = pos1+80
C      read(fbuff(pos1:pos2),'(a)') linetmp
C      pos2 = INDEX(fbuff(pos1:),NEW_LINE('a'))
C      pos1 = pos1+pos2
C      end
C
C      subroutine readline300( fbuff, pos1, pos2, linetmp )
C      integer*4 pos1, pos2
C      character*20000000 fbuff
C      character*300 linetmp
C      pos2 = pos1+300
C      read(fbuff(pos1:pos2),'(a)') linetmp
C      pos2 = INDEX(fbuff(pos1:),NEW_LINE('a'))
C      pos1 = pos1+pos2
C      end
