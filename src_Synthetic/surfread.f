       recursive subroutine surfread(feig_buff,eiglen,sigR,sigL,symb,nt,depth,
     +                     fr,cr,ur,wvr,cl,ul,wvl,v,dvdz,ampr,ampl,ratio,qR,qL,I0)
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
      character*100 bred,vzdor,mura
      character*300 linetmp
c      character*20000000 feig_buff
      character*99999999 feig_buff
      character*1 sigR,sigL
      character*2 symb
      real*4 depold, vaold, vdold, dfactor
      real*4 v(3,nt+10),dvdz(3,nt+10),ratio(nt+10)
      real*4 t(nt+10),fr(nt+10),ur(nt+10),cr(nt+10),wvr(nt+10),ampr(nt+10),qR(nt+10)
      real*4 ul(nt+10),cl(nt+10),wvl(nt+10),ampl(nt+10),qL(nt+10),I0(nt+10)
      data pi2/6.28318/,tlim/10000.0/,eps/1.0/
C-----------------------------------------------------------------
      if( eiglen > 99999999 ) STOP 'eigen file size exceeds the limit!'
      lsy=lnblnk(symb)
      wvl(1)=0.0
      qL(1)=0.0
      ul(1)=0.0
      t(1)=0.0
      t(nt+2)=10000.
      ires=0
      nd=10000
      nd_real=nd
C----------Reading Love stuff----------------------S
      if (sigL.eq.'+') THEN
c         open(1,file=infile,STATUS='OLD')
         pos1=1
         do m=1,2000000
             call readline100(feig_buff, pos1, pos2, bred)
             if(pos1.ge.eiglen) goto 11
c            read(feig_buff(pos1:pos2),'(a)',end =11),bred
c            pos2 = INDEX(feig_buff(pos1:),NEW_LINE('a'))
c            pos1 = pos1+pos2
c            pos2 = pos1+100
            do j=1,40
               if(bred(j:j+3).eq.'Love')then
                  if(lsy.eq.1.and.bred(j+11:j+11).eq.symb(1:1))go to 1
                  if(lsy.eq.2.and.bred(j+11:j+12).eq.symb(1:2))go to 1
               endif
            end do
         end do
         STOP 'NO LOVE'
C----------PERIOD LOOP----------------------S
1       ires=1
c      write(*,*) "Love found in line ",bred," with pos1 = ",pos1,
c     +           " and pos2 = ", pos2
        Do k=2,nt+1
          if(nd_real.eq.nd)  then
c           read line '@@@@@@'
               call readline100(feig_buff, pos1, pos2, vzdor)
c      write(*,*) "Line @@@ = ",vzdor
               if(pos1.ge.eiglen) goto 11
c               read(feig_buff(pos1:pos2),'(a)',end=11)vzdor
c               pos2 = INDEX(feig_buff(pos1:),NEW_LINE('a'))
c               pos1 = pos1+pos2
c               pos2 = pos1+100
            endif
c            read(1, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
            call readline300(feig_buff, pos1, pos2, linetmp)
c        write(*,*) "a new line (",pos1,",",pos2,") with k=",k," nt=",nt," : ",linetmp
            read(linetmp, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
            fr(k) = 1./t(k)
c            read(1,'(a)') ,vzdor
            call readline100(feig_buff, pos1, pos2, vzdor)
C            read(vzdor, '(E14.7)') I0(k)
            read(vzdor, *) I0(k)
c        write(*,*) "                second (",pos1,",",pos2,"): ",vzdor," ",I0(k)
c           Love component---S
            depold=0.
            vaold=0.
            vdold=0.
            found=0
            do l=1,nd
c               read(1,'(a)')mura
               call readline100(feig_buff, pos1, pos2, mura)
               if(pos1.ge.eiglen) goto 11
c      write(*,*) "   read line (",pos1,",",pos2,") ",mura
               if(mura(2:4).eq.'@@@')then
                  nd_real=l-1 
                  go to 5454
               endif
               if(found.eq.0) then
                  read(mura,'(3(E14.7,2X))',end=11),dep,va,vd
                  if(depold.le.depth.and.depth.lt.dep) then
c                 if(abs(dep-depth).lt.eps) then
                     dfactor=(depth-depold)/(dep-depold)
                     v(3,k)=vaold+(va-vaold)*dfactor
                     dvdz(3,k)=vdold+(vd-vdold)*dfactor
                     found = 1
                  end if
                  depold=dep
                  vaold=va
                  vdold=vd
               endif
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
         qL(nt+2)=0.0
         ul(nt+2)=0.0
      END if
C----------Reading Love stuff----------------------E
c      if(ires.eq.0) print *,'NO LOVE'
11    continue
C----------Reading Rayleigh stuff----------------------S
      wvr(1)=0.0
      qR(1)=0.0
      ur(1)=0.0
      if (sigR.eq.'+') Then
c         open(1,file=infile,STATUS='OLD')
C----------------------Search for Rayleigh--------S
         pos1=1
         do m=1,2000000
            call readline100(feig_buff, pos1, pos2, bred)
            if(pos1.ge.eiglen) goto 99
c            read(1,'(a)',end =99),bred
            do j=1,40
            if(bred(j:j+3).eq.'Rayl') then
c      write(*,*) "Rayl found in line ",bred," with pos1 = ",pos1," and
c     +                  eiglen = ", eiglen
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
c           read line '@@@@@@'
            call readline100(feig_buff, pos1, pos2, vzdor)
c        write(*,*) " line @@@ = ",vzdor
            if(pos1.ge.eiglen) goto 9797
c            read(1, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
c           read line1 with 7? params
            call readline300(feig_buff, pos1, pos2, linetmp)
c        write(*,*) "a new line (",pos1,",",pos2,") with k=",k," nt=",nt," : ",linetmp
            read(linetmp, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
            fr(k) = 1./t(k)
C           PRint*,k,t(k),ampr(k)
c            read(1,'(a)') ,vzdor
c           read line2 with 5? params
            call readline100(feig_buff, pos1, pos2, vzdor)
            read(vzdor, '(E14.7)') I0(k)
c        write(*,*) "                second (",pos1,",",pos2,"): ",vzdor," ",I0(k)
C----------Rayl. Horizontal component------S
            depold=0.
            vaold=0.
            vdold=0.
            do l=1,nd
c               read(1,'(a)')mura
               call readline100(feig_buff, pos1, pos2, mura)
c      write(*,*) "   read line 1 (",pos1,",",pos2,") ",mura
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
c5554        write(*,*) l,"  ",nd_real," lines read for part 1"
5554        if(nd_real.eq.nd) call readline100(feig_buff, pos1, pos2, vzdor)
C----------Rayl. Vertical   component------S
            depold=0.
            vaold=0.
            vdold=0.
            do l=1,nd_real
c               read(1,*),dep,va,vd
               call readline100(feig_buff, pos1, pos2, mura)
               read(mura,*),dep,va,vd
c      write(*,*) "   read line 2 (",pos1,",",pos2,") ",mura
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
c      write(*,*) l,"  ",nd_real," lines read for part 2"
C----------Rayl. Vertical   component------E
            wvr(1)=pi2/cr(2)/tlim
            wvr(nt+2)=0.0
            ampr(nt+2)=0.0
            qR(nt+2)=0.0
            ur(nt+2)=0.0
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

      subroutine readline100( fbuff, pos1, pos2, linetmp )
      integer*4 pos1, pos2
c      character*20000000 fbuff
      character*99999999 fbuff
      character*100 linetmp
c      pos2 = pos1+100
      pos2 = pos1 + INDEX(fbuff(pos1:),NEW_LINE('a'))
      read(fbuff(pos1:pos2),'(a)') linetmp
c      pos1 = pos1+pos2
      pos1 = pos2
      end

      subroutine readline300( fbuff, pos1, pos2, linetmp )
      integer*4 pos1, pos2
c      character*20000000 fbuff
      character*99999999 fbuff
      character*300 linetmp
c      pos2 = pos1+300
      pos2 = pos1 + INDEX(fbuff(pos1:),NEW_LINE('a'))
      read(fbuff(pos1:pos2),'(a)') linetmp
c      pos1 = pos1+pos2
      pos1 = pos2
      end
