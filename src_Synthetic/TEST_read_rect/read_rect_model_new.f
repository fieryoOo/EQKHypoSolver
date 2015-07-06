c ==========================================================
c read rectangular models into memory
c ==========================================================
c --- common /mdl/ ----------------------
      module mmodel
      type tmodel
         integer*4 ic,jc,n,nper,nfi,nla
         real*8 fi,sfi,la,sla,per,sper,bf,ef,bl,el,bp,ep,
     +          hfi(97),hla(201),hper(225),chfi(97)
         real*8 uw(97,201),cw(97,201)
         real*8 gw(97,201),aw(97,201)
      end type
      end module mmodel

      subroutine read_rect_model_new(namea,nmod,p,ierr)
      use mmodel
      implicit none
      integer*4 nmod,ierr
      real*8    p,geo,drad
      character*255 namea
C      character*90000000 fmod_buff
C      integer*4 buffsize
C      common /mdl/ic,jc,n,nper,nfi,nla,fi,sfi,la,sla,per,sper,
C     +        bf,ef,bl,el,bp,ep,
C     +        hfi,hla,hper,chfi,uw,cw,gw,aw

      integer*4 i

      type (tmodel) mdl

c ---
      ierr = 0
c --- read unformated file for period per--
      if(nmod.eq.0) then
        mdl%ic = 1
        mdl%jc = 1
        geo = 0.993277d0
        drad = datan(1.0d0)/45.0d0
        open(10,file=namea,form='unformatted',status='old')
        read(10) mdl%n,mdl%fi,mdl%nfi,mdl%sfi,mdl%la,mdl%nla,mdl%sla,mdl%per,mdl%nper,mdl%sper
C      write(*,*) " in read 1: ",mdl%n,mdl%fi,mdl%nfi,mdl%sfi,mdl%la,mdl%nla,mdl%sla,mdl%per,mdl%nper,mdl%sper
        mdl%n =0
        do i = 1,mdl%nfi
          mdl%hfi(i) = mdl%fi+(i-1)*mdl%sfi
          mdl%chfi(i) = datan(geo*dtan(drad*mdl%hfi(i)))/drad
        enddo
        mdl%bf = mdl%chfi(1)
        mdl%ef = mdl%chfi(mdl%nfi)
        do i = 1,mdl%nla
          mdl%hla(i) = mdl%la+(i-1)*mdl%sla
        enddo
        mdl%bl = mdl%hla(1)
        mdl%el = mdl%hla(mdl%nla)
        do i = 1,mdl%nper
          mdl%hper(i) = mdl%nper+(i-1)*mdl%sper
        enddo
        mdl%bp = mdl%hper(1)
        mdl%ep = mdl%hper(mdl%nper)
        return
      else if(nmod.eq.-1) then
        close(10)
        return
      endif
      p = mdl%per+mdl%n*mdl%sper
      read(10) mdl%uw,mdl%cw,mdl%gw,mdl%aw
C      write(*,*) " in read 2: ",uw,cw,gw,aw
      mdl%n = mdl%n+1
      end

