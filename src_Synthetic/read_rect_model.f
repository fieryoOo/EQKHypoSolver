c ==========================================================
c read rectangular models into memory
c ==========================================================

      recursive subroutine read_rect_model(model,nmod,p,ierr,mdl)
      use mmodel
      use omp_lib
      implicit none
      integer*4 i,nmod,ierr
      real*8    p,geo,drad
C      integer*4 TID

      type (tmdl) mdl
      type (tmodel) model

c ---
      ierr = 0
C      TID = OMP_GET_THREAD_NUM()
c --- read unformated file for period per--
      if(nmod.eq.0) then
        mdl%ic = 1
        mdl%jc = 1
        geo = 0.993277d0
        drad = datan(1.0d0)/45.0d0
C      open(fid,file=namea,form='unformatted')
      mdl%n = model%n; mdl%fi = model%fi; mdl%nfi = model%nfi; mdl%sfi = model%sfi
      mdl%la = model%la; mdl%nla = model%nla; mdl%sla = model%sla
      mdl%per = model%per; mdl%nper = model%nper; mdl%sper = model%sper
C      read(fid) mdl%n,mdl%fi,mdl%nfi,mdl%sfi,mdl%la,mdl%nla,mdl%sla,mdl%per,mdl%nper,mdl%sper
C      write(*,*) " in read 1: tid=",OMP_GET_THREAD_NUM()," params = ",
C     +           mdl%n,mdl%fi,mdl%nfi,mdl%sfi,mdl%la,mdl%nla,mdl%sla,mdl%per,mdl%nper,mdl%sper
        mdl%n = 0
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
        return
      endif
      p = mdl%per+mdl%n*mdl%sper
      mdl%n = mdl%n+1
C      read(fid) mdl%uw,mdl%cw,mdl%gw,mdl%aw
      mdl%uw(1:model%nfi,1:model%nla) = model%uw(mdl%n,:,:)
      mdl%cw(1:model%nfi,1:model%nla) = model%cw(mdl%n,:,:)
      mdl%gw(1:model%nfi,1:model%nla) = model%gw(mdl%n,:,:)
      mdl%aw(1:model%nfi,1:model%nla) = model%aw(mdl%n,:,:)
C      write(*,*) " in read 2: ",uw,cw,gw,aw
      end

