c****************************************************************
      recursive subroutine rbimod(e,m,vv,dbg,ierr, mdl)
      use mmodel
      implicit none

      type (tmdl) mdl
      integer*4 m,ierr,dbg
      real*8    e(3),vv(4)
C      common /mdl/ic,jc,n,nper,nfi,nla,fi,sfi,la,sla,per,sper,
C     +        bf,ef,bl,el,bp,ep,
C     +        hfi,hla,hper,chfi,uw,cw,gw,aw
c ---
      real*8 drad,dff,dfl
      real*8 hf,dhf,hl,dhl,v1,v2,v3,v4
c ---
      drad = datan(1.0d0)/45.0d0
      dfl = datan2(e(2),e(1))/drad
      ierr = 0
      if(dfl.lt.0.d0)dfl = dfl+360.d0
      dff = e(3)
      if(dabs(dff).gt.1.d0) dff = dsign(1.d0,e(3))
      dff = DASIN(dff)/drad
C      if(dbg.ne.0) write(*,*) 'Point: ',dff,dfl
c --- check that e is model's internal point
      if(dff.lt.mdl%bf.or.dff.gt.mdl%ef.or.dfl.lt.mdl%bl.or.dfl.gt.mdl%el) then
        goto 9
      endif
c --- check that e is cell's internal point
    2 if(dff.lt.mdl%chfi(mdl%ic)) then
        mdl%ic = mdl%ic-1
        if(mdl%ic.eq.0) then
          mdl%ic = 1
          goto 9
        endif
        goto 2
      endif
c ---
    3 if(dff.gt.mdl%chfi(mdl%ic+1)) then
        mdl%ic = mdl%ic+1
        if(mdl%ic.gt.mdl%nfi) then
          mdl%ic = mdl%nfi
          goto 9
        endif
        goto 3
      endif
c ---
    4 if(dfl.lt.mdl%hla(mdl%jc)) then
        mdl%jc = mdl%jc-1
        if(mdl%jc.eq.0) then
          mdl%jc = 1
          goto 9
        endif
        goto 4
      endif
c ---
    5 if(dfl.gt.mdl%hla(mdl%jc+1)) then
        mdl%jc = mdl%jc+1
        if(mdl%jc.gt.mdl%nla) then
          mdl%jc = mdl%nla
          goto 9
        endif
        goto 5
      endif
c --- four maps  interpolation ---
      m = 0
      hf = mdl%chfi(mdl%ic+1)-mdl%chfi(mdl%ic)
      dhf = dff-mdl%chfi(mdl%ic)
      hl = mdl%hla(mdl%jc+1)-mdl%hla(mdl%jc)
      dhl = dfl-mdl%hla(mdl%jc)
      v1 = mdl%uw(mdl%ic,mdl%jc)
      v2 = mdl%uw(mdl%ic+1,mdl%jc)
      v3 = mdl%uw(mdl%ic,mdl%jc+1)
      v4 = mdl%uw(mdl%ic+1,mdl%jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(1) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      v1 = mdl%cw(mdl%ic,mdl%jc)
      v2 = mdl%cw(mdl%ic+1,mdl%jc)
      v3 = mdl%cw(mdl%ic,mdl%jc+1)
      v4 = mdl%cw(mdl%ic+1,mdl%jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(2) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      v1 = mdl%gw(mdl%ic,mdl%jc)
      v2 = mdl%gw(mdl%ic+1,mdl%jc)
      v3 = mdl%gw(mdl%ic,mdl%jc+1)
      v4 = mdl%gw(mdl%ic+1,mdl%jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(3) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      v1 = mdl%aw(mdl%ic,mdl%jc)
      v2 = mdl%aw(mdl%ic+1,mdl%jc)
      v3 = mdl%aw(mdl%ic,mdl%jc+1)
      v4 = mdl%aw(mdl%ic+1,mdl%jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(4) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
C      write(*,*) "rbimod: ",vv(4)
      m = m+1
c ---
      return
   9  ierr = 1
      end
