c group velocity tracer for cross-correlation paths
c ---
      module mmodel
      integer*4 Cnper, Cnlat, Cnlon
      parameter (Cnper=225)
      parameter (Cnlat=97)
      parameter (Cnlon=201)
      type tmdl
C         integer, parameter :: Cnper = 225
C         integer, parameter :: Cnlat = 97
C         integer, parameter :: Cnlon = 201
         integer*4 ic,jc,n,nper,nfi,nla
         real*8 fi,sfi,la,sla,per,sper,bf,ef,bl,el,bp,ep,
     +          hfi(Cnlat),hla(Cnlon),hper(Cnper),chfi(Cnlat)
         real*8 uw(Cnlat,Cnlon),cw(Cnlat,Cnlon)
         real*8 gw(Cnlat,Cnlon),aw(Cnlat,Cnlon)
      end type

      type tmodel
         integer*4 n, nper, nfi, nla
         real*8 fi, sfi, la, sla, per, sper
         real*8, allocatable, dimension(:,:,:) :: uw,cw,gw,aw
C         real*8 uw(225,97,201), cw(225,97,201)
C         real*8 gw(225,97,201), aw(225,97,201)
      end type

      end module mmodel

      recursive subroutine atracer(fmodel,fsol,lsol,n,applyQ, nstai,fici,lami, cor)
      use mmodel
      use omp_lib
      implicit none

c ---
C      character*8 codi(2000)
      integer*4   nstai
      real*4      fici(2000),lami(2000)
C      common /stn/ nstai,codi,figi,fici,lami
C      common /mdl/ic,jc,nn,nper,nfi,nla,fi,sfi,la,sla,pper,sper,
C     +        bf,ef,bl,el,bp,ep,
C     +        hfi,hla,hper,chfi,uw,cw,gw,aw
c ---
      real*4 cor(500,2,2000)
C      common /trk/ cor
c ---
      logical*1 applyQ
      real*4 fsol,lsol
      real*8 afi,del,per
      real *8 GEO,rad,pi2,sol(3),dst(3),trres(4),sine,cosi
      integer*4 ierr,k,ntr,n,m,lnblnk
c     real*8     cor(500,2,2000)
      character *256 fmodel
C      data GEO/1.0/
      data GEO/0.993277d0/

      type (tmdl) mdl
      type (tmodel) model
      allocate (model%uw(Cnper,Cnlat,Cnlon))
      allocate (model%cw(Cnper,Cnlat,Cnlon))
      allocate (model%gw(Cnper,Cnlat,Cnlon))
      allocate (model%aw(Cnper,Cnlat,Cnlon))
c ---
      rad = datan(1.0d0)/45.0d0
      pi2 = datan(1.0d0)*8.0d0
C      write(*,*) fmodel(1:lnblnk(fmodel))
c --- get event coordinates ----
C      fsol = datan(GEO*dtan(rad*fsol))/rad
      afi = datan(GEO*dtan(rad*fsol))/rad
      sol(1)=DSIN((90.0d0-afi)*rad)*DCOS(lsol*rad)
      sol(2)=DSIN((90.0d0-afi)*rad)*DSIN(lsol*rad)
      sol(3)=DCOS((90.0d0-afi)*rad)
c --- MAIN LOOP ------------
      ntr = 4

!$OMP CRITICAL (IO)
      call read_model_file(fmodel, model)
!$OMP END CRITICAL (IO)

      call read_rect_model(model,0,per,ierr, mdl)
C      n = 0
      do k = 1,Cnper
        call read_rect_model(model,1,per,ierr, mdl)
        do m = 1,nstai
c     write(*,*) per,ierr
c         afi = datan(GEO*dtan(rad*figi(m)))/rad
          afi = fici(m)
          dst(1)=DSIN((90.0d0-afi)*rad)*DCOS(lami(m)*rad)
          dst(2)=DSIN((90.0d0-afi)*rad)*DSIN(lami(m)*rad)
          dst(3)=DCOS((90.0d0-afi)*rad)
          call tracer(sol,dst,.01d0,ntr,trres,del,0,ierr, mdl)
cMB    write(*,*) 'X ',per,trres,del
          if(ierr.eq.0) then  
cMB     trres(2) = trres(2)*pi2/per
cMB     cor(k) = dcmplx(dcos(trres(2)),dsin(trres(2)))*dexp(-trres(3))*trres(4)
cMB     write(*,*) 'A ', per,6371.0*del,trres(2)
            cor(k,1,m) = 6371.0*del/trres(2)
            if( applyQ ) then
               cor(k,2,m) = dexp(-trres(3))*trres(4)
            else
               cor(k,2,m) = trres(4)
            endif
CCC     write(*,*) m,per, 6371.0*del,cor(k,1,m)
c       write(*,*) per,trres,cor(k)
          else
            cor(k,1,m) = -1
            cor(k,2,m) = 0
          endif
        enddo
      enddo
C      n = 225
      call read_rect_model(model,-1,per,ierr, mdl)
      deallocate (model%uw)
      deallocate (model%cw)
      deallocate (model%gw)
      deallocate (model%aw)
C      write(*,*) "atracer done: ",cor(1,1,1)," ",cor(2,1,2)," ",cor(3,2,3)," ",cor(300,1,999)
      end
