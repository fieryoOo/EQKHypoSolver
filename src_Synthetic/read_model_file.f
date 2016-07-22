c ==========================================================
c read rectangular models into memory
c ==========================================================
      recursive subroutine read_model_file(fname,model)
      use mmodel
      implicit none
      integer*4 iper
      character*255 fname

      type (tmodel) model

c --- read unformated file for period per--
C      write(*,*) fname
      open(30,file=fname,form='unformatted')
      read(30) model%n,model%fi,model%nfi,model%sfi,model%la,model%nla,model%sla,model%per,model%nper,model%sper
      if( model%nfi.gt.nlatmax.OR.model%nla.gt.nlonmax.OR.model%nper.gt.npermax ) then
         STOP 'model dimension exceeds limits'
      endif
c      write(*,*) model%nfi, model%nla, model%nper

      allocate (model%uw(model%nper,model%nfi,model%nla))
      allocate (model%cw(model%nper,model%nfi,model%nla))
      allocate (model%gw(model%nper,model%nfi,model%nla))
      allocate (model%aw(model%nper,model%nfi,model%nla))

C      write(*,*) " in read 1: ",model%n,model%fi,model%nfi,model%sfi,model%la,model%nla,model%sla,model%per,model%nper,model%sper
      do iper=1,model%nper
         read(30) model%uw(iper,:,:),model%cw(iper,:,:),model%gw(iper,:,:),model%aw(iper,:,:)
C         write(*,*) model%uw(iper,:,:),model%cw(iper,:,:),model%gw(iper,:,:),model%aw(iper,:,:)
      enddo
      close(30)
      end

