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
      open(30,file=fname,form='unformatted')
      read(30) model%n,model%fi,model%nfi,model%sfi,model%la,model%nla,model%sla,model%per,model%nper,model%sper
      if( model%nfi.ne.Cnlat.OR.model%nla.ne.Cnlon.OR.model%nper.ne.Cnper ) then
         STOP 'incompatible model dimension'
      endif
C      write(*,*) " in read 1: ",model%n,model%fi,model%nfi,model%sfi,model%la,model%nla,model%sla,model%per,model%nper,model%sper
      do iper=1,225
         read(30) model%uw(iper,:,:),model%cw(iper,:,:),model%gw(iper,:,:),model%aw(iper,:,:)
C         write(*,*) model%uw(iper,:,:),model%cw(iper,:,:),model%gw(iper,:,:),model%aw(iper,:,:)
      enddo
      close(30)
      end

