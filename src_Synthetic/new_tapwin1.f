C###############################################################
        SUBROUTINE tapwin1(nf, nw1,nw2,n_beg,n_end,iq,s)
C       cosine unsymmetrical taper   
        parameter (nsize=8192)
        real*4 s(nsize)
        data pi/3.14159265358979323846/
C ???        do i=1,nl1
        do i=1,n_beg-1
        s(i)=0.0 
        enddo
        np=n_beg+nw1+1
        nq=n_end-nw2
        DO i=n_beg,np-1
        j=i+1-n_beg
        s(i)=(0.5*(1- cos(pi*float(j-1)/float(nw1))))**iq
        end do
        DO i=np,nq
        s(i)=1.0
        end do
        do i=nq+1,n_end
        j=i-nq
        s(i)=(0.5*(1+ cos(pi*float(j)/float(nw2))))**iq
        end do
        do i=n_end+1,nf
        s(i)=0.0
        end do
        return
        end
