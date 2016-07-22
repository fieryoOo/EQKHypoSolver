             subroutine syn(wv,dist,dt,f0,df,sre,sim,u_int,q_int,seism,tstart,
     + vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,f2,f3,k_spec)
C------------to calculate seismogram----------------------
C------------INPUT ARGUMENTS------------------------------
C---n2pow: power of 2; dist:distance; dt: time increment;vmax: max.gr.vel;
C---f0: starting frequency; df -frequency increment; npoints-n. of output points
C-----------wv: wavenumbers; u_int: gr.vel; q_int:apparent surface wave Q;
C-----------sre, sim: spectrum--------------------------------
C------------OUTPUT: seism------------------------------------
         real*4 sre(2048),sim(2048),seism(8192),wv(2048),q_int(2048),u_int(2048)
         real*4 asre(8192),asim(8192),amp(2048),T_curr(2048), ffactor
         logical*1 key_compr
         data pi2/6.2831854/
         nbase=2**n2pow
         n=nbase/2
         tstart=dist/vmax
         ffactor = 2./dt
         do i=2,n
             f_curr=f0+df*(i-1)
             om_curr=pi2*f_curr
             arg=wv(i)*dist
C      write(*,*) "wv = ",wv(i), " dis = ", dist, " arg = ", arg
             if(key_compr)arg=om_curr*dist/fix_vel
C-----dnom: geometrical spreading----
cYT             dnom=sqrt(arg)
             dnom=sqrt(dist)
C             dnom=sqrt(dist)
C---- aq: phase delay----------------
             aq=arg+pi2/8.-tstart*om_curr
C             if(mod(i,22).eq.0.and.i.le.100) write(*,*) i," ",aq," ",arg," ",tstart," ",om_curr," ",wv(i)," ",dist
cMB          aq=pi2/8.-tstart*om_curr
C      write(*,*) "debug2: ",aq," ",arg," ",pi2," ",tstart," ",om_curr," ",wv(i)," ",dist
C-----att: attenuation factor--------
             power=dist*om_curr/u_int(i)/q_int(i)/2.
             att=0.0
             if(power.lt.20.)  att=exp(-power)
C-----full spectrum=source spectrum*propagation factor----S
             att = 1.0
             cs=cos(aq)/dnom*att
             sc=sin(-aq)/dnom*att
cYT         correct spectrum amplitude (should probably be moved into FFT.f)
             sr=ffactor*sre(i)
             si=ffactor*sim(i)
C           write(*,*) i,sre(i), sim(i), dt, n2pow
C             if(mod(i,22).eq.0.and.i.le.100) write(*,*) i," ",dnom," ",att," ",dt," ",sre(i)," ",sim(i)
C             if(mod(i,22).eq.0.and.i.le.100) write(*,*) i," ",sre(i)," ",dt," ",sr," ",si
C            sr=sre(i)/dt
C            si=sim(i)/dt
             asre(i)=sr*cs-si*sc
             asim(i)=sr*sc+si*cs
C             if(mod(i,22).eq.0.and.i.le.100) write(*,*) i," ",asre(i)," ",asim(i)," ",sr," ",si," ",sc," ",cs
C            k=nbase-i+2
C            asre(k)=asre(i)
C            asim(k)=-asim(i)
             k=n+i
             asre(k)=0.0    
             asim(k)=0.0     
         enddo
         asre(1)=0.0
         asim(1)=0.0
         asre(n+1)=0.0
         asim(n+1)=0.0
C-----full spectrum=source spectrum*propagation factor----S
C----------seismogram outputting----------S
c         k_spec=0
ccYT         do kkk=2, n
ccYT            T_c=1./(f0+df*(kkk-1))
ccYT            k_spec=k_spec+1
c         do k_spec=2, n
c            T_c=1./(f0+df*(k_spec-1))
c            if(T_c.lt.(1.0/f3).or.T_c.gt.(1.0/f2))go to 2222 
c            T_curr(k_spec)=T_c
c            amp(k_spec)=sqrt(asre(k_spec)**2+asim(k_spec)**2)/ffactor
c      write(*,*) T_c,amp(k_spec)," dist =",dist," dnom =",dnom
c2222        continue
c         enddo
         call FFT(n2pow,asre,asim,1)
         do i=1,npoints
            seism(i)=asre(i)
         end do
C----------seismogram outputting----------E
         return
      end
