      subroutine nuspec(tg2_t,deltam,s2_2t13,sbins)
      implicit real*4 (a-h,o-z)
c.............. effective cross-sections and rates ..............
      dimension enu(0:500),specfunc(20,0:500),sbins(20)
      dimension spectmp1(20,0:500),spectmp2(0:500)
      common /crossec/enu,specfunc
c..................reactors ....................................
      dimension dist(22)
      dimension fr_iso(4)
      dimension a0(4),a1(4),a2(4)
      dimension contr(22) 
      common /flux/dist,fr_iso,a0,a1,a2,contr
c................................................................

      pi = acos(-1.0)
      gf=1.16632e-23
      ve=0.5*sqrt(2.)*gf*6.02e23/((5.068e4)**3.0) * 1.5 !em g/cm**3
      c2_2t13=1-s2_2t13
      s2_13=(1-sqrt(c2_2t13))/2
      c2_13=(1+sqrt(c2_2t13))/2
      s2_2t=4*tg2_t/((1.+tg2_t)**2.)
      s_2t=sqrt(s2_2t)
      c_2t=sqrt(1.-s2_2t)

      do ienu=0,500
         enn=enu(ienu)
c..
         flux=0.
         do j=1,4
            flux=flux+fr_iso(j)*exp(a0(j)+a1(j)*enn+a2(j)*enn**2.)
         enddo
c..
         dm4e=deltam/4.e6/enn
         t_2m=dm4e*s_2t/(dm4e*c_2t-ve*c2_13)
         s2_2m=t_2m**2./(1.+t_2m**2.)
         pmed=0.
         do ireactor=1,22
            pmed=pmed+1.e-2*contr(ireactor)*(c2_13**2.*(1.-s2_2m*
     %         (sin(1.267e3*deltam/enn*dist(ireactor)))**2.)+s2_13**2.)
         enddo
         do ispec=1,20
            spectmp1(ispec,ienu)=flux*pmed*specfunc(ispec,ienu)
         enddo         
      enddo
      do ispec=1,20
         do ienu=0,500
            spectmp2(ienu)=spectmp1(ispec,ienu)
         enddo
         sbins(ispec)=simps(spectmp2,enu(0),enu(500),500)
      enddo
      
      return
      end
