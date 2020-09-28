c..   calcula producao espectral de KamLAND, dados detalhes do sitio,
c..   e chi^2 para esta producao. Ver tambem cs/kamlanspec.f
c..
c..   utiliza "kamlspec.le.dat" 
c..
      implicit real*4 (a-h,o-z)
c.............. effective cross-sections and rates ..............
      dimension enu(0:500),specfunc(20,0:500),sbins(20)
      common /crossec/enu,specfunc
c.................. reactors ....................................
      dimension dist(22)
      data dist/160.,179.,191.,088.,138.,214.,146.,349.,351.,141.,295.,
     %          138.,401.,431.,561.,754.,830.,783.,712.,735.,709.,986./      
c.................. flux calculation ............................
      dimension fr_iso(4)
      dimension a0(4),a1(4),a2(4)
      data fr_iso/.570,.295,.078,.057/ !provided by kamland
      data a0/0.870,0.896,0.976,0.793/
      data a1/-0.160,-0.239,-0.162,-0.080/
      data a2/-0.0910,-0.0981,-0.0790,-0.1085/
      dimension contr(22) !contribuicoes percentuais de cada reator
      data contr/30.9, 13.8, 9.0, 7.9, 7.6, 7.5, 7.4, 3.8, 3.5, 1.3, 
     %   1.2, 0.9, 0.8, 0.7, 0.6, 0.4, 0.2, 0.2, 0.8, 0.6, 0.5, 0.5/
      common /flux/dist,fr_iso,a0,a1,a2,contr
c................................................................
      dimension chi2tgdm(0:200,0:80),chi2t12t13(0:200,0:40)
      dimension chi2dmt13(0:80,0:40)
      
      open(15,file='kamlspec.le.dat',status='old')
      do ienu=0,500
         read(15,*)enu(ienu),(specfunc(ispec,ienu),ispec=1,20)
      enddo
      close(15)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(21,file='sbins.tmp',status='unknown')
      open(22,file='sbins.dat',status='unknown')

      open(30,file='chi2pois.dat',status='unknown')
      chi2min=1.e6

      do is13=0,40
         do itg=0,200
            chi2t12t13(itg,is13)=1.e5
         enddo
      enddo
      do is13=0,40
         do idm=0,80
            chi2dmt13(idm,is13)=1.e5
         enddo
      enddo
      do itg=0,200
         do idm=0,80
            chi2tgdm(itg,idm)=1.e5
         enddo
      enddo
      
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do is13=0,40
      s2_2t13=0.4*float(is13)/20.!9.3e-2
      write(*,*)is13

c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do itg=0,200
      tg2_t=10.**(log10(0.1)+(log10(10.)-log10(0.1))*float(itg)/200.)
c      tg2_t=0.50
         
c----------------------------------------------------------------------
c----------------------------------------------------------------------
         do idm=0,80
            deltam=(6.+4.*float(idm)/80.)*1.e-5
c            deltam=8.e-5

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     goto 81  !comentar esta linha se sbins nao estiver calculado
            call nuspec(tg2_t,deltam,s2_2t13,sbins)
            write(21,121)tg2_t,deltam,(sbins(k),k=1,20)
 121        format(22e14.6)
            goto 82
 81         continue
            read(22,*)a,b,(sbins(k),k=1,20)
 82         continue
c----------------------------------------------------------------------
c----------------------------------------------------------------------
            call chi2pois(sbins,chi2,ff)
            write(30,*)s2_2t13,tg2_t,deltam,chi2,ff
            if(chi2.lt.chi2min)then
               chi2min=chi2
               dmmin=deltam
               tgmin=tg2_t
               fmin=ff
            endif
            if(chi2.lt.chi2tgdm(itg,idm))chi2tgdm(itg,idm)=chi2
            if(chi2.lt.chi2t12t13(itg,is13))chi2t12t13(itg,is13)=chi2
            if(chi2.lt.chi2dmt13(idm,is13))chi2dmt13(idm,is13)=chi2
c----------------------------------------------------------------------
c----------------------------------------------------------------------

 21         continue

         enddo
         write(30,*)
      enddo
      enddo
      
      open(31,file='chi2tgdm.dat',status='unknown')
      open(32,file='chi2t12t13.dat',status='unknown')
      open(33,file='chi2dmt13.dat',status='unknown')

      write(*,*)chi2min,dmmin,tgmin
      do is13=0,40
         s2_2t13=0.4*float(is13)/20.
         do itg=0,200
         tg2_t=10.**(log10(0.1)+(log10(10.)-log10(0.1))*float(itg)/200.)
            write(32,*)tg2_t,s2_2t13,chi2t12t13(itg,is13)
         enddo
         write(32,*)
      enddo
      do is13=0,40
         s2_2t13=0.4*float(is13)/20.
         do idm=0,80
            deltam=(6.+4.*float(idm)/80.)*1.e-5
            write(33,*)deltam,s2_2t13,chi2dmt13(idm,is13)
         enddo
         write(33,*)
      enddo
      do itg=0,200
         tg2_t=10.**(log10(0.1)+(log10(10.)-log10(0.1))*float(itg)/200.)
         do idm=0,80
            deltam=(6.+4.*float(idm)/80.)*1.e-5
            write(31,*)tg2_t,deltam,chi2tgdm(itg,idm)
         enddo
         write(31,*)
      enddo

      end

