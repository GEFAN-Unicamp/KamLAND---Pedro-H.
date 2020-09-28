c. utiliza arquivo crosskaml.dat, que fornece a secao de choque de KamLAND, 
c  criado pelo programa crosssection.f
c.
c. cria o arquivo kamlspec.dat, com a secao de choque efetiva para cada bin.
c.
      implicit real*4 (a-h,o-z)
      parameter (nenu=501)
      dimension t(101),sigel(101),sigmu(101)
      dimension tr(0:250),tt(0:250),ssigel(0:250)
      dimension ssig2el(0:250),ssig2mu(0:250)
      dimension eekl(nenu),ssklel(nenu,20)
      dimension sigenu(20)
      dimension crosstmp(0:100)
      dimension ef(501),fluxo(501)
      dimension rtt(20)
      dimension alovere(21)
      data alovere/15.6597,24.6363,27.6942,29.4698,31.3440,
     % 33.1196,35.0925,36.7694,38.8409,41.1097,43.9704,47.0284,
     % 50.3822,53.9334,57.6819,61.7263,66.1652,71.7879,
     % 78.8903,86.6831,110.000/
      dimension rtot1(0:500),rtot2(20)
      dimension effic(9),efficenu(9)
      data efficenu/.95,1.2,1.3,1.4,1.6,2.,2.8,3.,10./
      data effic/.45,.7,.76,.8,.79,.85,.85,.92,.92/

c      do i=1,20
c         write(*,*)alovere(i+1)-alovere(i)
c      enddo
c      stop

      open(10,file='crosskaml.dat',status='old')
      open(11,file='kamlspec.le.dat',status='unknown')

      open(15,file='fluxo.dat',status='old')
      do i=1,501
         read(15,*)ef(i),fluxo(i)
      enddo
      close(15)

      amn=939.56536  !1.0086649*931.494
      amp=938.27203  !1.0072765*931.494
      delta = amn-amp!1.29333
      pi=acos(-1.0)

      fator=283.91275

c      do ienu=1,nenu
c         read(11,*)eekl(ienu),(ssklel(ienu,ispec),ispec=1,20)
c      enddo
c      goto 14
      
c........................
      do 13 ienu=1,nenu   ! numero de valores em energia em crosskaml.dat, 

      write(*,*)ienu

      do i=1,101     !numero de pontos em cada arquivo
         read(10,*)enu,eel,sigel(i)
         t(i)=eel
         crosstmp(i-1)=sigel(i)
      enddo
      read(10,*)

      tmin = t(1)
      tmax = t(101)

c....................... mudando os limites de integracao (energia visivel)
      do 12 ispec=1,20
      enu1=180./alovere(22-ispec)
      enu2=180./alovere(21-ispec)
      evih=enu1-delta+0.511
      evim=enu2-delta+0.511

c..................... integracao em T medido
      tvih=max(evih-0.511,0.1*t(1))
      tvim=min(evim-0.511,5.*t(101))

      if(tvih.ge.tvim)then
         sigenuel=0.
         goto 21
      endif

c........................................
      do 11 i=0,250
         tt(i)=tvih+(tvim-tvih)*float(i)/250.

c..................... integracao em T real
            do j=0,250
            tr(j)=tmin+(tmax-tmin)*float(j)/250.
            sige=0.065*sqrt(tr(j))
            resol=exp(-0.5*((tt(i)-tr(j))/sige)**2.)
     %            /sqrt(2.*pi)/sige
            ttt = tr(j)
            ssigel(j)=resol*divdif(sigel,t,101,ttt,2)
            enddo
c........................................

         ssig2el(i)=simps(ssigel,tmin,tmax,250)
         energy=tt(i)+0.511
         efficiency=divdif(effic,efficenu,9,energy,1)
         ssig2el(i)=efficiency*ssig2el(i)

 11      continue
c.................. fim do loop na energia visivel

         sigenuel=simps(ssig2el,tvih,tvim,250)
 21      continue
c........................................


      eekl(ienu)=enu
      ssklel(ienu,21-ispec)=fator*sigenuel
 12   continue
c.................. fim do loop em bins espectrais
      write(11,111)enu,(ssklel(ienu,ispec),ispec=1,20)

 13   continue
c................... fim do loop na energia do neutrino


c.. teste de numero de eventos esperado
 14   continue
      stop
      
      rtot3=0.
      e1=eekl(1)
      e2=eekl(nenu)
      do ispec=1,20
      do ienu=1,501
         e = eekl(ienu)
         ff= divdif(fluxo,ef,501,e,2)
         ssklel(ienu,ispec)=ff*ssklel(ienu,ispec)
         rtot1(ienu-1)=ssklel(ienu,ispec)
      enddo
      rtot2(ispec)=simps(rtot1,e1,e2,500)
      rtot3=rtot3+rtot2(ispec)
      enddo
      write(*,*)
      write(*,*)rtot2
      write(*,*)
      write(*,*)rtot3,2879./rtot3

      do ispec=1,20
         write(98,*)alovere(ispec),rtot2(ispec)
         write(98,*)alovere(ispec+1),rtot2(ispec)
         e1=180./alovere(ispec+1)-(delta-0.511)
         e2=180./alovere(ispec)-(delta-0.511)
         write(97,*)e2,rtot2(ispec)*0.425/(e2-e1) !para comparar com grafico
         write(97,*)e1,rtot2(ispec)*0.425/(e2-e1) !para comparar com grafico
      enddo

 111  format(f6.3,1x,20(e12.5,1x))


      end



