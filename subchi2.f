c..   ultima modificacao: 06/2013
c..   calcula chi2 com espectro de KamLAND
c..   usando estatistica de Poisson
c..
      subroutine chi2pois(sbins,chi2,ff)
      implicit real*4 (a-h,o-z)
      dimension sbins(20),rkl(20),datakl(20)
      dimension datasm(20),datanu(20),alovere(21)
      common /alovere/alovere

      data datasm/59.547829, 105.21388, 91.706429, 119.44465, 131.14526, 
     % 161.02548, 145.14534, 184.27446, 201.80777, 245.68504, 243.55714, 
     % 230.59836, 200.76579, 176.71275, 153.93623, 131.16353, 117.47219, 
     % 92.720367, 57.012432, 30.065214/

      data datanu/0.395866,0.667727,0.629571,0.586645,0.804452,
     % 0.785374,0.732909,0.558029,0.558029,0.551669,0.489666,
     % 0.387917,0.467409,0.491256,0.585056,0.497615,0.575517,
     % 0.658188,0.612083,0.279809/

      do i=1,20
         datakl(i)=datanu(i)*datasm(i)
         rkl(i)=sbins(i)
      enddo

      chi2ff=1.e6
      do iff=1,100
         ff=0.75+0.5*float(iff)/100.
         chi2=0.

         do i=1,20
         chi2log=0.
         if(datakl(i).ne.0.)chi2log=datakl(i)*log(datakl(i)/(ff*rkl(i)))
         chi2=chi2+2.*(ff*rkl(i)-datakl(i)+chi2log)
         enddo
         chi2=chi2+((1.-ff)/0.041)**2.
         
         if(chi2.lt.chi2ff)then
         chi2ff=chi2
         fmin=ff
         endif
      enddo
      chi2=chi2ff
      ff=fmin

      return
      end


