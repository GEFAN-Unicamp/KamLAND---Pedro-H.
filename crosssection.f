c.. cross section from Vogel and Beacom, hep-ph/9903554
      implicit real*4 (a-h,o-z)
      dimension sigint(0:100)
      dimension enuenu(0:1000),Ee1ee1(0:1000,0:100),sigsig(0:1000,0:100)

      open(15,file='crosskaml.dat',status='unknown')
      open(16,file='crosskaml2.dat',status='unknown')

      amn=939.56536  !1.0086649*931.494
      amp=938.27203  !1.0072765*931.494
      ame=0.511

      delta = amn-amp

      f=1.
      g=1.26
      f2=3.706

      Gf=1.
      pi=acos(-1.0)
      cos2tc=0.974
      deltarin=0.024
      sig0=Gf*cos2tc/pi*(1.+deltarin)

      amm = (amn+amp)/2.

      do 22 ienu=0,500
      Enu=1.81+8.19*float(ienu)/500.

      Ee0=Enu-delta
      if(Ee0.lt.ame)then
         do jeel=0,100
            sigint(jeel)=0.
         enddo
         goto 22
      endif
      pe0=sqrt(Ee0**2.-ame**2.)
      ve0=pe0/Ee0

      y2 = (delta**2.-ame**2.)/2.

      Ee1min=Ee0*(1.-Enu/amm*(1.+ve0))-y2/amm
      Ee1max=Ee0*(1.-Enu/amm*(1.-ve0))-y2/amm

c      costmin = (1.+(amm/Enu)*((Ee1min + y2/amm)/Ee0-1.))/ve0
c      costmax = (1.+(amm/Enu)*((Ee1max + y2/amm)/Ee0-1.))/ve0

      do 21 jeel=0,100
c      Ee1=2.6+7.4*float(jeel)/100.
      Ee1=Ee1min+(Ee1max-Ee1min)*float(jeel)/100.
      pe1=sqrt(Ee1**2.-ame**2.)
      ve1=pe1/Ee1


c.. hep-ph/9903554 - eq.(13)
      y2 = (delta**2.-ame**2.)/2.
      cost = (1.+(amm/Enu)*((Ee1 + y2/amm)/Ee0-1.))/ve0


c.. hep-ph/9903554 - eq.(5) and (6)
c      cost=(amn**2.-amp**2.-2.*amp*(Enu-Ee1)-ame**2.+2.*Enu*Ee1)
c     %     /(2.*Enu*Ee1*ve1)

      gamma=2.*(f+f2)*g*((2.*Ee0+delta)*(1.-ve0*cost)-ame**2./Ee0)
     %      +(f**2.+g**2.)*(delta*(1.+ve0*cost)+ame**2./Ee0)
     %      +(f**2.+3.*g**2)*((Ee0+delta)*(1.-cost/ve0)-delta)
     %      +(f**2.-g**2.)*((Ee0+delta)*(1.-cost/ve0)-delta)*ve0*cost

      sig=sig0/2.*((f**2.+3.*g**2.)+(f**2.-g**2.)*ve1*cost)*Ee1*pe1
     %    -sig0/2.*gamma/amm*Ee0*pe0

      sigint(jeel)=amm/(Ee0*Enu*ve0)*sig

      write(15,*)Enu,Ee1,amm/(Ee0*Enu*ve0)*sig,cost
c      enuenu(ienu)=Enu
c      Ee1ee1(ienu,jeel)=Ee1
c      sigsig(ienu,jeel)=sigint(jeel)
 21   continue
      sig2=simps(sigint,Ee1min,Ee1max,100)

      write(16,*)Enu,sig2

      cost=1.
      sig1a=sig0/2.*((f**2.+3.*g**2.)+(f**2.-g**2.)*ve0*cost)*Ee0*pe0
      cost=-1.
      sig1b=sig0/2.*((f**2.+3.*g**2.)+(f**2.-g**2.)*ve0*cost)*Ee0*pe0
      sig1c=sig0*(f**2.+3.*g**2.)*Ee0*pe0
      write(97,*)Enu,sig1c,sig1a,sig1b

      write(15,*)
 22   continue

      end

