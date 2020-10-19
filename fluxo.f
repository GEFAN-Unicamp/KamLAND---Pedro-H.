      implicit real*4 (a-h,o-z)
      dimension fr_iso(4)
      dimension a0(4),a1(4),a2(4)
      dimension enu(0:500)
      data fr_iso/.570,.295,.078,.057/ !provided by kamland
      data a0/0.870,0.896,0.976,0.793/
      data a1/-0.160,-0.239,-0.162,-0.080/
      data a2/-0.0910,-0.0981,-0.0790,-0.1085/
      dimension contr(17) !contribuicoes percentuais de cada reator
      data contr/22.4400158,18.3808117,12.1540928,
     %  4.59471512,10.124464,6.19334984,10.1309271,3.85063148,
     %  3.57945991,1.21683121,0.711448312,0.611766875,0.377173662,
     %  0.711553276,0.172282934,0.250485122,4.500000/

      open(11,file='fluxo.dat',status='unknown')

      do ienu=0,500
      enu(ienu)=1.0+9.0*float(ienu)/500.
      stmp=0.
      do j=1,4
      stmp=stmp+fr_iso(j)*exp(a0(j)+a1(j)*enu(ienu)+a2(j)*enu(ienu)**2.)
      enddo

      events=0.
      do ireactor=1,17
         events=events+1.e-2*contr(ireactor)*stmp
      enddo

      write(11,*)enu(ienu),events 
      enddo

      end
