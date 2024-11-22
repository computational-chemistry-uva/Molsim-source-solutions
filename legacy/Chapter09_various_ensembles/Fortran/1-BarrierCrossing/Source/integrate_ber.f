      Subroutine Integrate_ber
      Implicit None
 
      Include 'system.inc'

C     Integrate The Equations Of Motion For An Nve System
C     Use Either Velocity Verlet Or Leap-Frog. You Do Not
C     Have To Declare Any New Variables
C
C     Hint: Use The Following Symbols:
C
C     Tstep : Timestep Integration
C     Xpos  : Old Position
C     Oldf  : Old Force
C     Cons  : Conserved Quantity
C
C     To Calculate The Force And Energy For A Given Position,
C     See Force.F

      Double Precision U,F,Vnew,lambda, Temperature
      Double Precision tau, Berendsen, Stoch, Scale, Ekin
      Double Precision RandomNumber

      tau = 0.01d0

c     Call this subroutine in mdloop.f with Choice = 5
c     For the Berendsen thermostat, use the Leap Frog

      Call Force(Xpos,U,F)

c     determine lambda
      Ekin = Vpos*Vpos
c     begin modification
      Berendsen = (Temperature - Ekin)*(Tstep/tau)
      Stoch = dsqrt((Ekin * Temperature)*(Tstep/tau))
      Stoch = Stoch * RandomNumber()*2.0d0
      Scale = Berendsen + Stoch
c     end modification
      lambda = dsqrt((Ekin + Scale)/Ekin)

      Vnew = lambda*(Vpos + Tstep*F)
      Xpos = Xpos + Tstep*Vnew
      Cons = U    + 0.125d0*((Vnew + Vpos)**2)
      Vpos = Vnew
      Oldf = F

 
      Return
      End
