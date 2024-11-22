      Subroutine Integrate_Nve
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

      Double Precision U,F,Vnew

C Start Modification
c     integration of equations of motion: Leap Frog

      Call Force(Xpos,U,F)

      Vnew = Vpos + Tstep*F
      Xpos = Xpos + Tstep*Vnew

c     for conserved quantity remember v(t) = (vpos + vnew)/2
c     this explains the factor 0.125d0
      Cons = U    + 0.125d0*((Vnew + Vpos)**2)
      Vpos = Vnew
      Oldf = F

C End Modification
 
      Return
      End
