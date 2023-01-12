      Subroutine Integrate_And
      Implicit None
 
      Include 'system.inc'
 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Integrate The Equations Of Motion For The Andersen Thermostat  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     RandomNumber : Uniform Random Number
C     Ran_Vel     : Generate A Velocity Form The Correct Distribution

      Double Precision U,F,RandomNumber,Ran_Vel

C Start Modification

c     integrate equations of motion: velocity verlet
      Xpos = Xpos + 0.5d0*Oldf*Tstep*Tstep + Tstep*Vpos
      Vpos = Vpos + 0.5d0*Oldf*Tstep
      
      Call Force(Xpos,U,F)

      Vpos = Vpos + 0.5d0*F*Tstep
      Cons = 0.0d0
      Oldf = F

Cccccccccccccccccccccccccc
C     Random Velocity    C
Cccccccccccccccccccccccccc

      If(RandomNumber().Lt.Nu*Tstep) Vpos = Ran_Vel()

C End Modification
 
      Return
      End
