      Subroutine Integrate_Mc
      Implicit None
 
      Include 'system.inc'
 
Cccccccccccccccccccccc
C     Mc Simulation  C
Cccccccccccccccccccccc

      Double Precision Uold,Unew,F,Xnew,RandomNumber
 
      Call Force(Xpos,Uold,F)

      Xnew = Xpos + 0.1d0*(RandomNumber()-0.5d0)

      Call Force(Xnew,Unew,F)

      If(RandomNumber().Lt.Dexp(-(Unew-Uold)/Temp)) Xpos = Xnew

      Oldf = 0.0d0
      Vpos = 0.0d0
      Cons = 0.0d0
 
      Return
      End
