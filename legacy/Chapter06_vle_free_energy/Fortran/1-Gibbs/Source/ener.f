      Subroutine Ener(En, Vir, R2, BoxID)

C  ---Calculates Energy (En) And Virial (Vir) For Given
C     Distance Squared Between (R2) Two Particles
 
      Implicit None
      Double Precision R2, R2i, R6i, En, Vir
      Integer BoxID
      
      Include 'potential.inc'
 
      If (R2.Le.Rc2(BoxID)) Then
         R2i = Sig2/R2
         R6i = R2i*R2i*R2i
         En  = Eps4*(R6i*R6i-R6i)
         Vir = Eps48*(R6i*R6i-0.5d0*R6i)
      Else
         En  = 0.D0
         Vir = 0.D0
      End If
      
      Return
      End
