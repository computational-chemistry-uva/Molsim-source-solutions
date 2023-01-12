      Subroutine Toterg(Ener, Vir, BoxID)

C     ---Calculates Total Energy

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
       
      Double Precision Xi, Yi, Zi, Ener, Eni, Viri, Vir
      Integer I, Jb, BoxID
 
      Ener = 0.0d0
      Vir = 0.0d0
      Do I = 1, Npart - 1
         If (Id(I).Eq.BoxID) Then
            Xi = X(I)
            Yi = Y(I)
            Zi = Z(I)
            Jb = I + 1
            Call Eneri(Xi, Yi, Zi, I, Jb, Eni, Viri, BoxID)
            Ener = Ener + Eni
            Vir = Vir + Viri
         End If
      End Do

      Return
      End
