      Subroutine Force(X,U,F)
      Implicit None
 
C     Calculate The Forces And Potential Energy
 
      Include 'system.inc'
 
      Double Precision X,U,F

      U = 0.0d0
      F = 0.0d0

      If(X.Ge.-1.0d0.And.X.Le.1.0d0) Then
         U =  0.5d0*(1.0d0 + Dcos(Onepi*X))
         F =  0.5d0*Onepi*Dsin(Onepi*X)
      Endif

      Return
      End
