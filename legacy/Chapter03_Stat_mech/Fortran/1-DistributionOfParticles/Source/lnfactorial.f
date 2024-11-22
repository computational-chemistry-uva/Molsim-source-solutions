      Function LnFactorial(N)
      Implicit None

      Double Precision LnFactorial
      Integer N,J

      LnFactorial = 0.0d0

      Do J=2,N
         LnFactorial = LnFactorial + Dlog(Dble(J))
      Enddo
      
      Return
      End
