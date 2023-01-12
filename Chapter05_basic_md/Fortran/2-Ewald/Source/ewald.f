      Program Ewald
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Compute The Ewald Sum Of A Cubic Lattice   C
C     (For Example Nacl)                         C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Double Precision Ureal,Ufourier,Uself

      Read(21,*) Ncell,Alpha,Kmax

      If(Kmax.gt.Maxkap) Then
         Write(6,*) Kmax, ' is larger than the maximum number of wavevectors allowed :', Maxkap 
         goto 10
      Endif

      Call Lattice
      Call Realspace(Ureal)
      Call Fourierspace(Ufourier,Uself)

      Ureal    = Ureal/Dble(Npart)
      Ufourier = Ufourier/Dble(Npart)
      Uself    = Uself/Dble(Npart)

      Write(6,*) 'Real Part            : ',Ureal
      Write(6,*) 'Fourier Part         : ',Ufourier
      Write(6,*) 'Self Part            : ',Uself
      Write(6,*) 'Total                : ',Ureal+Ufourier+Uself
      Write(6,*) 'Madelung Constant    : ',-2.0d0*(Ureal+Ufourier+Uself)
      
 10   Stop
      End
