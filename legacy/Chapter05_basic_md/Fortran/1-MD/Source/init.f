      Subroutine Init
      Implicit None
 
      Include 'system.inc'
 
C     Generates Initial Positions/Velocities
C     This Is Not To Easy; Do Not Look For Errors Here !!!!
 
      Integer I,J,K,Number,Nplace
 
      Double Precision RandomNumber,Rangauss,Fxo(MaxParticles),
     &                 Fyo(MaxParticles),Fzo(MaxParticles),Test,
     &                 Uold,Scale,Place,Size,MomentumX,MomentumY,MomentumZ
 
C     Generate Velocities From A Gaussian; Set Impulse To Zero
 
      MomentumX = 0.0d0
      MomentumY = 0.0d0
      MomentumZ = 0.0d0
      KineticEnergy = 0.0d0
 
      Do I = 1,NumberOfParticles
         Vxx(I) = Rangauss()
         Vyy(I) = Rangauss()
         Vzz(I) = Rangauss()
         MomentumX = MomentumX + Vxx(I)
         MomentumY = MomentumY + Vyy(I)
         MomentumZ = MomentumZ + Vzz(I)
      Enddo
 
      MomentumX = MomentumX/Dble(NumberOfParticles)
      MomentumY = MomentumY/Dble(NumberOfParticles)
      MomentumZ = MomentumZ/Dble(NumberOfParticles)

C     Calculate The Kinetic Energy
 
      Do I = 1,NumberOfParticles
         Vxx(I) = Vxx(I) - MomentumX
         Vyy(I) = Vyy(I) - MomentumY
         Vzz(I) = Vzz(I) - MomentumZ
 
         KineticEnergy = KineticEnergy + Vxx(I)*Vxx(I) + 
     &        Vyy(I)*Vyy(I) + Vzz(I)*Vzz(I)
      Enddo
 
C     Scale All Velocities To The Correct Temperature
 
      Scale = Dsqrt(Temperature*Dble(3*NumberOfParticles-3)/KineticEnergy)
 
      Do I = 1,NumberOfParticles
         Vxx(I) = Scale*Vxx(I)
         Vyy(I) = Scale*Vyy(I)
         Vzz(I) = Scale*Vzz(I)
      Enddo
 
C     Calculate Initial Positions On A Lattice
 
      Number = Int((Dble(NumberOfParticles)**(1.0d0/3.0d0)) + 1.5d0)
      Nplace = 0
 
      Size = Box/Dble(Number + 2)
 
      Place = 0.2d0*Size
 
      Do I = 1,Number
         Do J = 1,Number
            Do K = 1,Number
               Nplace = Nplace + 1
               If (Nplace.Le.NumberOfParticles) Then
                  Rxx(Nplace) = (Dble(I) + 
     &                 0.01d0*(RandomNumber()-0.5d0))*Size
                  Ryy(Nplace) = (Dble(J) + 
     &                 0.01d0*(RandomNumber()-0.5d0))*Size
                  Rzz(Nplace) = (Dble(K) + 
     &                 0.01d0*(RandomNumber()-0.5d0))*Size
               Endif
            Enddo
         Enddo
      Enddo
 
C     Calculate Better Positions Using A Steepest Decent Algorithm
 
      Do J = 1,50
 
         If (J.Eq.1) Then
            Call Force
            Uold = PotentialEnergy
            Write (6,*) 'Initial Energy        : ',Uold
         Endif
 
         Test = 0.0d0
 
C     Calculate Maximum Downhill Gradient
 
         Do I = 1,NumberOfParticles
            Rxf(I) = Rxx(I)
            Ryf(I) = Ryy(I)
            Rzf(I) = Rzz(I)
 
            Fxo(I) = Fxx(I)
            Fyo(I) = Fyy(I)
            Fzo(I) = Fzz(I)
 
            MomentumX = Dabs(Fxx(I))
            MomentumY = Dabs(Fyy(I))
            MomentumZ = Dabs(Fzz(I))
 
            If (MomentumX.Gt.Test) Test = MomentumX
            If (MomentumY.Gt.Test) Test = MomentumY
            If (MomentumZ.Gt.Test) Test = MomentumZ
         Enddo
 
         Test = Place/Test

C     Calculate Improved Positions
 
         Do I = 1,NumberOfParticles
 
            Rxx(I) = Rxx(I) + Test*Fxx(I)
            Ryy(I) = Ryy(I) + Test*Fyy(I)
            Rzz(I) = Rzz(I) + Test*Fzz(I)
 
C     Place Particles Back In The Box
 
            If (Rxx(I).Gt.Box) Then
               Rxx(I) = Rxx(I) - Box
            Elseif (Rxx(I).Lt.0.0d0) Then
               Rxx(I) = Rxx(I) + Box
            Endif
 
            If (Ryy(I).Gt.Box) Then
               Ryy(I) = Ryy(I) - Box
            Elseif (Ryy(I).Lt.0.0d0) Then
               Ryy(I) = Ryy(I) + Box
            Endif
 
            If (Rzz(I).Gt.Box) Then
               Rzz(I) = Rzz(I) - Box
            Elseif (Rzz(I).Lt.0.0d0) Then
               Rzz(I) = Rzz(I) + Box
            Endif
         Enddo
 
C     Calculate New Potential Energy
 
         Call Force
 
C     Check If The New Positions Are Acceptable
 
         If (PotentialEnergy.Lt.Uold) Then
            Uold  = PotentialEnergy
            Place = Place*1.2d0
 
            If (Place.Gt.Halfbox) Place = HalfBox
 
         Else
            Do I = 1,NumberOfParticles
               Fxx(I) = Fxo(I)
               Fyy(I) = Fyo(I)
               Fzz(I) = Fzo(I)
               Rxx(I) = Rxf(I)
               Ryy(I) = Ryf(I)
               Rzz(I) = Rzf(I)
            Enddo
            Place = Place*0.1d0
         Endif
      Enddo
 
C     Calculate Previous Position Using The Generated Velocity
 
      Do I = 1,NumberOfParticles
         Rxf(I) = Rxx(I) - Deltat*Vxx(I)
         Ryf(I) = Ryy(I) - Deltat*Vyy(I)
         Rzf(I) = Rzz(I) - Deltat*Vzz(I)
         Mxx(I) = Rxx(I)
         Myy(I) = Ryy(I)
         Mzz(I) = Rzz(I)
      Enddo
 
      Write (6,*) 'Final Energy          : ',Uold
 
      Return
      End
