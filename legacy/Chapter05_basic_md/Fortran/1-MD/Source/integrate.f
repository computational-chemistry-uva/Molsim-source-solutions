      Subroutine Integrate(Step,MomentumX,MomentumY,MomentumZ)
      Implicit None
 
      Include 'system.inc'
 
C     Integrate The Equations Of Motion And Calculate The Total Momentum 
 
      Integer I,Step
      Double Precision Xxx(MaxParticles),Yyy(MaxParticles),Zzz(MaxParticles),
     &                 Scale,MomentumX,MomentumY,MomentumZ

Ccccccccccccccccccccccccccccccccccccccccccccccccccc 
C     Set Kinetic Energy To Zero                  C
C     Verlet Integrator                           C
C                                                 C
C     Xxx/Yyy/Zzz = New Position (T + Delta T)    C
C     Rxx/Ryy/Rzz = Position     (T          )    C
C     Rxf/Ryf/Rzf = Old Position (T - Delta T)    C
C                                                 C
C     Vxx/Vyy/Vzz = Velocity     (T          )    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

      KineticEnergy = 0.0d0
 
      Do I = 1,NumberOfParticles
         Xxx(I) = 2.0d0*Rxx(I) - Rxf(I) + Fxx(I)*Deltat*Deltat
         Yyy(I) = 2.0d0*Ryy(I) - Ryf(I) + Fyy(I)*Deltat*Deltat
         Zzz(I) = 2.0d0*Rzz(I) - Rzf(I) + Fzz(I)*Deltat*Deltat
 
         Vxx(I) = (Xxx(I) - Rxf(I))/(2.0d0*Deltat)
         Vyy(I) = (Yyy(I) - Ryf(I))/(2.0d0*Deltat)
         Vzz(I) = (Zzz(I) - Rzf(I))/(2.0d0*Deltat)
 
         KineticEnergy = KineticEnergy + 0.5d0*(Vxx(I)*Vxx(I) + 
     &        Vyy(I)*Vyy(I) + Vzz(I)*Vzz(I))
      Enddo
 
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     For Nstep < Ninit; Use The Velocity Scaling To Get The Exact Temperature  C
C     Otherwise: Scale = 1. Beware That The Positions/Velocities Have To Be     C
C     Recalculated !!!!                                                         C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      If (Step.Le.NumberOfInitSteps) Then
         Scale = Dsqrt(Temperature*Dble(3*NumberOfParticles-3)/
     &       (2.0d0*KineticEnergy))
      Else
         Scale = 1.0d0
      Endif
 
      KineticEnergy = 0.0d0
      MomentumX = 0.0d0
      MomentumY = 0.0d0
      MomentumZ = 0.0d0
 
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Scale Velocities And Put Particles Back In The Box      C
C     Beware: The Old Positions Are Also Put Back In The Box  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      Do I = 1,NumberOfParticles
         Vxx(I) = Scale*Vxx(I)
         Vyy(I) = Scale*Vyy(I)
         Vzz(I) = Scale*Vzz(I)
 
         MomentumX = MomentumX + Vxx(I)
         MomentumY = MomentumY + Vyy(I)
         MomentumZ = MomentumZ + Vzz(I)
 
         Xxx(I) = Rxf(I) + 2.0d0*Vxx(I)*Deltat
         Yyy(I) = Ryf(I) + 2.0d0*Vyy(I)*Deltat
         Zzz(I) = Rzf(I) + 2.0d0*Vzz(I)*Deltat

         Mxx(I) = Mxx(I) + Xxx(I) - Rxx(I)
         Myy(I) = Myy(I) + Yyy(I) - Ryy(I)
         Mzz(I) = Mzz(I) + Zzz(I) - Rzz(I)

         KineticEnergy = KineticEnergy + 0.5d0*(Vxx(I)*Vxx(I) + Vyy(I)*Vyy(I) 
     &        + Vzz(I)*Vzz(I))
         
         Rxf(I) = Rxx(I)
         Ryf(I) = Ryy(I)
         Rzf(I) = Rzz(I)
 
         Rxx(I) = Xxx(I)
         Ryy(I) = Yyy(I)
         Rzz(I) = Zzz(I)

C     Put Particles Back In The Box
C     The Previous Position Has The Same
C     Transformation (Why ??)

         If (Rxx(I).Gt.Box) Then
            Rxx(I) = Rxx(I) - Box
            Rxf(I) = Rxf(I) - Box
         Elseif (Rxx(I).Lt.0.0d0) Then
            Rxx(I) = Rxx(I) + Box
            Rxf(I) = Rxf(I) + Box
         Endif
 
         If (Ryy(I).Gt.Box) Then
            Ryy(I) = Ryy(I) - Box
            Ryf(I) = Ryf(I) - Box
         Elseif (Ryy(I).Lt.0.0d0) Then
            Ryy(I) = Ryy(I) + Box
            Ryf(I) = Ryf(I) + Box
         Endif
 
         If (Rzz(I).Gt.Box) Then
            Rzz(I) = Rzz(I) - Box
            Rzf(I) = Rzf(I) - Box
         Elseif (Rzz(I).Lt.0.0d0) Then
            Rzz(I) = Rzz(I) + Box
            Rzf(I) = Rzf(I) + Box
         Endif
      Enddo
 
Ccccccccccccccccccccccccccccccccccccccccccccc
C     Add The Kinetic Part Of The Pressure  C
Ccccccccccccccccccccccccccccccccccccccccccccc
 
      Pressure = Pressure+ 2.0d0*KineticEnergy*Dble(NumberOfParticles)/
     &     (Box*Box*Box*Dble(3*NumberOfParticles-3))
 
      Return
      End





