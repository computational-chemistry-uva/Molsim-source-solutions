C     System.Inc
C
C     Information About The System
C
C     Maxpart     : Maximum Number Of Particles
C
C     Fxx/Fyy/Fzz : Forces
C     Rxx/Ryy/Rzz : Positions
C     Rxf/Ryf/Rzf : Old Positions
C     Mxx/Myy/Mzz : Positions That Are Not Put Back In The Box
C     Vxx/Vyy/Vzz : Velocities
C
C     Rcutsq      : Cut-Off Radius
C     Ecut        : Cut-Off Energy

      Integer MaxParticles

      Parameter (MaxParticles = 1000)

      Double Precision Fxx(MaxParticles),Fyy(MaxParticles),Fzz(MaxParticles),
     &                 Rxx(MaxParticles),Ryy(MaxParticles),Rzz(MaxParticles),
     &                 Rxf(MaxParticles),Ryf(MaxParticles),Rzf(MaxParticles),
     &                 Vxx(MaxParticles),Vyy(MaxParticles),Vzz(MaxParticles),
     &                 Mxx(MaxParticles),Myy(MaxParticles),Mzz(MaxParticles),
     &                 Deltat,Box,HalfBox,
     &                 Rcutsq,Ecut,KineticEnergy,PotentialEnergy,
     &                 TotalEnergy,Temperature,Pressure

      Integer          NumberOfSteps,NumberOfParticles,NumberOfInitSteps

      Common /Syss/    Fxx,Fyy,Fzz,Rxx,Ryy,Rzz,Rxf,Ryf,Rzf,
     &                 Temperature,Vxx,Vyy,Vzz,Deltat,
     &                 Box,HalfBox,KineticEnergy,PotentialEnergy,Mxx,Myy,
     &                 Mzz,TotalEnergy,Rcutsq,Ecut,Pressure,NumberOfSteps,
     &                 NumberOfParticles,NumberOfInitSteps
