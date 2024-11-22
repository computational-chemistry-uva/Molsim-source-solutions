      Subroutine Mdloop
      Implicit None
 
      Include 'system.inc'

      Integer I,J,CycleMultiplication
      Double Precision Ran_Vel,EnergyDrift,CumEnergyDrift,CumCycles

      Parameter (CycleMultiplication = 1000)

      Call Sample(1,0)

      CumEnergyDrift = 0.0d0
      CumCycles = 0.0d0

      Do I=1,Ncycle
         Do J=1,CycleMultiplication

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Initial Coordinates              C
C                                               C
C     Perform CycleMultiplication*Ncycle Md     C
C     Simulations With Different Initial        C
C     Conditions.                               C
C                                               C
C     Xpos   = Starting Position                C
C     Vpos   = Starting Velocity                C
C     Theta  = Storage For Initial Velocity     C
C     Qstar  = Place Of The Dividing Surface    C
C                                               C
C     To Program:                               C
C                                               C
C     -Generate Initial Position/Velocity       C
C     -Integrate The Equations Of Motion By A   C
C      Function Call To Subroutine Integrate    C
C     -Beware That In Integrate The Subroutine  C
C      Sample Is Called !!                      C
C                                               C
C     The Averaged Energy Drift (Over All       C
C     Simulations) Equals CumEnergyDrift/CumCycles                 C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

C Start Modification

        Xpos = Qstar
        Vpos = Ran_Vel()
        Theta = Vpos

        call integrate(EnergyDrift)

        CumEnergyDrift = CumEnergyDrift + EnergyDrift
        CumCycles = CumCycles + 1.0

C End Modification

         Enddo
      Enddo

      Call Sample(3,0)

      Write(6,*) 'Av. Energy Drift      :',CumEnergyDrift/CumCycles

      Return
      End
