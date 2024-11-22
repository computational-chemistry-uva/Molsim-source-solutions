      Subroutine Mcmove(RunningEnergy, RunningVirial, Attempts, Accepted, 
     &                  MaxDisplacement)

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Double Precision EnergyNew, EnergyOld, RunningEnergy, RandomNumber, 
     &     Xn, Yn, Zn, VirialOld, VirialNew, RunningVirial, MaxDisplacement
      Integer O, Attempts, Accepted
 
      Attempts = Attempts + 1

C     ---Select A Particle At Random

      O = Int(Dble(NumberOfParticles)*RandomNumber()) + 1

C     ---Calculate Energy Old Configuration

      Call EnergyParticle(X(O), Y(O), Z(O), O, 1, EnergyOld, VirialOld)

C     ---Give Particle A Random Displacement

      Xn = X(O) + (RandomNumber()-0.5d0)*MaxDisplacement
      Yn = Y(O) + (RandomNumber()-0.5d0)*MaxDisplacement
      Zn = Z(O) + (RandomNumber()-0.5d0)*MaxDisplacement

C     ---Calculate Energy New Configuration:

      Call EnergyParticle(Xn, Yn, Zn, O, 1, EnergyNew, VirialNew)

C     ---Acceptance Test

      If (RandomNumber().Lt.Exp(-Beta*(EnergyNew-EnergyOld))) Then

C     --Accepted

        Accepted = Accepted + 1
        RunningEnergy = RunningEnergy + (EnergyNew-EnergyOld)
        RunningVirial = RunningVirial + (VirialNew-VirialOld)
        
C        ---Put Particle In Simulation Box

        If (Xn.Lt.0.0d0) Then
          Xn = Xn + Box
        Elseif (Xn.Gt.Box) Then
          Xn = Xn - Box
        Endif
        
        If (Yn.Lt.0.0d0) Then
          Yn = Yn + Box
        Elseif (Yn.Gt.Box) Then
          Yn = Yn - Box
        Endif
        
        If (Zn.Lt.0.0d0) Then
          Zn = Zn + Box
        Elseif (Zn.Gt.Box) Then
          Zn = Zn - Box
        Endif
        
        X(O) = Xn
        Y(O) = Yn
        Z(O) = Zn
      End If
      Return
      End
