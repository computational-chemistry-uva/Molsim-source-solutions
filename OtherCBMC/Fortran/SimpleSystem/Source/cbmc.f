      Program Cbmc
      Implicit None

      Integer Sstmm,NumberOfTrials,Maxtrial,I,J,K,
     &        NumberOfCycles,Itrial,CycleMultiplication
      
      Parameter (Maxtrial = 100)
      Parameter (CycleMultiplication = 1000)

      Double Precision M1,Uold,CumW,Ws,Sumnew,
     &                 Sumold,RandomNumber,EnergySquaredSum,X1,
     &                 X2,X3,X1n,X2n,X3n,X1t(Maxtrial),
     &                 X2t(Maxtrial),X3t(Maxtrial),
     &                 Ut(Maxtrial),Deltax,EnergySum,CumCycle,
     &                 Accepted
 
Cccccccccccccccccccccccccccccccccccccccccccccccc
C     Written By Thijs Vlugt On 6-10-1998      C
C                                              C
C     Initialization random number generator   C
Cccccccccccccccccccccccccccccccccccccccccccccccc
 
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

Ccccccccccccccccccccccccccccccccccc
C     Read Input From Std Input   C
Ccccccccccccccccccccccccccccccccccc

      Write(*,*) 'How Many Cycles (x ',CycleMultiplication,')    ? '
      Read(*,*) NumberOfCycles

      Write(*,*) 'How Many Trial Directions ? '
      Read(*,*) NumberOfTrials

      If(NumberOfTrials.Gt.Maxtrial) then
        Write(*,*) 'NumberOfTrials > MaxTrials'
        Stop
      End If

      Write(*,*) 'Maximum Displacement      ? '
      Read(*,*) Deltax
 
      X1 = 0.0d0
      X2 = 0.0d0
      X3 = 0.0d0
      EnergySum = 0.0d0
      CumCycle = 0.0d0
      EnergySquaredSum = 0.0d0
      Accepted = 0.0d0

Cccccccccccccccccccccccccccccccccc
C     Start The Simulation       C
C     Calculate Initial Energy   C
Cccccccccccccccccccccccccccccccccc

      Uold = X1**2 + X2**2 + X3**2

      Do I=1,NumberOfCycles
         Do J=1,CycleMultiplication
                        
Ccccccccccccccccccccccccccccccccccccccccccccccccccccx
C     Cbmc; Generate NumberOfTrials Displacewments  C
C     Selectone According To Its Boltzmann          C
C     Factor...                                     C
C                                                   C
C     Old Coordinates : X1, X2, X3                  C
C     New Coordinates : X1t, Y1t, Z1t               C
C     Sumnew          : Weight New Config.          C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccc

            Sumnew = 0.0d0

            Do Itrial = 1,NumberOfTrials
               X1t(Itrial) = X1 + 
     &              (RandomNumber()-0.5d0)*Deltax
               X2t(Itrial) = X2 + 
     &              (RandomNumber()-0.5d0)*Deltax
               X3t(Itrial) = X3 + 
     &              (RandomNumber()-0.5d0)*Deltax
 
               Ut(Itrial)  = Dexp(-(X1t(Itrial)**2 + 
     &                             X2t(Itrial)**2 +
     &                             X3t(Itrial)**2))

               Sumnew = Sumnew + Ut(Itrial)
            Enddo

Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Select One Of The Trial Directions...   C
Ccccccccccccccccccccccccccccccccccccccccccccccc

            Ws = RandomNumber()*Sumnew
            CumW = 0 
            Itrial = 0

            Do While(CumW.Lt.Ws)
               Itrial = Itrial + 1
               CumW  = CumW  + Ut(Itrial)
            Enddo

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Old Configuration                          C
C     Generate NumberOfTrials-1 Positions Around C
C     The New Config; The First One Is The       C
C     Old Configuration...                       C
C                                                C
C     Sumold = Weight Of Old Configuration       C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

            Sumold = 0.0d0

            Do K=1,NumberOfTrials
               If(K.Eq.1) Then
                  X1n = X1
                  X2n = X2
                  X3n = X3
               Else
C     Start Modification
                  X1n = X1t(Itrial) + 
     &                 (RandomNumber()-0.5d0)*Deltax
                  X2n = X2t(Itrial) + 
     &                 (RandomNumber()-0.5d0)*Deltax
                  X3n = X3t(Itrial) + 
     &                 (RandomNumber()-0.5d0)*Deltax
               Endif
                  
               Sumold = Sumold + Dexp(-(X1n**2 + X2n**2 + X3n**2))

C     End   Modification
            Enddo

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Accept Or Reject This Configuration    C
Cccccccccccccccccccccccccccccccccccccccccccccc

            If(RandomNumber().Lt.(Sumnew/Sumold)) Then
               X1 = X1t(Itrial)
               X2 = X2t(Itrial)
               X3 = X3t(Itrial)
               
               Uold = X1**2 + X2**2 + X3**2
               Accepted  = Accepted + 1.0d0
            Endif

            EnergySum = EnergySum + Uold
            EnergySquaredSum = EnergySquaredSum + Uold**2
            CumCycle = CumCycle + 1.0d0

         Enddo
      Enddo
 
Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Write(6,*) 'Average Energy          : ',EnergySum/CumCycle
      Write(6,*) 'Sigma <E>               : ',
     &     Dsqrt((EnergySquaredSum/CumCycle)-(EnergySum/CumCycle)**2)
      Write(6,*) 'Fraction Accepted Moves : ',Accepted/CumCycle

      Stop
      End
