      Program Harmonic
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sets Of Harmonic Oscillators; E Or T Is Constant   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Sstmm, NumberOfOscillators, MaxOscillators, NumberOfCycles,
     &     MaxEnergy, Ninit, TotalEnergy, Utot, I, J, K, A, B,
     &     OscA, OscB, Choice,CycleMultiplication
      
      Parameter (MaxOscillators = 100000)
      Parameter (MaxEnergy = 100000)
      Parameter (CycleMultiplication = 1000)
      
      Integer Oscillator(0:MaxOscillators)
      
      Double Precision M1, RandomNumber, Beta, EnergySum, Count,
     &     Distribution(0:MaxEnergy), Normalization
      
Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc
      
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1) 
 
Ccccccccccccccccccccccccccccccc
C     Read Input Parameters   C
Ccccccccccccccccccccccccccccccc

      TotalEnergy = 0
      Beta = 0.0d0

      Write(*,*) 'Number of Oscillators  ? '
      Read(*,*) NumberOfOscillators

      Write(*,*) 'Number of Cycles (x ',CycleMultiplication,')      ? '
      Read(*,*) NumberOfCycles

      If(NumberOfCycles.Le.3) then
        Write(*,*) 'Please increase the number of cycles'
        Stop
      end if
      
      Write(*,*) '1. NVE Ensemble'
      Write(*,*) '2. NVT Ensemble'
      Read(*,*) Choice
      
      If(Choice.Eq.1) Then
        Write(*,*) 'Total Energy           ? '
        Read(*,*) TotalEnergy
        
        If(NumberOfOscillators.Le.3.Or.
     &       NumberOfOscillators.Ge.MaxOscillators.Or.
     &       TotalEnergy.Le.3.Or.TotalEnergy.Ge.MaxEnergy) then
          Write(*,*) 'Error in NumberOfOscillators or TotalEnergy....'
          Stop
        End if
      Else
        Write(*,*) 'Beta            ? '
        Read(*,*) Beta
        
        If(Beta.Lt.0.0d0.Or.Beta.Gt.1.0d7) then
          Write(*,*) 'Error in Beta....'
          Stop
        End if
        NumberOfOscillators = 1
      Endif
      
Cccccccccccccccccccc
C     Initialize   C
Cccccccccccccccccccc
      
      Do I=1,MaxOscillators
        Oscillator(I) = 0
      Enddo
      
      Do I=1,MaxEnergy
        Distribution(I) = 0.0d0
      Enddo
      
      Ninit = NumberOfCycles/2
      EnergySum = 0.0d0
      Count = 0.0d0
      Normalization = 0.0d0
      Utot = 0
      
Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Make An Initial Distribution Of The Total   C
C     Energy Over The Levels                      C
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      If(Choice.Eq.1) Then
        
        I = 0
 10     I = I + 1
        
        If(I.Gt.NumberOfOscillators) I = 1
        
        Oscillator(I) = Oscillator(I) + 1
        Utot = Utot + 1
        
        If(Utot.Ne.TotalEnergy) Goto 10
      Endif
      
Cccccccccccccccccccccccccccccccccccccccc
C     Calculate Total Initial Energy   C
Cccccccccccccccccccccccccccccccccccccccc
      
      Utot = 0
      
      Do I=1,NumberOfOscillators
        Utot = Utot + Oscillator(I)
      Enddo

      Write(6,*)
      Write(6,*) 'Initial Energy : ',Utot     

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do I=1,NumberOfCycles
        Do J=1,CycleMultiplication
          Do K=1,NumberOfOscillators

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Choose 2 Different Levels At Random    C
C     Exchange Energy At Random              C
Cccccccccccccccccccccccccccccccccccccccccccccc

            If(Choice.Eq.1) Then

Cccccccccccccccccccccc
C     Nve Ensemble   C
Cccccccccccccccccccccc

 20           OscA = 1 + Int(Dble(NumberOfOscillators)*RandomNumber())
              OscB = 1 + Int(Dble(NumberOfOscillators)*RandomNumber())

              If(OscA.Eq.OscB) Goto 20
              
              If(RandomNumber().Lt.0.5d0) Then
                A =  1
                B = -1
              Else
                A = -1
                B =  1
              Endif
                  
Cccccccccccccccccccccccc
C     Reject When E<0  C
Cccccccccccccccccccccccc

              If(Min(Oscillator(OscA)+A,Oscillator(OscB)+B).Ge.0) Then
                Oscillator(OscA) = Oscillator(OscA) + A
                Oscillator(OscB) = Oscillator(OscB) + B 
              Endif
              
            Else
                  
Ccccccccccccccccccccccc
C     Nvt Ensemble    C
Ccccccccccccccccccccccc

              If(RandomNumber().Lt.0.5d0) Then
                A =  1
              Else
                A = -1
              Endif

Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Accept/Reject Using The Metropolis Algorithm   C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc

              If((Oscillator(1)+A.Ge.0).And.
     &             (RandomNumber().Lt.Dexp(-Beta*Dble(A)))) Then
                Oscillator(1) = Oscillator(1) + A
              Endif

            Endif
          Enddo

Cccccccccccccccc
C     Sample   C
Cccccccccccccccc

          If(I.Ge.Ninit) Then
            Distribution(Oscillator(1)) = Distribution(Oscillator(1)) + 1.0d0
            Normalization = Normalization + 1.0d0
            EnergySum = EnergySum + Dble(Oscillator(1))
            Count = Count + 1.0d0
          Endif
        Enddo
      Enddo

Cccccccccccccccccccccccccccccc
C     Write Results          C
Cccccccccccccccccccccccccccccc

      Utot = 0

      Do I=1,NumberOfOscillators
        Utot = Utot + Oscillator(I)
      Enddo
      
      Write(6,*)
      Write(6,*) 'Final Energy             : ',Utot
      Write(6,*) 'Average Energy Level 1   : ',EnergySum/Count
      
      Open(21,File='results.dat')

      Do I=0,MaxEnergy
        If(Distribution(I).Gt.0.5d0) Then
          Write(21,'(I10,E20.10)') I,Distribution(I)/Normalization
        Endif
      Enddo
      
      Close(21)
      
      Stop
      End
