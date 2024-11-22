      Program Mc_Nvt
      Implicit None
C____________________________________________________C
C                                                    C
C                                                    C
C   Equation Of State Of The Lennard-Jones Fluid     C
C                                                    C
C____________________________________________________C

      Include 'potential.inc'
      Include 'parameter.inc'
      Include 'system.inc'
      Include 'conf.inc'

      Integer NumberOfEquilibrationCycles, NumberOfProductionCycles,
     &     SamplingFrequency, Mode, Icycl, NumberOfDisplacementsPerCycle,
     &     Attempts, Accepted, NumberOfCycles, Imove, J
      Double Precision RunningEnergy, TotalEnergy, RunningVirial,
     &     TotalVirial, MaxDisplacement, PressureSum, PressureCount, Pressure,
     &     ChemicalPotentialSum, ChemicalPotentialCount,Xi,Yi,Zi,
     &     EnergySum, EnergySquaredSum, EnergyCount, dummy1, dummy2, 
     &     RandomNumber
 
      Write (6, *) '**************** Mc_Nvt ***************'

C     ---Initialize Sysem

      Call Readdat(NumberOfEquilibrationCycles, NumberOfProductionCycles, 
     &     SamplingFrequency, NumberOfDisplacementsPerCycle, MaxDisplacement)

      PressureSum = 0.0d0
      PressureCount = 0.0d0
      EnergySum = 0.0d0
      EnergySquaredSum = 0.0d0
      ChemicalPotentialSum = 0.0d0
      ChemicalPotentialCount = 0.0d0

      Open (22,FILE='movie.pdb',FORM='FORMATTED',STATUS='REPLACE')

C     ---Total Energy Of The System

      Call EnergySystem(TotalEnergy, TotalVirial)
      Write (6, 99001) TotalEnergy, TotalVirial
      RunningEnergy = TotalEnergy
      RunningVirial = TotalVirial

C     ---Start Mc-Cycle

      Do Mode = 1, 2

C     --- Mode=1 Equilibration
C     --- Mode=2 Production
        
        If (Mode.Eq.1) Then
          NumberOfCycles = NumberOfEquilibrationCycles
          If (NumberOfCycles.Ne.0) Write (6, *) ' Start Equilibration '
        Else
          If (NumberOfCycles.Ne.0) Write (6, *) ' Start Production '
          NumberOfCycles = NumberOfProductionCycles
        End If
        
        Attempts = 0
        Accepted = 0
        
C     ---Intialize The Subroutine That Adjust The Maximum Displacement
        
        Call Adjust(Attempts, Accepted, MaxDisplacement)
        
        Do Icycl = 1, NumberOfCycles
          
          Do Imove = 1, NumberOfDisplacementsPerCycle
              
C     ---Attempt To Displace A Particle
              
            Call Mcmove(RunningEnergy, RunningVirial, Attempts, Accepted, 
     &           MaxDisplacement)
          End Do
            
          If (Mode.Eq.2) Then
              
C     ---Sample Averages
              
            If (Mod(Icycl,SamplingFrequency).Eq.0) Then
              Call Sample(Icycl, RunningEnergy, RunningVirial, Pressure)
              PressureSum = PressureSum + Pressure
              PressureCount = PressureCount + 1.0d0
              
              EnergySquaredSum=EnergySquaredSum+RunningEnergy*RunningEnergy
              EnergySum=EnergySum+RunningEnergy
              EnergyCount=EnergyCount+1.0
                
                
Ccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Chemical Potential: C
C     Do 10 Trial Chains                C
C     Calculate The Average Of          C
C     [Exp(-Beta*Energy)]               C
C     You Can Use The Subroutine        C
C     EnergyParticle For This.          C
Ccccccccccccccccccccccccccccccccccccccccc
                
              Do J=1,10

C     Start Modification

                Xi = RandomNumber()*Box
                Yi = RandomNumber()*Box
                Zi = RandomNumber()*Box
                Call EnergyParticle(Xi,Yi,Zi,NumberOfParticles+1,1,dummy1,dummy2)
                ChemicalPotentialSum = ChemicalPotentialSum + exp(-Beta*dummy1)
                ChemicalPotentialCount = ChemicalPotentialCount + 1.0
                  
C     End   Modification
                  
              Enddo
            Endif
          Endif
            
          If(Mod(Icycl,20).Eq.0) Call Writepdb(22)

          If (Mod(Icycl,NumberOfCycles/5).Eq.0) Then
            Write (6, *) '======>> Done ', Icycl, ' Out Of ', NumberOfCycles
              
C     ---Write Intermediate Configuration To File
              
            Call Store(MaxDisplacement)
              
C     ---Adjust Maximum Displacements
              
            Call Adjust(Attempts, Accepted, MaxDisplacement)
          End If
        End Do
        If (NumberOfCycles.Ne.0) Then
          If (Attempts.Ne.0) Write (6, 99003) Attempts, Accepted, 
     &         100.0d0*Dble(Accepted)/Dble(Attempts)

C     ---Test Total Energy

          Call EnergySystem(TotalEnergy, TotalVirial)
          If (Abs(TotalEnergy-RunningEnergy).Gt.1.D-6) Then
            Write (6, *) 
     &             ' ######### Problems Energy ################ '
          End If
          If (Abs(TotalVirial-RunningVirial).Gt.1.D-6) Then
            Write (6, *) 
     &             ' ######### Problems Virial ################ '
          End If
          Write (6, 99002) TotalEnergy, RunningEnergy, 
     &         TotalEnergy - RunningEnergy, TotalVirial, 
     &         RunningVirial, TotalVirial - RunningVirial
          Write (6,*)
            
Ccccccccccccccccccccccccccccccccccccccccccccccc
Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Print Chemical Potential And Pressure   C
Ccccccccccccccccccccccccccccccccccccccccccccccc
Ccccccccccccccccccccccccccccccccccccccccccccccc
            
          If(Mode.Eq.2) Then
              
            Write(6,*) 'Heat Capacity                      :',0.0
            Write (6,*) 'Average Pressure                  : ',
     &           PressureSum/PressureCount
            Write (6,*) 'Chemical Potential                : ',
     &           -Log((ChemicalPotentialSum/ChemicalPotentialCount)*
     &           (Box*Box*Box/Dble(NumberOfParticles)))/Beta
          End If
        End If
      End Do
      Call Store(MaxDisplacement)
      Close (22)
      Stop
        
99001 Format (' Total Energy Initial Configuration: ', F12.5, /, 
     &        ' Total Virial Initial Configuration: ', F12.5)
99002 Format (' Total Energy End Of Simulation    : ', F12.5, /, 
     &        '       Running Energy              : ', F12.5, /, 
     &        '       Difference                  :  ', E12.5, /, 
     &        ' Total Virial End Of Simulation    : ', F12.5, /, 
     &        '       Running Virial              : ', F12.5, /, 
     &        '       Difference                  :  ', E12.5)
99003 Format (' Number Of Att. To Displ. A Part.  : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
      End
