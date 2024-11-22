      Program Boltzmann
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Boltzmann Distribution       C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer MaxEnergyLevels,I,N

      Parameter (MaxEnergyLevels = 10000)

      Double Precision Distribution(0:MaxEnergyLevels),Normalization,
     &     Temperature,Beta,Tmp

      Write(*,*) 'Number Of Energy Levels (2-10000) ? '
      Read(*,*) N

      Write(*,*) 'Temperature ? '
      Read(*,*) Temperature

      If(N.Lt.2.Or.N.Ge.MaxEnergyLevels.Or.
     &     Temperature.Lt.1.0d-7.Or.Temperature.Gt.1.0d7) then
        Write(*,*) 'Input parameter error...'
        Stop
      End if

      Beta = 1.0d0/Temperature
      Normalization = 0.0d0

      Do I=0,(N-1)
         Tmp=Dexp(-Beta*Dble(I))

c     Start Modification
c     question 2: Degeneracy of energy level i equals i+1
c         Tmp= dble(i+1)*Dexp(-Beta*Dble(I))

c     question 3: Linear rotor. 
c         Tmp= dble(2*i+1)*Dexp(-Beta*0.5*Dble(I*(I+1)))
c     End Modification

         Normalization = Normalization + Tmp
         Distribution(I) = Tmp
      Enddo

      Open(21,File="results.dat")

      Do I=0,(N-1)
         Write(21,'(I8,F20.10)') I,Distribution(I)/Normalization
      Enddo

      Close(21)

      Stop
      End
