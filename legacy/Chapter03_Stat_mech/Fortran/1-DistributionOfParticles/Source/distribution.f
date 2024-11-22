      Program Distributions
      Implicit None

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Divide N Particles Among p Compartments    C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Sstmm,P,N,MaxParticles,MaxCompartments,NumberOfCycles,
     &     CycleMultiplication,I,J,K,Index

      Parameter (MaxParticles = 10000)
      Parameter (MaxCompartments = 100)
      Parameter (CycleMultiplication = 1000)

      Double Precision M1,Distribution(MaxCompartments,0:MaxParticles),
     &     RandomNumber,LnFactorial

      Integer Compartment(MaxCompartments)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1) 

Ccccccccccccccccccccccccccccccccccccccc
C     Read Info From Standard Input   C
Ccccccccccccccccccccccccccccccccccccccc
 
      Write(*,*) 'Number of Particles N? '
      Read(*,*) N

      Write(*,*) 'Number of Compartments p? '
      Read(*,*) P

      Write(*,*) 'Number of Cycles (x',CycleMultiplication,')? '
      Read(*,*) NumberOfCycles

      If(P.Lt.2.Or.P.Gt.MaxCompartments.Or.N.Lt.2.Or.N.Gt.MaxParticles) then
        Write(*,*) 'Error in input parameters....'
        Stop
      End If


Cccccccccccccccccccc
C     Initialize   C
Cccccccccccccccccccc
      
      Do I=0,N
        Do J=1,P
          Distribution(J,I) = 0.0d0
          Compartment(J)    = 0
        Enddo
      Enddo
      
Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do I=1,NumberOfCycles
        Do K=1,CycleMultiplication

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Distribute particles over the compartments                 C
C                                                                C   
C     Loop over all particles, pick a random compartment,        C
C     and add a particle to it.                                  C
C     A random number in the interval [0,1> can be               C
C     generated using  RandomNumber().                           C
C                                                                C
C     Compartment(index) = Nr Of Particles In Compartment Index  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

          Do J=1,N
            Index = 1 + Int(Dble(P)*RandomNumber())
            Compartment(Index) = Compartment(Index) + 1 
          Enddo 

C     End   Modification

Cccccccccccccccccccccccccccccc
C     Make Histogram         C
Cccccccccccccccccccccccccccccc

          Do J=1,P
            Distribution(J,Compartment(J)) = Distribution(J,Compartment(J)) + 1.0d0
            Compartment(J) = 0
          Enddo
        Enddo
      Enddo

Cccccccccccccccccccccccccccccc
C     Write Results          C
Cccccccccccccccccccccccccccccc

      Open(21,File='results.dat')

      Do J=1,P
        Do I=0,N
          If(Distribution(J,I).Gt.0.5d0)
     &         Distribution(J,I) = Distribution(J,I)/Dble(NumberOfCycles*CycleMultiplication)
          Write(21,*) I, Distribution(J,I)
        Enddo
        Write(21,*)
      Enddo

      Close(21)

Ccccccccccccccccccccccccccccccccccc
C     Write Analytical Dist.      C
Ccccccccccccccccccccccccccccccccccc

      Open(21,File='analytical.dat')

      Do J=0,N
         Write(21,*) J, Dexp(lnfactorial(N)-lnfactorial(J)-
     &       lnfactorial(N-J)-J*Dlog(dble(P))-(N-J)*Dlog(P/Dble(P-1)))
      Enddo

      Close(21)

Cccccccccccccccccccccccccccc
C     End Of The Program   C
Cccccccccccccccccccccccccccc

      Stop
      End
