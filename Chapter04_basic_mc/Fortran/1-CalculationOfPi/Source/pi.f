      Program CalcPi
      Implicit None

      Double Precision RandomNumber,Ratio,X,Y,Hits,Trials,Pi,M1
      Integer I,J,NumberOfCycles,Sstmm,CycleMultiplication
      
      Parameter (CycleMultiplication=1000)

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Calculates Pi Using The Circle/Square Problem   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

      Write(*,*) 'Number Of Cycles (x ',CycleMultiplication,') ?'
      Read(*,*)  NumberOfCycles

      Write(*,*) 'Ratio L/D        ? (Always >= 1 !)'
      Read(*,*)  Ratio

      If(Ratio.Lt.1.0d0) Then
         Write(6,*) 'Ratio Must Be At Least 1 !!!'
         Stop
      Endif

      Hits = 0.0d0
      Trials = 0.0d0

      Open(21,File='results.dat',Form='Formatted')

Ccccccccccccccccccccccccccccc
C     Loop Over All Cycles  C
Ccccccccccccccccccccccccccccc

      Do I=1,NumberOfCycles
         Do J=1,CycleMultiplication

Ccccccccccccccccccccccccccccccccccccc
C     Generate A Uniform Point      C
C     Check If It Is In The Circle  C
C     Trials = Number Of Points     C
C     Hits = Number In The Circle   C
Ccccccccccccccccccccccccccccccccccccc

C Start Modifications

c     generation of random x and y between -0.5 and +0.5
            X = (RandomNumber() - 0.5d0)*Ratio
            Y = (RandomNumber() - 0.5d0)*Ratio

c     if point is in the circle, Hits is increased by 1
            If(sqrt(X*X+Y*Y).Lt.0.5d0) Hits = Hits + 1.0d0

c     increase number of points
            Trials = Trials + 1.0d0

C End Modifications

         Enddo
      Enddo

      Close(21)

      Pi = 4*Ratio*Ratio*Hits/Trials

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     The Real Value Of Pi Can Be Calculated Using       C
C     Pi = 4.0 * Arctan (1.0)                            C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Write(6,*) 'Estimate Of Pi : ',Pi
      Write(6,*) 'Real Pi        : ',4.0d0*Datan(1.0d0)
      Write(6,*) 'Relative Error : ',
     &           Dabs(Pi-4.0d0*Datan(1.0d0))/(4.0d0*Datan(1.0d0))
      
      Stop
      End
