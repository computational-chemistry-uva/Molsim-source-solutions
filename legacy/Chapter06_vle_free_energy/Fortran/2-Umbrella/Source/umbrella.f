      Program Umbrella
      Implicit None

      Integer i,j,k,Maxx,CycleMultiplication,Sstmm
      Integer NumberOfCycles,NumberOfAcceptedMoves,NumberOfAttempts
      Double Precision Width,M1,RandomNumber
      Double Precision LeftBoundary,RightBoundary,Xold,Xnew,Uold,Unew

      Parameter (Maxx = 100000)
      Parameter (Width = 0.001d0)
      Parameter (CycleMultiplication = 1000)

      Double Precision MaximumDisplacement,PositionCounter(2*Maxx+1)

C   initialize the random number generator with the system time
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

      Write(6,*) 'How many cycles (x ',CycleMultiplication,')? '
      Read(*,*) NumberOfcycles

      Write(6,*) 'Maximum Displacement? '
      Read(*,*) MaximumDisplacement 

      Write(6,*) 'Minimum Value of X  ? '
      Read(*,*) LeftBoundary        

      Write(6,*) 'Maximum Value of X  ? '
      Read(*,*) RightBoundary 

      If(LeftBoundary.Ge.RightBoundary) Then
        Write(*,*) 'LeftBoundary should be smaller than RightBoundary' 
        Stop
      Endif
     
C  generate initial position and energy
      Xold=0.5d0*(RightBoundary+LeftBoundary)
      Uold=Xold*Xold
      Xnew=0.0d0
      Unew=0.0d0
      NumberOfAcceptedMoves=0
      NumberOfAttempts=0
      
      Do I=-Maxx,Maxx
         PositionCounter(I+Maxx+1)=0.0d0
      Enddo

      Write(*,*) 'Calculating.....'

C  loop over all cycles

      Do I=1,NumberOfCycles
         Do J=1,CycleMultiplication
            Xnew=Xold+(RandomNumber()-0.5d0)*MaximumDisplacement
            Unew=Xnew*Xnew
            NumberOfAttempts=NumberOfAttempts+1
            k=Int(Xold/Width)

            If((k.Lt.-Maxx).OR.(k.Gt.Maxx)) Then
               Write(*,*) 'Out of range : j !!!'
               Stop
            Else 
               If(k.eq.0) Then
                  PositionCounter(k+Maxx+1)=PositionCounter(k+Maxx+1)+0.5d0
               Else
                  PositionCounter(k+Maxx+1)=PositionCounter(k+Maxx+1)+1.0d0
               Endif
            Endif

C reject when outside slice
            If(Xnew.Gt.LeftBoundary.AND.Xnew.Lt.RightBoundary) Then
               If(RandomNumber().Lt.Dexp(Uold-Unew)) Then
                  Xold=Xnew
                  Uold=Unew
                  NumberOfAcceptedMoves=NumberOfAcceptedMoves+1
               Endif
            Endif
         Enddo
      Enddo

C write results

      Open(21,File="Umbrella.dat")
      Do I=-Maxx,Maxx
         If(PositionCounter(I+Maxx+1).Gt.0.5d0) Then
            Write(21,*) I*Width,PositionCounter(I+Maxx+1)/Dble(NumberOfAttempts)
         Endif
      Enddo
  
      Close(21)

      Write(*,*)
      Write(*,*) 'Fraction of accepted displacements: ',(Dble(NumberOfAcceptedMoves)/Dble(NumberOfAttempts))
      Write(*,*)
      Write(*,*) 'Histogram Written To Umbrella.dat'

      Stop
      End
