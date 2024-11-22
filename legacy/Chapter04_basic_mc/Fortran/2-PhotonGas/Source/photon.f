      Program Photon
      Implicit None

      Double Precision RandomNumber,Beta,Sum,Count,M1
      Integer NumberOfCycles,NumberOfInitCycles,New,Old,I,J,Sstmm
      Integer CycleMultiplication

      Parameter (CycleMultiplication=1000)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

Ccccccccccccccccccc
C     Read Info   C
Ccccccccccccccccccc

      Write(*,*) 'How Many Cycles (x ',CycleMultiplication,
     &')                       '
      Read(*,*)  NumberOfCycles

      Write(*,*) 'How Many Initialization Cycles ? (x '
     &        ,CycleMultiplication, ')'
      Read(*,*)  NumberOfInitCycles

      If(NumberOfInitCycles.Ge.NumberOfCycles) Then
        Write(*,*) 'Initialisation must be shorter than the run!'
        Stop
      Endif

      Write(*,*) 'Beta*Epsilon ? (Example: 1.0d0)'
      Read(*,*)  Beta

      New = 1
      Old = 1
      Sum = 0.0d0
      Count = 0.0d0

CCcccccccccccccccccccccccccccccccccccccccccc
C     Loop Over All Cycles                 C
C                                          C
C     Old = Old Position (Integer !!)      C
C     New = New Position (Integer !!)      C
Cccccccccccccccccccccccccccccccccccccccccccc

      Do I=1,NumberOfCycles
         Do J=1,CycleMultiplication

C Start Modification
            If(RandomNumber().Lt.0.5d0) Then
               New = Old + 1
            Else
               New = Old - 1
            Endif

            If(New.Lt.0) New = 0
C End   Modification

Cccccccccccccccccccccccccccccc
C     Check For Acceptance   C
Cccccccccccccccccccccccccccccc
            

            If(RandomNumber().Lt.Dexp(-Beta*(Dble(New-Old)))) 
     &           Old = New

cCccccccccccccccccccccccccccccccccccccccccccc
C     Calculate Average Occupancy Result    C
Ccccccccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.NumberOfInitCycles) Then
               Sum = Sum + Dble(Old)
               Count = Count + 1.0d0
            Endif

c     MODIFICATION
c     question 4: remove the above line "If(New.Lt.0) New = 0" and replace the above code by:
c            If((New>=0).and.(RandomNumber().Lt.Dexp(-Beta*(Dble(New-Old))))) then 
c              If((I.Gt.NumberOfInitCycles).and.(Old.ne.New)) Then
c                 Sum = Sum + Dble(New)
c                 Count = Count + 1.0d0
c              Endif
c             Old = New
c            Endif

         Enddo
      Enddo

Cccccccccccccccccccccccccccccccccccccccccccc
C     Write The Final Result               C
Cccccccccccccccccccccccccccccccccccccccccccc

      Write(6,*) 'Average Value     : ',Sum/Count
      Write(6,*) 'Theoretical Value : ',
     &     1.0d0/(Dexp(Beta)-1.0d0)
      Write(6,*) 'Relative Error    : ',
     &     Dabs((Dexp(Beta)-1.0d0)*((Sum/Count) - 
     &     (1.0d0/(Dexp(Beta)-1.0d0))))

      Stop
      End
