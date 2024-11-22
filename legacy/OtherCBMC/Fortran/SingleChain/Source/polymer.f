      Program Polymer
C     Samples configurations of a single polymer chain; the polymer chain is generated either using a "normal" MC  scheme or CBMC
      Implicit None

      Include 'system.inc'

      Integer Sstmm
      Double Precision M1
      
C     prepare random number generator  
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

C     Read input file (see file ../Source/run)
      Open(21, FILE="input")
 
      Read(21,*) NumberOfSteps,NumberOfInitializationSteps,ChainLength,
     & NumberOfTrialPositions
      Read(21,*) Lcbmc
      Read(21,*) Beta,Rcut,A,kt,Theta0

      Close(21)

      Write(6,*) 'Cbmc Program Of A Single Chain'
      Write(6,*)
      
C     Output input parameters 
      If(Lcbmc) Then
         Write(6,*) 'Cbmc Is Used'
      Else
         Write(6,*) 'Random Chains Are Grown'
      Endif

      Write(6,*)
      Write(6,*) 'Number Of Cycles  : ',NumberOfSteps
      Write(6,*) 'Number Of Init.   : ',NumberOfInitializationSteps
      Write(6,*) 'Temperature       : ',Beta
      Write(6,*) 'Chain Length      : ',ChainLength
      Write(6,*) 'Number Of Trials  : ',NumberOfTrialPositions
      Write(6,*)
      Write(6,*) 'Rcut              : ',Rcut
      Write(6,*) 'A                 : ',A
      Write(6,*) 'kt                : ',kt
      Write(6,*) 'Theta0 [Radians]  : ',Theta0
      Write(6,*)

      Beta   = 1.0d0/Beta
      A = A/(Rcut*Rcut)

C     If certain input parameters are no good, warn the user and kill the program 
      If(ChainLength.Gt.Maxchain) Then
         Write(6,*) 'Chain Is Too Long !!!'
         Stop
      Endif

      If(NumberOfTrialPositions.Gt.Maxtrial) Then
         Write(6,*) 'Too Many Trial Positions !!!'
         Stop
      Endif

C     Let the action begin!
      Call Mcloop
 
      Stop
      End
