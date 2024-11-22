      Program Ising
      Implicit None

Ccccccccccccccccccccccccccc
C     2d Ising Model      C
C     No Pbc              C
C     4 Neighbors         C
Ccccccccccccccccccccccccccc

      Integer Maxlat,LatticeSize,Iup(4),Jup(4),Magnet,Energy,
     &     I,J,Inew,Jnew,Lnew,Lold,Diff,Sstmm,Ilat,Jlat,
     &     K,L,NumberOfCycles,Ninit,Maxmag,Mnew,Mold
      Integer CycleMultiplication

      Parameter (Maxlat = 32   )
      Parameter (Maxmag = 1024 )
      Parameter (CycleMultiplication = 100 )

      Integer Lattice(0:Maxlat+1,0:Maxlat+1)

      Double Precision M1,RandomNumber,Beta,EnergySum,
     &     Dist2,Weight,Accepted,Attempts,Dist3,EnergyCount,
     &     W(0:Maxmag),Distribution(-Maxmag:Maxmag),
     &     Ttime,Tstart

Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1) 
 
Ccccccccccccccccccccccccc
C     Read Info         C
Ccccccccccccccccccccccccc

      Write(6,*) '2d Ising Model, 4 Neightbors, No Pbc'
      Write(6,*)
      Write(6,*) 'Written By Thijs J.H. Vlugt Nov. 1999'
      Write(6,*)
      Write(6,*)

      Write(6,*) 'Lattice Size                     ? '
      Read(*,*) LatticeSize

      If(Mod(LatticeSize,2).Ne.0.Or.Mod(Maxmag,2).Ne.0.Or.
     &     Mod(Maxlat,2).Ne.0) Stop
      
      Write(6,*) 'Beta                             ? '
      Read(*,*) Beta

      Write(6,*) 'Number Of Cycles                 ? '
      Read(*,*) NumberOfCycles

      If(LatticeSize.Lt.3.Or.LatticeSize.Gt.Maxlat) Stop
      If(NumberOfCycles.Lt.10) Stop

      Do I=0,Maxmag
         W(I) = 1.0d0
      Enddo
      
Cccccccccccccccccccccccccccccccccc
C     Multicanonical Algorithm   C
Cccccccccccccccccccccccccccccccccc

      Open(21,File="w.dat")

      Do I=0,Maxmag,2
         Read(21,*) K,EnergySum
         
         If(Abs(K).Gt.Maxmag.Or.Mod(K,2).Ne.0) Then
            Write(6,*) 'I Out Of Range !!'
            Stop
         Endif

         W(K) = EnergySum
      Enddo

      Close(21)
      
      Ninit = NumberOfCycles/3
      EnergySum = 0.0d0
      EnergyCount = 0.0d0
      Accepted = 0.0d0
      Attempts = 0.0d0
      Tstart = Ttime()
      
      If(Ninit.Gt.100) Ninit = 100

      Do I=-Maxmag,Maxmag
         Distribution(I) = 0.0d0
      Enddo

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Lattice                     C
C     Random Lattice                         C
C                                            C
C     Add One Boundary Layer With Spin=0     C
C     This Is A Way To Avoid If-Statements   C
C     In The Computation Of The Energy       C
C                                            C
C     LatticeSize  = Size Of The Lattice     C
C     Lattice(I,J) = Spin Of Site (I,J)      C
Cccccccccccccccccccccccccccccccccccccccccccccc

      Do I=0,Maxlat+1
         Do J=0,Maxlat+1
            Lattice(I,J) = 0
         Enddo
      Enddo

      Do I=1,LatticeSize
         Do J=1,LatticeSize
            If(RandomNumber().Lt.0.5d0) Then
               Lattice(I,J) =  1
            Else
               Lattice(I,J) = -1
            Endif
         Enddo
      Enddo

Ccccccccccccccccccccccccccccccc
C     Initialize Neighbours   C
Ccccccccccccccccccccccccccccccc

      Iup(1) =  1
      Jup(1) =  0
      
      Iup(2) = -1
      Jup(2) =  0

      Iup(3) =  0
      Jup(3) =  1
      
      Iup(4) =  0
      Jup(4) = -1

Ccccccccccccccccccccccccccccccccccccccc
C     Calculate Initial Energy        C
C                                     C
C     Magnet = Total Magnetisation    C
Ccccccccccccccccccccccccccccccccccccccc

      Energy = 0
      Magnet = 0

      Do I=1,LatticeSize
         Do J=1,LatticeSize

            Magnet = Magnet + Lattice(I,J)

            Do K=1,4
               Inew   = I + Iup(K)
               Jnew   = J + Jup(K)
               Energy = Energy - 
     &              Lattice(I,J)*Lattice(Inew,Jnew)
            Enddo
         Enddo
      Enddo

      Energy = Energy/2

      Write(6,*)
      Write(6,*) 'Lattice Size               : ',LatticeSize
      Write(6,*) 'Beta                       : ',Beta
      Write(6,*) 'Number Of Cycles           : ',NumberOfCycles
      Write(6,*) 'Number Of Init             : ',Ninit
      Write(6,*)
      Write(6,*)
      Write(6,*) 'Initial Energy             : ',Energy
      Write(6,*) 'Initial Magnetization      : ',Magnet
      Write(6,*)
      
      If(Abs(Magnet).Gt.Maxmag) Stop

Cccccccccccccccccccccccccccccc
C     Loop Over All Cycles   C
Cccccccccccccccccccccccccccccc

      Do K=1,NumberOfCycles
         Do L=1,CycleMultiplication*LatticeSize*LatticeSize
            
Ccccccccccccccccccccccccccccccccc
C     Flip A Single Spin        C
C     Metropolis Algorithm      C
Ccccccccccccccccccccccccccccccccc

            Ilat = 1 + Int(RandomNumber()*Dble(LatticeSize))
            Jlat = 1 + Int(RandomNumber()*Dble(LatticeSize))
            
            Mold  = Magnet
            Lold  = Lattice(Ilat,Jlat)
            Lnew  = -Lold
            Mnew  = Mold + Lnew - Lold
            Diff  = 0
            Attempts = Attempts + 1.0d0
            
Cccccccccccccccccccccccccccccccccccccccccc
C     Calculate The Energy Difference    C
C     Between New And Old                C
Cccccccccccccccccccccccccccccccccccccccccc

C     Start Modification
            Do I=1,4
               Inew = Ilat + Iup(I)
               Jnew = Jlat + Jup(I)

c     diff: factor 1/2 is not needed here: interaction
c     from particle i with j and j with i gives twice the
c     interaction.
               Diff = Diff + 
     &              Lattice(Inew,Jnew)*(Lold-Lnew)
            Enddo

C     End Modification

CCcccccccccccccccccccccccccccccccccc
C     Acceptance/Rejection Rule    C
C     Use Weight Function W        C
Cccccccccccccccccccccccccccccccccccc

            If(RandomNumber().Lt.
     &           Dexp((-Beta*Dble(Diff)) + 
     &           W(Abs(Mnew)) - W(Abs(Mold)))) Then

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Update The Lattice/Energy/Magnetisation    C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

               Accepted = Accepted + 1.0d0

C     Start Modification
               Lattice(Ilat,Jlat) = Lnew
               Energy             = Energy + Diff
               Magnet             = Mnew

C     End Modification
               
               If(Abs(Magnet).Gt.Maxmag) Then
                  Write(6,*) 'Increase The Value Of Maxmag !!'
                  Stop
               Endif
            Endif
                  
            If(K.Gt.Ninit) Then
               Weight = Dexp(-W(Abs(Magnet)))
               EnergySum = EnergySum + Weight*Dble(Energy)
               EnergyCount = EnergyCount + Weight
               
               Distribution(Magnet) = 
     &              Distribution(Magnet) + 1.0d0
            Endif
         Enddo
      Enddo

      Write(6,*)
      Write(6,*) 'Average Energy             : ',
     &     EnergySum/EnergyCount

      Write(6,*) 'Fraction Accepted Swaps    : ',
     &     Accepted/Attempts

      Write(6,*) 'Elapsed Time [S]           : ',
     &     Ttime()-Tstart

      Write(6,*)
      Write(6,*) 'Final Energy        (Simu) : ',Energy
      Write(6,*) 'Final Magnetization (Simu) : ',Magnet

Cccccccccccccccccccccccccccccccccc
C     Calculate Distributions    C
C     Magnetisation              C
Cccccccccccccccccccccccccccccccccc

      Dist2 = 0.0d0
      Dist3 = 0.0d0
            
      Do I=-Maxmag,Maxmag,2
         Weight = Dexp(-W(Abs(I)))

         Dist2  = Dist2 + Distribution(I)
         Dist3  = Dist3 + Distribution(I)*Weight
      Enddo

      Open(21,File="magnetic.dat")

      Do I=-Maxmag,Maxmag,2
         If(Distribution(I).Gt.0.5d0) Then
            Weight = Dexp(-W(Abs(I)))
            
            Write(21,*) I,Distribution(I)/Dist2,
     &           Distribution(I)/Dist3*Weight
         Endif
      Enddo   

      Close(21)

Ccccccccccccccccccccccccccccc
C     Weight Distribution   C
Ccccccccccccccccccccccccccccc

      Dist2 = 0.0d0

      Do I=0,Maxmag,2
         Distribution(I) = Max(1.0d0,Distribution(I) + Distribution(-I))*
     &        Dexp(-W(Abs(I)))

         Distribution(I) = -Dlog(Distribution(I))
         Dist2 = Min(Dist2,Distribution(I))
      Enddo

      Do I=0,Maxmag,2
         Distribution(I) = Distribution(I) - Dist2
         Distribution(I) = Min(10.0d0,1.05d0*Distribution(I))
      Enddo

      Open(21,File="w.dat")

      Do I=0,Maxmag,2
         Write(21,*) I,Distribution(I)
      Enddo

      Close(21)
      
Cccccccccccccccccccccccccccccccccc
C     Calculate Final Energy     C
C     Used To Check The Code     C
Cccccccccccccccccccccccccccccccccc

      Energy = 0
      Magnet = 0

      Do I=1,LatticeSize
         Do J=1,LatticeSize

            Magnet = Magnet + Lattice(I,J)

            Do K=1,4
               Inew   = I + Iup(K)
               Jnew   = J + Jup(K)
               Energy = Energy - 
     &              Lattice(I,J)*Lattice(Inew,Jnew)
            Enddo
         Enddo
      Enddo

      Energy = Energy/2

      Write(6,*) 'Final Energy        (Calc) : ',Energy
      Write(6,*) 'Final Magnetization (Calc) : ',Magnet
      Write(6,*)     

      Stop
      End
