      Program Mcnpt
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccc
C     Npt Simulation Of Hard-Spheres      C
Ccccccccccccccccccccccccccccccccccccccccccc

      Integer NumberOfCycles,I,J,K,N,M,CycleMultiplication
      Integer NumberOfInitCycles,Sstmm,Selected
      Logical Overlap

      Double Precision RandomNumber,R2,TrialDisplacements,
     &     AcceptedDisplacements,Xnew,Ynew,Znew,MaximumDisplacement,
     &     MaximumVolumeChange,TrialVolumeChanges,AcceptedVolumeChanges,
     &     Xold(1000),Yold(1000),Zold(1000),Beta,Pressure,Vnew,Vold,Boxold,
     &     Dx,Dy,Dz,VolumeCount,VolumeSum,M1

      Parameter (CycleMultiplication=1000)
      
Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Read Info From Disk                  C
Cccccccccccccccccccccccccccccccccccccccccccc

      Open(21,File="input")
      
      Read(21,*) NumberOfCycles,NumberOfInitCycles,NumberOfParticles,
     &     MaximumDisplacement,MaximumVolumeChange,Beta,Pressure

      Write(6,*) 'Number of cycles (x ',CycleMultiplication,') : ',NumberOfCycles
      Write(6,*) 'Number of init cycles : ',NumberOfInitCycles
      Write(6,*) 'Number of particles   : ',NumberOfParticles
      Write(6,*) 'Maximum displacement  : ',MaximumDisplacement
      Write(6,*) 'Maximum volume change : ',MaximumVolumeChange
      Write(6,*) 'Beta                  : ',Beta
      Write(6,*) 'Pressure              : ',Pressure
      
      If(NumberOfParticles.Gt.700.Or.NumberOfParticles.Lt.100) then
        Write(6,*) 'Error : Number of particles must be inbetween 100 
     & and 700'
        Stop
      End If
      
Cccccccccccccccccccccccccc
C     Initialize Stuff   C
Cccccccccccccccccccccccccc
      
      TrialDisplacements = 0.0d0
      AcceptedDisplacements = 0.0d0
      TrialVolumeChanges = 0.0d0
      AcceptedVolumeChanges = 0.0d0   
      VolumeCount = 0.0d0
      VolumeSum = 0.0d0
      Box  = 10.0d0
      
Cccccccccccccccccccccccccccccccccccccccccc
C     Put Particles On A Lattice         C
C     Number Of Particles Always Have    C
C     To Be Smaller Than 725             C
Cccccccccccccccccccccccccccccccccccccccccc
      
      N = 0
      Do I = 1,9
        Do J = 1,9
          Do K = 1,9
            N    = N + 1
            X(N) = Dble(I)*1.1d0
            Y(N) = Dble(J)*1.1d0
            Z(N) = Dble(K)*1.1d0
          Enddo
        Enddo
      Enddo
      
Ccccccccccccccccccccccccccccccccc
C     Start Of The Simulation   C
Ccccccccccccccccccccccccccccccccc

      Call Writepdb

      Do M = 1,NumberOfCycles
        Do N = 1,CycleMultiplication
            
Cccccccccccccccccccccccccccccccccccccccccccc
C     Select Particle And Move At Random   C
Cccccccccccccccccccccccccccccccccccccccccccc
          
          Selected = 1 + Int(RandomNumber()*Dble(NumberOfParticles+1))
          
          If(Selected.Gt.NumberOfParticles) Then

Ccccccccccccccccccccccccccccccc
C     Volume Change           C
C     Random Walk In Ln(V)    C
Ccccccccccccccccccccccccccccccc
            
            TrialVolumeChanges = TrialVolumeChanges + 1.0d0
            Vold = Box**3

c     Random walk in V instead of ln(V)

            Vnew = Dexp(Dlog(Vold) + 
     &           (RandomNumber()-0.5d0)*MaximumVolumeChange)

C     MODIFICATION
c     the above lines are replaced by:
c            Vnew = Vold + 
c     &           (RandomNumber()-0.5d0)*MaximumVolumeChange
c     Note that an additional change needs to be made in the acceptance
c     ratio, a bit below.

            Boxold = Box
            Box = Vnew**(1.0d0/3.0d0)

Ccccccccccccccccccccccccccccccc
C     Transform Coordinates   C
Ccccccccccccccccccccccccccccccc

            Do I = 1,NumberOfParticles
              Xold(I) = X(I)
              Yold(I) = Y(I)
              Zold(I) = Z(I)
              
              X(I) = X(I)*(Box/Boxold)
              Y(I) = Y(I)*(Box/Boxold)
              Z(I) = Z(I)*(Box/Boxold) 
            Enddo

Ccccccccccccccccccccccccccccc
C     Check For Overlaps    C
Ccccccccccccccccccccccccccccc

            Overlap = .False.
            
            Do I = 1,(NumberOfParticles-1)
              Do J = (I+1),NumberOfParticles
                
                Dx = X(I) - X(J)
                Dy = Y(I) - Y(J)
                Dz = Z(I) - Z(J)

Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc
 
                Dx = Dx - Box*
     &               Dble(Int(Dx/Box + 9999.5d0) - 9999)
                Dy = Dy - Box*
     &               Dble(Int(Dy/Box + 9999.5d0) - 9999)
                Dz = Dz - Box*
     &               Dble(Int(Dz/Box + 9999.5d0) - 9999)
                
                R2 = Dx**2 + Dy**2 + Dz**2 
                
                If(R2.Lt.1.0d0) Overlap = .True.
              Enddo
            Enddo

Ccccccccccccccccccccccccccccccccccccccccccc
C     No Overlap... Use Acceptance Rule   C
Ccccccccccccccccccccccccccccccccccccccccccc
 
C     MODIFICATION:
c     Here is the other change: NumberOfParticles + 1 is replaced by NumberOfParticles
            If(.Not.Overlap.And.RandomNumber().Lt.Dexp(-Beta*Pressure*
     &           (Vnew-Vold)+Dlog(Vnew/Vold)*Dble(NumberOfParticles+1))) Then
               
c     If(.Not.Overlap.And.RandomNumber().Lt.Dexp(-Beta*Pressure*(Vnew-Vold) + 
c     &           Dlog(Vnew/Vold)*Dble(NumberOfParticles))) Then

Cccccccccccccccccccccc
C     Accepted !!!!  C
Cccccccccccccccccccccc
                 
              AcceptedVolumeChanges = AcceptedVolumeChanges + 1.0d0
            Else
                 
Ccccccccccccccccccccccccccccc
C     Rejected !!!          C
C     Restore Coordinates   C
Ccccccccccccccccccccccccccccc

              Do I=1,NumberOfParticles
                X(I) = Xold(I)
                Y(I) = Yold(I)
                Z(I) = Zold(I)
              Enddo
              
              Box  = Boxold
              
            Endif
            
          Else
            
Cccccccccccccccccccccc
C     Displacement   C
Cccccccccccccccccccccc

            TrialDisplacements = TrialDisplacements + 1.0d0
            
            Xnew = X(Selected) + (RandomNumber()-0.5d0)*MaximumDisplacement
            Ynew = Y(Selected) + (RandomNumber()-0.5d0)*MaximumDisplacement
            Znew = Z(Selected) + (RandomNumber()-0.5d0)*MaximumDisplacement
            
Cccccccccccccccccccccccccccccccccccccc
C     Put Particle Back In The Box   C
Cccccccccccccccccccccccccccccccccccccc

            Xnew = Xnew - Box*
     &           (Dble(Int(Xnew/Box + 9999.0d0) - 9999))
            Ynew = Ynew - Box*
     &           (Dble(Int(Ynew/Box + 9999.0d0) - 9999))
            Znew = Znew - Box*
     &           (Dble(Int(Znew/Box + 9999.0d0) - 9999))
               
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     See If There Is An Overlap With Any Other Particle  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            Overlap = .False.
            
            Do I=1,NumberOfParticles
              If(I.Ne.Selected) Then
                
Ccccccccccccccccccccccc
C     Nearest Image   C
Ccccccccccccccccccccccc

                Dx = Xnew - X(I)
                Dy = Ynew - Y(I)
                Dz = Znew - Z(I)
                
                Dx = Dx - Box*
     &               Dble(Int(Dx/Box + 9999.5d0) - 9999)
                Dy = Dy - Box*
     &               Dble(Int(Dy/Box + 9999.5d0) - 9999)
                Dz = Dz - Box*
     &               Dble(Int(Dz/Box + 9999.5d0) - 9999)
                
                R2 = Dx**2 + Dy**2 + Dz**2 
                
                If(R2.Lt.1.0d0) Overlap = .True.
              Endif
            Enddo
            
Cccccccccccccccccccccccccccccccccc
C     No Overlaps, So Accepted   C
Cccccccccccccccccccccccccccccccccc
            
            If(.Not.Overlap) then
              AcceptedDisplacements = AcceptedDisplacements + 1.0d0
              X(Selected) = Xnew
              Y(Selected) = Ynew
              Z(Selected) = Znew
            Endif
          Endif  
        Enddo
          
Ccccccccccccccccccccccccccccccc
C     Make A Movie File       C
C     Print Current Volume    C
Ccccccccccccccccccccccccccccccc

        If(Mod(M,10).Eq.0) Then
          Call Writepdb
          
          Write(6,*) 'Volume / Box             : ',Box**3,Box
        Endif            
          
Cccccccccccccccccccccccc
C     Sample Volume    C
Cccccccccccccccccccccccc
          
        If(M.Gt.NumberOfInitCycles) Then
          VolumeCount = VolumeCount + 1.0d0
          VolumeSum = VolumeSum + (Box**3)
        Endif
      Enddo
      
Ccccccccccccccccccccccc
C     Print Results   C
Ccccccccccccccccccccccc

      Write(6,*)
      Write(6,*)
      Write(6,*) 'Fraction Succes (Displ.) : ',
     &     AcceptedDisplacements/TrialDisplacements
      Write(6,*) 'Fraction Succes (Volume) : ',
     &     AcceptedVolumeChanges/TrialVolumeChanges
      Write(6,*) 'Average Volume           : ',VolumeSum/VolumeCount
      Write(6,*) 'Average Density          : ',
     &     Dble(NumberOfParticles)/(VolumeSum/VolumeCount)
      Write(6,*)

      Stop
      End


