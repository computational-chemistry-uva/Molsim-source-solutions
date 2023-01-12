      Subroutine Mdloop
      Implicit None
 
      Include 'system.inc'
 
      Integer I
      Double Precision Av(6),MomentumX,MomentumY,MomentumZ,Tempz,utot0,dutot
 
      Do I = 1,6
         Av(I) = 0.0d0
      Enddo
      dutot=0.0d0

C     Initialize Radial Distribution Function

      Call Sample_Gyra(1)

C     Initialize Diffusion Coefficient

      Call Sample_Diff(1)

C     Loop Over All Cycles

      Do I = 1,NumberOfSteps
 
C     Calculate The Force
 
         Call Force
 
C     Integrate The Equations Of Motion
 
         Call Integrate(I,MomentumX,MomentumY,MomentumZ)
 
         If (I.Eq.1) Then
            Write (6,*)
            Write (6,*) 'Momentum X-Dir.        : ',MomentumX
            Write (6,*) 'Momentum Y-Dir.        : ',MomentumY
            Write (6,*) 'Momentum Z-Dir.        : ',MomentumZ
            Write (6,*)
            Write (6,*)
            Write (6,'(3a)') ' Step  TotalEnergy  KineticEnergy',
     &           '  PotentialEnergy  Temperature  Pressure'
            Write (6,*)
         Endif
 
         TotalEnergy = KineticEnergy + PotentialEnergy

c     choose a starting point to calculate the energy drift
c     declare utot0

c     Start Modification

         if(I.eq.NumberOfInitSteps) utot0=TotalEnergy
 
c     End Modification

         Tempz = 2.0d0*KineticEnergy/Dble(3*NumberOfParticles - 3)
 
         If (I.Gt.NumberOfInitSteps.And.Mod(I,50).Eq.0) 
     &        Write (6,'(I5,8(2x,E12.5))') I,TotalEnergy,KineticEnergy,
     &        PotentialEnergy,Tempz,Pressure

         If (Mod(I,100).Eq.0) Call Writepdb
 
         If (I.Eq.NumberOfSteps) Then
            Write (6,*)
            Write (6,*) 'Momentum X-Dir.        : ',MomentumX
            Write (6,*) 'Momentum Y-Dir.        : ',MomentumY
            Write (6,*) 'Momentum Z-Dir.        : ',MomentumZ
            Write (6,*)
         Endif
 
 
C     Calculate Averages
 
         If (I.Gt.NumberOfInitSteps) Then
            Av(1) = Av(1) + Tempz
            Av(2) = Av(2) + Pressure
            Av(3) = Av(3) + KineticEnergy
            Av(4) = Av(4) + PotentialEnergy
            Av(5) = Av(5) + TotalEnergy
            Av(6) = Av(6) + 1.0d0

c     update the energy drift
c     declare dutot and initialise it

c     Start Modification

            dutot = dutot + (utot0-TotalEnergy)/utot0

c     End Modification

C     Sample Radial Distribution Function

            If(Mod(I,100).Eq.0) Call Sample_Gyra(2)

            Call Sample_Diff(2)
         Endif
      Enddo
 
C     Print Averages To Screen
 
      If (Av(6).Gt.0.5d0) Then
         Av(6) = 1.0d0/Av(6)
         Do I = 1,5
            Av(I) = Av(I)*Av(6)
         Enddo
 
         Write (6,*)
         Write (6,*) 'Average Temperature     : ',Av(1)
         Write (6,*) 'Average Pressure        : ',Av(2)
         Write (6,*) 'Average KineticEnergy   : ',Av(3)
         Write (6,*) 'Average PotentialEnergy : ',Av(4)
         Write (6,*) 'Average TotalEnergy     : ',Av(5)

c     print the average energy drift
         Write (6,*) 'Average deltaU        : ',dutot*av(6)
      Endif

      Call Sample_Gyra(3)
      Call Sample_Diff(3)

      Return
      End
