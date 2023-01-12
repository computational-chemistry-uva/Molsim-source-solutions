      Program Hard2
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Hard Disks On A Square; Random Placement  C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer Nstep,I,J,n,k,Ipart,Pp,
     &     Ninit,Sstmm,CycleMultiplication,Overlap
      Double Precision Ran_Uniform,R2,NAttempted,NAccepted,Xn,Yn,
     &     Dispmax,M1,DX,DY,BoxSize

      Parameter (BoxSize = 10.0d0)
      Parameter (CycleMultiplication = 1000)
      
Cccccccccccccccccccccccccccccccccccccccccccc
C     Initialize Random Number Generator   C
Cccccccccccccccccccccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call Genrand(M1)

Cccccccccccccccccccccccccccccccccccccccccccc
C     Read Info From Disk                  C
Cccccccccccccccccccccccccccccccccccccccccccc

      Write(*,*) 'How Many Cycles (x ',CycleMultiplication,') ?)'
      Read(*,*)  Nstep

      Write(*,*) 'How Many Init. Cycles ?        (Example: 100      )'
      Read(*,*)  Ninit
      
      Write(*,*) 'How Many Particles    ?        (Always 2 < I < 80 )'
      Read(*,*)  Npart

      Write(*,*) 'Maximum Displacement  ?        (Disp > 0          )'
      Read(*,*)  Dispmax

      Write(*,*) 'Periodic boundary conditions ? (0 or 1)'
      Read(*,*) PBC

      If(Npart.Gt.80.Or.Npart.Lt.2) Stop

      Call Sample(1)

      NAttempted = 0.0d0
      NAccepted = 0.0d0
      Pp   = 0

Ccccccccccccccccccccccccccccccccccc
C     Put Particles On A Lattice  C
Ccccccccccccccccccccccccccccccccccc

      k = 0

      Do I=1,9
         Do J=1,9

            k    = k + 1
            X(k) = Dble(I)*1.1d0
            Y(k) = Dble(J)*1.1d0
         Enddo
      Enddo

      Do n=1,Nstep
         Do k=1,CycleMultiplication

            NAttempted = NAttempted + 1.0d0 

Ccccccccccccccccccccccccccccccccccc
C     Select Particle At Random   C
Ccccccccccccccccccccccccccccccccccc

            Ipart = 1        + Int(Ran_Uniform()*Dble(Npart))
            Xn    = X(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax
            Yn    = Y(Ipart) + (Ran_Uniform()-0.5d0)*Dispmax

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Always Reject When Particle Is Out Of The Box   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
            If((PBC.eq.1).or.((PBC.eq.0).and.(Min(Xn,Yn).Ge.0.0d0)
     &         .and.(Max(Xn,Yn).Lt.BoxSize))) then
          
               Overlap=0

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     See If There Is An Overlap With Any Other Particle  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               Do I=1,Npart
                  If(I.Ne.Ipart) Then
  
                     if(PBC.eq.0) then
                       R2 = (Xn-X(I))**2 + (Yn-Y(I))**2
                     else
                       DX=Xn-X(I)
                       DX=DX-BoxSize*DNINT(DX/BoxSize)
                       DY=Yn-Y(I)
                       DY=DY-BoxSize*DNINT(DY/BoxSize)
                       R2 = DX**2 + DY**2
                     endif

                     If(R2.Lt.1.0d0) Overlap=1
                  Endif
               Enddo

Cccccccccccccccccccccccccccccccccc
C     No Overlaps, So Accepted   C
Cccccccccccccccccccccccccccccccccc

               If(Overlap.eq.0) then
                  NAccepted    = NAccepted + 1.0d0
                  X(Ipart) = Xn
                  Y(Ipart) = Yn
               Endif
            Endif

            If(Mod(k,5).Eq.0.And.n.Gt.Ninit) Then
               Call Sample(2)
               Pp = Pp + 1
            Endif
         Enddo
         If(Mod(n,5).Eq.0) Call Writepdb
      Enddo

      If(Pp.Gt.4) Call Sample(3)
      Write(6,*) 'Fraction Success : ',NAccepted/NAttempted

      Stop
      End
