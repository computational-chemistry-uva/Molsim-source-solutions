      Program Hard1
      Implicit None

      Include 'system.inc'

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Hard Disks On A Square; Random Placement  C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer          Nstep,I,J,n,k,Sstmm,Overlap,CycleMultiplication
      Integer          BoxSize
      Double Precision Ran_Uniform,R2,NAttempted,NAccepted,M1,DX,DY
      
      Parameter (BoxSize=10.0d0)
      Parameter (CycleMultiplication=1000)

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

      Write(*,*) 'How Many Cycles (x ',CycleMultiplication,')?'
      Read(*,*)  Nstep
      
      Write(*,*) 'How Many Particles ?          (Always 2 < I < 80 )'
      Read(*,*)  Npart

      Write(*,*) 'Periodic boundary conditions ? (0 or 1)' 
      Read(*,*) PBC

      If(Npart.Gt.80.Or.Npart.Lt.2) Stop

      Call Sample(1)

      NAttempted = 0.0d0
      NAccepted = 0.0d0

      Do n=1,Nstep
         Do k=1,CycleMultiplication

            NAttempted = NAttempted + 1.0d0 

            Overlap = 0

            Do I=1,Npart
               If(Overlap.eq.0) then

Ccccccccccccccccccccc
C     New Particle  C
Ccccccccccccccccccccc

                  X(I) = Ran_Uniform()*BoxSize
                  Y(I) = Ran_Uniform()*BoxSize

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Check Distances With Previous Placed Particles  C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                  Do J=1,(I-1)
                     If(Overlap.eq.0) then
                        if(PBC.eq.0) then
                          R2 = (X(I)-X(J))**2 + (Y(I)-Y(J))**2
                        else
                          DX=X(I)-X(J)
                          DX=DX-BoxSize*DNINT(DX/BoxSize)
                          DY=Y(I)-Y(J)
                          DY=DY-BoxSize*DNINT(DY/BoxSize)
                          R2 = DX**2 + DY**2
                        endif

                        If(R2.Lt.1.0d0) Overlap = 1
                     Endif
                  Enddo
               Endif
            Enddo

Ccccccccccccccccccccccccccccccccc
C     No Overlaps Occured !!!   C
Ccccccccccccccccccccccccccccccccc

            If(Overlap.eq.0) then
               NAccepted = NAccepted + 1.0d0

               Call Sample(2)
            Endif

         Enddo
      Enddo

      If(NAccepted.Ge.0.5d0) Call Sample(3)
      Write(6,*) 'Fraction Succes : ',NAccepted/NAttempted

      Stop
      End
