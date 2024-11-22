      Subroutine EnergyParticle(Xi, Yi, Zi, I, Jb, Energy, Virial)
C
C    Calculates The Energy Of Particle I With Particles J=Jb,NumberOfParticles
C
C     Xi (Input)        X Coordinate Particle I
C     Yi (Input)        Y Coordinate Particle I
C     Zi (Input)        Z Coordinate Particle I
C     I  (Input)        Particle Number (Excluded !!!)
C     Energy  (Output)  Energy Particle I
C     Virial (Output)   Virial Particle I
C
C
      Implicit None
      
      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'
      Include 'potential.inc'

      Double Precision Xi,Yi,Zi,Energy,Dx,Dy,Dz,R2,R2i,R6i,Virial
      Integer I,J,Jb

      Energy = 0.0d0
      Virial = 0.0d0

      Do J = Jb, NumberOfParticles

Ccccccccccccccccccccccccccc
C     Excluse Particle I  C
Ccccccccccccccccccccccccccc

        If (J.Ne.I) Then
          
          Dx = Xi - X(J)
          Dy = Yi - Y(J)
          Dz = Zi - Z(J)
          
Ccccccccccccccccccccccccccccccccccc
C     Nearest Image Convention    C
Ccccccccccccccccccccccccccccccccccc

          If (Dx.Gt.HalfBox) Then
            Dx = Dx - Box
          Elseif (Dx.Lt.-HalfBox) Then
            Dx = Dx + Box
          End If
          
          If (Dy.Gt.HalfBox) Then
            Dy = Dy - Box
          Elseif (Dy.Lt.-HalfBox) Then
            Dy = Dy + Box
          End If
          
          If (Dz.Gt.HalfBox) Then
            Dz = Dz - Box
          Elseif (Dz.Lt.-HalfBox) Then
            Dz = Dz + Box
          End If
          
          R2 = Dx*Dx + Dy*Dy + Dz*Dz
          
Cccccccccccccccccccccccccccccc
C     Calculate The Energy   C
Cccccccccccccccccccccccccccccc

          If (R2.Lt.Rc2) Then
            R2i = Sig2/R2
            R6i = R2i*R2i*R2i
            Energy = Energy + Eps4*(R6i*R6i-R6i)
            Virial = Virial + 12.0*Eps4*(R6i*R6i-0.5*R6i)
          End If
        Endif
      Enddo
      
      Return
      End
