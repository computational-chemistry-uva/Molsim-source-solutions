      Subroutine Force
      Implicit None
 
C     Calculate The Forces And Potential Energy
 
      Include 'system.inc'
 
      Integer I,J
      Double Precision Dx,Dy,Dz,F,R2i,R6i
 
C     Set Forces, Potential Energy And Pressure To Zero
 
      Do I = 1,NumberOfParticles
         Fxx(I) = 0.0d0
         Fyy(I) = 0.0d0
         Fzz(I) = 0.0d0
      Enddo
 
      PotentialEnergy = 0.0d0
      Pressure = 0.0d0
 
C     Loop Over All Particle Pairs
 
      Do I = 1,NumberOfParticles - 1
         Do J = I + 1,NumberOfParticles
 
C     Calculate Distance And Perform Periodic
C     Boundary Conditions
 
            Dx = Rxx(I) - Rxx(J)
            Dy = Ryy(I) - Ryy(J)
            Dz = Rzz(I) - Rzz(J)
 
            If (Dx.Gt.HalfBox) Then
               Dx = Dx - Box
            Elseif (Dx.Lt. - HalfBox) Then
               Dx = Dx + Box
            Endif
 
            If (Dy.Gt.HalfBox) Then
               Dy = Dy - Box
            Elseif (Dy.Lt. - HalfBox) Then
               Dy = Dy + Box
            Endif
 
            If (Dz.Gt.HalfBox) Then
               Dz = Dz - Box
            Elseif (Dz.Lt. - HalfBox) Then
               Dz = Dz + Box
            Endif
 
            R2i = Dx*Dx + Dy*Dy + Dz*Dz

C     Check If The Distance Is Within The Cutoff Radius
 
            If (R2i.Lt.Rcutsq) Then
               R2i = 1.0d0/R2i
               R6i = R2i*R2i*R2i
 
               PotentialEnergy = PotentialEnergy + 4.0d0*R6i*(R6i-1.0d0) - Ecut
               F = 48.0d0*R6i*(R6i - 0.5d0)
               Pressure = Pressure + F
               F = F*R2i
 
               Fxx(I) = Fxx(I) + F*Dx
               Fyy(I) = Fyy(I) + F*Dy
               Fzz(I) = Fzz(I) + F*Dz
 
               Fxx(J) = Fxx(J) - F*Dx
               Fyy(J) = Fyy(J) - F*Dy
               Fzz(J) = Fzz(J) - F*Dz
 
            Endif
         Enddo
      Enddo
 
C     Scale The Pressure
 
      Pressure = Pressure/(3.0d0*Box*Box*Box)
 
      Return
      End
