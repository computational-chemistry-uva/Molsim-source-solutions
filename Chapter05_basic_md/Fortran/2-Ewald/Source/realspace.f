      Subroutine Realspace(Ureal)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccc
C     Real Path; Also Direct Calculation   C
C     Loop Over All Particle Pairs         C
Cccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J
      Double Precision Dx,Dy,Dz,Ureal,R,R2,Ir,Dderfc

      Ureal = 0.0d0

C     Start Modification
c    1. For all particle pairs calculate the distance in x, y, and z.
c    2. Apply periodic boundary conditions where necessary.
c    3. Calculate the real-space contribution to the energy. 
c    (Use the "Dderfc(x)" function to calculate the Error-Function Complement.) 

c     loop over pairs of particles
      Do I=1,(Npart-1)
         Do J=(I+1),Npart

c     calculate distances between the two particles 
c     in the three directions
            Dx = Rx(I) - Rx(J)
            Dy = Ry(I) - Ry(J)
            Dz = Rz(I) - Rz(J)

c     Apply periodic boundary conditions
            If(Dx.Gt.Hbox) Then
               Dx = Dx - Box
            Elseif(Dx.Lt.-Hbox) Then
               Dx = Dx + Box
            Endif

            If(Dy.Gt.Hbox) Then
               Dy = Dy - Box
            Elseif(Dy.Lt.-Hbox) Then
               Dy = Dy + Box
            Endif

            If(Dz.Gt.Hbox) Then
               Dz = Dz - Box
            Elseif(Dz.Lt.-Hbox) Then
               Dz = Dz + Box
            Endif

c     calculate distance
            R2 = Dx*Dx + Dy*Dy + Dz*Dz

c     real space part
            If(R2.Lt.Hbox2) Then
               R     = Dsqrt(R2)
               Ir    = 1.0d0/R
               Ureal = Ureal + Z(I)*Z(J)*Ir*Dderfc(Alpha*R)
            Endif
         Enddo
      Enddo
C     End   Modification

      Return
      End
