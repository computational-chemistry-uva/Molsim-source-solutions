      Subroutine Angle(Dx1,Dy1,Dz1,Dx2,Dy2,Dz2,U)
C     Calculate angular potential energy  associated with 3 consecutive beads:  I-2 .....  I-1 ..... I
      Implicit None
      Include 'system.inc'

C     Global variables
      Double Precision U,Dx1,Dx2,Dy1,Dy2,Dz1,Dz2
C       U:            angular potential energy
C       Dx1,Dy1,Dz1:  vector from I-1 to I
C       Dx2,Dy2,Dz2:  vector from I-1 to I-2

C     Local variables
      Double Precision R
C       R:            cosine of the angle formed by the three beads

C     Calculate the cosine R of the angle formed by the  vectors D1 and D2 by calculating their dot product: 
C          D1.D2 = |D1| |D2| cos(D1^D2) <=> 
C     <=> cos(D1^D2) = D1.D2/(|D1|*|D2|)
      R =        (Dx1*Dx2 + Dy1*Dy2 + Dz1*Dz2)/
     &     Dsqrt((Dx1*Dx1 + Dy1*Dy1 + Dz1*Dz1)*
     &           (Dx2*Dx2 + Dy2*Dy2 + Dz2*Dz2))

C    Calculate the potential energy associated with that angle
C      Note: the function Dacos(R) returns the angle between 0 and Pi associated with a cosine value of R
      U = 0.5d0*Kt*((Dacos(R)-Theta0)**2)

      Return
      End
