      Subroutine Repulsion(Dx,Dy,Dz,U)
C     Calculate repulsive potential energy between two beads
      Implicit None
      Include 'system.inc'

C     Global Variables   
      Double Precision Dx,Dy,Dz,U
C       Dx,Dy,Dz: vector connecting the two beads
C       U:        repulsive energy between the two beads   

C     Local variables
      Double Precision R 
C       R: distance between two beads


      R = Dsqrt(Dx*Dx + Dy*Dy + Dz*Dz)
      U = 0.0d0

C     If the distance between the two beads is less than cutoff, calculate the repulsive energy.  
C     Otherwise, the repulsive energy is 0.
      If(R.Lt.Rcut) U = A*((R-Rcut)**2)

      Return
      End
