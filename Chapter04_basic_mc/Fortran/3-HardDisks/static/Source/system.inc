Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     X,Y   = Coordinates Of Particles            C
C     Npart = Number Of Particles (Maximum 100)   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

      Double Precision X(100),Y(100), BoxSize
      Integer Npart,PBC

      Parameter (BoxSize = 10.0d0)

      Common /Sstt/ X,Y,Npart,PBC
