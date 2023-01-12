      Subroutine Store(MaxDisplacement)
C
C     Writes Configuration To Disk
C

      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer I
      Double Precision MaxDisplacement

      Open (43,FILE='restart.dat',FORM='FORMATTED')
 
      Write (43, *) Box, HalfBox
      Write (43, *) NumberOfParticles
      Write (43, *) MaxDisplacement

      Do I = 1, NumberOfParticles
         Write (43, *) X(I), Y(I), Z(I)
      End Do
      
      Close (43)

      Return
      End
