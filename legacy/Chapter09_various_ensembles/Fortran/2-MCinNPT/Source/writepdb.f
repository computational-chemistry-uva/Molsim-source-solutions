      Subroutine Writepdb
      Implicit None

      Include 'system.inc'
      
      Integer I,Countmodel
      
      Data Countmodel/ 0/

      Save Countmodel

      Countmodel = Countmodel + 1

      Write(22,'(A,I9)') 'MODEL',Countmodel
      Write(22,'(A6, 3f9.3, 3f7.2, 1x, A)') 'CRYST1',
     &   2.0d0 * Box, 2.0d0 * Box, 2.0d0 * Box,
     &   90.0, 90.0, 90.00, 'P 1         1'

      Do I=1,NumberOfParticles

         Write(22,'(A6,I5,1x,A4,1x,3x,1x,4x,1x,4x,3f8.3,12x,10x,A2)')
     &        'ATOM  ',I,' H  ',
     &        2.0d0 * X(I),2.0d0 * Y(I), 2.0d0 * Z(I),' H' 

      Enddo

      Write(22,'(A)') 'ENDMDL'

      Return
      End
