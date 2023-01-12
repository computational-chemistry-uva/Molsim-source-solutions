      Subroutine Writepdb
      Implicit None

      Include 'system.inc'

C     Make A Movie File Of The Simulation Box
C     Use Molmol To View It..

      Integer I,Countmodel,Countatom

      Data Countmodel/ 0/
      Data Countatom/ 0/

      Save Countmodel,Countatom

      Countmodel = Countmodel + 1

      Write(22,'(A,I9)') 'MODEL',Countmodel

      Write(22,'(A)') 'CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1         1'
      Do I=1,NumberOfParticles

         Write(22,'(A,I7,A,I12,4x,3f8.3)') 'ATOM',I,'  H',
     &        I,2.0d0*Rxx(I),2.0d0*Ryy(I),2.0d0*Rzz(I)
      Enddo

      Write(22,'(A)') 'ENDMDL'

      Return
      End
