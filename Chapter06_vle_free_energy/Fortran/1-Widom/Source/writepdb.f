      Subroutine Writepdb(FileNum)
      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'

      Integer I,Countmodel,Countatom,FileNum

      Data Countmodel/ 0/
      Data Countatom/ 0/

      Save Countmodel,Countatom

      Countmodel = Countmodel + 1

      Write(FileNum,'(A,I9)') 'MODEL',Countmodel

      Do I=1,NumberOfParticles

         Countatom = Countatom + 1

         Write(FileNum,'(A,I7,A,I12,4x,3f8.3)') 'ATOM',Countatom,'  O',
     &        Countatom,X(I),Y(I),Z(I)
      Enddo

      Write(FileNum,'(A)') 'ENDMDL'

      Return
      End
