      Subroutine Sample_Energy(Switch)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Sample The Total Energy Distribution   C
Cccccccccccccccccccccccccccccccccccccccccccccc


      Integer          Switch,I,J,NUMBEROFBINS
      Parameter (NUMBEROFBINS=1000)
      Double Precision Gg1(NUMBEROFBINS,MAXNUMBEROFSYSTEMS),Gg2,Emax,Deltax
      Character*5      Name

      Parameter (Emax = 10.0d0)

      Save Gg1,Gg2,Deltax

      If(Switch.Eq.1) Then
         
         Do I=1,NUMBEROFBINS
            Do J=1,MAXNUMBEROFSYSTEMS
               Gg1(I,J) = 0.0d0
            Enddo
         Enddo

         Gg2    = 0.0d0
         Deltax = Dble(NUMBEROFBINS)/Emax 

      Elseif(Switch.Eq.2) Then

         Gg2 = Gg2 + 1.0d0

         Do J=1,Ntemp
            I = 1 + Int(Deltax*Uold(J))

            If(I.Le.NUMBEROFBINS) Gg1(I,J) = Gg1(I,J) + 1.0d0
         Enddo

      Else

         Gg2    = 1.0d0/Gg2
         Deltax = 1.0d0/Deltax

         Do J=1,Ntemp
            Name = 'Ener'//Char(J+48)

            Open(32,File=Name)

            Do I=1,NUMBEROFBINS
               If(Gg1(I,J).Gt.0.5d0) 
     &              Write(32,*) Dble(I-1)*0.01d0,Gg1(I,J)*Gg2
            Enddo

            Close(32)
         Enddo
      Endif

      Return
      End
