      Subroutine Sample_Prof(Switch)
      Implicit None

      Include 'system.inc'

      Integer          Maxx,I,Switch
      Parameter(Maxx=500)

      Double Precision Dtot,Dens(Maxx),Delta
      Save Dtot,Dens,Delta

      If (Switch.Eq.1) Then

         Delta = Dble(Maxx-1)/4.0d0
         Dtot  = 0.0d0

         Do I=1,Maxx
            Dens(I) = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then
         Dtot = Dtot + 1.0d0
         I       = 1 + Int(Delta*(Xpos + 2.0))
         If(I.Ge.1.And.I.Lt.Maxx) Then
            Dens(I) = Dens(I) + 1.0d0
         Endif
      Else
         Do I=1,Maxx
            Write(27,*) ((Dble(I-1)/Delta) - 2.0d0),Dens(I)/Dtot
         Enddo
      Endif

      Return
      End
