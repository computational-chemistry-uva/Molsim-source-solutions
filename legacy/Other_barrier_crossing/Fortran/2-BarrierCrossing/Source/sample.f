      Subroutine Sample(Switch,Ttt)
      Implicit None

      Include 'system.inc'

      Integer          Maxx,I,Ttt,Switch
      Parameter(Maxx=50000)

      Double Precision Kt(Maxx),Ks(Maxx)
      Save Kt,Ks

      If (Switch.Eq.1) Then

C     Initialize Everything

         Do I=1,Maxx
            Kt(I) = 0.0d0
            Ks(I) = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Compute the Transmission Coeffient                  C
C     see also Case Study 23 on page 440/441 (equation c) C
C                                                         C
C     The Maximum Time Equals Maxx                        C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Begin Modification
         If(Ttt.Le.Maxx) Then

            If(Xpos.Gt.Qstar) Then
               Kt(Ttt) = Kt(Ttt) + Theta
            Endif

            Ks(Ttt) = Ks(Ttt) + Dabs(Theta)
         Endif

C     End Modification

      Else

C     Write Results

         Do I=1,Maxx
            If(I.Le.Nstep) Then

               If(Ks(I).Gt.0.5d0) Write(27,*) Tstep*Dble(I-1),
     &             2.0d0*Kt(I)/Ks(I)

            Endif
         Enddo

      Endif

      Return
      End
