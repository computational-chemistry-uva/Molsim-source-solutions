      Subroutine Adjust(Attempts, Accepted, MaxDisplacement)
      Implicit None
C
C     Adjusts Maximum Displacement Such That 50% Of The
C     Movels Will Be Accepted
C

      Include 'system.inc'
      
      Integer Attempts, Accepted, SavedAttempts, SavedAccepted
      Double Precision OldDisplacement, Frac, MaxDisplacement
      Save SavedAccepted, SavedAttempts
      Data SavedAccepted/0/
      Data SavedAttempts/0/
 
      If (Attempts.Eq.0.Or.SavedAttempts.Ge.Attempts) Then
         SavedAccepted = Accepted
         SavedAttempts = Attempts
      Else
         Frac = Dble(Accepted-SavedAccepted)/Dble(Attempts-SavedAttempts)
         OldDisplacement  = MaxDisplacement
         MaxDisplacement   = MaxDisplacement*Abs(Frac/0.5d0)

C        ---Limit The Change:

         If (MaxDisplacement/OldDisplacement.Gt.1.5d0) MaxDisplacement = OldDisplacement*1.5d0
         If (MaxDisplacement/OldDisplacement.Lt.0.5d0) MaxDisplacement = OldDisplacement*0.5d0
         If (MaxDisplacement.Gt.HalfBox/2.D0) MaxDisplacement = HalfBox/2.D0
         Write (6, 99001) MaxDisplacement, OldDisplacement, Frac, 
     &        Attempts - SavedAttempts, Accepted - SavedAccepted

C        ---Store Accepted And Attemp For Next Use

         SavedAccepted = Accepted
         SavedAttempts = Attempts
      End If
      Return
99001 Format (' Max. Displ. Set To : ', F6.3, ' (Old : ', F6.3, ')', /, 
     &        ' Frac. Acc.: ', F4.2, ' Attempts: ', I7, ' Succes: ', I7)
      End
