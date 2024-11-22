      Subroutine Integrate(DeltaUtotal)
      Implicit None
 
      Include 'system.inc'
 
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Integrate The Equations Of Motion For An Nve System  C
C     Velocity Verlet Integrator                           C
C                                                          C
C     DeltaUtotal Is The Energy Drift                             C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I
      Double Precision U,F,Consv,Utotal0,DeltaUtotal,DeltaUcounter
 
      Call Force(Xpos,U,F)

      DeltaUtotal = 0.0d0
      DeltaUcounter = 0.0d0
      Utotal0 = 0.0d0

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Integrate The E.O.M. For Nstep Steps   C
Cccccccccccccccccccccccccccccccccccccccccccccc

      Do I=1,Nstep
         Xpos = Xpos             + Vpos*Tstep      + 0.5d0*F*Tstep*Tstep
         Vpos = Vpos             + 0.5d0*Tstep*F

         Call Force(Xpos,U,F)

         Vpos  = Vpos            + 0.5d0*Tstep*F
         Consv = 0.5d0*Vpos*Vpos + U

Ccccccccccccccccccccccccccccccccccccccccccccccc
C     Used For Calculating The Energy Drift   C
C     Not Very Important Really...            C
Ccccccccccccccccccccccccccccccccccccccccccccccc

         If(I.Eq.1) Then
            Utotal0 = Consv
         Else
            If(.Not.(I.Gt.5.And.
     &               ((Xpos.Lt.-1.0d0.And.Vpos.Lt.0.0d0).Or.
     &                (Xpos.Gt.1.0d0.And.Vpos.Gt.0.0d0)))) Then

               DeltaUtotal = DeltaUtotal + Dabs((Consv-Utotal0)/Utotal0)
               DeltaUcounter = DeltaUcounter + 1.0d0
            Endif
         Endif

         Call Sample(2,I)

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Check If The Particle Passed The Boundary Or If It Can Never Recross The   C
C     Boundary Again...                                                          C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      Enddo

      If(DeltaUcounter.Gt.0.5d0) Then
         DeltaUtotal = DeltaUtotal/DeltaUcounter
      Else
         DeltaUtotal = 0.0d0
      Endif

      Return
      End
