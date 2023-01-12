      Subroutine Eneri(Xi, Yi, Zi, I, Jb, En, Vir, BoxID)
      Implicit None

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'
 
      Double Precision Xi, Yi, Zi, En, Dx, Dy, Dz, R2, Vir, Virij, Enij
      Integer I, J, Jb, BoxID

      En = 0.D0
      Vir = 0.D0
      Do J = Jb, Npart
         If (Id(J).Eq.BoxID) Then
            If (J.Ne.I) Then
               Dx = Xi - X(J)
               Dy = Yi - Y(J)
               Dz = Zi - Z(J)
               If (Dx.Gt.Hbox(BoxID)) Then
                  Dx = Dx - Box(BoxID)
               Else
                  If (Dx.Lt.-Hbox(BoxID)) Dx = Dx + Box(BoxID)
               End If
               If (Dy.Gt.Hbox(BoxID)) Then
                  Dy = Dy - Box(BoxID)
               Else
                  If (Dy.Lt.-Hbox(BoxID)) Dy = Dy + Box(BoxID)
               End If
               If (Dz.Gt.Hbox(BoxID)) Then
                  Dz = Dz - Box(BoxID)
               Else
                  If (Dz.Lt.-Hbox(BoxID)) Dz = Dz + Box(BoxID)
               End If
               R2 = Dx*Dx + Dy*Dy + Dz*Dz
               Call Ener(Enij, Virij, R2, BoxID)
               En = En + Enij
               Vir = Vir + Virij
            End If
         End If
      End Do

      Return
      End
