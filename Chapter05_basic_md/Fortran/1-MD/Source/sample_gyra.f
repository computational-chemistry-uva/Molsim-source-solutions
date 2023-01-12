      Subroutine Sample_Gyra(Ichoise)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Samples The Radial Distribution Function   C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J,Maxx,Ichoise,A

      Parameter(Maxx = 500)

      Double Precision Ggt,Gg(Maxx),Delta,R2,Dx,Dy,Dz

      Save Ggt,Gg,Delta

      If(Ichoise.Eq.1) Then
         Do I=1,Maxx
            Gg(I) = 0.0d0
         Enddo

         Ggt   = 0.0d0
         Delta = Dble(Maxx-1)/HalfBox
      
      Elseif(Ichoise.Eq.2) Then

         Ggt = Ggt + 1.0d0

C     Loop Over All Pairs

         Do I=1,NumberOfParticles-1
            Do J=I+1,NumberOfParticles

               Dx = Rxx(I) - Rxx(J)
               Dy = Ryy(I) - Ryy(J)
               Dz = Rzz(I) - Rzz(J)
 
               If (Dx.Gt.HalfBox) Then
                  Dx = Dx - Box
               Elseif (Dx.Lt. - HalfBox) Then
                  Dx = Dx + Box
               Endif
 
               If (Dy.Gt.HalfBox) Then
                  Dy = Dy - Box
               Elseif (Dy.Lt. - HalfBox) Then
                  Dy = Dy + Box
               Endif
 
               If (Dz.Gt.HalfBox) Then
                  Dz = Dz - Box
               Elseif (Dz.Lt. - HalfBox) Then
                  Dz = Dz + Box
               Endif

               R2 = Dsqrt(Dx*Dx + Dy*Dy + Dz*Dz)

C     Calculate In Which Bin This Interaction Is In
C
C     Delta = 1/Binsize
C     Maxx  = Number Of Bins

               If(R2.Lt.HalfBox) Then
                  
                  A     = Int(R2*Delta) + 1
                  Gg(A) = Gg(A)           + 2.0d0
               
               Endif
            Enddo
         Enddo

      Else

C     Write Results To Disk

         Ggt   = 1.0d0/(Ggt*Dble(NumberOfParticles))
         Delta = 1.0d0/Delta

         Do I=1,Maxx-1
            R2 = (16.0d0*Datan(1.0d0)/3.0d0)*
     &           (Dble(NumberOfParticles)/(Box**3))*(Delta**3)*
     &           ((Dble(I+1))**3 - (Dble(I)**3))
            
            Write(23,*) ((Dble(I)-0.5d0)*Delta),Gg(I)*Ggt/R2
         Enddo
      Endif

      Return
      End
