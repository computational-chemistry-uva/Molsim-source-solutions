      Subroutine Sample(Option)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Samples The Radial Distribution Function   C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J,Maxx,Option,A

      Parameter(Maxx = 500)

      Double Precision Ngr,Gg(Maxx),Delta,R2,DX,DY,R,ab,nid

      Save Ngr,Gg,Delta

      If(Option.Eq.1) Then
         Do I=1,Maxx
            Gg(I) = 0.0d0
         Enddo

         Ngr   = 0.0d0
         Delta = 5.0d0/Dble(Maxx-1)
      
      Elseif(Option.Eq.2) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sample The Radial Distribution Function   C
C     Loop Over All Particle Pairs              C
C     See Frenkel/Smit P. 86                    C
C                                               C
C     Delta  = Binsize                          C
C     Ggt/Gg = Counters                         C
C     Maxx   = Maximum Number Of Bins           C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification
         Ngr = Ngr + 1.0d0

         Do I=1,Npart-1
            Do J=I+1,Npart

c     calculate distance between two particles
               if(PBC.eq.0) then
                  R2 = (X(I)-X(J))**2 + (Y(I)-Y(J))**2
               else
                  DX=X(I)-X(J)
                  DX=DX-BoxSize*DNINT(DX/BoxSize)
                  DY=Y(I)-Y(J)
                  DY=DY-BoxSize*DNINT(DY/BoxSize)
                  R2 = DX**2 + DY**2
               endif

               If(R2.Lt.((BoxSize/2)*(BoxSize/2))) Then

                  R2    = Dsqrt(R2)
c     calculates bin
                  A     = Int(R2/Delta) + 1
c     add the TWO particles (both particle i and j)
                  Gg(A) = Gg(A)           + 2.0d0
               Endif
            Enddo
         Enddo
               
C     End   Modification

      Else

Ccccccccccccccccccccccccccccccc
C     Write Results To Disk   C
Ccccccccccccccccccccccccccccccc

         Do I=1,Maxx-1

            r=Delta*(I+0.5d0)
            ab=((i+1)**2-i**2)*Delta**2
            nid= 4.0d0*Datan(1.0d0)*Dble(Npart-1)/(BoxSize*BoxSize)*ab
            Gg(I)=Gg(I)/(Ngr*Dble(Npart)*nid)
         
            Write(21,*) r,Gg(I)
         Enddo
      Endif

      Return
      End
