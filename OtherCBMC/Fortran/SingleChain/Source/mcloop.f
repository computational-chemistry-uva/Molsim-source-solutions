      Subroutine Mcloop
C     Generate multiple configurations of a polymer chain and calculate its averaged properties (energy, end-to-end distance) 
      Implicit None

      Include 'system.inc'

C     Local variables
      Integer I,J,K

      Double Precision Xst(Maxchain),Yst(Maxchain),
     &     Zst(Maxchain),Weight,Weiold,Ubnew,
     &     Ubold,Unnew,Unold,R,RandomNumber,
     &     Bn,Bs,AvUbond,AvUnbond,NoConfs,Eb,En

C     Bn:              number of attempts to create a chain
C     Bs:              number of accepted chains 
C     AvUbond:         average bonded energy
C     AvUnbond:        average non bonded energy
C     NoConfs:         number of configurations used to compute averages 
C     Ubold, Ubnew:    bonded energy of the old and trial chains
C     Unoldn Unnew:    non bonded energy of the old and trial chains
C     Weiold, Weight:  statistical weight of the old and trial chains.  If random chains are generated, both values are 1.
C     R:               end-to-end distance
C     Eb,En:           Bonded and non bonded energy of the  chain 

C     initiallize variables
      Weight = 1.0d0
      Ubnew  = 0.0d0
      Unnew  = 0.0d0
      Weiold = 1.0d0
      Ubold  = 0.0d0
      Unold  = 0.0d0
      Bn     = 0.0d0
      Bs     = 0.0d0
      AvUbond    = 0.0d0
      AvUnbond    = 0.0d0
      NoConfs    = 0.0d0
      
C     Generate the very first Chain  
      Call Grow(.False.,Weiold,Ubold,Unold,Xst,Yst,Zst)

      Do J=1,ChainLength
         Xpos(J) = Xst(J)
         Ypos(J) = Yst(J)
         Zpos(J) = Zst(J)
      Enddo

      Eb = Ubold
      En = Unold

C     Initialize the subroutine "Sample" (the subroutine that will sample the average end-to-end distance)
      Call Sample(1,1.0d0,1.0d0)

C     Generate the Markov chain of states 
      Do I=1,NumberOfSteps
         if (MOD(I,20) .EQ. 0) THEN
           print*,"cycle",I
         ENDIF
         Do K = 1,1000

               Bn = Bn + 1.0d0

C              retrace old chain 
               Call Grow(.True. ,Weiold,Ubold,Unold,Xst,Yst,Zst)

C              create a new chain
               Call Grow(.False.,Weight,Ubnew,Unnew,Xst,Yst,Zst)

C              If the new chain is accepted, then Copy Coordinates   
               If(RandomNumber().Lt.(Weight/Weiold)) Then

                  Bs = Bs + 1.0d0
                  En = Unnew
                  Eb = Ubnew

                  Do J=1,ChainLength
                     Xpos(J) = Xst(J)
                     Ypos(J) = Yst(J)
                     Zpos(J) = Zst(J)
                  Enddo
               Endif

               Weight = 1.0d0

C           If the simulation is beyond the set number of equilibration steps,
C           sample End-To-End Distance And Average Energies     
            If(I.Gt.NumberOfInitializationSteps) Then
               R = Dsqrt((Xpos(ChainLength)-Xpos(1))**2 +
     &                 (Ypos(ChainLength)-Ypos(1))**2 +
     &                 (Zpos(ChainLength)-Zpos(1))**2)
               Call Sample(2,R,Weight)
               AvUbond = AvUbond + Eb
               AvUnbond = AvUnbond + En

C              update number of configurations sampled  
               NoConfs = NoConfs + 1.0d0

C              Periodically write a pdb file containing the particle coordinates 
               If(K.Eq.1.And.Mod(I,5).Eq.0) Call Writepdb
            Endif
         Enddo
      Enddo

C     Call subroutine "Sample" to output the histogram of the end-to-end distance
      Call Sample(3,1.0d0,1.0d0)

C    Output average energies and fraction of accepted chains 
      Write(6,*) 'Fraction Accepted : ',Bs/Bn
      Write(6,*) 'Average U-Bonded  : ',AvUbond/NoConfs
      Write(6,*) 'Average U-Nonb    : ',AvUnbond/NoConfs

      Return
      End
