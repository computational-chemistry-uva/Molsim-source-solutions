      Subroutine Grow(LOldChain,Weight,Ubonded,Unonb,Xst,Yst,Zst)
      Implicit None
Cccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Grow A Chain Using Cbmc Or Random Insertion  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccc

      Include 'system.inc'

      Logical LOldChain,LValidTrial
      Integer I,J,K,Ntrial,Iwalk

      Double Precision Ubonded,Unonb,Weight,Xst(Maxchain),
     &     Yst(Maxchain),Zst(Maxchain),Xtrial(Maxtrial),
     &     Ytrial(Maxtrial),Ztrial(Maxtrial),Dx1,Dx2,Dy1,
     &     Dy2,Dz1,Dz2,Ubend(Maxtrial),Uext(Maxtrial),X,Y,
     &     Z,Ws,RandomNumber,Cumw,Cmmm,U,Bfac(Maxtrial)

      Ubonded = 0.0d0
      Unonb   = 0.0d0
      Weight  = 1.0d0
      Ws      = 0.0d0

      If(Lcbmc) Then
         Ntrial = NumberOfTrialPositions
      Else
         Ntrial = 1
      Endif

C     the first bead of the trial chain is always at (0,0,0)
      Xst(1) = 0.0d0
      Yst(1) = 0.0d0
      Zst(1) = 0.0d0

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Grow The Chain; Loop Over All Segments Exept The First One  C
C                                                                 C
C     Xpos/Ypos/Zpos       : Position Of The Old Chain            C
C     Xst/Yst/Zst          : Position Of The Trial Chain          C
C     Xtrial/Ytrial/Ztrial : Position Of A Trial Segment          C
C     Lcbmc                : Do We Use Cbmc (.True. Or .False.)   C
C     LOldChain                 : Do We Have The Old Configuration     C
C                            (Only For Dynamic Schemes)           C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Loop over all the particles in the polymer chain 
      Do I=2,ChainLength

C     Loop Over All Trial Positions   C
         Do J=1,Ntrial

 10         Ubend(J) = 0.0d0
            Uext(J)  = 0.0d0
            LValidTrial   = .False.
            Ws       = 0.0d0
            
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Generate Trial Position                              C
C     J=1 And LOldChain : Old Position                          C
C     Otherwise    : New Position                          C
C                                                          C
C     Old Coordinates Are Stored In Xpos/Ypos/Zpos         C
C     The New Chain Is Stored In Xst/Yst/Zst               C
C     Trial Positions Are Stored In Xtrial/Ytrial/Ztrial   C
C                                                          C
C     Subroutine Ran_Sphere(X,Y,Z): Creates A Random Unit  C
C     Vector On A Sphere centered at (0,0,0)               C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C    if an old chain exists, we must retrace it to obtain its statistical weight; in that case the position of bead I in the old chain is also used as the first trial position of bead I

C    Note: retracing the old chain is not indispensable in our example because
c    we only have one chain; instead we could have just saved its statistical
C    weight from the previous step.
C    However, in a real CBMC implementation with more than one chain, the
C    environment of the chain may have changed since the previous time its
C    statistical weight was calculated; in that case, retracing is necessary. To
C    keep the algorithm general, CBMC is always implemented with retracing of the old
C    chain.
 
            If(J.Eq.1.And.LOldChain) Then
               Xtrial(J) = Xpos(I)
               Ytrial(J) = Ypos(I)
               Ztrial(J) = Zpos(I)
            Else  
C     Otherwise, generate a new trial position for bead I
C     Start Modification
               Call Ran_Sphere(X,Y,Z)
               
               Xtrial(J) = Xst(I-1) + X
               Ytrial(J) = Yst(I-1) + Y
               Ztrial(J) = Zst(I-1) + Z

C     End Modification
            Endif

Cccccccccccccccccccccccccccccccccccccccccc
C     Calculate Bond-Bending Potential   C
Cccccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.2) Then
C              Vector from bead I-1 to trial position J
               Dx1 = Xtrial(J) - Xst(I-1)
               Dy1 = Ytrial(J) - Yst(I-1)
               Dz1 = Ztrial(J) - Zst(I-1)
                                         
C              Vector from bead I-1 to bead I-2 
               Dx2 = Xst(I-2)  - Xst(I-1)
               Dy2 = Yst(I-2)  - Yst(I-1)
               Dz2 = Zst(I-2)  - Zst(I-1)

C              Calculate angular potential energy  associated with 3 consecutive beads:  I-2 .....  I-1 ..... J
               Call Angle(Dx1,Dy1,Dz1,Dx2,Dy2,Dz2,U)

               Ubend(J) = U
            Endif

Cccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Accept Or Reject The Chosen Configuration      C
C     This Is For Cmbc Only, Because Only In Cbmc    C
C     A Trial Position According To A Bond-Bending   C
C     Potential Has To Be Generated...               C
C                                                    C
C     When J=1.And.LOldChain : Always Accepted !! (Why?)  C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccc

            If(Lcbmc) Then

c              old chain is always accepted
               If(J.Eq.1.And.LOldChain) Then
                  LValidTrial = .True.
               Else
C              Start Modification
C              Accept or reject the other trial positions according to the Metropolis algorithm
                  If(RandomNumber().Lt.Dexp(-Beta*Ubend(J))) 
     &                 LValidTrial = .True.
 
C              End Modification
               Endif

C           If CMBC is not used, always accept  the single trial position that was randomly generated.
            Else
               LValidTrial = .True.
            Endif

            If(.Not.LValidTrial) Goto 10

Ccccccccccccccccccccccccccccccccccccccc
C     Calculate The External Energy   C
C     Only When There At Least 4      C
C     Beads...                        C
Ccccccccccccccccccccccccccccccccccccccc

            If(I.Gt.3) Then
               Do K=1,(I-3)

                  Dx1 = Xtrial(J) - Xst(K)
                  Dy1 = Ytrial(J) - Yst(K)
                  Dz1 = Ztrial(J) - Zst(K)
                  
                  Call Repulsion(Dx1,Dy1,Dz1,U)

                  Uext(J) = Uext(J) + U
               Enddo
            Endif

C        End Loop Over All Trial Positions           C
         Enddo

Ccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Choose A Trial Position For Cbmc                 C
C     For Cbmc And An Old Config.: take 1st Config.    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccc

         Iwalk = 1

         If(Lcbmc) Then
            Do K=1,Ntrial
               Bfac(K) = Dexp(-Beta*Uext(K))
               Ws      = Ws + Bfac(K)
            Enddo

            Cumw  = Bfac(1)
            Cmmm  = Ws*RandomNumber()

            If(.Not.LOldChain) Then
               Do While(Cumw.Lt.Cmmm)
                  Iwalk = Iwalk + 1
                  Cumw  = Cumw  + Bfac(Iwalk)
               Enddo
            Endif
         Else
            Ws = Dexp(-Beta*(Uext(1)+Ubend(1)))
         Endif

Ccccccccccccccccccccccccccccccccccccccccccc
C     Store The Chosen Configuration      C
C     Update Energies/Rosenbluth Weight   C
Ccccccccccccccccccccccccccccccccccccccccccc

         Xst(I) = Xtrial(Iwalk)
         Yst(I) = Ytrial(Iwalk)
         Zst(I) = Ztrial(Iwalk)

         Ubonded     = Ubonded + Ubend(Iwalk)
         Unonb       = Unonb   + Uext(Iwalk)
         Weight      = Weight*Ws/Dble(Ntrial)

C     END loop over all particles in chain
      Enddo

      Return
      End
