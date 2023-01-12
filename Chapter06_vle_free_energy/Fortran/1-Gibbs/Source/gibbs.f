      Program Gibbs
      Implicit None

C     ---Gibbs-Ensemble Simulation Of The Lennard-Joned Fluid

      Integer Equil, Prod, Nsamp, I, Icycl, Ndispl, Attempt, 
     &        Nacc, Ncycl, Nmoves, Imove, Nvol, AcceptVolume, AttemptVolume, BoxID, Nswap, 
     &        AcceptSwap, AttemptSwap, Sstmm
      Double Precision En(2), Ent(2), Vir(2), Virt(2), Vmax, Dr, 
     &                 Ran, Succ,RandomNumber,M1

      Include 'parameter.inc'
      Include 'chem.inc'
      Include 'conf.inc'
      Include 'potential.inc'
      Include 'system.inc'
 
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)

      Write (6, *) '**************** GIBBS ***************'
      Write(66,'(7A20)')"#Time","En. box 1","En. box 2","Press. box 1",
     & "Press. box 2","Dens. (N/V) box 1", "Dens. (N/V) box 2"
      write(45,'(3A20)')"# Time","Dens. (N/V) box 1","Dens. (N/V) box 2"


C     ---Initialize Sysem

      Call Readdat(Equil, Prod, Nsamp, Ndispl, Dr, Nvol, Vmax, Nswap, 
     &             Succ)
      Nmoves = Ndispl + Nvol + Nswap

C     ---Total Energy Of The System

      Do BoxID = 1, 2
         Call Toterg(En(BoxID), Vir(BoxID), BoxID)
         Write (6, 99001) BoxID, En(BoxID), Vir(BoxID)
      End Do

C     ---Start Mc-Cycle

      Do I = 1, 2

C        --- I=1 Equilibration
C        --- I=2 Production

         If (I.Eq.1) Then
            Ncycl = Equil
            If (Ncycl.Ne.0) Write (6, *) ' Start Equilibration '
         Else
            If (Ncycl.Ne.0) Write (6, *) ' Start Production '
            Ncycl = Prod
         End If
         Attempt = 0
         Nacc = 0
         AttemptVolume = 0
         AcceptVolume = 0
         AttemptSwap = 0
         AcceptSwap = 0

C        ---Initialize Calculation Chemical Potential

         Call Init_Chem(0)

C        ---Intialize The Subroutine That Adjust The Maximum Displacement

         Call Adjust(Attempt, Nacc, Dr, AttemptVolume, AcceptVolume, Vmax, Succ)

         Do Icycl = 1, Ncycl
            Do Imove = 1, Nmoves
               Ran = RandomNumber()*Dble(Ndispl+Nvol+Nswap)
               If (Ran.Lt.Dble(Ndispl)) Then

C                 ---Attempt To Displace A Particle

                  Call Mcmove(En, Vir, Attempt, Nacc, Dr)
               Else If (Ran.Lt.Dble(Ndispl+Nvol)) Then

C                 ---Attempt To Change The Volume

                  Call Mcvol(En, Vir, AttemptVolume, AcceptVolume, Vmax)
               Else

C                 ---Attemp To Exchange Particles

                  Call Mcswap(En, Vir, AttemptSwap, AcceptSwap)
 
               End If
            End Do
            If (I.Eq.2) Then

C              ---Sample Averages

               If (Mod(Icycl,Nsamp).Eq.0) Call Sample(Icycl, En, Vir)
            End If
            If (Mod(Icycl,Ncycl/5).Eq.0) Then
               Write (6, *) '======>> Done ', Icycl, ' Out Of ', Ncycl, 
     &                      Npbox(1), Npbox(2)

C              ---Write Intermediate Configuration To File

               Call Store(8, Dr, Vmax)

C              ---Adjust Maximum Displacements

               Call Adjust(Attempt, Nacc, Dr, AttemptVolume, AcceptVolume, Vmax, Succ)
            End If
         End Do
         If (Ncycl.Ne.0) Then
            If (Attempt.Ne.0) Write (6, 99003) Attempt, Nacc, 
     &                               100.*Float(Nacc)/Float(Attempt)
            If (AttemptVolume.Ne.0) Write (6, 99004) AttemptVolume, AcceptVolume, 100.*Float(AcceptVolume)
     &                            /Float(AttemptVolume)
            If (AttemptSwap.Ne.0) Write (6, 99005) AttemptSwap, AcceptSwap, 
     &                             100.*Float(AcceptSwap)/Float(AttemptSwap)
            Do BoxID = 1, 2

C              ---Test Total Energy

               Call Toterg(Ent(BoxID), Virt(BoxID), BoxID)
               If (Abs(Ent(BoxID)-En(BoxID)).Gt.1.D-6) Then
                  Write (6, *) 
     &                    ' ######### Problems Energy ################ '
               End If
               If (Abs(Virt(BoxID)-Vir(BoxID)).Gt.1.D-6) Then
                  Write (6, *) 
     &                    ' ######### Problems Virial ################ '
               End If
               Write (6, 99002) BoxID, Ent(BoxID), En(BoxID), Ent(BoxID) - En(BoxID), 
     &                          Virt(BoxID), Vir(BoxID), Virt(BoxID) - Vir(BoxID)
            End Do
C           ---Calculation Chemical Potential
            Call Init_Chem(2)
         End If
      End Do
      Call Store(21, Dr, Vmax)

 
99001 Format (' Box : ', I3, /, 
     &        ' Total Energy Initial Configuration : ', F12.5, /, 
     &        ' Total Virial Initial Configuration : ', F12.5)
99002 Format (' Box : ', I3, /, ' Total Energy End Of Simulation   : ', 
     &        F12.5, /, '       Running Energy             : ', F12.5, 
     &        /, '       Difference                 :  ', E12.5, /, 
     &        ' Total Virial End Of Simulation   : ', F12.5, /, 
     &        '       Running Virial             : ', F12.5, /, 
     &        '       Difference                 :  ', E12.5)
99003 Format (' Number Of Att. To Displ. A Part.  : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
99004 Format (' Number Of Att. To Chan. Volume    : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
99005 Format (' Number Of Att. To Exchange Part.  : ', I10, /, 
     &        ' Success: ', I10, '(= ', F5.2, '%)')
      
      Stop
      End
