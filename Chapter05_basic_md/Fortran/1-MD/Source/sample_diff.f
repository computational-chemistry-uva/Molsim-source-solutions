      Subroutine Sample_Diff(Switch)
      Implicit None

      Include 'system.inc'

C     MAXT  = Maximum Timesteps For The Correlation Time
C     MAXT0 = Maximum Number Of Time Origins
C     FREQT0 = frequency with which a new time origin is selected

C     time = Current time in units of Deltat
C     t0time[TMAXT0] = Time of the stored time origin in units of Deltat
C     t0Counter = Number of time origins stored so far
C     t0index = The index of the current t0 in the array of stored time origins.
C     CorrelTime = Time difference between the current time and the time origin
C     SampleCounter[MAXT] = Number of samples taken at a given CorrelTime.
 

      Integer   MAXT,MAXT0,FREQT0

      Parameter (MAXT=7500)
      Parameter (MAXT0=250)
      Parameter (FREQT0=50)

      Integer t0index,t0time(MAXT0),t0Counter,Switch,I,J,time,CorrelTime

      Double Precision Vxt0(MaxParticles,MAXT0),Vyt0(MaxParticles,MAXT0),
     &                 Vzt0(MaxParticles,MAXT0),Vacf(MAXT),R2(MAXT),
     &                 Rx0(MaxParticles,MAXT0),Ry0(MaxParticles,MAXT0),
     &                 Rz0(MaxParticles,MAXT0),SampleCounter(MAXT),CumIntegration,CorrelTimeime
       
      Save Vacf,SampleCounter,Vxt0,Vyt0,Vzt0,time,R2,Rx0,
     &     Ry0,Rz0,t0time,t0Counter
       
      If (Switch.Eq.1) Then

C     Initialize Everything

         time = 0
         t0Counter    = 0

         Do I = 1, MAXT
            R2(I)   = 0.0d0
            SampleCounter(I) = 0.0d0
            Vacf(I)  = 0.0d0
         Enddo

      Elseif(Switch.Eq.2) Then

         time = time + 1

         If (Mod(time,FREQT0).Eq.0) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     New Time Origin                                             C
C     Store The Positions/Velocities; The Current Velocities      C
C     Are Vxx(I) And The Current Positions Are Mxx(I).            C
C     Question: Why Do You Have To Be Careful With Pbc ?          C
C                                                                 C
C     Make Sure To Study Algorithm 8 (Page 82) Of Frenkel/Smit    C
C     Before You Start To Make Modifications                      C
C                                                                 C
C     Note That Most Variable Names Are Different Here Then In    C
C     Frenkel/Smit. In This Way, You Will Have To Think More...   C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification
C     New Time Origin

            t0Counter         = t0Counter + 1
            t0index       = Mod(t0Counter-1, MAXT0) + 1
            t0time(t0index)  = time

C     Store Particle Positions/Velocities
C     Question: Why Do We Use Mxx Instead Of Xxx/Rxx ?????

            Do I = 1, NumberOfParticles
               Rx0(I,t0index)  = Mxx(I)
               Ry0(I,t0index)  = Myy(I)
               Rz0(I,t0index)  = Mzz(I)
               Vxt0(I,t0index) = Vxx(I)
               Vyt0(I,t0index) = Vyy(I)
               Vzt0(I,t0index) = Vzz(I)
            Enddo
         Endif

C     End   Modification

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     Loop Over All Time Origins That Have Been Stored    C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification
         Do J = 1, Min(t0Counter, MAXT0)
            CorrelTime = time - t0time(J) + 1

C     Only If The Time Difference Is Shorter Than The Maximum Correlation Time

            If (CorrelTime.Lt.MAXT) Then
               SampleCounter(CorrelTime) = SampleCounter(CorrelTime) + 1.0d0
                  
               Do I = 1, NumberOfParticles
                  Vacf(CorrelTime) = Vacf(CorrelTime) + Vxx(I)*Vxt0(I,J) + 
     &                                  Vyy(I)*Vyt0(I,J) + 
     &                                  Vzz(I)*Vzt0(I,J)

                  R2(CorrelTime)  = R2(CorrelTime)  + (Mxx(I)-Rx0(I,J))**2 + 
     &                                  (Myy(I)-Ry0(I,J))**2 + 
     &                                  (Mzz(I)-Rz0(I,J))**2
               Enddo
            Endif
         Enddo

C     End   Modification

      Else

Ccccccccccccccccccccccccccccccccccc
C     Write Everything To Disk    C
Ccccccccccccccccccccccccccccccccccc

         CumIntegration = 0.0d0

         Do I = 1, MAXT-1
                           
            CorrelTimeime = Dble(I-1)*Deltat

            If (SampleCounter(I).Ge.0.5d0) Then
               Vacf(I) = Vacf(I)/(Dble(NumberOfParticles)*SampleCounter(I))
               R2(I)  = R2(I) /(Dble(NumberOfParticles)*SampleCounter(I))
            Else
               Vacf(I) = 0.0d0
               R2(I)  = 0.0d0
            Endif

            CumIntegration = CumIntegration + Deltat*Vacf(I)/3.0d0

            Write(24,*) CorrelTimeime,Vacf(I),CumIntegration
            
            If(I.Ne.1) Then
               Write(25,*) CorrelTimeime,R2(I),R2(I)/(6.0d0*CorrelTimeime)
            Else
               Write(25,*) CorrelTimeime,R2(I),' 0.0d0'
            Endif
         Enddo
      Endif

      Return
      End
