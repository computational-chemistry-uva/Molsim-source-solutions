      Subroutine Readdat(NumberOfEquilibrationCycles,NumberOfProductionCycles,
     &     SamplingFrequency,NumberOfDisplacementsPerCycle,MaxDisplacement)

      Implicit None

C     ---Read Input Data And Model Parameters
 
      Include 'parameter.inc'
      Include 'system.inc'
      Include 'potential.inc'
      Include 'conf.inc'

      Integer InitType,NumberOfEquilibrationCycles,NumberOfProductionCycles,
     &     I,NumberOfDisplacementsPerCycle,SamplingFrequency,Sstmm
      Double Precision Eps, Sig, Boxf, Rhof, Rho, MaxDisplacement, M1
 
C     ---Read Simulation Data
C
C     InitType = 0 : Initialize From A Lattice
C                1 : Read Configuration From Disk
C     Rho      = Density

      Open (15,FILE='input',FORM='FORMATTED')

      Read (15, *)
      Read (15, *) InitType,NumberOfEquilibrationCycles,
     &     NumberOfProductionCycles,SamplingFrequency
      Read (15, *)
      Read (15, *) MaxDisplacement
      Read (15, *)
      Read (15, *) NumberOfDisplacementsPerCycle
      Read (15, *)
      Read (15, *) NumberOfParticles,Temperature,Rho
      Read (15, *)
      Read (15, *) Eps, Sig, Mass, Rc

C     ---Initialize Random Number Generator
      
      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)
 
      If (NumberOfParticles.Gt.MaxParticles) Then
         Write (6, *) ' Error: Number Of Particles Too Large'
         Stop
      End If

      Box  = (Dble(NumberOfParticles)/Rho)**(1.D0/3.D0)
      HalfBox = Box/2

C     ---Read/Generate Configuration

      If (InitType.Eq.0) Then

C        ---Generate Configuration Form Lattice

         Call Lattice
      Else
         Write (6, *) ' Read Conf From Disk '
         Open (11,FILE='restart.dat',FORM='FORMATTED')
         Read (11, *) Boxf
         Read (11, *) NumberOfParticles
         Read (11, *) MaxDisplacement
         Rhof = Dble(NumberOfParticles)/Boxf**3
         If (Abs(Boxf-Box).Gt.1d-6) Then
            Write (6, 99007) Rho, Rhof
         End If
         Do I = 1, NumberOfParticles
            Read (11, *) X(I), Y(I), Z(I)
            X(I) = X(I)*Box/Boxf
            Y(I) = Y(I)*Box/Boxf
            Z(I) = Z(I)*Box/Boxf
         End Do
         Close(11)
      End If

C     ---Write Input Data

      Write (6, 99001) NumberOfEquilibrationCycles, NumberOfProductionCycles, SamplingFrequency
      Write (6, 99002) NumberOfDisplacementsPerCycle, MaxDisplacement
      Write (6, 99003) NumberOfParticles, Temperature, Rho, Box
      Write (6, 99004) Eps, Sig, Mass

C     ---Calculate Parameters:

      Beta = 1/Temperature
      
C     ---Calculate Cut-Off Radius Potential

      Rc    = Min(Rc, HalfBox)
      Rc2   = Rc*Rc
      Eps4  = 4.0d0*Eps
      Eps48 = 48.0d0*Eps
      Sig2  = Sig*Sig
      
      Return

99001 Format ('  Number Of Equilibration Cycles             :', I10, /, 
     &        '  Number Of Production Cycles                :', I10, /, 
     &        '  Sample Frequency                           :', I10, /)
99002 Format ('  Number Of Att. To Displ. A Part. Per Cycle :', I10, /, 
     &        '  Maximum Displacement                       :', F10.3, 
     &        //)
99003 Format ('  Number Of Particles                        :', I10, /, 
     &        '  Temperature                                :', F10.3, 
     &        /, '  Density                                    :', 
     &        F10.3, /, '  Box Length                                 :'
     &        , F10.3, /)
99004 Format ('  Model Parameters: ', /, '     Epsilon: ', F5.3, /, 
     &        '     Sigma  : ', F5.3, /, '     Mass   : ', F5.3)
99007 Format (' Requested Density: ', F5.2, 
     &        ' Different From Density On Disk: ', F5.2, /, 
     &        ' Rescaling Of Coordinates!')
      End
