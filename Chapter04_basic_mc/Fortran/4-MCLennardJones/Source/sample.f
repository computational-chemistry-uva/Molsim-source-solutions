      Subroutine Sample(I, TotalEnergy, TotalVirial, Pressure)
      Implicit None
C
C      Write Quantities (Pressure And Energy) To File
C

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'
      Include 'potential.inc'


      Integer I
      Double Precision TotalEnergy,Energy,Eps,Sig,TotalVirial,Pressure,Volume
 
      If (NumberOfParticles.Ne.0) Then
        Energy   = TotalEnergy/Dble(NumberOfParticles)
        Volume   = Box**3
        Pressure = (Dble(NumberOfParticles)/Volume)/Beta + 
     &       TotalVirial/(3.0d0*Volume)
        Pressure = Pressure + (16.0d0/3.0d0)*4.0d0*Atan(1.0)*
     &             (Dble(NumberOfParticles)/Volume)**2*Eps*
     &             ((1.0d0/3.0d0)*(Sig/Rc)**9-(Sig/Rc)**3)
      Else
         Energy   = 0.0d0
         Pressure = 0.0d0
      End If

      Write (66, *) I, Energy
      Write (67, *) I, Pressure

      Return
      End
