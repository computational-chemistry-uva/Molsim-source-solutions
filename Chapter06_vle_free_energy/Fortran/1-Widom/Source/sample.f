      Subroutine Sample(I, TotalEnergy, TotalVirial, Pressure)
      Implicit None
C
C      Write Quantities (Pressure And Energy) To File
C

      Include 'parameter.inc'
      Include 'conf.inc'
      Include 'system.inc'

      Integer I
      Double Precision TotalEnergy,Energy,TotalVirial,Pressure,Volume
 
      If (NumberOfParticles.Ne.0) Then
        Energy   = TotalEnergy/Dble(NumberOfParticles)
        Volume   = Box**3
        Pressure = (Dble(NumberOfParticles)/Volume)/Beta + 
     &       TotalVirial/(3.0d0*Volume)
      Else
         Energy   = 0.0d0
         Pressure = 0.0d0
      End If

      Write (66, *) I, Energy
      Write (67, *) I, Pressure

      Return
      End
