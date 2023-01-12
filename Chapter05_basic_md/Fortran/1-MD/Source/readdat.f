      Subroutine Readdat
      Implicit None
 
      Include 'system.inc'
 
C     Read In System Information

      Open (21, File="input")
 
      Read (21,*) Box,NumberOfParticles,NumberOfSteps,Temperature,
     &     Deltat,NUmberOfInitSteps
 
      Close (21)

      If (NumberOfParticles.Gt.MaxParticles) Then
         Write (6,*) 'Maximum No. Of Particles Is : ',MaxParticles
         Stop
      Endif
 
C     Calculate Some Parameters
 
      HalfBox = 0.5d0*Box
 
      Rcutsq = (0.49999d0*Box)**2
      Ecut = 4.0d0*((Rcutsq**(-6.0d0)) - (Rcutsq**(-3.0d0)))
 
C     Print Information To The Screen
 
      Write (6,*) 'Molecular Dynamics Program'
      Write (6,*)
      Write (6,*) 'Number Of Particles   : ',NumberOfParticles
      Write (6,*) 'Boxlength             : ',Box
      Write (6,*) 'Density               : ',Dble(NumberOfParticles)/(Box*Box*Box)
      Write (6,*) 'Temperature           : ',Temperature
      Write (6,*) 'Cut-Off Radius        : ',Dsqrt(Rcutsq)
      Write (6,*) 'Cut-Off Energy        : ',Ecut
      Write (6,*) 'Number Of Steps       : ',NumberOfSteps
      Write (6,*) 'Number Of Init Steps  : ',NUmberOfInitSteps
      Write (6,*) 'Delta t               : ',Deltat  
 
      Return
      End
