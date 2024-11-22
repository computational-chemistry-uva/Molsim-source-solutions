C     System.Inc
C
C     Information About The System; Various Constants
C
C     Maxchain = Maximum Chain Length
C     Maxtrial = Maximum Number Of Trial Positions
C
C     Lcbmc    = Do We Use Cbmc ?
C     
C     ChainLength    = Chain Length
C     NumberOfTrialPositions    = Number Of Trial Positions
C
C     X/Y/Zpos = Coordinates Of The Chain

      Integer          NumberOfSteps,NumberOfInitializationSteps,
     & Maxchain,
     & Maxtrial,NumberOfTrialPositions,  ChainLength
      Logical          Lstatic,Lcbmc

      Parameter(Maxchain = 100)
      Parameter(Maxtrial = 20)

      Double Precision Xpos(Maxchain),Ypos(Maxchain),Zpos(Maxchain),
     &                 Rcut,A,kt,Theta0,Beta

      Common /Syss/    Xpos,Ypos,Zpos,Rcut,A,kt,Theta0,Beta,
     &                 NumberOfSteps,NumberOfTrialPositions,
     & ChainLength,NumberOfInitializationSteps,Lstatic,Lcbmc  
