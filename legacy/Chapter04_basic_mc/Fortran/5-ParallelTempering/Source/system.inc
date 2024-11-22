C     Rx/Ry : Particle Position
C     Uold  : Energy
C     Box   : Boxsize
C     Ntemp : Number Of Systems

      Integer NUMBEROFPARTICLES,MAXNUMBEROFSYSTEMS

      Parameter (NUMBEROFPARTICLES=9)
      Parameter (MAXNUMBEROFSYSTEMS=7)

      Double Precision Rx(NUMBEROFPARTICLES,MAXNUMBEROFSYSTEMS)
      Double Precision Ry(NUMBEROFPARTICLES,MAXNUMBEROFSYSTEMS)
      Double Precision Uold(MAXNUMBEROFSYSTEMS), Box

      Integer          Ntemp

      Common /Ptty/ Rx,Ry,Uold,Box,Ntemp
