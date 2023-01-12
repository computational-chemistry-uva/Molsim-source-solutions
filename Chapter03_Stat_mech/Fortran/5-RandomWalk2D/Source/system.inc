      Integer MaxLattice,MaxParticles,NumberOfParticles

      Parameter (MaxLattice = 65 )
      Parameter (MaxParticles = 2000 )

      Integer Lattice(MaxLattice,MaxLattice),ParticlePosX(MaxParticles),
     &        ParticlePosY(MaxParticles),Mxx(MaxParticles),Myy(MaxParticles)

      Common /Sstrr/ Lattice,ParticlePosX,ParticlePosY,NumberOfParticles,
     &               Mxx,Myy
