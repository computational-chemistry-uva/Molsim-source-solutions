from mc import MonteCarlo

a = MonteCarlo(
    numberOfParticles=600,
    temperature=1.0,
    boxSize=12.0,
    maxDisplacement=3.0,
    numberOfInitCycles=100000,
    numberOfProdCycles=100000,
    sampleFrequency=100,
)
a.run()
