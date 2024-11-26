from molsim import MonteCarlo

mc = MonteCarlo(
    numberOfParticles=500,
    temperature=1.0,
    boxSize=8.0,
    maxDisplacement=0.5,
    numberOfInitCycles=500000,
    numberOfProdCycles=500000,
    sampleFrequency=1000,
    sigma=1.0,
    epsilon=1.0,
    logLevel=0,
    seed=12,
    pressure=1.0,
    volumeProbability=0.001,
)

mc.run()
