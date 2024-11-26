from molsim import MonteCarlo

mc = MonteCarlo(
    numberOfParticles=200,
    temperature=1.0,
    boxSize=6.0,
    maxDisplacement=0.5,
    numberOfInitCycles=10000,
    numberOfProdCycles=100000,
    sampleFrequency=100,
    sigma=1.0,
    epsilon=1.0,
    logLevel=0,
    seed=12,
    pressure=0.1,
    volumeProbability=0.0,
    swapProbability=0.0,
)

mc.run()
