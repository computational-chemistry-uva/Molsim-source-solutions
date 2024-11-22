from molsim import MonteCarlo

mc = MonteCarlo(
    numberOfParticles=100,
    temperature=1.0,
    boxSize=12.0,
    maxDisplacement=0.5,
    numberOfInitCycles=10000,
    numberOfProdCycles=10000,
    sampleFrequency=100,
    sigma=1.0,
    epsilon=1.0,
    logLevel=0,
    seed=12,
)

mc.run()
