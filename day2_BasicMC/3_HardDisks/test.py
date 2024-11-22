import HardDisks
import numpy as np
import matplotlib.pyplot as plt

hd = HardDisks.HardDisks(
    numberOfInitCycles=int(1e4),
    numberOfProdCycles=int(1e4),
    numberOfParticles=16,
    maxDisplacement=1.0,
    sampleFrequency=1000,
    boxSize=96.0,
    rdfBins=100,
    method=HardDisks.Method.Static
)
hd.run()
rdf = hd.getRDF()
