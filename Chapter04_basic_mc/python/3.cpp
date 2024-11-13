#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

HardDisks::HardDisks(size_t numberOfInitCycles; size_t numberOfProdCycles; size_t numberOfParticles;
                     double maxDisplacement; size_t sampleFrequency; double boxSize, Method method)
    : numberOfInitCycles(numberOfInitCycles),
      numberOfProdCycles(numberOfProdCycles),
      numberOfParticles(numberOfParticles),
      maxDisplacement(maxDisplacement),
      sampleFrequency(sampleFrequency),
      boxSize(boxSize),
      method(method)
{
  initialize();
}

void HardDisks::initialize()
{
  size_t latticeSites = static_cast<size_t> std::ceil(std::sqrt(static_cast<double>(numberOfParticles)));

  std::vector<std::pair<double, double>> latticePositions(latticeSites * latticeSites);
  for (size_t ix = 0; ix < latticeSites; ix++)
  {
    for (size_t iy = 0; iy < latticeSites; iy++)
    {
      double x = (ix / latticeSites) * boxSize;
      double y = (iy / latticeSites) * boxSize;
      latticePositions[ix + latticeSites * iy] = std::make_pair(x, y);
    }
  }
  std::sample(latticePositions.begin(), latticePositions.end(), std::back_inserter(positions), positions.size(), gen);
}

void run()
{
  for (size_t cycle = 0; cycle < numberOfInitCycles; cycle++)
  {
    numberOfAttemptedMoves++;
  }
}