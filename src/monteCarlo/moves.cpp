#include <cmath>
#include <vector>

#include "mc.h"

void MonteCarlo::translationMove()
{
  translationsAttempted++;
  // Generate random displacement
  size_t particleIdx = static_cast<size_t>(uniform() * numberOfParticles);
  double3 displacement = maxDisplacement * (double3(uniform(), uniform(), uniform()) - 0.5);
  double3 trialPosition = positions[particleIdx] + displacement;

  // Calculate old and new energy and virial
  EnergyVirial oldEnergyVirial =
      particleEnergyVirial(positions, positions[particleIdx], particleIdx, boxSize, cutOff, sigma, epsilon);
  EnergyVirial newEnergyVirial =
      particleEnergyVirial(positions, trialPosition, particleIdx, boxSize, cutOff, sigma, epsilon);

  // Accept or reject the move based on Metropolis criterion
  if (uniform() < std::exp(-beta * (newEnergyVirial.energy - oldEnergyVirial.energy)))
  {
    translationsAccepted++;
    positions[particleIdx] = trialPosition;
    runningEnergyVirial += (newEnergyVirial - oldEnergyVirial);
  }
}

void MonteCarlo::optimizeMaxDisplacement()
{
  // Adjust max displacement to maintain optimal acceptance ratio
  if (translationsAttempted > 100)
  {
    double acceptance = translationsAccepted / translationsAttempted;
    double scaling = std::clamp(2.0 * acceptance, 0.5, 1.5);
    maxDisplacement = std::clamp(scaling * maxDisplacement, 0.0001, 0.25 * boxSize);
    translationsAccepted = 0;
    translationsAttempted = 0;
  }
}

void MonteCarlo::volumeMove()
{
  volumeAttempted++;

  double volumeChange = (uniform() - 0.5) * maxVolumeChange;
  double newVolume = volume + volumeChange;
  double newBoxSize = std::cbrt(newVolume);
  double scale = newBoxSize / boxSize;

  std::vector<double3> trialPositions(positions);
  for (size_t i = 0; i < numberOfParticles; i++)
  {
    trialPositions[i] *= scale;
  }
  EnergyVirial newEnergyVirial =
      systemEnergyVirial(trialPositions, newBoxSize, cutOff, sigma, epsilon, cutOffPrefactor);

  if (uniform() < std::exp((numberOfParticles + 1.0) * std::log(newVolume / volume) -
                           beta * (pressure * volumeChange + (newEnergyVirial.energy - runningEnergyVirial.energy))))
  {
    volumeAccepted++;
    positions = trialPositions;
    boxSize = newBoxSize;
    volume = boxSize * boxSize * boxSize;
    density = numberOfParticles / volume;
    runningEnergyVirial = newEnergyVirial;
  }
}

void MonteCarlo::optimizeVolumeChange()
{
  if (volumeAttempted > 100)
  {
    double acceptance = volumeAccepted / volumeAttempted;
    double scaling = std::clamp(2.0 * acceptance, 0.5, 1.5);
    maxVolumeChange = std::clamp(scaling * maxVolumeChange, 0.001 * volume, 0.5 * volume);
    volumeAccepted = 0;
    volumeAttempted = 0;
  }
}
