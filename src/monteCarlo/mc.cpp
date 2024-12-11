#include "mc.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "writePDB.h"

MonteCarlo::MonteCarlo(int numberOfParticles, int numberOfInitCycles, int numberOfProdCycles, double temperature,
                       double boxSize, double maxDisplacement, double translationProbability, double pressure,
                       double volumeProbability, double maxVolumeChange, int sampleFrequency, double sigma,
                       double epsilon, int logLevel, int seed)
    : numberOfParticles(numberOfParticles),
      numberOfInitCycles(numberOfInitCycles),
      numberOfProdCycles(numberOfProdCycles),
      temperature(temperature),
      boxSize(boxSize),
      maxDisplacement(maxDisplacement),
      translationProbability(translationProbability),
      pressure(pressure),
      volumeProbability(volumeProbability),
      maxVolumeChange(maxVolumeChange),
      sampleFrequency(sampleFrequency),
      sigma(sigma),
      epsilon(epsilon),
      positions(numberOfParticles),
      mt(seed),
      uniform_dist(0.0, 1.0),
      logger(logLevel)

{
  // Calculate cutoff radius, cutoff energy, and system properties
  cutOff = 3.0 * sigma;

  volume = boxSize * boxSize * boxSize;
  density = numberOfParticles / volume;
  beta = 1 / temperature;

  double sigmaOverCutoff = sigma / cutOff;
  double r3i = sigmaOverCutoff * sigmaOverCutoff * sigmaOverCutoff;
  cutOffPrefactor = EnergyVirial(epsilon * (8.0 / 3.0) * M_PI * ((1.0 / 3.0) * r3i * r3i * r3i - r3i),
                                 epsilon * (16.0 / 3.0) * M_PI * ((2.0 / 3.0) * r3i * r3i * r3i - r3i));

  // Empty movie.pdb file
  std::ofstream file("movie.pdb");
  file.close();

  // Initialize MD simulation
  logger.info("Class MC created.");
  logger.debug(repr());

  // Calculate number of grid points and grid spacing based on number of particles
  int numGrids = static_cast<int>(std::round(std::pow(numberOfParticles, 1.0 / 3.0) + 0.5));
  double gridSize = boxSize / (static_cast<double>(numGrids) + 2.0);
  logger.debug("numGrids " + std::to_string(numGrids) + " gridSize " + std::to_string(gridSize));

  // Initialize particle positions on a grid with slight random perturbations
  int counter = 0;
  for (int i = 0; i < numGrids; ++i)
  {
    for (int j = 0; j < numGrids; ++j)
    {
      for (int k = 0; k < numGrids; ++k)
      {
        if (counter < numberOfParticles)
        {
          positions[counter] =
              double3((i + 0.01 * (uniform() - 0.5)) * gridSize, (j + 0.01 * (uniform() - 0.5)) * gridSize,
                      (k + 0.01 * (uniform() - 0.5)) * gridSize);
          ++counter;
        }
      }
    }
  }

  double totalProbability = translationProbability + volumeProbability;
  volumeProbability = volumeProbability / totalProbability;

  writePDB("movie.pdb", positions, boxSize, frameNumber);
  ++frameNumber;

  totalEnergyVirial = systemEnergyVirial(positions, boxSize, cutOff, sigma, epsilon, cutOffPrefactor);
  runningEnergyVirial = totalEnergyVirial;

  logger.info("(Init) completed. Total energy: " + std::to_string(totalEnergyVirial.energy) +
              ", Total virial: " + std::to_string(totalEnergyVirial.virial));
  logger.info(repr());
};

void MonteCarlo::computePressure()
{
  // Compute pressure using virial theorem and record it
  double virialPressure = (density / beta) + (runningEnergyVirial.virial / (3.0 * volume));
  virialPressure += (cutOffPrefactor.virial * density * density);
  pressures.push_back(virialPressure);
}

void MonteCarlo::computeChemicalPotential()
{
  // Estimate chemical potential by inserting random particles
  EnergyVirial chemPotEV;
  for (int k = 0; k < 10; k++)
  {
    double3 randomPosition = double3(uniform(), uniform(), uniform());
    randomPosition *= boxSize;
    chemPotEV += particleEnergyVirial(positions, randomPosition, -1, boxSize, cutOff, sigma, epsilon);
  }
  double chemicalPotential = std::exp(-beta * chemPotEV.energy);
  chemicalPotentials.push_back(chemicalPotential);
}

void MonteCarlo::run()
{
  // Main simulation loop
  for (cycle = 0; cycle < numberOfInitCycles + numberOfProdCycles; ++cycle)
  {
    // Attempt translation moves
    for (int i = 0; i < numberOfParticles; ++i)
    {
      double rand = uniform();
      if (rand < volumeProbability)
      {
        volumeMove();
      }
      else
      {
        translationMove();
      }
    }

    // Every sampleFrequency cycles, compute and log properties
    if (cycle % sampleFrequency == 0)
    {
      totalEnergyVirial = systemEnergyVirial(positions, boxSize, cutOff, sigma, epsilon, cutOffPrefactor);
      logger.info(repr());
      writePDB("movie.pdb", positions, boxSize, frameNumber);
      frameNumber++;

      optimizeMaxDisplacement();
      optimizeVolumeChange();
      // If in production phase, compute observables
      if (cycle > numberOfInitCycles)
      {
        computePressure();
        computeChemicalPotential();
        energies.push_back(runningEnergyVirial.energy);
      }
    }
  }

  totalEnergyVirial = systemEnergyVirial(positions, boxSize, cutOff, sigma, epsilon, cutOffPrefactor);
  logger.info(repr());
  writePDB("movie.pdb", positions, boxSize, frameNumber);
}

std::string MonteCarlo::repr()
{
  std::string s;
  EnergyVirial drift = (runningEnergyVirial - totalEnergyVirial);

  double averageChemicalPotential = average(chemicalPotentials);
  double realChemPot =
      chemicalPotentials.size()
          ? -temperature * (std::log(averageChemicalPotential / density) + (cutOffPrefactor.energy * density))
          : 0.0;

  s += "Monte Carlo program\n";
  s += "----------------------------\n";
  s += "Number of particles  : " + std::to_string(numberOfParticles) + "\n";
  s += "Temperature          : " + std::to_string(temperature) + "\n";
  if (volumeProbability != 0.0)
  {
    s += "Pressure             : " + std::to_string(pressure) + "\n";
  }
  s += "Box length           : " + std::to_string(boxSize) + "\n";
  s += "Volume               : " + std::to_string(volume) + "\n";
  s += "Density              : " + std::to_string(density) + "\n";
  s += "CutOff radius        : " + std::to_string(cutOff) + "\n";
  s += "CutOff energy        : " + std::to_string(cutOffPrefactor.energy * density) + "\n";
  s += "Steps run            : " + std::to_string(cycle) + "\n";
  s += "Max displacement     : " + std::to_string(maxDisplacement) + "\n";
  s += "Translation acc.     : " + std::to_string(translationsAccepted / translationsAttempted) + "\n";
  if (volumeProbability != 0.0)
  {
    s += "Max volume change    : " + std::to_string(maxVolumeChange) + "\n";
    s += "Volume acc.          : " + std::to_string(volumeAccepted / volumeAttempted) + "\n";
  }
  s += "Running energy       : " + std::to_string(runningEnergyVirial.energy) + "\n";
  s += "Total energy         : " + std::to_string(totalEnergyVirial.energy) + "\n";
  s += "Drift energy         : " + std::to_string(drift.energy / totalEnergyVirial.energy) + "\n";
  s += "Running virial       : " + std::to_string(runningEnergyVirial.virial) + "\n";
  s += "Total virial         : " + std::to_string(totalEnergyVirial.virial) + "\n";
  s += "Drift virial         : " + std::to_string(drift.virial / totalEnergyVirial.virial) + "\n";
  s += "Average pressure     : " + std::to_string(average(pressures)) + "\n";
  s += "Average chem. pot.   : " + std::to_string(realChemPot) + "\n";
  s +=
      "Average Heat Cap.    : " + std::to_string(variance(energies) / (temperature * temperature * numberOfParticles)) +
      "\n";
  s += "\n";
  return s;
}
