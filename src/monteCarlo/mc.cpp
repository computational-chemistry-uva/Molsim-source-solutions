#include "mc.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "writePDB.h"

MonteCarlo::MonteCarlo(size_t numberOfParticles, size_t numberOfInitCycles, size_t numberOfProdCycles,
                       double temperature, double boxSize, double maxDisplacement, double translationProbability,
                       double pressure, double volumeProbability, double maxVolumeChange, size_t sampleFrequency,
                       double sigma, double epsilon, size_t logLevel, size_t seed)
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
  size_t numGrids = static_cast<size_t>(std::round(std::pow(numberOfParticles, 1.0 / 3.0) + 0.5));
  double gridSize = boxSize / (static_cast<double>(numGrids) + 2.0);
  logger.debug("numGrids " + std::to_string(numGrids) + " gridSize " + std::to_string(gridSize));

  // Initialize particle positions on a grid with slight random perturbations
  size_t counter = 0;
  for (size_t i = 0; i < numGrids; ++i)
  {
    for (size_t j = 0; j < numGrids; ++j)
    {
      for (size_t k = 0; k < numGrids; ++k)
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
  for (size_t k = 0; k < 10; k++)
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
    for (size_t i = 0; i < numberOfParticles; ++i)
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
  s += std::format("Number of particles  : {}\n", numberOfParticles);
  s += std::format("Temperature          : {:6.3f}\n", temperature);
  if (volumeProbability != 0.0)
  {
    s += std::format("Pressure             : {:6.3f}\n", pressure);
  }
  s += std::format("Box length           : {:6.3f}\n", boxSize);
  s += std::format("Volume               : {:6.3f}\n", volume);
  s += std::format("Density              : {:6.3f}\n", density);
  s += std::format("CutOff radius        : {:6.3f}\n", cutOff);
  s += std::format("CutOff energy        : {:6.3f}\n", cutOffPrefactor.energy * density);
  s += std::format("Steps run            : {}\n", cycle);
  s += std::format("Max displacement     : {:6.3f}\n", maxDisplacement);
  s += std::format("Translation acc.     : {:6.3f}\n", translationsAccepted / translationsAttempted);
  if (volumeProbability != 0.0)
  {
    s += std::format("Max volume change    : {:6.3f}\n", maxVolumeChange);
    s += std::format("Volume acc.          : {:6.3f}\n", volumeAccepted / volumeAttempted);
  }
  s += std::format("Running energy       : {:6.3f}\n", runningEnergyVirial.energy);
  s += std::format("Total energy         : {:6.3f}\n", totalEnergyVirial.energy);
  s += std::format("Drift energy         : {:6.3e}\n", drift.energy);
  s += std::format("Running virial       : {:6.3f}\n", runningEnergyVirial.virial);
  s += std::format("Total virial         : {:6.3f}\n", totalEnergyVirial.virial);
  s += std::format("Drift virial         : {:6.3e}\n", drift.virial);
  s += std::format("Average pressure     : {:6.3f}\n", average(pressures));
  s += std::format("Average chem. pot.   : {:6.3f}\n", realChemPot);
  s += std::format("Average Heat Cap.    : {:6.3f}\n",
                   variance(energies) / (temperature * temperature * numberOfParticles));
  s += "\n";
  return s;
}
