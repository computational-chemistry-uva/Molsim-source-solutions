#include "mc.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "writePDB.h"

MonteCarlo::MonteCarlo(size_t numberOfParticles, double temperature, double boxSize, double maxDisplacement,
                       size_t numberOfInitCycles, size_t numberOfProdCycles, size_t sampleFrequency, double sigma,
                       double epsilon, size_t logLevel, size_t seed)
    : numberOfParticles(numberOfParticles),
      temperature(temperature),
      boxSize(boxSize),
      maxDisplacement(maxDisplacement),
      sampleFrequency(sampleFrequency),
      sigma(sigma),
      epsilon(epsilon),
      numberOfInitCycles(numberOfInitCycles),
      numberOfProdCycles(numberOfProdCycles),
      positions(numberOfParticles),
      mt(seed),
      uniform_dist(0.0, 1.0),
      logger(logLevel)

{
  // Calculate cutoff radius, cutoff energy, and system properties
  cutOff = 0.4999 * boxSize;
  cutOffSquared = cutOff * cutOff;

  volume = boxSize * boxSize * boxSize;
  density = numberOfParticles / volume;
  beta = 1 / temperature;

  double sigmaOverCutoff = sigma / cutOff;
  double r3i = sigmaOverCutoff * sigmaOverCutoff * sigmaOverCutoff;
  cutOffEnergy = epsilon * numberOfParticles * (8.0 / 3.0) * M_PI * density * ((1.0 / 3.0) * r3i * r3i * r3i - r3i);
  cutOffVirial = epsilon * (16.0 / 3.0) * M_PI * density * density * ((2.0 / 3.0) * r3i * r3i * r3i - r3i);

  // empty movie.pdb file
  std::ofstream file("movie.pdb");
  file.close();

  // Initialize MD simulation
  logger.info("Class MC created.");
  logger.debug(repr());

  // Calculate number of grid points and grid spacing based on number of
  // particles
  size_t numGrids = static_cast<size_t>(std::round(std::pow(numberOfParticles, 1.0 / 3.0) + 1.5));
  double gridSize = boxSize / (static_cast<double>(numGrids) + 2.0);
  logger.debug("numGrids " + std::to_string(numGrids) + " gridSize " + std::to_string(gridSize));

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

  writePDB("movie.pdb", positions, boxSize, frameNumber);
  ++frameNumber;

  totalEnergyVirial = systemEnergyVirial();
  runningEnergyVirial = totalEnergyVirial;

  logger.info("(Init) completed. Total energy: " + std::to_string(totalEnergyVirial.energy) +
              ", Total virial: " + std::to_string(totalEnergyVirial.virial));
  logger.info(repr());
};

void MonteCarlo::translationMove(size_t particleIdx)
{
  numberOfAttemptedMoves++;
  double3 displacement = maxDisplacement * (double3(uniform(), uniform(), uniform()) - 0.5);
  double3 newPosition = positions[particleIdx] + displacement;

  EnergyVirial oldEnergyVirial = particleEnergyVirial(positions[particleIdx], particleIdx);
  EnergyVirial newEnergyVirial = particleEnergyVirial(newPosition, particleIdx);

  if (uniform() < std::exp(-beta * (newEnergyVirial.energy - oldEnergyVirial.energy)))
  {
    numberOfAcceptedMoves++;
    positions[particleIdx] = newPosition;
    runningEnergyVirial += (newEnergyVirial - oldEnergyVirial);
  }
}

EnergyVirial MonteCarlo::particleEnergyVirial(double3 position, size_t particleIdx)
{
  EnergyVirial particleEnergyVirial;
  for (size_t i = 0; i < numberOfParticles; i++)
  {
    if (i != particleIdx)
    {
      double3 dr = position - positions[i];
      dr = wrap(dr, boxSize);
      double r2 = double3::dot(dr, dr);
      if (r2 < cutOffSquared)
      {
        double r2i = (sigma * sigma) / r2;
        double r6i = r2i * r2i * r2i;

        particleEnergyVirial.energy += 4.0 * epsilon * (r6i * r6i - r6i);
        particleEnergyVirial.virial += 48.0 * epsilon * (r6i * r6i - 0.5 * r6i);
      }
    }
  }
  return particleEnergyVirial;
}

EnergyVirial MonteCarlo::systemEnergyVirial()
{
  EnergyVirial systemEnergyVirial;
  for (size_t i = 0; i < numberOfParticles; i++)
  {
    systemEnergyVirial += particleEnergyVirial(positions[i], i);
  }
  systemEnergyVirial.energy /= 2.0;
  systemEnergyVirial.virial /= 2.0;
  systemEnergyVirial.energy += cutOffEnergy;
  return systemEnergyVirial;
}

void MonteCarlo::computePressure()
{
  pressure = (density / beta) + (runningEnergyVirial.virial / (3.0 * volume));
  pressure += cutOffVirial;
  pressures.push_back(pressure);
}

void MonteCarlo::computeChemicalPotential()
{
  EnergyVirial chemPotEV;
  for (size_t k = 0; k < 10; k++)
  {
    double3 randomPosition = double3(uniform(), uniform(), uniform());
    randomPosition *= boxSize;
    chemPotEV += particleEnergyVirial(randomPosition, -1);
  }
  double chemicalPotential = std::exp(-beta * chemPotEV.energy);
  chemicalPotentials.push_back(chemicalPotential);
}

void MonteCarlo::computeHeatCapacity() {}

void MonteCarlo::optimizeMaxDisplacement()
{
  double acceptance = numberOfAcceptedMoves / numberOfAttemptedMoves;
  double scaling = std::clamp(2.0 * acceptance, 0.5, 1.5);
  maxDisplacement = std::clamp(scaling * maxDisplacement, 0.0001, 0.25 * boxSize);
}

void MonteCarlo::run()
{
  size_t numberOfAttemptedMoves = 0;
  size_t numberOfAcceptedMoves = 0;
  for (cycle = 0; cycle < numberOfInitCycles + numberOfProdCycles; ++cycle)
  {
    for (size_t i = 0; i < numberOfParticles; ++i)
    {
      size_t particleIdx = static_cast<size_t>(uniform() * numberOfParticles);
      translationMove(particleIdx);
    }

    if (cycle % sampleFrequency == 0)
    {
      totalEnergyVirial = systemEnergyVirial();
      logger.info(repr());
      writePDB("movie.pdb", positions, boxSize, frameNumber);
      frameNumber++;

      optimizeMaxDisplacement();
      if (cycle > numberOfInitCycles)
      {
        computePressure();
        computeChemicalPotential();
        // computeHeatCapacity()
      }
    }
  }


  totalEnergyVirial = systemEnergyVirial();
  logger.info(repr());
  writePDB("movie.pdb", positions, boxSize, frameNumber);
}

std::string MonteCarlo::repr()
{
  std::string s;
  EnergyVirial drift = (runningEnergyVirial - totalEnergyVirial);

  double averageChemicalPotential = average(chemicalPotentials);
  double realChemPot = chemicalPotentials.size() ? -temperature * (std::log(averageChemicalPotential / density) + cutOffEnergy) : 0.0;

  s += "Monte Carlo program\n";
  s += "----------------------------\n";
  s += "Number of particles  : " + std::to_string(numberOfParticles) + "\n";
  s += "Temperature          : " + std::to_string(temperature) + "\n";
  s += "Box length           : " + std::to_string(boxSize) + "\n";
  s += "Density              : " + std::to_string(density) + "\n";
  s += "CutOff radius        : " + std::to_string(cutOff) + "\n";
  s += "CutOff energy        : " + std::to_string(cutOffEnergy) + "\n";
  s += "Steps run            : " + std::to_string(cycle) + "\n";
  s += "Max displacement     : " + std::to_string(maxDisplacement) + "\n";
  s += "Acceptance fraction  : " + std::to_string(numberOfAcceptedMoves / numberOfAttemptedMoves) + "\n";
  s += "Running energy       : " + std::to_string(runningEnergyVirial.energy) + "\n";
  s += "Total energy         : " + std::to_string(totalEnergyVirial.energy) + "\n";
  s += std::format("Drift energy         : {:6.3e}\n", drift.energy);
  s += "Running virial       : " + std::to_string(runningEnergyVirial.virial) + "\n";
  s += "Total virial         : " + std::to_string(totalEnergyVirial.virial) + "\n";
  s += "Drift virial         : " + std::to_string(drift.virial) + "\n";
  s += "Current Pressure     : " + std::to_string(pressure) + "\n";
  s += "Average pressure     : " + std::to_string(average(pressures)) + "\n";
  s += "Average chem. pot.   : " + std::to_string(realChemPot) + "\n";
  s += "Average Heat Cap.    : " + std::to_string(average(heatCapacities)) + "\n";
  s += "\n";
  return s;
}

PYBIND11_MODULE(mc, m)
{
  pybind11::class_<MonteCarlo>(m, "MonteCarlo")
      .def(pybind11::init<size_t, double, double, double, size_t, size_t, size_t, double, double, size_t, size_t>(),
           pybind11::arg("numberOfParticles"), pybind11::arg("temperature"), pybind11::arg("boxSize"),
           pybind11::arg("maxDisplacement"), pybind11::arg("numberOfInitCycles"), pybind11::arg("numberOfProdCycles"),
           pybind11::arg("sampleFrequency") = 100, pybind11::arg("sigma") = 1.0, pybind11::arg("epsilon") = 1.0,
           pybind11::arg("logLevel") = 0, pybind11::arg("seed") = 12)
      .def("__repr__", &MonteCarlo::repr)
      .def("run", &MonteCarlo::run);
}