#pragma once

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "double3.h"
#include "utils.h"

struct EnergyVirial
{
  double energy;
  double virial;

  EnergyVirial(double energy = 0.0, double virial = 0.0) : energy(energy), virial(virial) {};

  EnergyVirial& operator+=(const EnergyVirial& other)
  {
    this->energy += other.energy;
    this->virial += other.virial;
    return *this;
  };
};

EnergyVirial operator-(const EnergyVirial& a, const EnergyVirial& b)
{
  return {a.energy - b.energy, a.virial - b.virial};
};

EnergyVirial operator/(const EnergyVirial& a, const EnergyVirial& b)
{
  return {a.energy / b.energy, a.virial / b.virial};
};

struct MonteCarlo
{
  // system settings
  size_t numberOfParticles;
  double temperature;
  double boxSize;
  double maxDisplacement;
  double sigma;
  double epsilon;
  size_t numberOfInitCycles;
  size_t numberOfProdCycles;
  size_t sampleFrequency;

  double cutOff;
  double cutOffSquared;
  double cutOffEnergy;
  double cutOffVirial;
  double chemPotConversion;
  double volume;
  double density;
  double beta;

  size_t cycle;
  double numberOfAttemptedMoves{0};
  double numberOfAcceptedMoves{0};

  std::vector<double3> positions;

  // random generator
  std::mt19937 mt;
  std::uniform_real_distribution<double> uniform_dist;

  EnergyVirial runningEnergyVirial;
  EnergyVirial totalEnergyVirial;

  // observables
  EnergyVirial drift;
  double pressure = 0.0;
  std::vector<double> pressures;
  std::vector<double> chemicalPotentials;
  std::vector<double> heatCapacities;

  // initialize logger
  Logger logger;
  size_t frameNumber = 1;

  MonteCarlo(size_t numberOfParticles, double temperature, double boxSize, double maxDisplacement,
             size_t numberOfInitCycles, size_t numberOfProdCycles, size_t sampleFrequency = 100, double sigma = 1.0,
             double epsilon = 1.0, size_t logLevel = 0, size_t seed = 12);

  double uniform() { return uniform_dist(mt); }

  void translationMove(size_t particleIdx);
  void optimizeMaxDisplacement();
  void computePressure();
  void computeChemicalPotential();
  void computeHeatCapacity();

  EnergyVirial particleEnergyVirial(double3 position, size_t particleIdx);
  EnergyVirial systemEnergyVirial();

  void run();
  std::string repr();
};