#pragma once

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "double3.h"
#include "sample.h"
#include "utils.h"

struct MolecularDynamics
{
  // system settings
  size_t numberOfParticles;
  double temperature;
  double dt;
  double boxSize;
  size_t sampleFrequency;
  size_t totalSteps = 0;

  double cutOff;
  double cutOffEnergy;
  double volume;
  double density;
  size_t degreesOfFreedom;
  double3 totalMomentum;
  double gridSize;
  bool equilibration = false;

  std::vector<double3> positions;
  std::vector<double3> oldPositions;
  std::vector<double3> unwrappedPositions;
  std::vector<double3> velocities;
  std::vector<double3> forces;

  // random generator
  std::mt19937 mt;
  std::uniform_real_distribution<double> uniform_dist;
  std::normal_distribution<double> normal_dist;

  // observables
  double pressure = 0.0;
  double kineticEnergy = 0.0;
  double potentialEnergy = 0.0;
  double totalEnergy = 0.0;
  double observedTemperature = 0.0;

  double baselineEnergy = 0.0;
  double driftEnergy = 0.0;

  SampleThermodynamicalAverages thermoSampler;
  SampleRDF rdfSampler;
  SampleMSD msdSampler;

  // initialize logger
  Logger logger;
  size_t frameNumber = 1;

  MolecularDynamics(size_t numberOfParticles, double temperature, double dt, double boxSize,
                    size_t sampleFrequency = 100, size_t logLevel = 0, size_t seed = 12);

  double uniform() { return uniform_dist(mt); }
  double normal() { return normal_dist(mt); }

  void init();
  void latticeInitialization();
  void gradientDescent();
  void calculateForce();
  void integrate();
  void run(size_t& numberOfCycles, bool equilibrate, bool outputPDB = true);
  std::string repr();
};
