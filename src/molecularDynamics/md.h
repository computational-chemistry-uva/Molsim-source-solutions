#pragma once

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "double3.h"
#include "sample.h"
#include "thermostats.h"
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

  double cutOff = 3.0;
  double cutOffEnergy;
  double volume;
  double density;
  size_t degreesOfFreedom;
  double3 totalMomentum;
  double gridSize;
  bool equilibration = false;

  bool useNoseHoover;
  VelocityScaling velocityScaling;
  NoseHooverNVT noseHoover;

  std::vector<double3> positions;
  std::vector<double3> momenta;
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
                    size_t sampleFrequency = 100, size_t logLevel = 0, size_t seed = 12, bool useNoseHoover = false);

  double uniform() { return uniform_dist(mt); }
  double normal() { return normal_dist(mt); }

  /**
   * @brief Initializes particle positions on a cubic lattice with slight random perturbations.
   *
   * Distributes particles evenly on a grid within the simulation box, adding small random displacements to avoid
   * overlap.
   */
  void latticeInitialization();

  /**
   * @brief Minimizes potential energy using steepest descent method.
   *
   * Iteratively adjusts particle positions in the direction of the force to reduce potential energy.
   * Accepts or rejects steps based on energy reduction and adjusts step size accordingly.
   */
  void gradientDescent();

  /**
   * @brief Calculates forces between particles and updates potential energy and pressure.
   *
   * Implements the Lennard-Jones potential with a cutoff radius.
   * Uses pairwise interactions and accumulates forces, potential energy, and pressure.
   */
  void calculateForce();

  void thermostat();

  /**
   * @brief Integrates particle positions and momenta using the Verlet algorithm.
   *
   * Updates positions and momenta, applies velocity scaling during equilibration,
   * and recalculates kinetic energy and total momentum.
   */
  void integrate();

  /**
   * @brief Runs the molecular dynamics simulation for a specified number of steps.
   *
   * Performs force calculations and integrations in each step, handles equilibration,
   * logs information, and samples properties.
   *
   * @param numberOfSteps The number of integration steps to perform.
   * @param equilibrate Whether the simulation is in the equilibration phase.
   */
  void run(size_t& numberOfCycles, bool equilibrate, bool outputPDB = true);

  /**
   * @brief Returns a string representation of the current state of the simulation.
   *
   * Provides detailed information about simulation parameters and current measurements.
   *
   * @return A formatted string containing simulation data.
   */
  std::string repr();
};
