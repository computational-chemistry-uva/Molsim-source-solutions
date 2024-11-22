#pragma once

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "double3.h"
#include "energyVirial.h"
#include "utils.h"

/**
 * \brief Monte Carlo simulation class for particle systems.
 *
 * The MonteCarlo class encapsulates the Monte Carlo simulation methods and properties
 * for simulating particle interactions in a defined system. It includes system settings,
 * particle positions, energy calculations, and simulation control functions.
 */
struct MonteCarlo
{
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

  std::mt19937 mt;
  std::uniform_real_distribution<double> uniform_dist;

  EnergyVirial runningEnergyVirial;
  EnergyVirial totalEnergyVirial;

  EnergyVirial drift;
  double pressure = 0.0;
  std::vector<double> pressures;
  std::vector<double> chemicalPotentials;
  std::vector<double> heatCapacities;

  Logger logger;
  size_t frameNumber = 1;

  /**
   * \brief Constructs a MonteCarlo simulation object.
   *
   * Initializes the Monte Carlo simulation with the specified parameters.
   *
   * \param numberOfParticles Number of particles in the system.
   * \param temperature Temperature of the system.
   * \param boxSize Size of the simulation box.
   * \param maxDisplacement Maximum displacement in a move.
   * \param numberOfInitCycles Number of initialization cycles.
   * \param numberOfProdCycles Number of production cycles.
   * \param sampleFrequency Frequency of sampling the system.
   * \param sigma Lennard-Jones sigma parameter.
   * \param epsilon Lennard-Jones epsilon parameter.
   * \param logLevel Logging level.
   * \param seed Seed for random number generator.
   */
  MonteCarlo(size_t numberOfParticles, double temperature, double boxSize, double maxDisplacement,
             size_t numberOfInitCycles, size_t numberOfProdCycles, size_t sampleFrequency = 100, double sigma = 1.0,
             double epsilon = 1.0, size_t logLevel = 0, size_t seed = 12);

  /**
   * \brief Generates a uniform random number between 0 and 1.
   *
   * \return A random double between 0 and 1.
   */
  double uniform() { return uniform_dist(mt); }

  /**
   * \brief Performs a translation move on a particle.
   *
   * Attempts to move a particle to a new position and accepts or rejects based on energy change.
   *
   * \param particleIdx Index of the particle to move.
   */
  void translationMove(size_t particleIdx);

  /**
   * \brief Optimizes the maximum displacement based on acceptance rate.
   *
   * Adjusts the maximum displacement to achieve an optimal acceptance rate.
   */
  void optimizeMaxDisplacement();

  /**
   * \brief Computes the pressure of the system.
   *
   * Calculates the current pressure and records it.
   */
  void computePressure();

  /**
   * \brief Computes the chemical potential of the system.
   *
   * Estimates the chemical potential and records it.
   */
  void computeChemicalPotential();

  /**
   * \brief Computes the heat capacity of the system.
   *
   * Calculates the heat capacity (functionality to be implemented).
   */
  void computeHeatCapacity();

  /**
   * \brief Calculates the energy and virial for a particle.
   *
   * Computes the energy and virial contributions of a particle at a given position.
   *
   * \param position Position of the particle.
   * \param particleIdx Index of the particle.
   * \param startIndex Start loop over other particles from this index, allowing lower triangular computation
   * \return EnergyVirial object containing energy and virial.
   */
  EnergyVirial particleEnergyVirial(double3 position, size_t particleIdx);

  /**
   * \brief Calculates the total energy and virial of the system.
   *
   * Computes the total energy and virial contributions of all particles.
   *
   * \return EnergyVirial object containing total energy and virial.
   */
  EnergyVirial systemEnergyVirial();

  /**
   * \brief Runs the Monte Carlo simulation.
   *
   * Executes the simulation for the specified number of cycles.
   */
  void run();

  /**
   * \brief Returns a string representation of the Monte Carlo simulation.
   *
   * Provides detailed information about the current state of the simulation.
   *
   * \return A string containing the representation.
   */
  std::string repr();
};