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
  int numberOfParticles;
  double temperature;
  double boxSize;
  double sigma;
  double epsilon;
  int numberOfInitCycles;
  int numberOfProdCycles;
  int sampleFrequency;
  double maxDisplacement;
  double translationProbability;
  double pressure;
  double volumeProbability;
  double maxVolumeChange;

  double cutOff;
  EnergyVirial cutOffPrefactor;
  double volume;
  double density;
  double beta;

  int cycle;
  double translationsAttempted{0};
  double translationsAccepted{0};
  double volumeAttempted{0};
  double volumeAccepted{0};

  std::vector<double3> positions;

  std::mt19937 mt;
  std::uniform_real_distribution<double> uniform_dist;

  EnergyVirial runningEnergyVirial;
  EnergyVirial totalEnergyVirial;

  EnergyVirial drift;
  std::vector<double> pressures;
  std::vector<double> chemicalPotentials;
  std::vector<double> energies;

  Logger logger;
  int frameNumber = 1;

  /**
   * \brief Constructs a MonteCarlo simulation object.
   *
   * Initializes the Monte Carlo simulation with the specified parameters.
   *
   * \param numberOfParticles Number of particles in the system.
   * \param numberOfInitCycles Number of initialization cycles.
   * \param numberOfProdCycles Number of production cycles.
   * \param temperature Temperature of the system.
   * \param boxSize Size of the simulation box.
   * \param maxDisplacement Maximum displacement in a move.
   * \param translationProbability
   * \param pressure Pressure to couple to for NPT.
   * \param volumeProbability Probability of performing a volume move between (0, 1)
   * \param maxVolumeChange Maximum volume change of a volume move
   * \param swapProbability
   * \param sampleFrequency Frequency of sampling the system.
   * \param sigma Lennard-Jones sigma parameter.
   * \param epsilon Lennard-Jones epsilon parameter.
   * \param logLevel Logging level.
   * \param seed Seed for random number generator.
   */
  MonteCarlo(int numberOfParticles, int numberOfInitCycles, int numberOfProdCycles, double temperature, double boxSize,
             double maxDisplacement, double translationProbability = 1.0, double pressure = 0.0,
             double volumeProbability = 0.0, double maxVolumeChange = 1.0, int sampleFrequency = 100,
             double sigma = 1.0, double epsilon = 1.0, int logLevel = 0, int seed = 12);

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
   */
  void translationMove();

  /**
   * \brief Optimizes the maximum displacement based on acceptance rate.
   *
   * Adjusts the maximum displacement to achieve an optimal acceptance rate.
   */
  void optimizeMaxDisplacement();

  /**
   * \brief Performs a volume move on the box.
   *
   * Attempts to change the boxSize and accepts or rejects based on energy change.
   */
  void volumeMove();

  void optimizeVolumeChange();

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