#include "md.h"

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "writePDB.h"

/**
 * @brief Constructs a MolecularDynamics simulation with the given parameters.
 *
 * Initializes positions, velocities, and forces for the specified number of particles.
 * Also sets up random number generators and samplers.
 *
 * @param numberOfParticles The number of particles in the simulation.
 * @param temperature The initial temperature of the system.
 * @param dt The time step size for integration.
 * @param boxSize The size of the simulation box (assuming cubic box).
 * @param logLevel The logging level for the logger.
 * @param seed The seed for the random number generator.
 */
MolecularDynamics::MolecularDynamics(size_t numberOfParticles, double temperature, double dt, double boxSize,
                                     size_t sampleFrequency, size_t logLevel, size_t seed)
    : numberOfParticles(numberOfParticles),
      temperature(temperature),
      dt(dt),
      boxSize(boxSize),
      sampleFrequency(sampleFrequency),
      positions(numberOfParticles),
      oldPositions(numberOfParticles),
      unwrappedPositions(numberOfParticles),
      velocities(numberOfParticles),
      forces(numberOfParticles),
      mt(seed),
      uniform_dist(0.0, 1.0),
      normal_dist(0.0, 1.0),
      rdfSampler(numberOfParticles, boxSize),
      msdSampler(numberOfParticles, boxSize, dt * sampleFrequency),
      logger(logLevel)

{
  // Calculate cutoff radius, cutoff energy, and system properties
  cutOff = 0.4999 * boxSize;
  double cutOffi6 = std::pow(cutOff, -6.0);
  cutOffEnergy = 4.0 * (cutOffi6 * cutOffi6 - cutOffi6);
  volume = boxSize * boxSize * boxSize;
  density = numberOfParticles / volume;
  degreesOfFreedom = 3 * numberOfParticles - 3;

  // empty movie.pdb file
  std::ofstream file("movie.pdb");
  file.close();

  // Initialize MD simulation
  logger.info("Class MD created.");
  logger.debug(repr());
  init();
  logger.debug(repr());
};

/**
 * @brief Initializes the simulation by setting velocities and positions.
 *
 * Randomizes velocities, zeroes total momentum, scales velocities to the desired temperature,
 * initializes particle positions on a lattice, and performs gradient descent to minimize energy.
 */
void MolecularDynamics::init()
{
  // Initialize velocities with random Gaussian distribution and compute total momentum
  for (double3& velocity : velocities)
  {
    velocity = double3(normal(), normal(), normal());
    totalMomentum += velocity;
  }
  totalMomentum /= static_cast<double>(numberOfParticles);
  logger.debug("(Init) initialized velocities.");

  // Zero total momentum and calculate initial kinetic energy
  kineticEnergy = 0.0;
  for (double3& velocity : velocities)
  {
    velocity -= totalMomentum;
    kineticEnergy += double3::dot(velocity, velocity);
  }
  observedTemperature = 2.0 * kineticEnergy / degreesOfFreedom;
  logger.debug("(Init) Zeroed momentum.");

  // Scale velocities to match the desired temperature
  double scale = std::sqrt(temperature * static_cast<double>(degreesOfFreedom) / kineticEnergy);
  for (double3& velocity : velocities)
  {
    velocity *= scale;
  }
  logger.debug("(Init) scaled velocities to temperature " + std::to_string(temperature));

  latticeInitialization();
  writePDB("movie.pdb", positions, boxSize, frameNumber);
  ++frameNumber;
  gradientDescent();

  // Recalculate old positions for Verlet integration
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    oldPositions[i] = positions[i] - dt * velocities[i];
  }
  logger.info("(Init) completed. Final energy: " + std::to_string(potentialEnergy));
}

/**
 * @brief Initializes particle positions on a cubic lattice with slight random perturbations.
 *
 * Distributes particles evenly on a grid within the simulation box, adding small random displacements to avoid overlap.
 */
void MolecularDynamics::latticeInitialization()
{
  // Calculate number of grid points and grid spacing based on number of particles
  size_t numGrids = static_cast<size_t>(std::round(std::pow(numberOfParticles, 1.0 / 3.0) + 1.5));
  gridSize = boxSize / (static_cast<double>(numGrids) + 2.0);
  logger.info("numGrids " + std::to_string(numGrids) + " gridSize " + std::to_string(gridSize));

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
}

/**
 * @brief Minimizes potential energy using steepest descent method.
 *
 * Iteratively adjusts particle positions in the direction of the force to reduce potential energy.
 * Accepts or rejects steps based on energy reduction and adjusts step size accordingly.
 */
void MolecularDynamics::gradientDescent()
{
  // Initialize step size and variables for gradient descent
  double dx = 0.2 * gridSize;
  double previousEnergy = 0.0;
  std::vector<double3> oldForces(numberOfParticles);
  for (size_t step = 0; step < 1000; ++step)
  {
    if (step == 0)
    {
      // On first step, calculate initial forces and energy
      calculateForce();
      previousEnergy = potentialEnergy;
      logger.debug("(Init) initial energy: " + std::to_string(previousEnergy));
    }

    // Save current positions and forces before attempting step
    std::copy(positions.begin(), positions.end(), oldPositions.begin());
    std::copy(forces.begin(), forces.end(), oldForces.begin());

    // Calculate maximum force component to determine scaling factor
    double maxForce;
    for (double3& force : forces)
    {
      maxForce = std::max(maxForce, abs(force.x));
      maxForce = std::max(maxForce, abs(force.y));
      maxForce = std::max(maxForce, abs(force.z));
    }
    double scaleGradient = dx / maxForce;

    // Update positions in the direction of forces (gradient descent step)
    for (size_t i = 0; i < numberOfParticles; ++i)
    {
      positions[i] += scaleGradient * forces[i];
      positions[i] = wrap(positions[i], boxSize);
    }

    // Recalculate forces and potential energy after position update
    calculateForce();

    // Accept or reject the new positions based on energy change, adjust step size accordingly
    if (potentialEnergy < previousEnergy)
    {
      logger.debug("(Init) accepted gradient descent at step " + std::to_string(step) +
                   " energy difference = " + std::to_string(potentialEnergy - previousEnergy));
      previousEnergy = potentialEnergy;
      dx = std::min(0.5 * boxSize, 1.2 * dx);
    }
    else
    {
      logger.debug("(Init) rejected gradient descent at step " + std::to_string(step) +
                   " energy difference = " + std::to_string(potentialEnergy - previousEnergy));
      std::copy(oldPositions.begin(), oldPositions.end(), positions.begin());
      std::copy(oldForces.begin(), oldForces.end(), forces.begin());
      dx *= 0.1;
    }
  }
  logger.debug(repr());
}

/**
 * @brief Calculates forces between particles and updates potential energy and pressure.
 *
 * Implements the Lennard-Jones potential with a cutoff radius.
 * Uses pairwise interactions and accumulates forces, potential energy, and pressure.
 */
void MolecularDynamics::calculateForce()
{
  // Initialize forces to zero and reset potential energy and pressure
  double cutOffSquared = cutOff * cutOff;

  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    forces[i] = double3();
  }

  potentialEnergy = 0.0;
  pressure = 0.0;

  // Loop over all unique particle pairs to compute interactions
  for (size_t i = 0; i < positions.size() - 1; ++i)
  {
    for (size_t j = i + 1; j < positions.size(); ++j)
    {
      // Calculate minimum image distance between particles
      double3 dr = positions[i] - positions[j];
      dr = wrap(dr, boxSize);

      // Get squared distance
      double r2 = double3::dot(dr, dr);

      // If within cutoff radius, compute Lennard-Jones potential and force
      if (r2 < cutOffSquared)
      {
        double r2i = 1.0 / r2;
        double r6i = r2i * r2i * r2i;

        // Update potential energy and pressure
        potentialEnergy += 4.0 * r6i * (r6i - 1.0) - cutOffEnergy;
        double virial = 48.0 * r6i * (r6i - 0.5);
        pressure += virial;

        // Add forces to particles i and j
        forces[i] += virial * r2i * dr;
        forces[j] -= virial * r2i * dr;
      }
    }
  }

  // Normalize pressure by dividing by 3 times the volume
  pressure /= 3.0 * volume;
}

/**
 * @brief Integrates particle positions and velocities using the Verlet algorithm.
 *
 * Updates positions and velocities, applies velocity scaling during equilibration,
 * and recalculates kinetic energy and total momentum.
 */
void MolecularDynamics::integrate()
{
  double dt2 = dt * dt;

  // Compute new positions and velocities using Verlet integration
  kineticEnergy = 0.0;
  std::vector<double3> newPositions(numberOfParticles);
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    newPositions[i] = 2.0 * positions[i] - oldPositions[i] + forces[i] * dt2;
    double3 dr = wrap(newPositions[i] - oldPositions[i], boxSize);
    velocities[i] = dr / (2.0 * dt);
    kineticEnergy += 0.5 * double3::dot(velocities[i], velocities[i]);
  }

  // Apply velocity scaling if in equilibration phase
  double velocityScaling = equilibration ? std::sqrt(temperature * degreesOfFreedom / (2.0 * kineticEnergy)) : 1.0;

  kineticEnergy = 0.0;
  totalMomentum = double3();

  // Recalculate kinetic energy and total momentum after scaling
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    velocities[i] *= velocityScaling;
    newPositions[i] = oldPositions[i] + 2.0 * velocities[i] * dt;
    kineticEnergy += 0.5 * double3::dot(velocities[i], velocities[i]);

    totalMomentum += velocities[i];

    unwrappedPositions[i] = newPositions[i];
    newPositions[i] = wrapFloor(newPositions[i], boxSize);
    positions[i] = wrapFloor(positions[i], boxSize);
  }

  // Update positions
  std::copy(positions.begin(), positions.end(), oldPositions.begin());
  std::copy(newPositions.begin(), newPositions.end(), positions.begin());

  pressure += 2.0 * kineticEnergy * numberOfParticles / (volume * 3.0 * numberOfParticles);
}

/**
 * @brief Runs the molecular dynamics simulation for a specified number of steps.
 *
 * Performs force calculations and integrations in each step, handles equilibration,
 * logs information, and samples properties.
 *
 * @param numberOfSteps The number of integration steps to perform.
 * @param equilibrate Whether the simulation is in the equilibration phase.
 */
void MolecularDynamics::run(size_t& numberOfSteps, bool equilibrate, bool outputPDB)
{
  if (!equilibrate && equilibration)
  {
    // Start modification
    // set baseline energy for drift
    baselineEnergy = totalEnergy;
    // End modification
  }

  equilibration = equilibrate;
  // Main simulation loop
  for (size_t step = 0; step < numberOfSteps; ++step)
  {
    // Update forces and integrate positions and velocities
    calculateForce();
    integrate();
    ++totalSteps;

    totalEnergy = kineticEnergy + potentialEnergy;
    if (step == 1) logger.info(repr());

    observedTemperature = 2.0 * kineticEnergy / degreesOfFreedom;

    // Start modification
    // Update drift energy when not in equilibration
    if (!equilibration)
    {
      driftEnergy += std::abs((totalEnergy - baselineEnergy) / (baselineEnergy * numberOfSteps));
    }
    // End modification

    // Log and sample data every 100 steps
    if (step % sampleFrequency == 0)
    {
      logger.info(repr());
      if (outputPDB) writePDB("movie.pdb", positions, boxSize, frameNumber);
      if (!equilibration)
      {
        ++frameNumber;
        thermoSampler.sample(observedTemperature, pressure, potentialEnergy, kineticEnergy, totalEnergy);
        rdfSampler.sample(positions);
        msdSampler.sample(unwrappedPositions, velocities);
      }
    }
  }

  if (!equilibrate)
  {
    logger.info(thermoSampler.repr());
  }
}

/**
 * @brief Returns a string representation of the current state of the simulation.
 *
 * Provides detailed information about simulation parameters and current measurements.
 *
 * @return A formatted string containing simulation data.
 */
std::string MolecularDynamics::repr()
{
  std::string s;
  s += "Molecular Dynamics program\n";
  s += "----------------------------\n";
  s += "Number of particles  : " + std::to_string(numberOfParticles) + "\n";
  s += "Temperature          : " + std::to_string(temperature) + "\n";
  s += "delta t              : " + std::to_string(dt) + "\n";
  s += "Box length           : " + std::to_string(boxSize) + "\n";
  s += "Density              : " + std::to_string(density) + "\n";
  s += "CutOff radius        : " + std::to_string(cutOff) + "\n";
  s += "CutOff energy        : " + std::to_string(cutOffEnergy) + "\n";
  s += "Steps run            : " + std::to_string(totalSteps) + "\n";
  s += "Equilibration        : " + std::string(equilibration ? "true" : "false") + "\n";
  s += "Observed temperature : " + std::to_string(observedTemperature) + "\n";
  s += "Pressure             : " + std::to_string(pressure) + "\n";
  s += "Potential energy     : " + std::to_string(potentialEnergy) + "\n";
  s += "Kinetic energy       : " + std::to_string(kineticEnergy) + "\n";
  s += "Total energy         : " + std::to_string(totalEnergy) + "\n";
  s += "Drift energy         : " +
       std::to_string(driftEnergy / std::max(1.0, static_cast<double>(rdfSampler.numberOfSamples))) + "\n";
  s += "\n";
  return s;
}
