#include <random>
#include <vector>

#include "double3.h"

struct VelocityScaling
{
  double temperature;
  size_t degreesOfFreedom;

  VelocityScaling(double temperature, size_t degreesOfFreedom);
  void scale(std::vector<double3>& velocities, double& kineticEnergy);
};

struct NoseHooverNVT
{
  double temperature;
  double degreesOfFreedom;
  double timescaleParameter;
  double timeStep;
  double thermostatMass;
  double thermostatVelocity;
  double thermostatForce;
  double thermostatPosition;

  std::mt19937 mt;
  std::normal_distribution<double> normal_dist;

  NoseHooverNVT(double temperature, size_t degreesOfFreedom, double timescaleParameter, double timeStep,
                size_t seed = 12);

  void scale(std::vector<double3>& velocities, double& kineticEnergy);
  double getEnergy();
};