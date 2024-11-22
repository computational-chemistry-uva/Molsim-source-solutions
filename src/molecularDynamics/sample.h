#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vector>

#include "double3.h"

struct SampleThermodynamicalAverages
{
  size_t numberOfSamples = 0;

  std::vector<double> vTemperature;
  std::vector<double> vPressure;
  std::vector<double> vPotentialEnergy;
  std::vector<double> vKineticEnergy;
  std::vector<double> vTotalEnergy;

  SampleThermodynamicalAverages() {};
  void sample(double temperature, double pressure, double potentialEnergy, double kineticEnergy, double totalEnergy);
  std::string repr();
};

struct SampleRDF
{
  size_t numberOfSamples = 0;
  size_t numberOfBins = 500;
  std::vector<double> histogram;

  size_t numberOfParticles;
  double boxSize;

  double cutOff;
  double delta;

  std::vector<double> r;

  SampleRDF(size_t numberOfParticles, double boxSize);

  void sample(std::vector<double3>& positions);
  pybind11::array_t<double> getResults();
};

struct SampleMSD
{
  // Time tracking variables
  size_t time = 0;
  size_t originTimeCounter = 0;
  size_t numberOfCorrelationTimes = 7500;
  size_t maxOriginTimes = 250;
  size_t resetOriginInterval = 50;

  size_t numberOfParticles;
  double boxSize;
  double sampleTime;

  // Data storage vectors
  std::vector<size_t> sampleCounts;
  std::vector<size_t> originTimes;
  std::vector<double> meanSquareDisplacements;
  std::vector<double> velocityAutocorrelation;
  std::vector<std::vector<double3>> velocityAtOrigin;
  std::vector<std::vector<double3>> positionAtOrigin;

  SampleMSD(size_t numberOfParticles, double boxSize, double sampleTime);

  void sample(std::vector<double3>& unwrappedPositions, std::vector<double3>& velocities);
  pybind11::array_t<double> getResults();
};
