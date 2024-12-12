#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <vector>

#include "double3.h"

struct SampleRDF
{
  int numberOfSamples = 0;
  int numberOfBins = 100;
  std::vector<double> histogram;

  int numberOfParticles;
  double boxSize;

  double cutOff;
  double delta;

  std::vector<double> r;

  SampleRDF(int numberOfParticles, double boxSize, double cutOff);

  void sample(std::vector<double3>& positions);
  pybind11::array_t<double> getResults();
};

struct SampleMSD
{
  // Time tracking variables
  int time = 0;
  int originTimeCounter = 0;
  int numberOfCorrelationTimes = 7500;
  int maxOriginTimes = 250;
  int resetOriginInterval = 50;

  int numberOfParticles;
  double boxSize;
  double sampleTime;

  // Data storage vectors
  std::vector<int> sampleCounts;
  std::vector<int> originTimes;
  std::vector<double> meanSquareDisplacements;
  std::vector<double> velocityAutocorrelation;
  std::vector<std::vector<double3>> velocityAtOrigin;
  std::vector<std::vector<double3>> positionAtOrigin;

  SampleMSD(int numberOfParticles, double boxSize, double sampleTime);

  void sample(std::vector<double3>& unwrappedPositions, std::vector<double3>& velocities);
  pybind11::array_t<double> getResults();
};
