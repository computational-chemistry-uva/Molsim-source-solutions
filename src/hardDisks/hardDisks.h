#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <random>
#include <string>
#include <vector>

#include "double2.h"

struct HardDisks
{
  enum class Method
  {
    Dynamic,
    Static
  };

  size_t numberOfInitCycles;
  size_t numberOfProdCycles;
  size_t numberOfParticles;
  double maxDisplacement;
  size_t sampleFrequency;
  double boxSize;

  Method method;

  std::vector<double2> positions;
  std::vector<std::vector<double2>> samplePositions;

  size_t numberOfAttemptedMoves{0};
  size_t numberOfAcceptedMoves{0};
  double acceptanceRatio;
  size_t numberOfSamples{0};

  std::vector<double> rdf;
  size_t rdfBins;
  double delta;
  double halfBoxSizeSq;

  std::random_device rd;
  std::mt19937 gen{rd()};

  HardDisks(size_t numberOfInitCycles, size_t numberOfProdCycles, size_t numberOfParticles, double maxDisplacement,
            size_t sampleFrequency, double boxSize, size_t rdfBins, bool runStatic);

  void initialize();
  void run();
  void dynamicRun();
  void staticRun();
  void sampleRDF();
  void reset();
  pybind11::array_t<double> getRDF();
};
