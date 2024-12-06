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

  int numberOfInitCycles;
  int numberOfProdCycles;
  int numberOfParticles;
  double maxDisplacement;
  int sampleFrequency;
  double boxSize;

  Method method;

  std::vector<double2> positions;
  std::vector<std::vector<double2>> samplePositions;

  int numberOfAttemptedMoves{0};
  int numberOfAcceptedMoves{0};
  double acceptanceRatio;
  int numberOfSamples{0};

  std::vector<double> rdf;
  int rdfBins;
  double delta;
  double halfBoxSizeSq;

  std::random_device rd;
  std::mt19937 gen{rd()};

  HardDisks(int numberOfInitCycles, int numberOfProdCycles, int numberOfParticles, double maxDisplacement,
            int sampleFrequency, double boxSize, int rdfBins, bool runStatic);

  void initialize();
  void run();
  void dynamicRun();
  void staticRun();
  void sampleRDF();
  void reset();
  pybind11::array_t<double> getRDF();
};
