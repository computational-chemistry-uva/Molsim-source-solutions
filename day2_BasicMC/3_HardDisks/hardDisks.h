#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <iostream>
#include <random>
#include <string>
#include <vector>

struct double2
{
  double x;
  double y;

  double2(double x = 0.0, double y = 0.0) : x(x), y(y) {};
};
double dot(const double2& a, const double2& b) { return a.x * b.x + a.y * b.y; };
double2 operator+(const double2& a, const double2& b) { return double2(a.x + b.x, a.y + b.y); };
double2 operator-(const double2& a, const double2& b) { return double2(a.x - b.x, a.y - b.y); };

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
            size_t sampleFrequency, double boxSize, size_t rdfBins, Method method);

  void initialize();
  void run();
  void dynamicRun();
  void staticRun();
  void sampleRDF();
  void reset();
  pybind11::array_t<double> getRDF();
};
