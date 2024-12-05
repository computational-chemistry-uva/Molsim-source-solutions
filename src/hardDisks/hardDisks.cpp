#include "hardDisks.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

HardDisks::HardDisks(size_t numberOfInitCycles, size_t numberOfProdCycles, size_t numberOfParticles,
                     double maxDisplacement, size_t sampleFrequency, double boxSize, size_t rdfBins, bool runStatic)
    : numberOfInitCycles(numberOfInitCycles),
      numberOfProdCycles(numberOfProdCycles),
      numberOfParticles(numberOfParticles),
      maxDisplacement(maxDisplacement),
      sampleFrequency(sampleFrequency),
      boxSize(boxSize),
      positions(numberOfParticles),
      rdf(rdfBins),
      rdfBins(rdfBins)
{
  initialize();
  delta = 0.5 * boxSize / (static_cast<double>(rdfBins));
  halfBoxSizeSq = 0.25 * boxSize * boxSize;

  method = runStatic ? Method::Static : Method::Dynamic;
}

void HardDisks::initialize()
{
  size_t latticeSites = static_cast<size_t>(boxSize);
  std::vector<double2> latticePositions(latticeSites * latticeSites);
  for (size_t ix = 0; ix < latticeSites; ix++)
  {
    for (size_t iy = 0; iy < latticeSites; iy++)
    {
      latticePositions[ix + latticeSites * iy] = double2(ix, iy);
    }
  }
  std::sample(latticePositions.begin(), latticePositions.end(), positions.begin(), positions.size(), gen);
}

void HardDisks::run()
{
  reset();
  switch (method)
  {
    case Method::Dynamic:
      dynamicRun();
    case Method::Static:
      staticRun();
    default:
      break;
  }

  acceptanceRatio = static_cast<double>(numberOfAcceptedMoves) / static_cast<double>(numberOfAttemptedMoves);
  std::cout << "Acceptance fraction: " << acceptanceRatio << std::endl;
}

void HardDisks::reset()
{
  if (numberOfAttemptedMoves > 0)
  {
    std::vector<double> rdf(rdfBins);
    numberOfAttemptedMoves = 0;
    numberOfAcceptedMoves = 0;
    numberOfSamples = 0;
    initialize();
  }
}

void HardDisks::dynamicRun()
{
  std::uniform_int_distribution<> indexDist(0, numberOfParticles - 1);
  std::uniform_real_distribution<> uniform(-0.5, 0.5);
  for (size_t cycle = 0; cycle < numberOfInitCycles + numberOfProdCycles; cycle++)
  {
    numberOfAttemptedMoves++;

    int particleIdx = indexDist(gen);
    double2 newPosition = positions[particleIdx];
    newPosition.x += maxDisplacement * uniform(gen);
    newPosition.y += maxDisplacement * uniform(gen);

    bool overlap = false;
    for (size_t i = 0; i < numberOfParticles; i++)
    {
      if (i != particleIdx)
      {
        double2 dr = newPosition - positions[i];

        // first wrap to [-L/2, L/2]
        dr.x -= boxSize * std::rint(dr.x / boxSize);
        dr.y -= boxSize * std::rint(dr.y / boxSize);

        double r2 = dot(dr, dr);
        if (r2 < 1.0)
        {
          overlap = true;
        }
      }
    }

    if (!overlap)
    {
      positions[particleIdx] = newPosition;
      numberOfAcceptedMoves++;
    }
    if (cycle % sampleFrequency == 0 && cycle > numberOfInitCycles)
    {
      sampleRDF();
    }
  }
}

void HardDisks::staticRun()

{
  std::uniform_real_distribution<> uniform(0.0, boxSize);

  for (size_t cycle = 0; cycle < numberOfInitCycles + numberOfProdCycles; cycle++)
  {
    numberOfAttemptedMoves++;
    std::vector<double2> newPositions(positions.size());
    for (size_t i = 0; i < numberOfParticles; i++)
    {
      newPositions[i] = double2(uniform(gen), uniform(gen));
    }

    bool overlap = false;
    for (size_t i = 0; i < numberOfParticles - 1; i++)
    {
      for (size_t j = i + 1; j < numberOfParticles; j++)
      {
        double2 dr = newPositions[i] - newPositions[j];

        // first wrap to [-L/2, L/2]
        dr.x -= boxSize * std::rint(dr.x / boxSize);
        dr.y -= boxSize * std::rint(dr.y / boxSize);

        double r2 = dot(dr, dr);
        if (r2 < 1.0)
        {
          overlap = true;
          break;
        }
      }
      if (overlap)
      {
        break;
      }
    }
    if (!overlap)
    {
      positions = newPositions;
      numberOfAcceptedMoves++;
    }
    if (cycle % sampleFrequency == 0 && cycle > numberOfInitCycles)
    {
      sampleRDF();
    }
  }
}

void HardDisks::sampleRDF()
{
  numberOfSamples++;
  for (size_t i = 0; i < numberOfParticles - 1; i++)
  {
    for (size_t j = i + 1; j < numberOfParticles; j++)
    {
      double2 dr = positions[i] - positions[j];
      dr.x -= boxSize * std::floor(dr.x / boxSize);
      dr.y -= boxSize * std::floor(dr.y / boxSize);

      double r2 = dot(dr, dr);
      if (r2 < halfBoxSizeSq)
      {
        size_t idx = static_cast<size_t>(std::sqrt(r2) / delta);
        rdf[idx] += 2.0;
      }
    }
  }
}

pybind11::array_t<double> HardDisks::getRDF()
{
  std::vector<double> normalizedRDF(rdf.size());
  for (size_t i = 0; i < rdfBins; i++)
  {
    double areaDiff = M_PI * delta * delta * ((i + 1) * (i + 1) - i * i);
    double invDensitySq = boxSize * boxSize / (numberOfParticles * (numberOfParticles - 1));
    normalizedRDF[i] = 4.0 * rdf[i] * invDensitySq / (areaDiff * static_cast<double>(numberOfSamples));
  }
  return pybind11::array_t<double>(normalizedRDF.size(), normalizedRDF.data());
}
