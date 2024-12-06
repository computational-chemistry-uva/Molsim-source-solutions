#pragma once

#include <pybind11/iostream.h>

#include "double3.h"

class Logger
{
 public:
  Logger(int logLevel) : logLevel(logLevel) {}

  void debug(const std::string &message)
  {
    if (logLevel < 1)
    {
      pybind11::print("[DEBUG]: ", message);
    }
  }

  void info(const std::string &message)
  {
    if (logLevel < 2)
    {
      pybind11::print("[INFO]: ", message);
    }
  }

  void error(const std::string &message) { pybind11::print("[ERROR]: ", message); }

 private:
  int logLevel{0};
};

/**
 * @brief Wraps a vector into the simulation box using periodic boundary
 * conditions.
 *
 * Adjusts each component of the vector by subtracting the box size times the
 * nearest integer of the component divided by the box size.
 *
 * @param v The vector to be wrapped.
 * @param boxSize The size of the simulation box.
 * @return The wrapped vector within the simulation box.
 */
static double3 wrap(const double3 &v, const double &boxSize)
{
  // Apply periodic boundary conditions using round
  return double3(v.x - (boxSize * std::round(v.x / boxSize)), v.y - (boxSize * std::round(v.y / boxSize)),
                 v.z - (boxSize * std::round(v.z / boxSize)));
}

/**
 * @brief Wraps a vector into the simulation box using floor-based periodic
 * boundary conditions.
 *
 * Adjusts each component of the vector by subtracting the box size times the
 * floor of the component divided by the box size.
 *
 * @param v The vector to be wrapped.
 * @param boxSize The size of the simulation box.
 * @return The wrapped vector within the simulation box.
 */
static double3 wrapFloor(const double3 &v, const double &boxSize)
{
  // Apply periodic boundary conditions using floor
  return double3(v.x - (boxSize * std::floor(v.x / boxSize)), v.y - (boxSize * std::floor(v.y / boxSize)),
                 v.z - (boxSize * std::floor(v.z / boxSize)));
}

static double average(std::vector<double> &data)
{
  return (data.size()) ? std::accumulate(data.begin(), data.end(), 0.0) / static_cast<double>(data.size()) : 0.0;
}

static double variance(std::vector<double> &data)
{
  double size = static_cast<double>(data.size());
  if (size == 0.0)
  {
    return 0.0;
  }
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / size;
  double meanOfSquares = std::inner_product(data.begin(), data.end(), data.begin(), 0.0) / size;
  return meanOfSquares - mean * mean;
}

/**
 * @brief Calculates the block average and variance of a dataset.
 *
 * Divides the data into 5 blocks, computes the average for each block, and then
 * calculates the overall average and variance from these block averages.
 *
 * @param data The dataset to be averaged.
 * @return A pair containing the total average and the variance between block
 * averages.
 */
static std::pair<double, double> blockAverage(std::vector<double> &data)
{
  std::vector<double> averages(5);
  std::vector<int> counts(5);
  int samplesPerBin = std::round(data.size() / 5);

  // Sum data into blocks
  for (int i = 0; i < data.size(); ++i)
  {
    int bin = std::floor(i * 5 / data.size());
    averages[bin] += data[i];
    ++counts[bin];
  }

  double totalAverage = 0.0;
  // Calculate average for each block
  for (int i = 0; i < 5; ++i)
  {
    averages[i] /= counts[i];
    totalAverage += averages[i] / 5;
  }

  double totalVariance = 0.0;
  // Compute variance from block averages
  for (int i = 0; i < 5; ++i)
  {
    totalVariance += (averages[i] - totalAverage) * (averages[i] - totalAverage);
  }
  totalVariance /= 4;  // Using N-1 for variance calculation

  return std::make_pair(totalAverage, totalVariance);
};
