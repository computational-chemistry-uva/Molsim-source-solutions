#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "double3.h"

void writePDB(const std::string &fileName, std::vector<double3> &positions, double &boxSize, size_t &frameNumber)
{
  std::ofstream file(fileName, std::ios::app);

  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << fileName << std::endl;
    return;
  }

  file << "MODEL " << std::setw(9) << frameNumber << "\n";
  file << "CRYST1   " << std::setw(6) << std::fixed << std::setprecision(3) << boxSize << "   " << std::setw(6)
       << boxSize << "   " << std::setw(6) << boxSize << "  90.00  90.00  90.00 P 1         1\n";

  for (size_t i = 0; i < positions.size(); ++i)
  {
    double3 wrapped = wrapFloor(positions[i], boxSize);
    file << std::left << "ATOM" << std::right << std::setw(7) << i << "  H" << std::setw(16) << i << std::fixed
         << std::setprecision(3) << std::setw(8) << wrapped.x << std::setw(8) << wrapped.y << std::setw(8) << wrapped.z
         << std::setw(21) << " "  // To maintain space between the coordinates and "H"
         << "H" << std::endl;
  }

  file << "ENDMDL\n\n";
  file.close();
}
