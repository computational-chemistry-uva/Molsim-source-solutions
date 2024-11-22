#pragma once

/**
 * \brief Represents the energy and virial of a system or particle.
 *
 * The EnergyVirial struct encapsulates the total energy and virial contributions
 * which are used in the calculations of system properties in Monte Carlo simulations.
 */
struct EnergyVirial
{
  double energy;
  double virial;

  /**
   * \brief Constructs an EnergyVirial object with specified energy and virial.
   *
   * \param energy Initial energy value.
   * \param virial Initial virial value.
   */
  EnergyVirial(double energy = 0.0, double virial = 0.0) : energy(energy), virial(virial) {};

  EnergyVirial &operator+=(const EnergyVirial &other)
  {
    this->energy += other.energy;
    this->virial += other.virial;
    return *this;
  };
};
static inline EnergyVirial operator-(const EnergyVirial &a, const EnergyVirial &b)
{
  return {a.energy - b.energy, a.virial - b.virial};
};

static inline EnergyVirial operator/(const EnergyVirial &a, const EnergyVirial &b)
{
  return {a.energy / b.energy, a.virial / b.virial};
};