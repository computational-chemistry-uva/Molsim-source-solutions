#pragma once

#include <algorithm>
#include <cmath>
#include <string>

union double3
{
  double v[4];
  struct
  {
    double x, y, z, w;
  };

  double3() : x(0.0), y(0.0), z(0.0), w(0.0) {}
  double3(double v) : x(v), y(v), z(v), w(0.0) {}
  double3(double x, double y, double z) : x(x), y(y), z(z), w(0.0) {}
  double3(double a[3]) : x(double(a[0])), y(double(a[1])), z(double(a[2])), w(0.0) {}

  bool operator==(double3 const& rhs) const { return (x == rhs.x) && (y == rhs.y) && (z == rhs.z); }

  inline double& operator[](size_t i) { return v[i]; }
  inline const double& operator[](size_t i) const { return v[i]; }

  double3 operator-() const { return double3(-this->x, -this->y, -this->z); }
  double3& operator+=(const double3& b)
  {
    this->x += b.x, this->y += b.y, this->z += b.z;
    return *this;
  }
  double3& operator-=(const double3& b)
  {
    this->x -= b.x, this->y -= b.y, this->z -= b.z;
    return *this;
  }
  double3& operator/=(const double& s)
  {
    this->x /= s, this->y /= s, this->z /= s;
    return *this;
  }
  double3& operator*=(const double& s)
  {
    this->x *= s, this->y *= s, this->z *= s;
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const double3& vec);

  struct KeyHash
  {
    size_t operator()(const double3& k) const
    {
      size_t h1 = std::hash<double>()(k.x);
      size_t h2 = std::hash<double>()(k.y);
      size_t h3 = std::hash<double>()(k.z);
      return (h1 ^ (h2 << 1)) ^ h3;
    }
  };
  struct KeyEqual
  {
    bool operator()(const double3& lhs, const double3& rhs) const
    {
      return (std::fabs(lhs.x - rhs.x) < 1e-3) && std::fabs(lhs.y - rhs.y) < 1e-3 && std::fabs(lhs.z == rhs.z) < 1e-3;
    }
  };

  inline static double3 abs(double3 v1) { return double3(std::abs(v1.x), std::abs(v1.y), std::abs(v1.z)); }
  inline static double dot(const double3& v1, const double3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }

  inline static double3 max(const double3& v1, const double3& v2)
  {
    return double3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
  }
  inline static double3 min(const double3& v1, const double3& v2)
  {
    return double3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
  }
  inline static double3 cross(const double3& v1, const double3& v2)
  {
    return double3(v1.y * v2.z - v2.y * v1.z, v1.z * v2.x - v2.z * v1.x, v1.x * v2.y - v2.x * v1.y);
  }
  inline static double3 normalize(const double3& v)
  {
    double f = 1.0 / sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
    return double3(f * v.x, f * v.y, f * v.z);
  }

  std::string to_string()
  {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
  }
};

inline double3 operator*(const double3& a, const double3& b) { return double3(a.x * b.x, a.y * b.y, a.z * b.z); }

inline double3 operator-(const double3& a, const double3& b) { return double3(a.x - b.x, a.y - b.y, a.z - b.z); }

inline double3 operator+(const double3& a, const double3& b) { return double3(a.x + b.x, a.y + b.y, a.z + b.z); }

inline double3 operator/(const double3& a, const double3& b) { return double3(a.x / b.x, a.y / b.y, a.z / b.z); }

inline double3 operator*(const double3& a, double b) { return double3(a.x * b, a.y * b, a.z * b); }

inline double3 operator*(const double& a, const double3& b) { return double3(a * b.x, a * b.y, a * b.z); }

inline double3 operator/(const double3& a, double b) { return double3(a.x / b, a.y / b, a.z / b); }

inline double3 sqrt(const double3& a) { return double3(std::sqrt(a.x), std::sqrt(a.y), std::sqrt(a.z)); }
