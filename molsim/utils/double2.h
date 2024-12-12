#pragma once

struct double2
{
  double x;
  double y;

  double2(double x = 0.0, double y = 0.0) : x(x), y(y) {};
};
static inline double2 operator+(const double2& a, const double2& b) { return double2(a.x + b.x, a.y + b.y); };
static inline double2 operator-(const double2& a, const double2& b) { return double2(a.x - b.x, a.y - b.y); };
static inline double dot(const double2& a, const double2& b) { return a.x * b.x + a.y * b.y; };
