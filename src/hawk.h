#pragma once

#include <cmath>
#include <limits>
#include <memory>
#include <random>

// Constants

constexpr double infinity = std::numeric_limits<double>::infinity();
constexpr double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
  return degrees * pi / 180.0;
}

inline double clamp(double x, double min, double max) {
  if (x < min)
    return min;
  if (x > max)
    return max;
  return x;
}

// Common Headers

#include "ray.hpp"
#include "vec3.hpp"
