#include "vec3.hpp"
#include "hawk.h"

vec3 vec3::random() {
  return vec3(random_double(), random_double(), random_double());
}

vec3 vec3::random(double min, double max) {
  return vec3(random_double(min, max), random_double(min, max),
              random_double(min, max));
}

vec3 random_in_unit_disk() {
  while (true) {
    auto p = vec3(random_double(-1, 1), random_double(-1, 1), 0);
    if (p.length_squared() >= 1)
      continue;
    return p;
  }
}
