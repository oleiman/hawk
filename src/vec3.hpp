#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

class vec3 {
public:
  vec3() : e{0, 0, 0} {}
  vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

  auto x() const { return e[0]; }
  auto y() const { return e[1]; }
  auto z() const { return e[2]; }

  vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
  double operator[](int i) const {
    assert(i < e.size());
    return e[i];
  }
  double &operator[](int i) {
    assert(i < e.size());
    return e[i];
  }

  vec3 &operator+=(const vec3 &v) {
    e[0] += v[0];
    e[1] += v[1];
    e[2] += v[2];
    return *this;
  }

  vec3 &operator*=(const double t) {
    std::transform(
        e.begin(), e.end(), e.begin(),
        std::bind(std::multiplies<double>(), std::placeholders::_1, t));
    return *this;
  }

  vec3 &operator/=(const double t) { return *this *= 1 / t; }

  double length() const { return std::sqrt(length_squared()); }

  double length_squared() const {
    return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
  }

  bool near_zero() const {
    // Return true if the vector is close to zero in all dimensions
    const auto s = 1e-8;
    return (std::fabs(e[0]) < s) && //
           (std::fabs(e[1]) < s) && //
           (std::fabs(e[2]) < s);
  }

  static vec3 random();
  static vec3 random(double min, double max);

private:
  std::array<double, 3> e;
};

using point3 = vec3; // 3D point
using color = vec3;  // RGB color

// vec3 Utility Functions

inline std::ostream &operator<<(std::ostream &out, const vec3 &v) {
  return out << v[0] << ' ' << v[1] << ' ' << v[2];
}

inline vec3 operator+(const vec3 &u, const vec3 &v) {
  return vec3(u[0] + v[0], u[1] + v[1], u[2] + v[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
  return vec3(u[0] - v[0], u[1] - v[1], u[2] - v[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
  return vec3(u[0] * v[0], u[1] * v[1], u[2] * v[2]);
}

inline vec3 operator*(double t, const vec3 &v) {
  return vec3(t * v[0], t * v[1], t * v[2]);
}

inline vec3 operator*(const vec3 &v, double t) { return t * v; }

inline vec3 operator/(vec3 v, double t) { return (1 / t) * v; }

inline double dot(const vec3 &u, const vec3 &v) {
  return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
  return vec3(u[1] * v[2] - u[2] * v[1], //
              u[2] * v[0] - u[0] * v[2], //
              u[0] * v[1] - u[1] * v[0]);
}

inline vec3 unit_vector(vec3 v) { return v / v.length(); }

inline vec3 random_in_unit_sphere() {
  while (true) {
    auto p = vec3::random(-1, 1);
    if (p.length_squared() >= 1)
      continue;
    return p;
  }
}

inline vec3 random_unit_vector() {
  return unit_vector(random_in_unit_sphere());
}

inline vec3 random_in_hemisphere(const vec3 &normal) {
  vec3 in_unit_sphere = random_in_unit_sphere();
  if (dot(in_unit_sphere, normal) >
      0.0) { // In the same hemisphere as the normal
    return in_unit_sphere;
  } else {
    return -in_unit_sphere;
  }
}

vec3 random_in_unit_disk();

inline vec3 reflect(const vec3 &v, const vec3 &n) {
  return v - (2 * dot(v, n) * n);
}

inline vec3 refract(const vec3 &uv, const vec3 &n, double etai_over_etat) {
  auto cos_theta = std::fmin(dot(-uv, n), 1.0);
  vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
  vec3 r_out_parallel =
      -std::sqrt(std::fabs(1.0 - r_out_perp.length_squared())) * n;
  return r_out_perp + r_out_parallel;
}
