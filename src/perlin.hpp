#pragma once

#include "hawk.h"

#include <array>

class perlin {
public:
  perlin() {
    std::transform(ranvec.begin(), ranvec.end(), ranvec.begin(),
                   [](auto &f) { return unit_vector(vec3::random(-1, 1)); });
    perlin_generate_perm(perm_x);
    perlin_generate_perm(perm_y);
    perlin_generate_perm(perm_z);
  }

  double noise(const point3 &p) const {
    auto u = p.x() - std::floor(p.x());
    auto v = p.y() - std::floor(p.y());
    auto w = p.z() - std::floor(p.z());

    auto i = static_cast<int>(std::floor(p.x()));
    auto j = static_cast<int>(std::floor(p.y()));
    auto k = static_cast<int>(std::floor(p.z()));
    vec3 c[2][2][2];

    for (uint8_t di = 0; di < 2; ++di) {
      for (uint8_t dj = 0; dj < 2; ++dj) {
        for (uint8_t dk = 0; dk < 2; ++dk) {
          c[di][dj][dk] = ranvec[       //
              perm_x[(i + di) & 0xFF] ^ //
              perm_y[(j + dj) & 0xFF] ^ //
              perm_z[(k + dk) & 0xFF]   //
          ];
        }
      }
    }

    return perlin_interp(c, u, v, w);
  }

  double turb(const point3 &p, int depth = 7) const {
    auto accum = 0.0;
    auto temp_p = p;
    auto weight = 1.0;
    for (int i = 0; i < depth; ++i) {
      accum += weight * noise(temp_p);
      weight *= 0.5;
      temp_p *= 2;
    }
    return std::fabs(accum);
  }

private:
  static constexpr int point_count = 256;
  std::array<vec3, point_count> ranvec = {};
  std::array<uint8_t, point_count> perm_x = {};
  std::array<uint8_t, point_count> perm_y = {};
  std::array<uint8_t, point_count> perm_z = {};

  static void perlin_generate_perm(std::array<uint8_t, point_count> &arr) {
    for (int i = 0; i < arr.size(); ++i) {
      arr[i] = i;
    }
    permute(arr);
  }

  static void permute(std::array<uint8_t, point_count> &arr) {
    for (int i = arr.size() - 1; i > 0; --i) {
      int target = random_int(0, i);
      std::swap(arr[i], arr[target]);
    }
  }

  static double trilinear_interp(double c[2][2][2], double u, double v,
                                 double w) {
    auto accum = 0.0;
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          accum += (i * u + (1 - i) * (1 - u)) * //
                   (j * v + (1 - j) * (1 - v)) * //
                   (k * w + (1 - k) * (1 - w)) * //
                   c[i][j][k];
        }
      }
    }
    return accum;
  }

  static double perlin_interp(vec3 c[2][2][2], double u, double v, double w) {
    auto uu = u * u * (3 - 2 * u);
    auto vv = v * v * (3 - 2 * v);
    auto ww = w * w * (3 - 2 * w);
    auto accum = 0.0;

    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          vec3 weight_v(u - i, v - j, w - k);
          accum += (i * uu + (1 - i) * (1 - uu)) * //
                   (j * vv + (1 - j) * (1 - vv)) * //
                   (k * ww + (1 - k) * (1 - ww)) * //
                   dot(c[i][j][k], weight_v);
        }
      }
    }
    return accum;
  }
};
