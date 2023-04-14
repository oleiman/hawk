#pragma once

#include "hawk.h"

class aabb {
public:
  aabb() {}
  aabb(const point3 &a, const point3 &b) : minimum(a), maximum(b) {}

  point3 min() const { return minimum; }
  point3 max() const { return maximum; }

  bool hit_slow(const ray &r, double t_min, double t_max) const {
    for (int a = 0; a < 3; ++a) {
      auto t0 = std::fmin((minimum[a] - r.origin()[a]) / r.direction()[a],
                          (maximum[a] - r.origin()[a]) / r.direction()[a]);
      auto t1 = std::fmax((minimum[a] - r.origin()[a]) / r.direction()[a],
                          (maximum[a] - r.origin()[a]) / r.direction()[a]);

      t_min = std::fmax(t0, t_min);
      t_max = std::fmin(t1, t_max);
      if (t_max <= t_min) {
        return false;
      }
    }
    return true;
  }

  bool hit(const ray &r, double t_min, double t_max) const {
    for (int a = 0; a < 3; ++a) {
      auto invD = 1.0f / r.direction()[a];
      auto t0 = (min()[a] - r.origin()[a]) * invD;
      auto t1 = (max()[a] - r.origin()[a]) * invD;
      if (invD < 0.0f) {
        std::swap(t0, t1);
      }
      t_min = t0 > t_min ? t0 : t_min;
      t_max = t1 < t_max ? t1 : t_max;
      if (t_max <= t_min) {
        return false;
      }
    }
    return true;
  }

private:
  point3 minimum;
  point3 maximum;
};

aabb surrounding_box(const aabb &box0, const aabb &box1);
