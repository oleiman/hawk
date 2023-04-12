#pragma once

#include "hittable.hpp"
#include "vec3.hpp"

class sphere : public hittable {
public:
  sphere() = default;
  sphere(point3 cen, double r, std::shared_ptr<material> m)
      : center(cen), radius(r), mat_ptr(m) {}
  virtual ~sphere() = default;

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

private:
  point3 center;
  double radius;
  std::shared_ptr<material> mat_ptr;
};
