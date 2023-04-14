#pragma once

#include "hawk.h"

#include "hittable.hpp"

class moving_sphere : public hittable {
public:
  moving_sphere() = default;
  moving_sphere(point3 cen0, point3 cen1, double time0, double time1, double r,
                std::shared_ptr<material> m)
      : center0(cen0), center1(cen1), time0(time0), time1(time1), radius(r),
        mat_ptr(m) {}

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;
  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

  inline point3 center(double time) const {
    return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
  }

private:
  point3 center0, center1;
  double time0, time1;
  double radius;
  std::shared_ptr<material> mat_ptr;
};
