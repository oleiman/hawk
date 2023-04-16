#pragma once

#include "hawk.h"

#include "hittable.hpp"
#include "material.hpp"
#include "texture.hpp"

class constant_medium : public hittable {
public:
  constant_medium(std::shared_ptr<hittable> b, double d,
                  std::shared_ptr<texture> a)
      : boundary(b), neg_inv_density(-1 / d),
        phase_function(std::make_shared<isotropic>(a)) {}

  constant_medium(std::shared_ptr<hittable> b, double d, color c)
      : boundary(b), neg_inv_density(-1 / d),
        phase_function(std::make_shared<isotropic>(c)) {}

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

private:
  std::shared_ptr<hittable> boundary;
  double neg_inv_density;
  std::shared_ptr<material> phase_function;
};
