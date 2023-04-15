#pragma once

#include "hittable.hpp"
#include "vec3.hpp"

#include <cassert>

class sphere : public hittable {
public:
  sphere() = default;
  sphere(point3 cen, double r, std::shared_ptr<material> m)
      : center(cen), radius(r), mat_ptr(m) {}
  virtual ~sphere() = default;

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

private:
  point3 center;
  double radius;
  std::shared_ptr<material> mat_ptr;

  static void get_sphere_uv(const point3 &p, double &u, double &v) {
    // p: point on a unit sphere, centered at the origin
    // u: returned value [0,1] of angle around the Y axis from X=-1 (azimuth
    // angle)
    // v: return value [0,1] of angle from Y=-1 to Y+1 (polar angle)
    //    <1  0  0> yields <0.50, 0.50>    <-1  0  0> yields <0.00 0.50>
    //    <0  1  0> yields <0.50, 1.00>    < 0 -1  0> yields <0.50 0.00>
    //    <0  0  1> yields <0.25, 0.50>    < 0  0 -1> yields <0.75 0.50>

    auto theta = std::acos(-p.y());
    auto phi = std::atan2(-p.z(), p.x()) + pi;

    // normalize aximuth and polar angles to give an index into a texture
    u = phi / (2 * pi);
    v = theta / pi;
  }
};
