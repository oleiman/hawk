#pragma once

#include "hawk.h"

#include "hittable.hpp"

class xy_rect : public hittable {
public:
  xy_rect() {}

  xy_rect(double x0, double x1, double y0, double y1, double k,
          std::shared_ptr<material> mat)
      : x0(x0), x1(x1), y0(y0), y1(y1), k(k), mp(mat) {}

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

private:
  double x0, x1, y0, y1, k;
  std::shared_ptr<material> mp;
};

class xz_rect : public hittable {
public:
  xz_rect() {}

  xz_rect(double x0, double x1, double z0, double z1, double k,
          std::shared_ptr<material> mat)
      : x0(x0), x1(x1), z0(z0), z1(z1), k(k), mp(mat) {}

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

private:
  double x0, x1, z0, z1, k;
  std::shared_ptr<material> mp;
};

class yz_rect : public hittable {
public:
  yz_rect() {}

  yz_rect(double y0, double y1, double z0, double z1, double k,
          std::shared_ptr<material> mat)
      : y0(y0), y1(y1), z0(z0), z1(z1), k(k), mp(mat) {}

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

private:
  double y0, y1, z0, z1, k;
  std::shared_ptr<material> mp;
};
