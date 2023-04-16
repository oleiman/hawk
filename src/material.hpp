#pragma once

#include "hawk.h"

#include "texture.hpp"

struct hit_record;

class material {
public:
  virtual ~material() = default;
  virtual color emitted(double u, double v, const point3 &p) const {
    return color(0, 0, 0);
  }
  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       color &attenuation, ray &scattered) const = 0;
};

class lambertian : public material {
public:
  lambertian(const color &a) : albedo(std::make_shared<solid_color>(a)) {}
  lambertian(std::shared_ptr<texture> a) : albedo(a) {}
  virtual ~lambertian() = default;

  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       color &attenuation, ray &scattered) const override;

private:
  std::shared_ptr<texture> albedo;
};

class metal : public material {
public:
  metal(const color &a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
  virtual ~metal() = default;

  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       color &attenuation, ray &scattered) const override;

private:
  color albedo;
  double fuzz;
};

class dielectric : public material {
public:
  dielectric(double index_of_refraction) : ir(index_of_refraction) {}
  virtual ~dielectric() = default;

  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       color &attenuation, ray &scattered) const override;

private:
  double ir;

  static double reflectance(double cosine, double ref_idx) {
    // Use Schlick's approximation for reflectance
    auto r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * std::pow((1 - cosine), 5);
  }
};

class diffuse_light : public material {
public:
  diffuse_light(std::shared_ptr<texture> a) : emit(a) {}
  diffuse_light(color c) : emit(std::make_shared<solid_color>(c)) {}

  virtual bool scatter(const ray &ray_in, const hit_record &rec,
                       color &attenuation, ray &scattered) const override {
    return false;
  }

  virtual color emitted(double u, double v, const point3 &p) const override {
    return emit->value(u, v, p);
  }

public:
  std::shared_ptr<texture> emit;
};

class isotropic : public material {
public:
  isotropic(color c) : albedo(std::make_shared<solid_color>(c)) {}
  isotropic(std::shared_ptr<texture> a) : albedo(a) {}

  virtual bool scatter(const ray &r_in, const hit_record &rec,
                       color &attenuation, ray &scattered) const override;

private:
  std::shared_ptr<texture> albedo;
};
