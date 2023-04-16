#include "material.hpp"
#include "hittable.hpp"

bool lambertian::scatter(const ray &r_in, const hit_record &rec,
                         color &attenuation, ray &scattered) const {
  auto scatter_direction = rec.normal + random_unit_vector();

  // Catch degenerate scatter direction
  if (scatter_direction.near_zero()) {
    scatter_direction = rec.normal;
  }
  scattered = ray(rec.p, scatter_direction, r_in.time());
  attenuation = albedo->value(rec.u, rec.v, rec.p);
  return true;
}

bool metal::scatter(const ray &r_in, const hit_record &rec, color &attenuation,
                    ray &scattered) const {
  vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
  scattered =
      ray(rec.p, reflected + fuzz * random_in_unit_sphere(), r_in.time());
  attenuation = albedo;
  return (dot(scattered.direction(), rec.normal) > 0);
}

bool dielectric::scatter(const ray &r_in, const hit_record &rec,
                         color &attenuation, ray &scattered) const {
  attenuation = color(1.0, 1.0, 1.0);
  double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;
  vec3 unit_direction = unit_vector(r_in.direction());

  double cos_theta = std::fmin(dot(-unit_direction, rec.normal), 1.0);
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

  bool cannot_refract = refraction_ratio * sin_theta > 1.0;
  vec3 direction;

  if (cannot_refract ||
      reflectance(cos_theta, refraction_ratio) > random_double()) {
    direction = reflect(unit_direction, rec.normal);
  } else {
    direction = refract(unit_direction, rec.normal, refraction_ratio);
  }

  scattered = ray(rec.p, direction, r_in.time());
  return true;
}

bool isotropic::scatter(const ray &r_in, const hit_record &rec,
                        color &attenuation, ray &scattered) const {
  scattered = ray(rec.p, random_in_unit_sphere(), r_in.time());
  attenuation = albedo->value(rec.u, rec.v, rec.p);
  return true;
}
