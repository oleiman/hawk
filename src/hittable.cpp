#include "aarect.hpp"
#include "hittable_list.hpp"
#include "moving_sphere.hpp"
#include "sphere.hpp"

aabb surrounding_box(const aabb &box0, const aabb &box1) {
  point3 small(std::fmin(box0.min().x(), box1.min().x()),
               std::fmin(box0.min().y(), box1.min().y()),
               std::fmin(box0.min().z(), box1.min().z()));
  point3 big(std::fmax(box0.max().x(), box1.max().x()),
             std::fmax(box0.max().y(), box1.max().y()),
             std::fmax(box0.max().z(), box1.max().z()));
  return aabb(small, big);
}

bool hittable_list::hit(const ray &r, double t_min, double t_max,
                        hit_record &rec) const {
  hit_record temp_rec;
  bool hit_anything = false;
  auto closest_so_far = t_max;

  for (const auto &object : objects()) {
    if (object->hit(r, t_min, closest_so_far, temp_rec)) {
      hit_anything = true;
      closest_so_far = temp_rec.t;
      rec = temp_rec;
    }
  }
  return hit_anything;
}

bool hittable_list::bounding_box(double time0, double time1,
                                 aabb &output_box) const {
  if (objects().empty())
    return false;

  aabb temp_box;
  bool first_box = true;

  for (const auto &object : objects()) {
    if (!object->bounding_box(time0, time1, temp_box)) {
      return false;
    }
    output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
    first_box = false;
  }

  return true;
}

bool sphere::hit(const ray &r, double t_min, double t_max,
                 hit_record &rec) const {
  vec3 oc = r.origin() - center;
  auto a = r.direction().length_squared();
  auto half_b = dot(oc, r.direction());
  auto c = oc.length_squared() - radius * radius;

  auto discriminant = half_b * half_b - a * c;
  if (discriminant < 0)
    return false;

  auto sqrtd = std::sqrt(discriminant);

  auto root = (-half_b - sqrtd) / a;
  if (root < t_min || t_max < root) {
    root = (-half_b + sqrtd) / a;
    if (root < t_min || t_max < root) {
      return false;
    }
  }

  rec.t = root;
  rec.p = r.at(rec.t);

  vec3 outward_normal = (rec.p - center) / radius;
  rec.set_face_normal(r, outward_normal);

  // treating the unit normal vector as though it were a point on a unit sphere
  // centered at the origin allows us to calculate the spherical coordinates of
  // that point and use those to index into a texture.
  sphere::get_sphere_uv(outward_normal, rec.u, rec.v);
  rec.mat_ptr = mat_ptr;

  return true;
}

bool sphere::bounding_box(double time0, double time1, aabb &output_box) const {
  output_box = aabb(center - vec3(radius, radius, radius),
                    center + vec3(radius, radius, radius));
  return true;
}

// TODO(oren): remove duplicated code (only diff from sphere::hit is how center
// is obtained)
bool moving_sphere::hit(const ray &r, double t_min, double t_max,
                        hit_record &rec) const {
  vec3 oc = r.origin() - center(r.time());
  auto a = r.direction().length_squared();
  auto half_b = dot(oc, r.direction());
  auto c = oc.length_squared() - radius * radius;

  auto discriminant = half_b * half_b - a * c;
  if (discriminant < 0)
    return false;

  auto sqrtd = std::sqrt(discriminant);

  auto root = (-half_b - sqrtd) / a;
  if (root < t_min || t_max < root) {
    root = (-half_b + sqrtd) / a;
    if (root < t_min || t_max < root) {
      return false;
    }
  }

  rec.t = root;
  rec.p = r.at(rec.t);

  vec3 outward_normal = (rec.p - center(r.time())) / radius;
  rec.set_face_normal(r, outward_normal);
  rec.mat_ptr = mat_ptr;

  return true;
}

bool moving_sphere::bounding_box(double time0, double time1,
                                 aabb &output_box) const {
  aabb box0(center(time0) - vec3(radius, radius, radius),
            center(time0) + vec3(radius, radius, radius));
  aabb box1(center(time1) - vec3(radius, radius, radius),
            center(time1) + vec3(radius, radius, radius));
  output_box = surrounding_box(box0, box1);
  return true;
}

bool xy_rect::hit(const ray &r, double t_min, double t_max,
                  hit_record &rec) const {

  // check whether the ray intersects z = k
  auto t = (k - r.origin().z()) / r.direction().z();
  if (t < t_min || t > t_max) {
    return false;
  }

  auto x = r.origin().x() + t * r.direction().x();
  auto y = r.origin().y() + t * r.direction().y();

  // does the ray intersect z = k inside the rectangle?
  if (x < x0 || x > x1 || y < y0 || y > y1) {
    return false;
  }

  // calculate texture coords (fractional index into the rectangle)
  rec.u = (x - x0) / (x1 - x0);
  rec.v = (y - y0) / (y1 - y0);
  rec.t = t;

  // normal is in z direction because we're axis alligned
  auto outward_normal = vec3(0, 0, 1);

  rec.set_face_normal(r, outward_normal);
  rec.mat_ptr = mp;
  rec.p = r.at(t);
  return true;
}

bool xy_rect::bounding_box(double time0, double time1, aabb &output_box) const {
  // bounding box must have non-zero width in each dimension, so pad the Z
  // dimension a small amount (rectangle in XY plane will have an infinitely
  // thin depth).
  output_box = aabb(point3(x0, y0, k - 0.0001), point3(x1, y1, k + 0.0001));
  return true;
}

bool xz_rect::hit(const ray &r, double t_min, double t_max,
                  hit_record &rec) const {

  // check whether the ray intersects z = k
  auto t = (k - r.origin().y()) / r.direction().y();
  if (t < t_min || t > t_max) {
    return false;
  }

  auto x = r.origin().x() + t * r.direction().x();
  auto z = r.origin().z() + t * r.direction().z();

  // does the ray intersect z = k inside the rectangle?
  if (x < x0 || x > x1 || z < z0 || z > z1) {
    return false;
  }

  // calculate texture coords (fractional index into the rectangle)
  rec.u = (x - x0) / (x1 - x0);
  rec.v = (z - z0) / (z1 - z0);
  rec.t = t;

  // normal is in y direction because we're axis alligned
  auto outward_normal = vec3(0, 1, 0);

  rec.set_face_normal(r, outward_normal);
  rec.mat_ptr = mp;
  rec.p = r.at(t);
  return true;
}

bool xz_rect::bounding_box(double time0, double time1, aabb &output_box) const {
  // bounding box must have non-zero width in each dimension, so pad the Y
  // dimension a small amount (rectangle in XZ plane will have an infinitely
  // thin depth).
  output_box = aabb(point3(x0, k - 0.0001, z0), point3(x1, k + 0.0001, z1));
  return true;
}

bool yz_rect::hit(const ray &r, double t_min, double t_max,
                  hit_record &rec) const {

  // check whether the ray intersects z = k
  auto t = (k - r.origin().x()) / r.direction().x();
  if (t < t_min || t > t_max) {
    return false;
  }

  auto y = r.origin().y() + t * r.direction().y();
  auto z = r.origin().z() + t * r.direction().z();

  // does the ray intersect z = k inside the rectangle?
  if (y < y0 || y > y1 || z < z0 || z > z1) {
    return false;
  }

  // calculate texture coords (fractional index into the rectangle)
  rec.u = (y - y0) / (y1 - y0);
  rec.v = (z - z0) / (z1 - z0);
  rec.t = t;

  // normal is in y direction because we're axis alligned
  auto outward_normal = vec3(1, 0, 0);

  rec.set_face_normal(r, outward_normal);
  rec.mat_ptr = mp;
  rec.p = r.at(t);
  return true;
}

bool yz_rect::bounding_box(double time0, double time1, aabb &output_box) const {
  // bounding box must have non-zero width in each dimension, so pad the Y
  // dimension a small amount (rectangle in YZ plane will have an infinitely
  // thin depth).
  output_box = aabb(point3(k - 0.0001, y0, z0), point3(k + 0.0001, y1, z1));
  return true;
}
