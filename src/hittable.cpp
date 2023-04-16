#include "aarect.hpp"
#include "box.hpp"
#include "constant_medium.hpp"
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

box::box(const point3 &p0, const point3 &p1, std::shared_ptr<material> ptr)
    : box_min(p0), box_max(p1) {

  // xy sides at zmax and zmin
  sides.add(
      std::make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
  sides.add(
      std::make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));

  // xz sides at ymax and ymin
  sides.add(
      std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
  sides.add(
      std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));

  // yz sides at xmax and xmin
  sides.add(
      std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
  sides.add(
      std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
}

bool box::hit(const ray &r, double t_min, double t_max, hit_record &rec) const {
  return sides.hit(r, t_min, t_max, rec);
}

bool box::bounding_box(double time0, double time1, aabb &output_box) const {
  output_box = aabb(box_min, box_max);
  return true;
}

bool translate::hit(const ray &r, double t_min, double t_max,
                    hit_record &rec) const {
  ray moved_r(r.origin() - offset, r.direction(), r.time());

  if (!ptr->hit(moved_r, t_min, t_max, rec)) {
    return false;
  }

  rec.p += offset;
  rec.set_face_normal(moved_r, rec.normal);

  return true;
}

bool translate::bounding_box(double time0, double time1,
                             aabb &output_box) const {
  if (!ptr->bounding_box(time0, time1, output_box)) {
    return false;
  }
  output_box = aabb(output_box.min() + offset, output_box.max() + offset);
  return true;
}

rotate_y::rotate_y(std::shared_ptr<hittable> p, double angle) : ptr(p) {
  auto radians = degrees_to_radians(angle);
  sin_theta = std::sin(radians);
  cos_theta = std::cos(radians);
  hasbox = ptr->bounding_box(0, 1, bbox);

  point3 min(infinity, infinity, infinity);
  point3 max(-infinity, -infinity, -infinity);

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        auto x = i * bbox.max().x() + (1 - i) * bbox.min().x();
        auto y = j * bbox.max().y() + (1 - j) * bbox.min().y();
        auto z = k * bbox.max().z() + (1 - k) * bbox.min().z();

        auto newx = cos_theta * x + sin_theta * z;
        auto newz = -sin_theta * x + cos_theta * z;

        vec3 tester(newx, y, newz);

        for (int c = 0; c < 3; c++) {
          min[c] = std::fmin(min[c], tester[c]);
          max[c] = std::fmax(max[c], tester[c]);
        }
      }
    }
  }
  bbox = aabb(min, max);
}

bool rotate_y::hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const {
  auto origin = r.origin();
  auto direction = r.direction();

  origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
  origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];

  direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
  direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];

  ray rotated_r(origin, direction, r.time());

  if (!ptr->hit(rotated_r, t_min, t_max, rec)) {
    return false;
  }

  auto p = rec.p;
  auto normal = rec.normal;

  p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
  p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];

  normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
  normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];

  rec.p = p;
  rec.set_face_normal(rotated_r, normal);

  return true;
}

bool rotate_y::bounding_box(double time0, double time1,
                            aabb &output_box) const {
  output_box = bbox;
  return hasbox;
}

bool constant_medium::hit(const ray &r, double t_min, double t_max,
                          hit_record &rec) const {
  // Print occasional samples when debugging as needed
  const bool enableDebug = false;
  const bool debugging = enableDebug && random_double() < 0.00001;

  hit_record rec1, rec2;

  if (!boundary->hit(r, -infinity, infinity, rec1)) {
    return false;
  }

  if (!boundary->hit(r, rec1.t + 0.0001, infinity, rec2)) {
    return false;
  }

  if (debugging) {
    std::cerr << "\nt_min=" << rec1.t << ", t_max=" << rec2.t << "\n";
  }

  if (rec1.t < t_min) {
    rec1.t = t_min;
  }

  if (rec2.t > t_max) {
    rec2.t = t_max;
  }

  if (rec1.t >= rec2.t) {
    return false;
  }

  if (rec1.t < 0) {
    rec1.t = 0;
  }

  const auto ray_length = r.direction().length();
  const auto distance_inside_boundary = (rec2.t - rec1.t) * ray_length;
  const auto hit_distance = neg_inv_density * std::log(random_double());

  if (hit_distance > distance_inside_boundary) {
    return false;
  }

  rec.t = rec1.t + hit_distance / ray_length;
  rec.p = r.at(rec.t);

  if (debugging) {
    std::cerr << "hit_distance = " << hit_distance << "\n" //
              << "rec.t = " << rec.t << "\n"               //
              << "rec.p = " << rec.p << "\n";
  }

  // normal is arbitrary
  rec.normal = vec3(1, 0, 0);
  rec.front_face = true;
  rec.mat_ptr = phase_function;
  return true;
}

bool constant_medium::bounding_box(double time0, double time1,
                                   aabb &output_box) const {
  return boundary->bounding_box(time0, time1, output_box);
}
