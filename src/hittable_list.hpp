#pragma once

#include "hittable.hpp"

#include <memory>
#include <vector>

class hittable_list : public hittable {
public:
  hittable_list() = default;
  hittable_list(std::shared_ptr<hittable> object) { add(object); }
  virtual ~hittable_list() = default;

  void clear() { _objects.clear(); }
  void add(std::shared_ptr<hittable> object) { _objects.push_back(object); }
  const std::vector<std::shared_ptr<hittable>> &objects() const {
    return _objects;
  }

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

  virtual bool bounding_box(double time0, double time1,
                            aabb &output_box) const override;

private:
  std::vector<std::shared_ptr<hittable>> _objects;
};
