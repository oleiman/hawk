#pragma once

#include "hittable.hpp"

#include <memory>
#include <vector>

class hittable_list : public hittable {
public:
  hittable_list() = default;
  hittable_list(std::shared_ptr<hittable> object) { add(object); }
  virtual ~hittable_list() = default;

  void clear() { objects.clear(); }
  void add(std::shared_ptr<hittable> object) { objects.push_back(object); }

  virtual bool hit(const ray &r, double t_min, double t_max,
                   hit_record &rec) const override;

private:
  std::vector<std::shared_ptr<hittable>> objects;
};
