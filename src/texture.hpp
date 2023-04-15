#pragma once

#include "hawk.h"

#include "perlin.hpp"

class texture {
public:
  virtual ~texture() = default;
  virtual color value(double u, double v, const point3 &p) const = 0;
};

class solid_color : public texture {
public:
  solid_color() {}
  solid_color(color c) : color_value(c) {}

  solid_color(double red, double green, double blue)
      : solid_color(color(red, green, blue)) {}

  virtual color value(double u, double v, const point3 &p) const override {
    return color_value;
  }

private:
  color color_value;
};

class checker_texture : public texture {
public:
  checker_texture() {}

  checker_texture(std::shared_ptr<texture> even, std::shared_ptr<texture> odd)
      : even(even), odd(odd) {}

  checker_texture(color c1, color c2)
      : even(std::make_shared<solid_color>(c1)),
        odd(std::make_shared<solid_color>(c2)) {}

  virtual color value(double u, double v, const point3 &p) const override {
    auto sines =
        std::sin(10 * p.x()) * std::sin(10 * p.y()) * std::sin(10 * p.z());
    if (sines < 0) {
      return odd->value(u, v, p);
    } else {
      return even->value(u, v, p);
    }
  }

private:
  std::shared_ptr<texture> even;
  std::shared_ptr<texture> odd;
};

class noise_texture : public texture {
public:
  noise_texture() {}
  noise_texture(double sc) : scale(sc) {}

  virtual color value(double u, double v, const point3 &p) const override {
    // return color(1, 1, 1) * noise.turb(scale * p);
    return color(1, 1, 1) * 0.5 *
           (1 + std::sin(scale * p.z() + 10 * noise.turb(p)));
  }

private:
  perlin noise;
  double scale;
};