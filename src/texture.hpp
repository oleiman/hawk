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

class image_texture : public texture {
  constexpr static int B_PER_PIXEL = 3;

public:
  image_texture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
  image_texture(const std::string &filename);
  ~image_texture() { delete data; }

  virtual color value(double u, double v, const vec3 &p) const override {
    // If we have no texture data, then return solid cyan
    if (data == nullptr) {
      return color(0, 1, 1);
    }

    u = clamp(u, 0.0, 1.0);
    v = 1.0 - clamp(v, 0.0, 1.0); // Flip V to image coordinates

    auto i = static_cast<int>(u * width);
    auto j = static_cast<int>(v * height);

    if (i >= width) {
      i = width - 1;
    }

    if (j >= height) {
      j = height - 1;
    }

    const auto color_scale = 1.0 / 255.0;
    auto pixel = data + j * bytes_per_scanline + i * B_PER_PIXEL;

    return color(color_scale * pixel[0], color_scale * pixel[1],
                 color_scale * pixel[2]);
  }

private:
  unsigned char *data;
  int width, height;
  int bytes_per_scanline;
};
