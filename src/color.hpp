#pragma once

#include "hawk.h"

#include <iostream>

void write_color(std::ostream &out, color pixel_color) {
  // Write the translated [0,255] value of each color component.
  out << static_cast<int>(255.999 * pixel_color.x()) << ' ' //
      << static_cast<int>(255.999 * pixel_color.y()) << ' ' //
      << static_cast<int>(255.999 * pixel_color.z()) << '\n';
}

void write_color(std::ostream &out, color pixel_color, int samples_per_pixel) {
  auto r = pixel_color.x();
  auto g = pixel_color.y();
  auto b = pixel_color.z();

  // Divide the color by number of samples and gamma-correct for gamma=2.0
  auto scale = 1.0 / samples_per_pixel;
  r = std::sqrt(scale * r);
  g = std::sqrt(scale * g);
  b = std::sqrt(scale * b);

  out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' ' //
      << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' ' //
      << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}
