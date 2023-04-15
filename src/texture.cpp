#include "texture.hpp"

#include "hawk_stb.h"

#include <string>

image_texture::image_texture(const std::string &filename) {
  auto components_per_pixel = B_PER_PIXEL;

  data = stbi_load(filename.c_str(), &width, &height, &components_per_pixel,
                   components_per_pixel);

  if (!data) {
    std::cerr << "ERROR: Could not load texture image file '" << filename
              << "'.\n";
    width = height = 0;
  }

  bytes_per_scanline = B_PER_PIXEL * width;
}
