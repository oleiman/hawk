#include "hawk.h"

#include "aarect.hpp"
#include "box.hpp"
#include "bvh.hpp"
#include "camera.hpp"
#include "color.hpp"
#include "hittable_list.hpp"
#include "material.hpp"
#include "moving_sphere.hpp"
#include "sphere.hpp"
#include "texture.hpp"

#include <iostream>

color ray_color(const ray &r, const color &background, const hittable &world,
                int depth) {
  hit_record rec;
  if (depth <= 0) {
    return color(0, 0, 0);
  }

  if (!world.hit(r, 0.001, infinity, rec)) {
    return background;
  }

  ray scattered;
  color attenuation;
  color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

  if (rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
    return emitted +
           attenuation * ray_color(scattered, background, world, depth - 1);
  } else {
    return emitted;
  }

  // vec3 unit_direction = unit_vector(r.direction());
  // auto t = 0.5 * (unit_direction.y() + 1.0);
  // return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
}

hittable_list random_scene() {
  hittable_list world;

  auto checker = std::make_shared<checker_texture>(color(0.2, 0.3, 0.1),
                                                   color(0.9, 0.9, 0.9));
  // auto checker =
  //     std::make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0, 0,
  //     0));
  auto ground_material = std::make_shared<lambertian>(checker);
  world.add(
      std::make_shared<sphere>(point3(0, -1000, 0), 1000, ground_material));

  for (int a = -11; a < 11; ++a) {
    for (int b = -11; b < 11; ++b) {
      auto choose_mat = random_double();
      point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

      if ((center - point3(4, 0.2, 0)).length() > 0.9) {
        std::shared_ptr<material> sphere_material;

        if (choose_mat < 0.8) {
          // diffuse
          auto albedo = color::random() * color::random();
          sphere_material = std::make_shared<lambertian>(albedo);
          auto center2 = center + vec3(0, random_double(0, 0.5), 0);
          world.add(std::make_shared<moving_sphere>(center, center2, 0.0, 1.0,
                                                    0.2, sphere_material));
        } else if (choose_mat < 0.95) {
          // metal
          auto albedo = color::random(0.5, 1);
          auto fuzz = random_double(0, 0.5);
          sphere_material = std::make_shared<metal>(albedo, fuzz);
          world.add(std::make_shared<sphere>(center, 0.2, sphere_material));
        } else {
          // glass
          sphere_material = std::make_shared<dielectric>(1.5);
          world.add(std::make_shared<sphere>(center, 0.2, sphere_material));
        }
      }
    }
  }

  auto material1 = std::make_shared<dielectric>(1.5);
  world.add(std::make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

  auto material2 = std::make_shared<lambertian>(color(0.4, 0.2, 0.1));
  world.add(std::make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

  auto material3 = std::make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
  world.add(std::make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

  return world;
}

hittable_list two_spheres() {
  hittable_list objects;

  auto checker = std::make_shared<checker_texture>(color(0.2, 0.3, 0.1),
                                                   color(0.9, 0.9, 0.9));

  objects.add(std::make_shared<sphere>(point3(0, -10, 0), 10,
                                       std::make_shared<lambertian>(checker)));
  objects.add(std::make_shared<sphere>(point3(0, 10, 0), 10,
                                       std::make_shared<lambertian>(checker)));
  return objects;
}

hittable_list two_perlin_spheres() {
  hittable_list objects;

  auto pertext = std::make_shared<noise_texture>(4);
  objects.add(std::make_shared<sphere>(point3(0, -1000, 0), 1000,
                                       std::make_shared<lambertian>(pertext)));
  objects.add(std::make_shared<sphere>(point3(0, 2, 0), 2,
                                       std::make_shared<lambertian>(pertext)));

  return objects;
}

hittable_list earth() {
  auto earth_texture = std::make_shared<image_texture>("../bin/earthmap.jpg");
  auto earth_surface = std::make_shared<lambertian>(earth_texture);
  auto globe = std::make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);

  return hittable_list(globe);
}

hittable_list simple_light() {
  hittable_list objects;

  auto pertext = std::make_shared<noise_texture>(4);
  objects.add(std::make_shared<sphere>(point3(0, -1000, 0), 1000,
                                       std::make_shared<lambertian>(pertext)));
  objects.add(std::make_shared<sphere>(point3(0, 2, 0), 2,
                                       std::make_shared<dielectric>(2)));

  auto difflight = std::make_shared<diffuse_light>(color(4, 4, 4));
  objects.add(std::make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));

  auto soft_light = std::make_shared<diffuse_light>(color(2, 2, 2));
  objects.add(std::make_shared<sphere>(point3(0, 2, 0), 0.5, difflight));

  return objects;
}

hittable_list cornell_box() {
  hittable_list objects;

  auto red = std::make_shared<lambertian>(color(0.65, 0.05, 0.05));
  auto white = std::make_shared<lambertian>(color(0.73, 0.73, 0.73));
  auto green = std::make_shared<lambertian>(color(0.12, 0.45, 0.15));
  auto light = std::make_shared<diffuse_light>(color(15, 15, 15));

  objects.add(std::make_shared<yz_rect>(0, 555, 0, 555, 555, green));
  objects.add(std::make_shared<yz_rect>(0, 555, 0, 555, 0, red));
  objects.add(std::make_shared<xz_rect>(213, 343, 227, 332, 554, light));
  objects.add(std::make_shared<xz_rect>(0, 555, 0, 555, 0, white));
  objects.add(std::make_shared<xz_rect>(0, 555, 0, 555, 555, white));
  objects.add(std::make_shared<xy_rect>(0, 555, 0, 555, 555, white));

  std::shared_ptr<hittable> box1 =
      std::make_shared<box>(point3(0, 0, 0), point3(165, 330, 165), white);
  box1 = std::make_shared<rotate_y>(box1, 15);
  box1 = std::make_shared<translate>(box1, vec3(265, 0, 295));
  objects.add(box1);

  std::shared_ptr<hittable> box2 =
      std::make_shared<box>(point3(0, 0, 0), point3(165, 165, 165), white);
  box2 = std::make_shared<rotate_y>(box2, -18);
  box2 = std::make_shared<translate>(box2, vec3(130, 0, 65));
  objects.add(box2);

  return objects;
}

int main() {

  // Image

  auto aspect_ratio = 16.0 / 9.0;
  int image_width = 400;
  int image_height = static_cast<int>(image_width / aspect_ratio);
  int samples_per_pixel = 400;
  constexpr int max_depth = 50;

  // World
  hittable_list objects;

  point3 lookfrom;
  point3 lookat;
  auto vfov = 40.0;
  auto aperture = 0.0;
  color background(0, 0, 0);

  switch (0) {
  case 1:
    objects = random_scene();
    background = color(0.70, 0.80, 1.00);
    lookfrom = point3(13, 2, 3);
    lookat = point3(0, 0, 0);
    vfov = 20.0;
    aperture = 0.1;
    break;
  case 2:
    objects = two_spheres();
    background = color(0.70, 0.80, 1.00);
    lookfrom = point3(13, 2, 3);
    lookat = point3(0, 0, 0);
    vfov = 20.0;
    break;
  case 3:
    objects = two_perlin_spheres();
    background = color(0.70, 0.80, 1.00);
    lookfrom = point3(13, 2, 3);
    lookat = point3(0, 0, 0);
    vfov = 20.0;
    break;
  case 4:
    objects = earth();
    background = color(0.70, 0.80, 1.00);
    lookfrom = point3(13, 2, 3);
    lookat = point3(0, 0, 0);
    vfov = 20.0;
    break;
  case 5:
    objects = simple_light();
    background = color(0, 0, 0);
    lookfrom = point3(26, 3, 6);
    lookat = point3(0, 2, 0);
    vfov = 20.0;
    break;
  default:
  case 6:
    objects = cornell_box();
    aspect_ratio = 1.0;
    image_width = 600;
    image_height = static_cast<int>(image_width / aspect_ratio);
    samples_per_pixel = 200;
    background = color(0, 0, 0);
    lookfrom = point3(278, 278, -800);
    lookat = point3(278, 278, 0);
    vfov = 40.0;
    break;
  }

  bvh_node world(objects, 0.0, 1.0);
  // auto world = objects;

  // Camera

  vec3 vup(0, 1, 0);
  auto dist_to_focus = 10.0;

  camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus,
             0.0, 1.0);

  // Render

  std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

  for (int j = image_height - 1; j >= 0; --j) {
    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
    for (int i = 0; i < image_width; ++i) {
      color pixel_color(0, 0, 0);
      for (int s = 0; s < samples_per_pixel; ++s) {
        auto u = (i + random_double()) / (image_width - 1);
        auto v = (j + random_double()) / (image_height - 1);
        ray r = cam.get_ray(u, v);
        pixel_color += ray_color(r, background, world, max_depth);
      }
      write_color(std::cout, pixel_color, samples_per_pixel);
    }
  }
  std::cerr << "\nDone.\n";
}
