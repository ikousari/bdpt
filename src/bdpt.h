#ifndef BDPT_H
#define BDPT_H
#include "rtweekend.h"
#include "hittable.h"
#include "material.h"
#include <vector>
#include <random>

class camera;

namespace bdpt {

struct Vertex {
    Vertex(){}
    hit_record rec;

    vec3 wi;
    vec3 wo;
    bool isDelta = false;
    // pdf of sampling the *next* direction from this vertex, measured per solid angle
    double pdfFwd = 0.0;
    // pdf of sampling the *previous* direction (reverse) per solid angle
    double pdfRev = 0.0;
    // Cumulative throughput up to (and excluding) this vertex. Follows typical path tracer convention.
    color throughput = color(1,1,1);
};

struct Path {
    std::vector<Vertex> vertices;
    color throughput = color(1,1,1);
};

Path generate_light_subpath(const hittable& world, const hittable& lights,  int maxDepth);
Path generate_camera_subpath(const camera& cam, const hittable& world, int maxDepth, int px, int py, int sx, int sy);

color ConnectPaths(const hittable& world, const hittable& lights, Path cameraPath, Path lightPath);
color EvaluatePath(const hittable& scene, const hittable& lights, Path& cameraPath, Path& lightPath, int s, int t);
color ConnectVertices(const hittable& scene, Path& cameraPath, Path& lightPath, int s, int t);
double MISWeight(const hittable& scene, Path& cameraPath, Path& lightPath, int s, int t);
double PathPdf(const hittable& scene, Path& cameraPath, Path& lightPath, int s, int t);
color EvaluateDirectLightHit(Path& cameraPath);
color EvaluateDirectIllumination(const hittable& scene, Path& cameraPath, Path& lightPath);
color ConnectToLight(const hittable& scene, const hittable& lights, Path& cameraPath, int t);

// // Top-level pixel sample using BDPT (one light-subpath per call, simple interface for drop-in)
// color bdpt_sample_pixel(const camera& cam, const hittable& world, int image_w, int image_h, int px, int py, int maxDepth);

} // namespace bdpt

#endif /* BDPT_H */