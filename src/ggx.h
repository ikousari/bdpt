#ifndef GGX_H
#define GGX_H
#include "vec3.h"
#include "rtweekend.h"
#include "onb.h"

color schlick_fresnel(double cos_theta, const color& F0)
{
    double m = std::max(0.0, std::min(1.0, 1.0 - cos_theta));
    double m2 = m * m;
    double m5 = m2 * m2 * m;
    
    return F0 + (color(1, 1, 1) - F0) * m5;
}

double GGX_D(const vec3& h, const vec3& n, double alpha)
{
    double a2 = alpha * alpha;
    double cos_h = dot(h, n);
    
    if (cos_h <= 0) return 0;
    
    double cos2_h = cos_h * cos_h;
    double tan2_h = (1.0 - cos2_h) / cos2_h;
    
    double denom = pi * cos2_h * cos2_h * (a2 + tan2_h) * (a2 + tan2_h);
    
    return a2 / (denom + 1e-8);
}

double G1(const vec3& v, const vec3& n, double alpha)
{
    double cos_v = dot(v, n);
    if (cos_v <= 0) return 0;

    double a2 = alpha * alpha;
    double cos2_v = cos_v * cos_v;
    double tan2_v = (1.0 - cos2_v) / cos2_v;

    return 2.0 / (1.0 + sqrt(1.0 + a2 * tan2_v));
}

double GGX_G(const vec3& wi, const vec3& wo, const vec3& n, double alpha)
{
    return G1(wi, n, alpha) * G1(wo, n, alpha);
}



vec3 GGX_VNDF(const vec3& wi, const vec3& n, double alpha, double r1, double r2)
{
    double phi = 2.0 * pi * r1;
    double cos_theta = sqrt((1.0 - r2) / (1.0 + (alpha * alpha - 1.0) * r2));
    double sin_theta = sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
    
    double x = cos(phi) * sin_theta;
    double y = sin(phi) * sin_theta;
    double z = cos_theta;

    onb basis(n);
    return unit_vector(basis.transform(vec3(x,y,z)));
}
#endif