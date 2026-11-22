#pragma once
// Background "celestial sphere" lookups: given the direction an escaped ray
// ends up traveling, return a linear-RGB color.

#include "vec3.h"
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

struct Color {
    float r = 0, g = 0, b = 0;
    Color operator+(const Color& o) const { return {r + o.r, g + o.g, b + o.b}; }
    Color operator*(float s) const { return {r * s, g * s, b * s}; }
    Color& operator+=(const Color& o) { r += o.r; g += o.g; b += o.b; return *this; }
};

// Direction -> equirectangular (u, v) in [0,1)^2. World is z-up:
//   u = longitude (atan2(y, x)), v = polar angle from +z.
inline void dir_to_equirect(const Vec3& d, double& u, double& v) {
    const double PI = 3.14159265358979323846;
    u = (std::atan2(d.y, d.x) + PI) / (2 * PI);
    v = std::acos(std::clamp(d.z, -1.0, 1.0)) / PI;
}

// An equirectangular image in linear RGB, sampled bilinearly
// (with longitude wrapping around).
struct EquirectMap {
    int w = 0, h = 0;
    std::vector<float> rgb;  // linear, 3 floats per pixel

    Color texel(int x, int y) const {
        x = ((x % w) + w) % w;
        y = std::clamp(y, 0, h - 1);
        const float* p = &rgb[3 * (size_t(y) * w + x)];
        return {p[0], p[1], p[2]};
    }
    Color sample(const Vec3& dir) const {
        double u, v;
        dir_to_equirect(dir, u, v);
        double fx = u * w - 0.5, fy = v * h - 0.5;
        int x0 = (int)std::floor(fx), y0 = (int)std::floor(fy);
        float tx = float(fx - x0), ty = float(fy - y0);
        Color a = texel(x0, y0) * (1 - tx) + texel(x0 + 1, y0) * tx;
        Color b = texel(x0, y0 + 1) * (1 - tx) + texel(x0 + 1, y0 + 1) * tx;
        return a * (1 - ty) + b * ty;
    }
};

// Procedural starfield used when no star map is given: random points on the
// sphere splatted into an equirect texture, so it goes through the same
// sampling path as a real star map.
inline EquirectMap make_procedural_stars(int w = 4096, int h = 2048, int n = 30000,
                                         unsigned seed = 184) {
    EquirectMap m;
    m.w = w; m.h = h;
    m.rgb.assign(size_t(w) * h * 3, 0.0f);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int i = 0; i < n; ++i) {
        // Uniform on the sphere: z uniform in [-1,1], longitude uniform.
        double z = 2 * U(rng) - 1, phi = 2 * M_PI * U(rng);
        double s = std::sqrt(1 - z * z);
        Vec3 d(s * std::cos(phi), s * std::sin(phi), z);
        double u, v;
        dir_to_equirect(d, u, v);
        int x = std::min(int(u * w), w - 1), y = std::min(int(v * h), h - 1);
        // Mostly faint stars, a few bright ones (rough power law), slight color tint.
        float bright = float(0.15 * std::pow(U(rng), -1.3));
        bright = std::min(bright, 8.0f);
        float t = float(U(rng));
        Color c = {bright * (0.8f + 0.4f * t), bright * 0.9f, bright * (1.2f - 0.4f * t)};
        float* p = &m.rgb[3 * (size_t(y) * w + x)];
        p[0] += c.r; p[1] += c.g; p[2] += c.b;
    }
    return m;
}

// Checkerboard in (longitude, latitude), for checking that the lensing looks
// right. Thin colored lines mark the equator (red) and longitude 0 (green),
// so you can follow how the sky gets warped and flipped around the hole.
inline Color checker_sky(const Vec3& dir) {
    double u, v;
    dir_to_equirect(dir, u, v);
    const int NU = 36, NV = 18;  // 10-degree cells
    int iu = int(std::floor(u * NU)), iv = int(std::floor(v * NV));
    bool odd = (iu + iv) & 1;
    Color c = odd ? Color{0.9f, 0.9f, 0.9f} : Color{0.08f, 0.08f, 0.12f};
    if (std::fabs(v - 0.5) < 0.002) c = {0.9f, 0.1f, 0.1f};                       // equator
    if (std::fabs(u - 0.5) < 0.002 || u < 0.001 || u > 0.999) c = {0.1f, 0.8f, 0.2f};  // lon 0 / 180
    return c;
}
