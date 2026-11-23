#pragma once

#include "geodesic.h"
#include "sky.h"  // Color stores linear RGB.
#include <cmath>
#include <algorithm>

enum class DiskShading { Stripes, Blackbody };

// A thin, opaque ring in the z = 0 plane. Distances use r_s = 1.
struct DiskParams {
    double inner_radius = 3.0;
    double outer_radius = 15.0;
    DiskShading shading = DiskShading::Blackbody;
    double peak_temperature = 6500.0;  // Kelvin; chosen for visible warm/white colors.
};

// Approximate one curved integration step by a straight segment and find
// its intersection with the disk. A miss leaves `hit` unchanged, allowing
// the caller to keep tracing and test later plane crossings.
inline bool intersect_disk(const RayState& prev, const RayState& cur,
                           Vec3& hit, const DiskParams& disk = DiskParams{}) {
    const double z0 = prev.x.z;
    const double z1 = cur.x.z;

    // Include a step ending on the plane, but exclude one starting there:
    // this avoids counting the same crossing twice. Coplanar segments have
    // no unique crossing and are ignored by this zero-thickness model.
    const bool crossed = (z0 > 0.0 && z1 <= 0.0) ||
                         (z0 < 0.0 && z1 >= 0.0);
    if (!crossed) return false;

    const double t = z0 / (z0 - z1);
    Vec3 candidate = prev.x + (cur.x - prev.x) * t;
    candidate.z = 0.0;
    const double r = std::hypot(candidate.x, candidate.y);
    if (!(r >= disk.inner_radius && r <= disk.outer_radius)) return false;

    hit = candidate;
    return true;
}

// Diagnostic concentric stripes, each one Schwarzschild radius wide.
// These are deliberately simple colors for checking lensing and disk hits.
inline Color shade_disk_stripes(const Vec3& hit) {
    const double r = std::hypot(hit.x, hit.y);
    const bool bright = std::fmod(std::floor(r), 2.0) >= 1.0;
    return bright ? Color{1.0f, 0.6f, 0.2f}
                  : Color{0.2f, 0.05f, 0.01f};
}

// Thin-disk temperature profile, normalized so its maximum is peak_temperature.
// T^4 is proportional to r^-3 * (1 - sqrt(r_in/r)); its peak is at
// r = (49/36) * r_in. The exact inner edge has zero temperature.
inline double disk_temperature(double r, const DiskParams& disk = DiskParams{}) {
    if (r <= disk.inner_radius || r > disk.outer_radius) return 0.0;
    const double q = r / disk.inner_radius;
    const double peak_q = 49.0 / 36.0;
    const double profile = (1.0 - 1.0 / std::sqrt(q)) / (q * q * q);
    const double peak_profile = (1.0 - 1.0 / std::sqrt(peak_q)) /
                                (peak_q * peak_q * peak_q);
    return disk.peak_temperature * std::pow(profile / peak_profile, 0.25);
}

// Tanner Helland's approximate temperature-to-sRGB fit (1000-40000 K):
// https://tannerhelland.com/2012/09/18/convert-temperature-rgb-algorithm-code.html
// This is an appearance approximation, not spectral radiometry. Convert the
// fitted sRGB channels to linear RGB before the renderer averages samples.
inline Color blackbody_color(double kelvin) {
    if (kelvin <= 0.0) return {};
    const double t = std::clamp(kelvin, 1000.0, 40000.0) / 100.0;
    const double red = t <= 66.0 ? 255.0
        : 329.698727446 * std::pow(t - 60.0, -0.1332047592);
    const double green = t <= 66.0
        ? 99.4708025861 * std::log(t) - 161.1195681661
        : 288.1221695283 * std::pow(t - 60.0, -0.0755148492);
    const double blue = t >= 66.0 ? 255.0 : (t <= 19.0 ? 0.0
        : 138.5177312231 * std::log(t - 10.0) - 305.0447927307);
    auto linear = [](double channel) -> float {
        const double s = std::clamp(channel / 255.0, 0.0, 1.0);
        return float(s <= 0.04045 ? s / 12.92 : std::pow((s + 0.055) / 1.055, 2.4));
    };
    return {linear(red), linear(green), linear(blue)};
}

inline Color shade_disk(const Vec3& hit, const DiskParams& disk = DiskParams{}) {
    if (disk.shading == DiskShading::Stripes) return shade_disk_stripes(hit);
    const double temperature = disk_temperature(std::hypot(hit.x, hit.y), disk);
    if (temperature <= 0.0 || disk.peak_temperature <= 0.0) return {};
    // Relative T^4 emission makes the inner region bright and the outer disk
    // dimmer. Multiplying the approximate RGB hue by this is a display model,
    // not an integration of visible spectral radiance.
    const double relative_temperature = temperature / disk.peak_temperature;
    return blackbody_color(temperature) * float(std::pow(relative_temperature, 4.0));
}
