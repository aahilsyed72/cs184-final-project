#pragma once
// Null-geodesic integrator for a Schwarzschild black hole, units r_s = 1.
//
// Trick: for a photon, the Schwarzschild orbit equation is
//     d²u/dφ² + u = (3/2) r_s u²,    u = 1/r
// The Newtonian-looking 3D ODE
//     d²x/dλ² = -1.5 h² x / |x|^5,   h = |x × v|
// is a central force, so it keeps angular momentum h and stays in a plane,
// and its orbit equation in that plane is exactly the one above. So it traces
// the same spatial paths as real light rays, just with a non-affine
// parameter λ. That's all we need: only *where* a ray goes matters for
// rendering, not *when* it gets there.

#include "vec3.h"
#include <algorithm>

struct TraceParams {
    bool adaptive = true;      // true: dλ = step_scale * r;  false: dλ = fixed_dl everywhere
    double fixed_dl = 0.03;    // fixed-step size (0.03 = adaptive's step at the photon sphere, r = 1.5)
    double step_scale = 0.02;  // dλ = step_scale * r   (smaller steps near the hole)
    double dl_min = 1e-4;      // clamp so we never take a zero-size step
    double dl_max = 1.0;       // clamp so we never jump over the hole from far away
    double r_horizon = 1.0;    // r < 1  -> captured (event horizon, r_s = 1)
    double r_escape = 50.0;    // r > 50 and moving outward -> escaped
    int max_steps = 20000;     // rays that orbit near r = 1.5 can loop for a long time
};

enum class RayResult { Captured, Escaped, MaxSteps };

struct RayState {
    Vec3 x;  // position
    Vec3 v;  // direction (dx/dλ)
};

// Acceleration from the geodesic equation above. h2 = h² is constant along a ray.
inline Vec3 geodesic_accel(const Vec3& x, double h2) {
    double r2 = dot(x, x);
    double r5 = r2 * r2 * std::sqrt(r2);
    return x * (-1.5 * h2 / r5);
}

// One classic RK4 step for the 2nd-order system, written as the 1st-order
// system x' = v, v' = a(x).
inline void rk4_step(RayState& s, double h2, double dl) {
    const Vec3 &x = s.x, &v = s.v;
    Vec3 k1x = v;                       Vec3 k1v = geodesic_accel(x, h2);
    Vec3 k2x = v + k1v * (0.5 * dl);    Vec3 k2v = geodesic_accel(x + k1x * (0.5 * dl), h2);
    Vec3 k3x = v + k2v * (0.5 * dl);    Vec3 k3v = geodesic_accel(x + k2x * (0.5 * dl), h2);
    Vec3 k4x = v + k3v * dl;            Vec3 k4v = geodesic_accel(x + k3x * dl, h2);
    s.x += (k1x + 2.0 * k2x + 2.0 * k3x + k4x) * (dl / 6.0);
    s.v += (k1v + 2.0 * k2v + 2.0 * k3v + k4v) * (dl / 6.0);
}

// Integrate one ray until it falls in, escapes, or runs out of steps.
//
// `on_step(prev, cur)` is called after every step with the state before and
// after it. If it returns true, tracing stops (and trace() returns MaxSteps,
// which the caller ignores because it set its own flag). Stage 1 passes a
// no-op; Stage 2 uses it to find accretion-disk plane crossings in between
// steps.
template <class OnStep>
RayResult trace_ray(RayState s, const TraceParams& p, OnStep&& on_step,
                    RayState* final_state = nullptr, int* steps_taken = nullptr) {
    const double h2 = dot(cross(s.x, s.v), cross(s.x, s.v));  // conserved
    RayResult result = RayResult::MaxSteps;
    int i = 0;
    for (; i < p.max_steps; ++i) {
        double r = length(s.x);
        if (r < p.r_horizon) { result = RayResult::Captured; break; }
        if (r > p.r_escape && dot(s.x, s.v) > 0) { result = RayResult::Escaped; break; }

        double dl = p.adaptive ? std::clamp(p.step_scale * r, p.dl_min, p.dl_max) : p.fixed_dl;
        RayState prev = s;
        rk4_step(s, h2, dl);
        if (on_step(prev, s)) break;
    }
    if (final_state) *final_state = s;
    if (steps_taken) *steps_taken = i;
    return result;
}

inline RayResult trace_ray(const RayState& s, const TraceParams& p) {
    return trace_ray(s, p, [](const RayState&, const RayState&) { return false; });
}
