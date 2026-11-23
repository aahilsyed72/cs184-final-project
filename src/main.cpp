// Schwarzschild black hole renderer — CS 184 final project.
//
// For every pixel: build a camera ray, integrate it as a null geodesic
// (geodesic.h), and color it by what happens:
//   captured (r < 1)  -> black (this is the shadow)
//   escaped  (r > 50) -> look up the sky in the ray's final direction
//   max steps         -> black (ray stuck orbiting near the photon sphere)

#include "disk.h"
#include "geodesic.h"
#include "heatmap.h"
#include "sky.h"
#include "vec3.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// ---------------------------------------------------------------- options
struct Options {
    Vec3 pos{0, -20, 2};
    Vec3 lookat{0, 0, 0};
    Vec3 up{0, 0, 1};
    double fov_deg = 40;  // vertical field of view
    double inclination_deg = -1;  // Unset: keep the original camera position.
    int width = 1280, height = 720;
    int spp = 1;          // samples per pixel (stratified sqrt(spp) x sqrt(spp) grid)
    std::string sky = "stars";  // "stars" | "checker" | path to an equirect image
    std::string out = "renders/out.png";
    double exposure = 1.0;
    int threads = 0;      // 0 = OpenMP default (all cores)
    bool flatcam = false; // skip the static-observer correction (for comparison)
    std::string heatmap_out;          // if set, also write a step-count heatmap here
    double heat_lo = 0, heat_hi = 0;  // heatmap color range (0 = auto from data)
    std::string csv_out;              // if set, append a timing row to this CSV
    TraceParams trace;
    DiskParams disk;
};

static void usage() {
    std::printf(
        "usage: blackhole [options]\n"
        "  --pos x,y,z         camera position          (default 0,-20,2)\n"
        "  --lookat x,y,z      look-at point            (default 0,0,0)\n"
        "  --up x,y,z          up vector                (default 0,0,1)\n"
        "  --inclination deg   angle from +z, 0 face-on, 90 edge-on; aims at origin\n"
        "                      uses --pos distance, places camera in the negative-y half-plane\n"
        "  --fov deg           vertical FOV             (default 40)\n"
        "  --res WxH           resolution               (default 1280x720)\n"
        "  --spp n             samples per pixel        (default 1; use 4 or 9 to antialias)\n"
        "  --sky S             stars | checker | path/to/equirect.(png|jpg|hdr)\n"
        "  --disk MODE         blackbody | stripes (default blackbody)\n"
        "  --exposure e        brightness multiplier    (default 1)\n"
        "  --step k            adaptive step factor, dλ = k*r (default 0.02)\n"
        "  --fixed dl          use a FIXED step dλ = dl instead of adaptive (e.g. 0.03)\n"
        "  --heatmap file.png  also write a per-pixel step-count heatmap (log scale, with legend)\n"
        "  --heat-range lo,hi  fix the heatmap color range (default: auto from data)\n"
        "  --csv file.csv      append timing/step stats as a CSV row\n"
        "  --rmax R            escape radius            (default 50)\n"
        "  --maxsteps N        max integration steps    (default 20000)\n"
        "  --threads N         OpenMP threads (1 = single-threaded)\n"
        "  --flatcam           launch rays in coordinate directions (no local-frame correction)\n"
        "  -o file.png         output path              (default renders/out.png)\n");
}

static Vec3 parse_vec(const char* s) {
    Vec3 v;
    if (std::sscanf(s, "%lf,%lf,%lf", &v.x, &v.y, &v.z) != 3) {
        std::fprintf(stderr, "bad vector '%s' (expected x,y,z)\n", s);
        std::exit(1);
    }
    return v;
}

static Options parse_args(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", a.c_str()); std::exit(1); }
            return argv[++i];
        };
        if (a == "--pos") o.pos = parse_vec(next());
        else if (a == "--lookat") o.lookat = parse_vec(next());
        else if (a == "--up") o.up = parse_vec(next());
        else if (a == "--fov") o.fov_deg = std::atof(next());
        else if (a == "--inclination") {
            const char* value = next();
            char* end = nullptr;
            o.inclination_deg = std::strtod(value, &end);
            if (end == value || *end != '\0' || !std::isfinite(o.inclination_deg) ||
                o.inclination_deg < 0.0 || o.inclination_deg > 90.0) {
                std::fprintf(stderr, "--inclination must be a number from 0 to 90 degrees\n");
                std::exit(1);
            }
        }
        else if (a == "--res") {
            if (std::sscanf(next(), "%dx%d", &o.width, &o.height) != 2) { std::fprintf(stderr, "bad --res\n"); std::exit(1); }
        }
        else if (a == "--spp") o.spp = std::max(1, std::atoi(next()));
        else if (a == "--sky") o.sky = next();
        else if (a == "--disk") {
            const std::string mode = next();
            if (mode == "blackbody") o.disk.shading = DiskShading::Blackbody;
            else if (mode == "stripes") o.disk.shading = DiskShading::Stripes;
            else {
                std::fprintf(stderr, "bad --disk mode '%s' (expected blackbody or stripes)\n", mode.c_str());
                std::exit(1);
            }
        }
        else if (a == "--exposure") o.exposure = std::atof(next());
        else if (a == "--step") o.trace.step_scale = std::atof(next());
        else if (a == "--rmax") o.trace.r_escape = std::atof(next());
        else if (a == "--maxsteps") o.trace.max_steps = std::atoi(next());
        else if (a == "--threads") o.threads = std::atoi(next());
        else if (a == "--flatcam") o.flatcam = true;
        else if (a == "--fixed") { o.trace.adaptive = false; o.trace.fixed_dl = std::atof(next()); }
        else if (a == "--heatmap") o.heatmap_out = next();
        else if (a == "--heat-range") {
            if (std::sscanf(next(), "%lf,%lf", &o.heat_lo, &o.heat_hi) != 2) { std::fprintf(stderr, "bad --heat-range\n"); std::exit(1); }
        }
        else if (a == "--csv") o.csv_out = next();
        else if (a == "-o") o.out = next();
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); usage(); std::exit(1); }
    }
    // Apply after parsing so --pos/--lookat ordering does not change the result.
    if (o.inclination_deg >= 0.0) {
        const double radius = length(o.pos);
        if (!std::isfinite(radius) || radius <= o.trace.r_horizon) {
            std::fprintf(stderr, "--inclination requires a camera distance outside the horizon\n");
            std::exit(1);
        }
        const double angle = o.inclination_deg * M_PI / 180.0;
        o.pos = Vec3(0, -radius * std::sin(angle), radius * std::cos(angle));
        o.lookat = Vec3(0, 0, 0);
        o.up = Vec3(0, 0, 1);
    }
    return o;
}

// ---------------------------------------------------------------- camera
// Pinhole camera: an orthonormal basis (right, cam_up, forward) built from
// pos / lookat / up. Pixel (px, py) maps to a point on an image plane one
// unit in front of the eye.
struct Camera {
    Vec3 pos, forward, right, cam_up;
    double tan_half_fov, aspect;

    Camera(const Options& o) {
        pos = o.pos;
        forward = normalize(o.lookat - o.pos);
        Vec3 up = o.up;
        if (length(cross(forward, up)) < 1e-6) up = Vec3(0, 1, 0);  // looking straight along up
        right = normalize(cross(forward, up));
        cam_up = cross(right, forward);
        tan_half_fov = std::tan(o.fov_deg * M_PI / 360.0);
        aspect = double(o.width) / o.height;
    }

    bool static_observer = true;

    // (sx, sy) are continuous pixel coordinates, (0,0) = top-left corner.
    Vec3 ray_dir(double sx, double sy, int W, int H) const {
        double ndc_x = (2.0 * sx / W - 1.0) * aspect * tan_half_fov;
        double ndc_y = (1.0 - 2.0 * sy / H) * tan_half_fov;
        Vec3 n = normalize(forward + right * ndc_x + cam_up * ndc_y);
        if (!static_observer) return n;
        // The pinhole gives a direction n in the *local frame* of an observer
        // hovering at radius r. Near the hole, local radial lengths are
        // stretched by 1/sqrt(1 - 1/r) compared with the coordinate r, so the
        // matching coordinate direction has its radial part scaled by
        // sqrt(1 - 1/r). With this, the shadow's angular radius comes out to
        // exactly sin(a) = (3*sqrt(3)/2 / r) * sqrt(1 - 1/r)  (Synge 1966).
        double r = length(pos);
        Vec3 rhat = pos / r;
        double nr = dot(n, rhat);
        Vec3 nt = n - rhat * nr;
        return normalize(rhat * (nr * std::sqrt(1.0 - 1.0 / r)) + nt);
    }
};

// ---------------------------------------------------------------- color output
static unsigned char to_srgb8(float lin) {
    lin = std::max(0.0f, lin);
    float s = lin <= 0.0031308f ? 12.92f * lin : 1.055f * std::pow(lin, 1.0f / 2.4f) - 0.055f;
    return (unsigned char)std::lround(std::clamp(s, 0.0f, 1.0f) * 255.0f);
}

static bool load_equirect(const std::string& path, EquirectMap& m) {
    int n;
    // stbi_loadf returns linear floats: LDR files (png/jpg) get their sRGB
    // gamma (2.2) removed automatically, and .hdr files are already linear.
    float* data = stbi_loadf(path.c_str(), &m.w, &m.h, &n, 3);
    if (!data) return false;
    m.rgb.assign(data, data + size_t(m.w) * m.h * 3);
    stbi_image_free(data);
    return true;
}

// ---------------------------------------------------------------- main
int main(int argc, char** argv) {
    Options opt = parse_args(argc, argv);

#ifdef _OPENMP
    if (opt.threads > 0) omp_set_num_threads(opt.threads);
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif

    // Sky
    bool use_checker = (opt.sky == "checker");
    EquirectMap sky;
    if (!use_checker) {
        if (opt.sky == "stars") sky = make_procedural_stars();
        else if (!load_equirect(opt.sky, sky)) {
            std::fprintf(stderr, "could not load sky image '%s'\n", opt.sky.c_str());
            return 1;
        }
    }
    auto sample_sky = [&](const Vec3& d) { return use_checker ? checker_sky(d) : sky.sample(d); };

    Camera cam(opt);
    cam.static_observer = !opt.flatcam;
    const int W = opt.width, H = opt.height;
    const int grid = std::max(1, (int)std::lround(std::sqrt((double)opt.spp)));
    const int spp = grid * grid;
    std::vector<Color> img(size_t(W) * H);
    std::vector<float> steps_px(size_t(W) * H);  // average RK4 steps per sample, per pixel

    const DiskParams& disk = opt.disk;
    long n_captured = 0, n_escaped = 0, n_maxsteps = 0, n_disk = 0;
    long long total_steps = 0;
    auto t0 = std::chrono::steady_clock::now();

    // Rows are independent, so we parallelize over them. schedule(dynamic)
    // matters: rows through the shadow and photon ring take many more steps
    // than rows of empty sky, so a static split leaves threads idle.
    // The counters are combined with an OpenMP reduction: each thread keeps a
    // private copy and they are summed at the end (no locks or atomics).
#pragma omp parallel for schedule(dynamic, 1) \
    reduction(+ : n_captured, n_escaped, n_maxsteps, n_disk, total_steps)
    for (int py = 0; py < H; ++py) {
        for (int px = 0; px < W; ++px) {
            Color acc;
            long long pixel_steps = 0;
            for (int sy = 0; sy < grid; ++sy)
                for (int sx = 0; sx < grid; ++sx) {
                    double fx = px + (sx + 0.5) / grid, fy = py + (sy + 0.5) / grid;
                    RayState s{cam.pos, cam.ray_dir(fx, fy, W, H)};
                    RayState end;
                    bool disk_hit = false;
                    Vec3 hit_position;
                    int nsteps = 0;
                    RayResult res = trace_ray(
                        s, opt.trace,
                        [&](const RayState& prev, const RayState& cur) {
                            disk_hit = intersect_disk(prev, cur, hit_position, disk);
                            return disk_hit;  // Opaque disk: stop at the first valid hit.
                        }, &end,
                        &nsteps);
                    pixel_steps += nsteps;
                    // Callback termination returns MaxSteps, so check the disk flag first.
                    if (disk_hit) {
                        acc += shade_disk(hit_position, disk);
                        ++n_disk;
                    } else if (res == RayResult::Escaped) {
                        acc += sample_sky(normalize(end.v));
                        ++n_escaped;
                    } else if (res == RayResult::Captured) {
                        ++n_captured;  // black
                    } else {
                        ++n_maxsteps;  // black
                    }
                }
            size_t idx = size_t(py) * W + px;
            img[idx] = acc * float(opt.exposure / spp);
            steps_px[idx] = float(double(pixel_steps) / spp);
            total_steps += pixel_steps;
        }
    }

    // Timing covers only the ray-tracing loop: not sky generation/loading or PNG writing.
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::vector<unsigned char> out(size_t(W) * H * 3);
    for (size_t i = 0; i < img.size(); ++i) {
        out[3 * i + 0] = to_srgb8(img[i].r);
        out[3 * i + 1] = to_srgb8(img[i].g);
        out[3 * i + 2] = to_srgb8(img[i].b);
    }
    if (!stbi_write_png(opt.out.c_str(), W, H, 3, out.data(), W * 3)) {
        std::fprintf(stderr, "failed to write %s\n", opt.out.c_str());
        return 1;
    }

    const std::string mode = opt.trace.adaptive ? "adaptive" : "fixed";
    const double step_param = opt.trace.adaptive ? opt.trace.step_scale : opt.trace.fixed_dl;
    long total = n_captured + n_escaped + n_maxsteps + n_disk;
    double steps_per_sec = total_steps / secs;
    double mean_steps = double(total_steps) / total;
    float max_px_steps = *std::max_element(steps_px.begin(), steps_px.end());

    if (!opt.heatmap_out.empty()) {
        char title[64];
        std::snprintf(title, sizeof title, "%s %g", opt.trace.adaptive ? "ADAPT" : "FIXED", step_param);
        int hw;
        std::vector<uint8_t> hm = heatmap::render(steps_px, W, H, opt.heat_lo, opt.heat_hi, title, hw);
        if (!stbi_write_png(opt.heatmap_out.c_str(), hw, H, 3, hm.data(), hw * 3)) {
            std::fprintf(stderr, "failed to write %s\n", opt.heatmap_out.c_str());
            return 1;
        }
        std::printf("wrote heatmap %s\n", opt.heatmap_out.c_str());
    }

    if (!opt.csv_out.empty()) {
        std::FILE* probe = std::fopen(opt.csv_out.c_str(), "r");
        bool need_header = true;
        if (probe) { need_header = std::fgetc(probe) == EOF; std::fclose(probe); }
        std::FILE* f = std::fopen(opt.csv_out.c_str(), "a");
        if (!f) { std::fprintf(stderr, "cannot open %s\n", opt.csv_out.c_str()); return 1; }
        if (need_header)
            std::fprintf(f, "width,height,spp,threads,step_mode,step_param,render_s,total_steps,"
                            "steps_per_s,mean_steps_per_ray,max_steps_pixel,captured_pct,maxstep_rays\n");
        std::fprintf(f, "%d,%d,%d,%d,%s,%g,%.4f,%lld,%.4e,%.1f,%.0f,%.3f,%ld\n", W, H, spp, nthreads,
                     mode.c_str(), step_param, secs, total_steps, steps_per_sec, mean_steps,
                     max_px_steps, 100.0 * n_captured / total, n_maxsteps);
        std::fclose(f);
    }

    std::printf("wrote %s  (%dx%d, %d spp, %d threads, %s step %g)  %.3f s\n", opt.out.c_str(), W, H,
                spp, nthreads, mode.c_str(), step_param, secs);
    std::printf("rays: %ld captured (%.2f%%), %ld escaped, %ld hit disk, %ld hit max steps\n", n_captured,
                100.0 * n_captured / total, n_escaped, n_disk, n_maxsteps);
    std::printf("steps: %lld total, %.1f per ray, %.3g steps/s\n", total_steps, mean_steps,
                steps_per_sec);
    return 0;
}
