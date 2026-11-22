# Black Hole Renderer — Development Notes

Schwarzschild (non-spinning) black hole, units r_s = 1. C++17 / CMake / OpenMP, PNG via stb_image_write.

Build:
```
cmake -S . -B build && cmake --build build -j
./build/blackhole --help
```

Machine: Apple Silicon (arm64), 8 cores, Apple clang 15 + Homebrew libomp.

---

## Stage 1 — Shadow + lensed starfield

### What was implemented
- `src/geodesic.h` — ray integrator.
  - Equation of motion: `d²x/dλ² = -1.5 h² x / |x|^5`, with `h = |x × v|` computed once per ray.
    This is a central force, so motion stays in a plane and h is conserved. Its Binet orbit
    equation is `u'' + u = 1.5 u²` (u = 1/r), the Schwarzschild photon orbit equation. So it
    traces the correct light paths, just with a non-affine parameter λ.
  - Classic RK4 on the first-order system (x' = v, v' = a(x)).
  - Adaptive step `dλ = 0.02 · r`, clamped to [1e-4, 1]. Steps are small near the hole and large far away.
  - Termination: `r < 1` → captured; `r > 50` *and* moving outward → escaped; 20000 steps → give up (black).
  - Generic `on_step(prev, cur)` callback for detecting disk crossings in Stage 2.
- `src/sky.h` — equirectangular sky lookup (bilinear, wraps in longitude), procedural starfield,
  and a 10° checkerboard with a red equator and green longitude-0/180 lines for debugging.
- `src/main.cpp` — pinhole camera (`--pos --lookat --up --fov --res`), stratified supersampling (`--spp`),
  OpenMP over rows with `schedule(dynamic)`, linear → sRGB output.
- Star maps load through `stbi_loadf`, which returns linear floats. It removes gamma from PNG/JPG and leaves .hdr as is.

### Problems hit and fixes
1. **Shadow ~3% too large.** The first render's shadow was 260 px across at r = 20; theory says 252.5 px.
   The camera launched rays in *coordinate* directions. A real observer hovering at radius r measures
   angles in their local frame, where radial lengths are stretched by 1/√(1−1/r). So
   tan α_coord = tan α_local / √(1−1/r).
   Substituting into the orbit invariant (u')² + u² − u³ = 1/b² gives b = r·sin α_local / √(1−1/r),
   which is Synge's shadow formula.
   **Fix:** scale the radial component of each launch direction by √(1−1/r) (`Camera::ray_dir`).
   `--flatcam` turns the correction off for comparison.

   | camera r | predicted diameter | flat camera | corrected camera |
   |---|---|---|---|
   | 20 | 252.5 px | 260 px | 252 px |
   | 8  | 630.8 px | 674 px | 630 px |

   (1280x720, 40° vertical FOV, camera in the equatorial plane looking at the hole.)
2. Small LSP errors in the editor before the first CMake configure. These came from the language
   server parsing the headers as pre-C++11 without a compile database. `CMAKE_EXPORT_COMPILE_COMMANDS`
   fixes it; there was no real compile error.
3. Apple clang has no OpenMP runtime → installed `libomp` with Homebrew, set `OpenMP_ROOT` in CMake.

### Known approximations
- Escaped rays are colored by their direction at r = 50, not the true asymptotic direction. The
  bending left after r = 50 is of order r_s/r ≈ 0.02 rad (about 1°) at most, and smooth across the image. Raise it with `--rmax 200` if needed.

- **Why stars near the Einstein ring look like streaks but stars by the shadow look like normal dots.**
  Both are lensing, in opposite directions. Near the Einstein ring, images are stretched along the ring
  and barely changed toward the center, so stars become arcs. Just outside the shadow, rays skim the
  photon sphere and swing ≥180°, so the whole sky (including behind the camera) is squeezed into a thin
  band. Stars there shrink to sub-pixel dots, which look like a normal dense starfield.
  The placeholder starfield adds to the streaking: each star is one pixel in the 4096x2048 sky image
  (0.09°), about 3 screen pixels at FOV 40°, so any stretching is easy to see.

### Renders
| file | command | time (8 threads) |
|---|---|---|
| `renders/stage1_checker.png` | `--sky checker --spp 4` | 5.5 s |
| `renders/stage1_stars.png`   | `--sky stars --spp 4`   | 6.5 s |

Default camera: pos (0,−20,2), look-at origin, FOV 40°, 1280x720. 5.4% of rays are captured; none hit the step limit.

---

## Instrumentation — step counts, heatmaps, timing, fixed vs adaptive steps

### What was implemented
- **Fixed-step mode** (`--fixed dl`): `TraceParams::adaptive = false` makes `trace_ray` use a constant
  `dλ = dl` instead of `0.02·r`. The default comparison uses `--fixed 0.03`, which equals the adaptive
  step at the photon sphere (0.02 × 1.5), so both modes resolve the hardest region equally well.
- **Per-pixel step counting:** `trace_ray` already reported steps taken. The render loop now saves the
  average steps per sample for every pixel (`steps_px`) and adds them into `total_steps` with an OpenMP
  `reduction(+:...)`, so threads don't share a counter.
- **Heatmap** (`--heatmap file.png`, `--heat-range lo,hi`): `src/heatmap.h` maps step counts on a
  **log scale** through matplotlib's *inferno* colormap (`src/colormap.h`, 33 sampled entries,
  interpolated) and adds a legend strip: color bar, 1-2-5 ticks per decade, labels, and the step mode.
  Text uses a built-in 5x7 bitmap font, so no font library is needed. Log scale, because counts range
  from ~150 (empty sky) to 20000 (rays stuck orbiting) — on a linear scale everything except the
  photon ring would be one color.
- **CSV timing** (`--csv file.csv`): appends one row per render with resolution, spp, threads, step mode,
  render time, total steps, steps/s, mean steps/ray, max steps in a pixel, captured %, and max-step rays.
  The timer covers only the ray-tracing loop (not sky creation or PNG writing).
- `scripts/bench.sh` runs {640x360, 1280x720} × {adaptive 0.02, fixed 0.03} × {1, 4, 8 threads}, 3 repeats
  each, into `results/timing.csv`. It then renders both heatmaps on the same color range and runs an
  accuracy sweep. `scripts/summarize_bench.py` takes medians and writes `results/bench_summary.md`.

### Results (Apple M3: 4 performance + 4 efficiency cores)

#### Render timing (median of repeats, 1 spp, checker sky, camera at r = 20)

| resolution | step mode | threads | render time (s) | speedup vs 1 thread | total RK4 steps | mean steps/ray | steps/s |
|---|---|---|---|---|---|---|---|
| 640x360 | adaptive 0.02 | 1 | 1.704 | 1.00x | 52,895,924 | 230 | 3.1e+07 |
| 640x360 | adaptive 0.02 | 4 | 0.477 | 3.57x | 52,895,924 | 230 | 1.11e+08 |
| 640x360 | adaptive 0.02 | 8 | 0.314 | 5.42x | 52,895,924 | 230 | 1.68e+08 |
| 640x360 | fixed 0.03 | 1 | 16.415 | 1.00x | 512,400,312 | 2224 | 3.12e+07 |
| 640x360 | fixed 0.03 | 4 | 4.408 | 3.72x | 512,400,312 | 2224 | 1.16e+08 |
| 640x360 | fixed 0.03 | 8 | 2.861 | 5.74x | 512,400,312 | 2224 | 1.79e+08 |
| 1280x720 | adaptive 0.02 | 1 | 7.713 | 1.00x | 211,571,772 | 230 | 2.74e+07 |
| 1280x720 | adaptive 0.02 | 4 | 2.170 | 3.55x | 211,571,772 | 230 | 9.75e+07 |
| 1280x720 | adaptive 0.02 | 8 | 1.484 | 5.20x | 211,571,772 | 230 | 1.43e+08 |
| 1280x720 | fixed 0.03 | 1 | 70.411 | 1.00x | 2,049,486,832 | 2224 | 2.91e+07 |
| 1280x720 | fixed 0.03 | 4 | 20.383 | 3.45x | 2,049,486,832 | 2224 | 1.01e+08 |
| 1280x720 | fixed 0.03 | 8 | 15.169 | 4.64x | 2,049,486,832 | 2224 | 1.35e+08 |

#### Accuracy vs. cost (640x360, vs. adaptive step 0.002 reference)

| step mode | mean steps/ray | render time (s) | captured % | pixels differing from reference |
|---|---|---|---|---|
| adaptive 0.002 (reference) | 2288 | 4.247 | 5.375 | 0.000% |
| adaptive 0.01 | 458 | 0.768 | 5.375 | 0.002% |
| adaptive 0.02 | 230 | 0.387 | 5.375 | 0.001% |
| adaptive 0.05 | 103 | 0.170 | 5.375 | 0.001% |
| adaptive 0.1 | 74 | 0.143 | 5.375 | 0.004% |
| fixed 0.03 | 2224 | 3.543 | 5.375 | 0.000% |
| fixed 0.1 | 668 | 1.004 | 5.375 | 0.001% |
| fixed 0.3 | 223 | 0.343 | 5.375 | 0.012% |
| fixed 1 | 67 | 0.106 | 4.823 | 1.123% |

Heatmaps (1280x720, same 100–20000 color range): `renders/heatmap_adaptive.png`, `renders/heatmap_fixed.png`.

### Observations
- **OpenMP scaling:** ~3.5–3.7x on 4 threads, ~5.2–5.7x on 8. The 4→8 gain is smaller because threads
  5–8 run on the M3's slower efficiency cores. The longest runs (70 s single-thread at 1280x720 fixed)
  scaled worst (4.6x on 8). I suspect thermal throttling during those long runs, but I didn't measure it.
- **Raw speed is identical** in both modes (~3e7 steps/s single-thread, ~33 ns per RK4 step = 4 force
  evaluations). All of the time difference comes from step *count*.
- **Fixed 0.03 vs adaptive 0.02:** same image (both match the reference to ≤0.001% of pixels), but fixed
  takes **9.7x more steps and time**. It uses the tiny photon-sphere step all the way out to r = 50.
  The adaptive heatmap shows work concentrated in a thin ring at the shadow edge (rays that skim
  r = 1.5); the fixed heatmap is ~2200 steps almost everywhere.
- **Fair comparison at equal cost:** fixed 0.3 uses about as many steps as adaptive 0.02 (223 vs 230)
  and is still nearly correct (0.012% of pixels wrong vs 0.001%) — 10x more errors, but both are small,
  because RK4 is 4th order and in this scene the camera is only at r = 20. Fixed steps fail quickly once
  they get too coarse: at `--fixed 1` the shadow shrinks (4.82% captured vs 5.375%) and 1.1% of pixels
  are wrong. The likely cause is that a step of 1 is as big as the horizon radius, so near-grazing
  rays can step past the r < 1 test or drift off their true path (not checked ray by ray). Adaptive steps stay correct even at 0.1 (74 steps/ray).
  Adaptive should pull further ahead for distant cameras, where fixed steps waste even more work
  crossing empty space (not measured).
- Caveat: the "pixels differing" metric uses the checker sky at 640x360, so it under-counts errors in
  the very thin photon ring. The Stage 3 critical-impact-parameter test is a sharper accuracy check.

### Problems hit
- The accuracy summary first crashed: the renderer writes step 1.0 as `1` (`%g`), but the script saved
  `fixed_1.0.png`. Fixed by using `1` in `bench.sh`.
