# render_image

A high-performance N-body simulation volume renderer. Reads particle snapshots in Gadget binary or HDF5 format and produces PNG images with full control over colourmap, opacity, density scaling, and SPH kernel smoothing. Supports single-frame renders, zoom sequences, and rotation animations. Parallelised with OpenMP and optionally MPI.

---

## Features

- **CIC and SPH deposit modes** — cloud-in-cell for speed, or kernel-smoothed SPH projection for publication-quality images
- **Tiled CIC deposit** — processes large images (4K, 8K) in horizontal strips, keeping peak RAM under 300 MB regardless of output resolution
- **O(N) smoothing length estimator** (`-fast_smooth`) — replaces ~160s exact KNN search with a sub-second density-grid estimate; smoothing lengths are cached to disk
- **Direct 2D SPH projection** — projects particles straight onto the image plane using the analytically integrated cubic spline kernel, bypassing the 3D grid entirely
- **10 built-in colourmap palettes** — viridis, magma, inferno, plasma, hot, fire, ice, grayscale, coolwarm, custom
- **6 opacity transfer functions** — flat, linear, sqrt, power, log, threshold
- **Auto-levelling** — percentile-based vmin/vmax with optional frame-to-frame locking for animations
- **Scene presets** — cluster, scattered, filament
- **Animation support** — zoom sequences and rotation sequences with Rodrigues rotation
- **Multi-particle-type support** — gas, dark matter, stars, or any bitmask combination
- **OpenMP parallelism** — deposit and neighbour search scale across all available cores
- **MPI parallelism** — x-slab decomposition for distributed-memory clusters

---

## Dependencies

| Library | Purpose |
|---------|---------|
| libpng | PNG output |
| libhdf5 | HDF5 snapshot reading |
| libz | zlib compression (HDF5 dependency) |
| OpenMP | Shared-memory parallelism (optional) |
| MPI | Distributed-memory parallelism (optional) |

On macOS with Homebrew: `brew install gcc hdf5 libpng`

---

## Building

Edit the `COMPILE_ON_SYSTEM` variable at the top of `makefile` to match your platform, then:

```bash
make clean && make
```

Supported platform targets: `MacBook`, `MacPro`, `Magnus`, `OzSTAR`.

### Compile-time options

Set these in the `OPT` lines in the makefile or pass on the command line:

| Flag | Description |
|------|-------------|
| `-DCIC` | Cloud-in-cell deposit (default, fast) |
| `-DKERNEL_SMOOTHING` | SPH kernel deposit (slower, smoother) |
| `-DNONPERIODIC` | Non-periodic boundary conditions |
| `-DENABLE_OPENMP` | Enable OpenMP threading |
| `-DENABLE_MPI` | Enable MPI for distributed runs |
| `-DLONG_IDS` | 64-bit particle IDs |
| `-DMAX_H_VOXELS=N` | Max kernel radius in pixels at 768px reference (default 8) |
| `-DMAX_H_DEPOSIT_PX=N` | Hard pixel cap applied during SPH deposit (default 8) |
| `-DMIN_RHO_FRAC=f` | Skip particles below this fraction of mean density (default 0.1) |
| `-DIMAGE_DIMENSIONX=N` | Output image width in pixels (default 768) |
| `-DIMAGE_DIMENSIONY=N` | Output image height in pixels (default 768) |
| `-DCIC_STRIP_HEIGHT=N` | CIC strip height in rows (default 64); set to 0 for original monolithic allocation |

### Recommended compiler flags

Add these to `OPTS` in the makefile for best performance:

```makefile
OPTS += -O3 -ffast-math
```

`-O3 -ffast-math` gives 20–40% speedup on the deposit inner loops by enabling vectorisation and floating-point reassociation. Safe for rendering — small floating-point differences are visually imperceptible.

---

## Usage

```
./render_image -input <snapshot> -output <prefix> [options]
```

### Required

| Flag | Description |
|------|-------------|
| `-input <file>` | Snapshot root path (without `.hdf5` extension) |
| `-output <file>` | Output image filename root |

### Particle selection

| Flag | Description |
|------|-------------|
| `-isHDF5` | HDF5 input format (default: Gadget binary) |
| `-dark_matter` | Render type 1 particles (default) |
| `-gas` | Render type 0 particles |
| `-stars` | Render type 4 particles |
| `-all_types` | Render all particle types |
| `-ptype <N>` | Keep particle type N (repeatable) |

### View volume

| Flag | Description |
|------|-------------|
| `-xc <val>` | X centre of view volume |
| `-yc <val>` | Y centre of view volume |
| `-zc <val>` | Z centre of view volume |
| `-lbox <val>` | Side length of view volume (same units as snapshot) |

### Colour and opacity

| Flag | Description |
|------|-------------|
| `-colormap <n>` | Palette: `viridis` `magma` `inferno` `plasma` `hot` `fire` `ice` `grayscale` `coolwarm` `custom` |
| `-reverse_colormap` | Reverse the palette |
| `-bg_color R,G,B[,A]` | Background colour in [0,1] (default `0,0,0,1`) |
| `-opacity <val>` | Global opacity multiplier 0–1 (default 1) |
| `-opacity_func <n>` | Transfer function: `flat` `linear` `sqrt` `power` `log` `threshold` |
| `-opacity_gamma <val>` | Exponent for `power` function |
| `-opacity_threshold <val>` | Cutoff for `threshold` function |

### Density scaling

| Flag | Description |
|------|-------------|
| `-auto_pct_lo <val>` | Low clip percentile 0–1 (default 0.001) |
| `-auto_pct_hi <val>` | High clip percentile 0–1 (default 0.999) |
| `-no_auto_levels` | Use fixed vmin/vmax instead of auto-levelling |
| `-vmin <val>` | Lower density clip in log10 units (with `-no_auto_levels`) |
| `-vmax <val>` | Upper density clip in log10 units (with `-no_auto_levels`) |
| `-linear_scale` | Use linear density scaling instead of log10 |

### Scene presets

| Flag | Description |
|------|-------------|
| `-scene cluster` | Magma colormap, sqrt opacity, 5% low clip — good for galaxy clusters |
| `-scene scattered` | Plasma colormap, flat opacity, 10% low clip — good for sparse distributions |
| `-scene filament` | Inferno colormap, power γ=0.5 opacity, 2% low clip — good for cosmic web |

### Animation

| Flag | Description |
|------|-------------|
| `-itmax <N>` | Render N frames with fixed view volume |
| `-zoom <N>` | Render N frames, shrinking the box each frame |
| `-zoom_factor <f>` | Multiply box side by f per zoom frame, 0 < f < 1 (default 0.5) |
| `-rot_dangle <deg>` | Rotate view by this many degrees per frame |
| `-rot_axis x,y,z` | Rotation axis (default `0,0,1`) |
| `-lock_levels` | Auto-level on frame 0, lock vmin/vmax for all subsequent frames |

Output frames are named `<prefix>.NNNN.png`.

### SPH kernel smoothing (requires `-DKERNEL_SMOOTHING`)

| Flag | Description |
|------|-------------|
| `-fast_smooth` | O(N) grid-based h estimator (sub-second) instead of exact KNN (~160s) |
| `-sph_cache <file>` | Cache smoothing lengths; reloaded if particle count and parameters match |
| `-num_ngb <N>` | Target neighbour count for smoothing length (default 32) |
| `-sph_eta <val>` | h = eta × dist_to_Nth_neighbour (default 1.2; larger = smoother) |
| `-ngrid_z <N>` | Depth of 3D density grid for CIC mode (default = image width) |

---

## Examples

### Single cluster image

```bash
./render_image \
    -input snapshot_122 -output cluster \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 0.2 \
    -dark_matter -scene cluster \
    -fast_smooth -sph_cache snap122.hcache -itmax 1
```

### Large volume render

```bash
./render_image \
    -input snapshot_122 -output volume \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 1.0 \
    -dark_matter -colormap plasma \
    -auto_pct_lo 0.1 -auto_pct_hi 0.9999 -bg_color 0,0,0,1 \
    -fast_smooth -sph_cache snap122.hcache -itmax 1
```

### 4K or 8K movie frame (CIC mode)

Compile with the desired resolution and tiled strip support:

```bash
make clean && make OPTS="-O3 -ffast-math -DCIC -DNONPERIODIC -DENABLE_OPENMP \
    -DIMAGE_DIMENSIONX=4096 -DIMAGE_DIMENSIONY=4096 -DCIC_STRIP_HEIGHT=64"
```

Then render:

```bash
./render_image \
    -input snapshot_122 -output frames_4k \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 1.0 \
    -dark_matter -colormap magma -ngrid_z 128 \
    -auto_pct_lo 0.1 -auto_pct_hi 0.9999 -bg_color 0,0,0,1 \
    -itmax 1
```

Peak RAM with `CIC_STRIP_HEIGHT=64`: 134 MB at 4K×4K×128, 268 MB at 8K×8K×128. Without tiling the equivalent monolithic grids would be 8 GB and 32 GB respectively.

### 36-frame rotation animation

```bash
./render_image \
    -input snapshot_122 -output frames \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 0.2 \
    -dark_matter -scene cluster \
    -fast_smooth -sph_cache snap122.hcache \
    -itmax 36 -rot_dangle 10 -rot_axis 0,1,0 -lock_levels

ffmpeg -framerate 24 -i frames.%04d.png -c:v libx264 -pix_fmt yuv420p rotation.mp4
```

### Zoom sequence

```bash
./render_image \
    -input snapshot_122 -output zoom \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 2.0 \
    -dark_matter -scene cluster \
    -fast_smooth -sph_cache snap122.hcache \
    -zoom 20 -zoom_factor 0.85 -lock_levels
```

---

## Source file overview

| File | Purpose |
|------|---------|
| `render_image.c` | Main entry point, CLI parsing, render loop, rotation/zoom logic |
| `io.c` | HDF5 and Gadget binary snapshot reader |
| `find_neighbours.c` | Smoothing length computation — O(N) density-grid estimator and exact KNN; disk caching |
| `flat_kd_tree.h` | Header-only flat-array kd-tree for exact KNN search |
| `deposit_sph_2d.c` | Direct 2D SPH kernel projection (separate translation unit for compiler inlining) |
| `smooth_to_mesh.c` | Density field construction — SPH 2D deposit or tiled CIC 3D deposit + projection |
| `kernels.c` | Cubic spline kernel and 2D projected kernel; lookup tables |
| `colormap.c` | Colour palettes, opacity transfer functions, `render_config_t` pipeline |
| `select_particles.c` | View-volume selection and particle type masking |
| `make_tree.c` | Octree construction for spatial indexing |
| `walk_tree.c` | Octree traversal functions |
| `split_across_tasks.c` | MPI x-slab decomposition |
| `header.c` | Global state definitions (ThisTask, NTask, slab boundaries) |
| `write_to_ppm.c` | PNG output with colormap, opacity, and background compositing |

### Header files

The headers are split into focused modules to minimise recompilation when individual interfaces change:

| Header | Contents |
|--------|---------|
| `types.h` | All POD structs and constants (`sim_info`, tree nodes, `ngb_buf_t`, `deposit_particle_t`, `pixel_t`, etc.) |
| `globals.h` | `extern` globals only (`ThisTask`, `NTask`, `slab_x_lo/hi`, `SnapFormat`) |
| `tree.h` | Octree build and walk function prototypes |
| `kernels.h` | SPH kernel function prototypes |
| `io.h` | Snapshot reader prototypes; `sim_info_alloc/free`; `ngb_buf_alloc/free` |
| `render.h` | Rendering pipeline prototypes (`smooth_to_mesh`, `deposit_sph_2d`, `write_to_png_ex`, etc.) |
| `colormap.h` | `render_config_t`, palette and opacity enums, colormap API |
| `header.h` | Compatibility umbrella — includes all of the above; existing code using `#include "header.h"` continues to work unchanged |

---

## Performance notes

Typical timings on a 4-core Apple M2 with `-O3 -ffast-math`:

| Step | 768×768 | 4096×4096 |
|------|---------|-----------|
| HDF5 read + particle selection | ~4s | ~4s |
| `-fast_smooth` smoothing lengths | ~0.9s | ~0.9s |
| Smoothing length cache load | <0.1s | <0.1s |
| SPH 2D deposit (`-DKERNEL_SMOOTHING`) | ~15s | ~300s |
| CIC tiled deposit | ~3s | ~50s |
| PNG write | <0.1s | ~1s |

The SPH deposit is memory-bandwidth limited. For 4K/8K movies CIC mode is strongly recommended — it is substantially faster and the difference in visual quality is imperceptible at large scales.

Smoothing lengths are stable between frames of a rotation animation, so the cache (`-sph_cache`) eliminates the neighbour search cost entirely after the first frame.

### CIC strip memory usage

| Resolution | Depth | Strip height | Peak RAM per strip |
|-----------|-------|-------------|-------------------|
| 768×768 | 128 | 64 | 25 MB |
| 2048×2048 | 128 | 64 | 67 MB |
| 4096×4096 | 128 | 64 | 134 MB |
| 8192×8192 | 128 | 64 | 268 MB |
| 8192×8192 | 256 | 64 | 536 MB |

Use `-DCIC_STRIP_HEIGHT=128` for fewer passes at the cost of 2× more RAM per strip. Use `-DCIC_STRIP_HEIGHT=0` to disable tiling and restore the original monolithic allocation (only practical for small images).

---

## Smoothing length tuning (`-fast_smooth`)

The O(N) estimator deposits particles onto a 128³ density grid using CIC, applies a single box-smooth pass, then estimates `h` from the local density via trilinear interpolation. Particles whose local density falls below `MIN_RHO_FRAC` (default 0.1) times the volume mean are skipped — they contribute negligible signal and their exclusion eliminates isolated-dot and grid-boundary artefacts in void regions. Skipped pixels are filled by the zero-hole inpainting pass in the render loop.

To adjust the density threshold:
- `-DMIN_RHO_FRAC=0.05` — keep more sparse particles (more structure visible in voids)
- `-DMIN_RHO_FRAC=0.2` — skip more sparse particles (cleaner background, more inpainting)

The smoothing length cap is expressed in physical units relative to a 768px reference width, keeping it consistent regardless of output resolution. At `lbox=1`, `MAX_H_VOXELS=8` gives `h_max = 8/768 ≈ 1% of lbox` whether rendering at 768px or 8192px.

---

## Custom colourmap

Pass a text file with one `R G B` triplet per line (values 0–255, 256 lines) via `-colormap_file <path>` combined with `-colormap custom`.

---

## Time sequences

### Seamless transitions between snapshots

Two things cause jarring transitions in a time sequence: brightness flicker (different auto-levels per snapshot) and positional jumps (particles teleporting between outputs). Both are solved independently.

**Fixing brightness flicker — global colour scale**

Run a quick pass over all snapshots to find the global density range, then render with fixed levels:

```python
import subprocess, re

snapshots = [f"snapshot_{i:03d}" for i in range(100, 130)]
vmins, vmaxs = [], []

for snap in snapshots:
    result = subprocess.run([
        "./render_image", "-input", snap, "-output", "/tmp/probe",
        "-isHDF5", "-xc", "56.18", "-yc", "43.35", "-zc", "49.13",
        "-lbox", "0.5", "-dark_matter", "-itmax", "1"
    ], capture_output=True, text=True)
    for line in result.stdout.splitlines():
        m = re.search(r"vmin=([\d.]+)\s+vmax=([\d.]+)", line)
        if m:
            vmins.append(float(m.group(1)))
            vmaxs.append(float(m.group(2)))

# Use a robust global range
global_vmin = sorted(vmins)[len(vmins) // 20]    # 5th percentile
global_vmax = sorted(vmaxs)[-len(vmaxs) // 20]   # 95th percentile

# Render full sequence with fixed scale
for i, snap in enumerate(snapshots):
    subprocess.run([
        "./render_image", "-input", snap,
        "-output", f"frame_{i:04d}",
        "-isHDF5", "-xc", "56.18", "-yc", "43.35", "-zc", "49.13",
        "-lbox", "0.5", "-dark_matter", "-colormap", "magma",
        "-no_auto_levels", "-vmin", str(global_vmin), "-vmax", str(global_vmax),
        "-itmax", "1"
    ])
```

**Fixing positional jumps — velocity interpolation**

The `-interp_frac` flag shifts particle positions forward by a fraction of the snapshot interval using the velocity field stored in the HDF5 file. This generates smooth sub-frame transitions without reading a second snapshot.

```bash
# Single interpolated frame: half-way between snapshot_122 and snapshot_123
./render_image -input snapshot_122 -output frame_half \
    -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 0.5 \
    -dark_matter -colormap magma \
    -no_auto_levels -vmin 0.1 -vmax 3.5 \
    -interp_frac 0.5 -snap_dt 0.05 -itmax 1
```

| Flag | Description |
|------|-------------|
| `-interp_frac <f>` | Shift positions by `f × snap_dt × velocity`; f ∈ [0,1] |
| `-snap_dt <val>` | Time interval between snapshots in simulation units |

If `-snap_dt` is omitted, velocities are used as raw offsets scaled only by `interp_frac` — useful for exploring the velocity field visually.

**Generating a full interpolated sequence between two snapshots:**

```bash
N=10
for i in $(seq 0 $N); do
    FRAC=$(python3 -c "print($i / $N)")
    ./render_image -input snapshot_122 \
        -output "interp_${i}" \
        -isHDF5 -xc 56.18 -yc 43.35 -zc 49.13 -lbox 0.5 \
        -dark_matter -colormap magma \
        -no_auto_levels -vmin 0.1 -vmax 3.5 \
        -interp_frac "$FRAC" -snap_dt 0.05 -itmax 1
done
```

**Full pipeline — combining both fixes in the notebook:**

```python
snapshots = [f"snapshot_{i:03d}" for i in range(100, 130)]
n_interp   = 10       # sub-frames between each snapshot pair
global_vmin = 0.1     # determined from the probe pass above
global_vmax = 3.5

frame = 0
for snap in snapshots:
    for i in range(n_interp):
        frac = i / n_interp
        subprocess.run([
            "./render_image", "-input", snap,
            "-output", f"movie_{frame:05d}",
            "-isHDF5", "-xc", "56.18", "-yc", "43.35", "-zc", "49.13",
            "-lbox", "0.5", "-dark_matter", "-colormap", "magma",
            "-no_auto_levels",
            "-vmin", str(global_vmin), "-vmax", str(global_vmax),
            "-interp_frac", str(frac), "-snap_dt", "0.05",
            "-itmax", "1"
        ])
        frame += 1

# Assemble with ffmpeg
subprocess.run([
    "ffmpeg", "-y", "-framerate", "24",
    "-i", "movie_%05d.0000.png",
    "-c:v", "libx264", "-pix_fmt", "yuv420p",
    "simulation.mp4"
])
```

This produces a movie with `n_snapshots × n_interp` frames. Particles move smoothly between outputs; brightness is consistent across the full sequence.

**Limitations of velocity interpolation**

The interpolation is first-order (straight-line motion), which is accurate for small `interp_frac` values but breaks down near strong interactions or mergers where the velocity field changes rapidly within a snapshot interval. For cosmological dark matter simulations with `n_interp = 10` the result is visually seamless. Reduce `n_interp` or increase `snap_dt` if you see unphysical particle trajectories.
