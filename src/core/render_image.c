#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <hdf5.h>

#include "header.h"
#include "colormap.h"
#include "args.h"
#include "config.h"

#define BUFSIZE 500
#define rhocrit 27.755

#include <png.h>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

void find_neighbours(int num_part, float *smoothing_length, int num_ngb,
                     float *posx, float *posy, float *posz,
                     float xmin, float ymin, float zmin, float sph_eta);
void find_neighbours_fast(int num_part, float *smoothing_length, int num_ngb,
                           float *posx, float *posy, float *posz,
                           float xmin, float ymin, float zmin,
                           float sph_eta, const char *cache_file, int render_width);
void find_neighbours_cached(int num_part, float *smoothing_length, int num_ngb,
                             float *posx, float *posy, float *posz,
                             float xmin, float ymin, float zmin,
                             float sph_eta, const char *cache_file);

#ifndef IMAGE_DIMENSIONX
#define IMAGE_DIMENSIONX  768
#endif
#ifndef IMAGE_DIMENSIONY
#define IMAGE_DIMENSIONY  768
#endif
#define MAX_COLOUR_COMPONENT_VALUE 255

static void print_usage(const char *prog)
{
    fprintf(stdout,
        "Usage: %s -input <file>           Input snapshot root\n"
        "          -output <file>          Output image root\n"
        "          [-isHDF5]               HDF5 format (default: binary)\n"
        "          [-units 0|1]            0=raw, 1=box-normalised coords\n"
        "          [-xc <val>]             x centre\n"
        "          [-yc <val>]             y centre\n"
        "          [-zc <val>]             z centre\n"
        "          [-lbox <val>]           View box side length\n"
        "          [-itmax <N>]            N frames, fixed box (default 1)\n"
        "          [-zoom <N>]             N frames, box shrinks each iter\n"
        "          [-zoom_factor <f>]      Zoom factor per iter, 0<f<1 (default 0.5)\n"
        "          [-rot_dangle <deg>]     Degrees to rotate per frame\n"
        "          [-rot_axis x,y,z]       Rotation axis (default 0,0,1=z)\n"
        "          [-ptype <N>]            Keep particle type N (repeatable)\n"
        "          [-all_types]            Keep all particle types\n"
        "          [-gas]                  Shorthand: keep type 0\n"
        "          [-dark_matter]          Shorthand: keep type 1\n"
        "          [-stars]                Shorthand: keep type 4\n"
        "\n"
        "Colour / opacity options:\n"
        "          [-scene <n>]      Preset: cluster | scattered | filament\n"
        "          [-colormap <name>]      Palette name (see below)\n"
        "          [-reverse_colormap]     Reverse the palette\n"
        "          [-opacity <val>]        Global opacity [0,1] (default 1)\n"
        "          [-opacity_func <name>]  Transfer function (see below)\n"
        "          [-opacity_gamma <val>]  Exponent for 'power' function\n"
        "          [-opacity_threshold <v>] Cutoff value for 'threshold'\n"
        "          [-vmin <val>]           Lower density clip (log10 units)\n"
        "          [-vmax <val>]           Upper density clip (log10 units)\n"
        "          [-linear_scale]         Use linear (not log10) density\n"
        "          [-colormap_file <path>] Load custom palette from file\n"
        "          [-bg_color R,G,B[,A]]   Background colour in [0,1]\n"
        "\n"
        "Available palettes: viridis, magma, inferno, plasma, hot, fire,\n"
        "                    ice, grayscale, coolwarm, custom\n"
        "Opacity functions:  flat, linear, sqrt, power, log, threshold\n"
        "\nDensity scaling (auto-levels on by default):\n"
        "          [-no_auto_levels]       Use fixed vmin/vmax instead\n"
        "          [-lock_levels]          Auto-level on frame 0, lock for all frames\n"
        "          [-auto_pct_lo <val>]    Low  clip percentile 0-1 (default 0.001)\n"
        "          [-auto_pct_hi <val>]    High clip percentile 0-1 (default 0.999)\n"
        "  With -no_auto_levels, -vmin/-vmax set fixed log10 density bounds\n"
        "\nTime interpolation options:\n"
        "          [-interp_frac <f>]      Interpolate positions by fraction f in [0,1]\n"
        "                                  using velocities from the input snapshot.\n"
        "                                  f=0: snapshot positions; f=1: one full\n"
        "                                  snapshot interval ahead. Use with -itmax N\n"
        "                                  to generate sub-frame interpolated sequences.\n"
        "          [-snap_dt <val>]        Time interval between snapshots in simulation\n"
        "                                  units (used to scale velocity displacement).\n"
        "                                  If omitted, velocities are used as raw offsets\n"
        "                                  scaled only by interp_frac.\n"
        "\nKernel smoothing options (requires -DKERNEL_SMOOTHING):\n"
        "          [-num_ngb <N>]          Neighbours for h estimate (default 32)\n"
        "          [-sph_eta <val>]        h = eta*dist_Nth_ngb (default 1.2)\n"
        "          [-sph_cache <file>]     cache smoothing lengths to file\n"
        "                                  (reloaded on next run if n/ngb/eta match)\n"
        "          [-fast_smooth]          O(N) grid h estimator (~0.3s) vs KNN (~30s)\n"
        "          Larger eta: smoother; smaller eta: sharper but noisier\n"
        "\n3-D grid options:\n"
        "          [-ngrid_z <N>]          Depth of 3-D density grid (default=min(width,256))\n"
        "          Reduce for lower memory use, e.g. -ngrid_z 64\n"
        "          At 4K+ always set this explicitly to avoid huge allocations\n",
        prog);
    fflush(stdout);
}

int main(int argc, char *argv[])
{
    int i;
    char filename[256]        = "";
    char image_file[256]      = "";

    long long NumPart = 0, NumPartRead = 0;

    float *x = NULL, *y = NULL, *z = NULL;
    float *vx = NULL, *vy = NULL, *vz = NULL;
    int   *ptype = NULL;
    float *smoothing_length = NULL;

    int isDistributed = 0;

    struct sim_info header;
    memset(&header, 0, sizeof(header));

    int ptype_mask = (1 << 1);   /* default: dark matter */

    double xc= 0.0, yc= 0.0, zc = 0.0;
    clock_t tstart, tfinish;

    cli_args_t cfg;
    cli_args_default(&cfg);

#ifdef ENABLE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
    ThisTask = 0;
    NTask    = 1;
#endif

    if (argc < 2) {
        if (ThisTask == 0) print_usage(argv[0]);
#ifdef ENABLE_MPI
        MPI_Finalize();
#endif
        return 0;
    }
    
    for (int k = 1; k < argc - 1; k++) {
        if (strcmp(argv[k], "-config") == 0) {
            snprintf(cfg.config_path, sizeof(cfg.config_path),
                     "%s", argv[k + 1]);
            break;
        }
    }
    if (cfg.config_path[0])
        load_yaml_config(cfg.config_path, &cfg);   /* YAML baseline */
    
    parse_args(argc, argv, &cfg);                  /* CLI wins      */
    
    /* Derive simulation-unit centre from user input */
    xc = cfg.xcen; yc = cfg.ycen; zc = cfg.zcen;

    /* ---- Setup ---- */
    check_input_filenames(filename, cfg.file_root, cfg.isHDF5, &isDistributed);

    if (ThisTask == 0) {
        if (cfg.xcen > 0 || cfg.lbox > 0) {
            fprintf(stdout, "Centre: (%g|%g|%g)\n", cfg.xcen, cfg.ycen, cfg.zcen);
            fprintf(stdout, "Box Length: %g\n", cfg.lbox);
        }
        fprintf(stdout, "Render configuration:\n");
        render_config_print(&cfg.rcfg);
        fprintf(stdout, "Reading header...\n");
        fflush(stdout);
    }

    if (cfg.isHDF5)
        read_hdf5_header(filename, &header, &NumPart);
    else
        read_gadget_binary_header(filename, &header, &NumPart);

    if (ThisTask == 0) {
        fprintf(stdout, "Number of particle types in file: %d\n",
                header.num_types);
        for (i = 0; i < (int)(sizeof(int)*8); i++) {
            if (cfg.ptype_mask != -1 && !(cfg.ptype_mask & (1 << i))) continue;
            if (i >= header.num_types || header.nall[i] == 0)
                fprintf(stdout,
                    "Warning: type %d requested but not present\n", i);
        }
        fflush(stdout);
    }

    if (header.BoxSize == 0)
        header.BoxSize = 1.e6;
    else if (cfg.boxunits == 1) {
        cfg.xcen *= header.BoxSize; cfg.ycen *= header.BoxSize;
        cfg.zcen *= header.BoxSize; cfg.lbox *= header.BoxSize;
    }
    xc= cfg.xcen; yc= cfg.ycen; zc = cfg.zcen;

    if (ThisTask == 0) {
        fprintf(stdout, "Number of files: %d\n",       header.NumFiles);
        fprintf(stdout, "Number of particles: %lld\n", NumPart);
        fprintf(stdout, "Time/Expansion Factor: %g\n", header.time);
        fflush(stdout);
    }

    x     = (float *)malloc(sizeof(float) * NumPart);
    y     = (float *)malloc(sizeof(float) * NumPart);
    z     = (float *)malloc(sizeof(float) * NumPart);
    ptype = (int   *)malloc(sizeof(int)   * NumPart);
    
    int do_interp  = (cfg.interp_frac != 0.0f && cfg.isHDF5);
    
    if (do_interp) {
        vx = (float *)malloc(sizeof(float) * NumPart);
        vy = (float *)malloc(sizeof(float) * NumPart);
        vz = (float *)malloc(sizeof(float) * NumPart);
        
        if (!vx || !vy || !vz) {
            fprintf(stderr, "Failed to allocate velocity arrays\n");
            exit(1);
        }
    }
    
    if (ThisTask == 0) { fprintf(stdout, "Reading particles...\n"); fflush(stdout); }
    
    if (cfg.isHDF5)
        read_particles_from_hdf5(cfg.file_root, x, y, z, vx, vy, vz, ptype,
                                  header.NumFiles, &NumPartRead, do_interp);
    else
        read_particles_from_gadget_binary(cfg.file_root, x, y, z, ptype,
                                           header.NumFiles, &NumPartRead);
   
    fprintf(stdout, "NumPart: %llu\t NumPartRead: %llu\n", NumPart,NumPartRead);
    fflush(stdout);

    
    if (do_interp) {
        /*
         * Velocity interpolation: shift particle positions by frac * v * dt
         * to generate sub-snapshot frames between two outputs.
         *
         * This produces smooth transitions in a time sequence without needing
         * to store or read a second snapshot.  The approximation is first-order
         * (straight-line motion), which is accurate for small fractions of the
         * snapshot interval and breaks down near strong interactions.
         *
         * Usage:
         *   -interp_frac 0.5          half-step forward using snapshot velocities
         *   -snap_dt 0.05             snapshot interval in simulation time units
         *
         * If -snap_dt is omitted, velocities are used as raw offsets (useful when
         * you just want to explore the velocity field visually).
         *
         * To generate N interpolated frames between snapshots A and B:
         *   for f in $(seq 0 0.1 1.0); do
         *     ./render_image.exe -input snapA ... -interp_frac $f -snap_dt 0.05
         *   done
         */
        float scale = (cfg.snap_dt > 0.0f) ? cfg.interp_frac * cfg.snap_dt : cfg.interp_frac;
        
        fprintf(stdout,
                "Interpolating positions: frac=%.3f dt=%g scale=%g\n",
                cfg.interp_frac, cfg.snap_dt, scale);
        fflush(stdout);
        
        for (long long k = 0; k < NumPartRead; k++) {
            x[k] += scale * vx[k];
            y[k] += scale * vy[k];
            z[k] += scale * vz[k];
        }
        
        free(vx); free(vy); free(vz);
    }

    select_particles(x, y, z, ptype, header.BoxSize, NumPartRead,
                     xc, yc, zc, cfg.lbox, cfg.ptype_mask, &NumPart);

#ifdef ENABLE_MPI
    split_across_tasks_as_slabs(x, y, z, &NumPart, xc, yc, zc,
                             (float)cfg.lbox, (float)header.BoxSize);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (NumPart == 0) {
        if (ThisTask == 0)
            fprintf(stdout, "Error: no particles selected (mask=0x%x)\n",
                    (unsigned)cfg.ptype_mask);
#ifdef ENABLE_MPI
        MPI_Finalize();
#endif
        return 1;
    }
    
    if (ThisTask == 0) {
        fprintf(stdout, "Plotting %lld particles (type mask 0x%x)...\n",
                NumPart, (unsigned)ptype_mask);
        fflush(stdout);
    }

#ifdef KERNEL_SMOOTHING
    /*
     * Adaptive smoothing lengths for SPH kernel deposit.
     * Not needed for CIC mode — smooth_to_mesh ignores smoothing_length
     * when KERNEL_SMOOTHING is not defined.
     */
    smoothing_length = (float *)malloc(sizeof(float) * NumPart);
    if (!smoothing_length) {
        fprintf(stderr, "malloc failed for smoothing_length\n");
        return 1;
    }
    {
        float xmin_v = (float)(xc - 0.5 * cfg.lbox);
        float ymin_v = (float)(yc - 0.5 * cfg.lbox);
        float zmin_v = (float)(zc - 0.5 * cfg.lbox);
        tstart = clock();
        /* Pass physical voxel size as xmin_v so find_neighbours can apply
         * a resolution-based h cap.  ymin_v/zmin_v are unused by the new code. */
        float vox_size_v = (float)cfg.lbox / (float)IMAGE_DIMENSIONX;
        if (cfg.fast_smooth) {
            find_neighbours_fast((int)NumPart, smoothing_length, cfg.num_ngb,
                                 x, y, z, vox_size_v, 0.0f, 0.0f, cfg.sph_eta,
                                 cfg.sph_cache[0] ? cfg.sph_cache : NULL,
                                 IMAGE_DIMENSIONX);
        } else {
            find_neighbours_cached((int)NumPart, smoothing_length, cfg.num_ngb,
                                   x, y, z, vox_size_v, 0.0f, 0.0f, cfg.sph_eta,
                                   cfg.sph_cache[0] ? cfg.sph_cache : NULL);
        }
        tfinish = clock();
        fprintf(stdout, "find_neighbours - time: %.2f s\n",
                (double)(tfinish - tstart) / CLOCKS_PER_SEC);
        fflush(stdout);
    }
#else
    smoothing_length = NULL;   /* CIC deposit does not use smoothing lengths */
    fprintf(stdout, "CIC deposit mode (no kernel smoothing)\n");
    fflush(stdout);
#endif

    /* ---- Render loop ---- */
    /*
     * smooth_to_mesh fills data[][] with raw density values.
     * We pass the raw data and the full render_config_t to
     * write_to_png_ex, which applies normalisation, colourmap, and
     * opacity in one pass — avoiding an extra loop over the image.
     */
    float *data        = (float *)calloc(IMAGE_DIMENSIONX * IMAGE_DIMENSIONY,
                                          sizeof(float));
    float *global_data = (float *)calloc(IMAGE_DIMENSIONX * IMAGE_DIMENSIONY,
                                          sizeof(float));

    float BoundingBox = (float)cfg.lbox;
    float theta       = 0.0f;

    /*
     * Render loop — itmax frames, fixed bounding box by default.
     *
     * Modes:
     *   -itmax N              N frames, fixed box (rotation if -rot_dangle set)
     *   -zoom  N [-zoom_factor f]  N frames, box multiplied by f each iter (default 0.5)
     *   -itmax N -rot_dangle D  N frames rotating D degrees per frame around
     *                            rot_axis through (xc,yc,zc)
     *
     * Rotation is implemented by rotating a working copy of the particle
     * coordinates around the view centre before each deposit.  The original
     * arrays are untouched so frame 0 always matches a non-rotating render.
     *
     * Rodrigues rotation formula:
     *   v_rot = v*cos(a) + (k x v)*sin(a) + k*(k.v)*(1-cos(a))
     * where k is the unit rotation axis and a is the angle in radians.
     */
    int n_frames = (cfg.zoom > 0) ? cfg.zoom : cfg.itmax;

    /* Working coordinate arrays for rotation */
    float *rx = (float *)malloc(NumPart * sizeof(float));
    float *ry = (float *)malloc(NumPart * sizeof(float));
    float *rz = (float *)malloc(NumPart * sizeof(float));
    if (!rx || !ry || !rz) {
        fprintf(stderr, "malloc failed for rotation buffers");
        return 1;
    }

    for (int iter = 0; iter < n_frames; iter++) {
        memset(data,        0, IMAGE_DIMENSIONX * IMAGE_DIMENSIONY * sizeof(float));
        memset(global_data, 0, IMAGE_DIMENSIONX * IMAGE_DIMENSIONY * sizeof(float));

        /* Apply cumulative rotation for this frame.
         * Angle = iter * rot_dangle (degrees).  For iter=0 this is 0 so
         * rx/ry/rz == x/y/z and the first frame is always un-rotated. */
        {
            float angle_rad = (float)(iter) * cfg.rot_dangle * (float)M_PI / 180.0f;
            float ca = cosf(angle_rad), sa = sinf(angle_rad);
            float kx = cfg.rot_axis[0], ky = cfg.rot_axis[1], kz = cfg.rot_axis[2];
            for (long long pi = 0; pi < NumPart; pi++) {
                float vx = x[pi] - (float)xc;
                float vy = y[pi] - (float)yc;
                float vz = z[pi] - (float)zc;
                float kdotv = kx*vx + ky*vy + kz*vz;
                /* cross product k x v */
                float cx_ = ky*vz - kz*vy;
                float cy_ = kz*vx - kx*vz;
                float cz_ = kx*vy - ky*vx;
                rx[pi] = (float)xc + vx*ca + cx_*sa + kx*kdotv*(1.0f - ca);
                ry[pi] = (float)yc + vy*ca + cy_*sa + ky*kdotv*(1.0f - ca);
                rz[pi] = (float)zc + vz*ca + cz_*sa + kz*kdotv*(1.0f - ca);
            }
        }

        /* ngrid_z: depth of the 3D CIC grid in z-slices.
         * Defaulting to IMAGE_DIMENSIONX would give a cubic grid —
         * 4096^3 at 4K = 256 GB. Cap at 256 slices which is sufficient
         * for all projection renders and keeps memory sane.
         * Override explicitly with -ngrid_z if you need finer z-sampling. */
        int effective_ngrid_z = cfg.ngrid_z > 0 ? cfg.ngrid_z
                              : (IMAGE_DIMENSIONX < 256 ? IMAGE_DIMENSIONX : 256);

        tstart = clock();
        smooth_to_mesh(NumPart, smoothing_length, rx, ry, rz,
                       xc, yc, zc, BoundingBox,
                       (float)effective_ngrid_z,
                       IMAGE_DIMENSIONX, IMAGE_DIMENSIONY, data);
        tfinish = clock();
        fprintf(stdout, "smooth_to_mesh - time: %.2f s\n",
                (double)(tfinish - tstart) / CLOCKS_PER_SEC);
        fflush(stdout);

#ifdef ENABLE_MPI
        MPI_Reduce(data, global_data,
                   IMAGE_DIMENSIONX * IMAGE_DIMENSIONY,
                   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        memcpy(global_data, data,
               IMAGE_DIMENSIONX * IMAGE_DIMENSIONY * sizeof(float));
#endif

        if (ThisTask == 0) {
            /* Apply noise floor globally after MPI_Reduce(SUM).
             * Must be done here — not per-task — because each task has a
             * different local dmax, so per-task noise floors are inconsistent
             * and create visible stripe artefacts at slab boundaries.
             * A single global dmax gives a consistent threshold everywhere. */
            int npix_g = IMAGE_DIMENSIONX * IMAGE_DIMENSIONY;
            float gdmax = 0.0f;
            for (int p = 0; p < npix_g; p++)
                if (global_data[p] > gdmax) gdmax = global_data[p];
            /*
             * Noise floor: zero pixels below 0.01% of the peak.
             * 1e-4 (vs the old 1e-6) aggressively suppresses corner
             * speckle from boundary particles with poorly-sampled kernels,
             * while staying well below any real density structure.
             * At dmax~5600 (log10=3.75) this zeros anything below ~0.56,
             * i.e. log10 < -0.25 — safely above the corner noise level.
             */
            float noise_floor = gdmax * 1.0e-4f;
            if (noise_floor > 0.0f)
                for (int p = 0; p < npix_g; p++)
                    if (global_data[p] < noise_floor) global_data[p] = 0.0f;

            /*
             * CIC zero-hole inpainting.
             *
             * With ~0.03 particles/voxel on a 768³ grid the CIC 2×2×2
             * stencil leaves contiguous empty regions several pixels wide
             * after projection.  A fixed-radius blur needs many passes to
             * fill large holes and risks blurring real structure.
             *
             * Instead we use iterative 3×3 propagation: each pass expands
             * the filled frontier by one pixel in every direction.  We
             * repeat until no zero pixels remain inside the data extent, or
             * we hit MAX_FILL_PASSES.  Only zero pixels are ever written;
             * non-zero pixels are never modified, so dynamic range and
             * sharpness of real structure are fully preserved.
             *
             * MAX_FILL_PASSES = 64 handles holes up to ~64 px across, which
             * comfortably covers the worst-case CIC aliasing at this
             * resolution.  Each pass is O(npix) so total cost is negligible.
             */
            {
#define MAX_FILL_PASSES 64
                float *buf = (float *)malloc((size_t)npix_g * sizeof(float));
                if (buf) {
                    memcpy(buf, global_data, (size_t)npix_g * sizeof(float));
                    int changed = 1;
                    for (int pass = 0; pass < MAX_FILL_PASSES && changed; pass++) {
                        changed = 0;
                        for (int row = 1; row < IMAGE_DIMENSIONY - 1; row++) {
                            for (int col = 1; col < IMAGE_DIMENSIONX - 1; col++) {
                                int pidx = row * IMAGE_DIMENSIONX + col;
                                if (buf[pidx] > 0.0f) continue; /* already filled */
                                float sum = 0.0f; int nn = 0;
                                for (int dy = -1; dy <= 1; dy++)
                                    for (int dx = -1; dx <= 1; dx++) {
                                        float v = buf[(row+dy)*IMAGE_DIMENSIONX+(col+dx)];
                                        if (v > 0.0f) { sum += v; nn++; }
                                    }
                                if (nn > 0) {
                                    global_data[pidx] = sum / (float)nn;
                                    changed = 1;
                                }
                            }
                        }
                        /* Sync buf so next pass sees this pass's fills */
                        memcpy(buf, global_data, (size_t)npix_g * sizeof(float));
                    }
                    free(buf);
                }
#undef MAX_FILL_PASSES
            }

            sprintf(image_file, "%s.%04d.png", cfg.image_file_root, iter);

            /*
             * For -lock_levels: compute auto-levels NOW, before inpainting,
             * so inpainted void pixels don't drag vmin down and corrupt the
             * colour scale.  With SPH + MIN_RHO_FRAC, large void regions are
             * zero before inpainting; after inpainting they get small
             * interpolated values that pull vmin very low, compressing the
             * colour scale so everything meaningful saturates to the top colour.
             *
             * Computing levels from the pre-inpainting data means vmin/vmax
             * reflect only real deposit values — the same pixels that carry
             * actual density information.  Inpainting then fills voids with
             * visually consistent colours without affecting the stretch.
             */
            if (cfg.rcfg.auto_levels && cfg.lock_levels && iter == 0) {
                auto_levels_from_data(&cfg.rcfg, global_data, npix_g);
                fprintf(stdout, "Levels locked (pre-inpaint): vmin=%.3f vmax=%.3f\n",
                        cfg.rcfg.vmin, cfg.rcfg.vmax);
                fflush(stdout);
                /* Freeze levels now — write_to_png_ex won't recompute them */
                cfg.rcfg.auto_levels = 0;
            }

            write_to_png_ex(image_file,
                            IMAGE_DIMENSIONX, IMAGE_DIMENSIONY,
                            global_data, &cfg.rcfg);

            /*
             * For non-lock_levels animation: levels are still auto-computed
             * inside write_to_png_ex per frame (auto_levels remains 1).
             * For lock_levels: already frozen above; write_to_png_ex uses
             * the fixed vmin/vmax without recomputing.
             */
        }


        /* Shrink the view box only when in zoom mode */
        if (cfg.zoom > 0) BoundingBox *= cfg.zoom_factor;
    }
    free(rx); free(ry); free(rz);

    free(data); free(global_data);
    free(x); free(y); free(z); free(ptype);
    if (smoothing_length) free(smoothing_length);
    sim_info_free(&header);

    fprintf(stdout, "Finished...\n");
    fflush(stdout);

#ifdef ENABLE_MPI
    MPI_Finalize();
#endif
    return 0;
}
