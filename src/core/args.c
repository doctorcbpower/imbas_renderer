/*
 * args.c — CLI argument parsing for render_image.
 *
 * Contains:
 *   cli_args_default()  — canonical default values (previously in main())
 *   parse_args()        — the while(i < argc) loop (previously in main())
 *   args_unpack()       — copies cli_args_t back into main()'s flat locals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "args.h"
#include "colormap.h"   /* render_config_default(), render_config_parse_arg(),
                           OPACITY_* and CMAP_* constants                    */

/* ------------------------------------------------------------------ */
/* Defaults                                                             */
/* ------------------------------------------------------------------ */

void cli_args_default(cli_args_t *args)
{
    memset(args, 0, sizeof(*args));

    args->ptype_mask  = (1 << 1);   /* dark matter */
    args->xcen        = 0.5;
    args->ycen        = 0.5;
    args->zcen        = 0.5;
    args->itmax       = 1;
    args->zoom_factor = 0.5f;
    args->rot_axis[0] = 0.0f;
    args->rot_axis[1] = 0.0f;
    args->rot_axis[2] = 1.0f;       /* z-axis */
    args->num_ngb     = 32;
    args->sph_eta     = 1.2f;

    render_config_default(&args->rcfg);
}

/* ------------------------------------------------------------------ */
/* Scene preset helper (shared with config.c)                          */
/* ------------------------------------------------------------------ */

void apply_scene_preset(cli_args_t *args, const char *scene)
{
    if (strcmp(scene, "cluster") == 0) {
        /* One dominant halo: sqrt opacity compresses bright core
         * so outskirts remain visible. magma, tight lo clip. */
        args->rcfg.auto_pct_lo  = 0.05f;
        args->rcfg.auto_pct_hi  = 1.0f;
        args->rcfg.opacity_func = OPACITY_SQRT;
        args->rcfg.global_alpha = 1.0f;
        args->rcfg.colormap     = CMAP_MAGMA;
        args->rcfg.log_scale    = 1;
    } else if (strcmp(scene, "scattered") == 0) {
        /* Many halos at similar densities: flat opacity gives
         * equal visual weight to all. plasma, wider lo clip. */
        args->rcfg.auto_pct_lo  = 0.10f;
        args->rcfg.auto_pct_hi  = 1.0f;
        args->rcfg.opacity_func = OPACITY_FLAT;
        args->rcfg.global_alpha = 1.0f;
        args->rcfg.colormap     = CMAP_PLASMA;
        args->rcfg.log_scale    = 1;
    } else if (strcmp(scene, "filament") == 0) {
        /* Large-scale structure: power opacity (gamma=0.5)
         * compresses dense nodes so filaments stay visible.
         * inferno, very low lo clip to show faint sheets. */
        args->rcfg.auto_pct_lo   = 0.02f;
        args->rcfg.auto_pct_hi   = 1.0f;
        args->rcfg.opacity_func  = OPACITY_POWER;
        args->rcfg.opacity_gamma = 0.5f;
        args->rcfg.global_alpha  = 1.0f;
        args->rcfg.colormap      = CMAP_INFERNO;
        args->rcfg.log_scale     = 1;
    } else {
        fprintf(stderr,
                "Unknown -scene '%s'. Options: cluster, scattered, filament\n",
                scene);
    }
}

/* ------------------------------------------------------------------ */
/* Parser                                                               */
/* ------------------------------------------------------------------ */

int parse_args(int argc, char *argv[], cli_args_t *args)
{
    int i = 1;
    while (i < argc) {

        /* ---- Config file (recorded; loading is the caller's job) ---- */
        if (strcmp(argv[i], "-config") == 0)
            snprintf(args->config_path, sizeof(args->config_path),
                     "%s", argv[++i]);

        /* ---- I/O ---- */
        else if (strcmp(argv[i], "-input")  == 0)
            snprintf(args->file_root, sizeof(args->file_root),
                     "%s", argv[++i]);
        else if (strcmp(argv[i], "-output") == 0)
            snprintf(args->image_file_root, sizeof(args->image_file_root),
                     "%s", argv[++i]);
        else if (strcmp(argv[i], "-isHDF5") == 0)
            args->isHDF5 = 1;

        /* ---- View / simulation ---- */
        else if (strcmp(argv[i], "-units") == 0) args->boxunits = atoi(argv[++i]);
        else if (strcmp(argv[i], "-xc")    == 0) args->xcen     = atof(argv[++i]);
        else if (strcmp(argv[i], "-yc")    == 0) args->ycen     = atof(argv[++i]);
        else if (strcmp(argv[i], "-zc")    == 0) args->zcen     = atof(argv[++i]);
        else if (strcmp(argv[i], "-lbox")  == 0) args->lbox     = atof(argv[++i]);
        else if (strcmp(argv[i], "-itmax") == 0) args->itmax    = atoi(argv[++i]);
        else if (strcmp(argv[i], "-zoom")  == 0) args->zoom     = atoi(argv[++i]);
        else if (strcmp(argv[i], "-zoom_factor") == 0) {
            float zf = (float)atof(argv[++i]);
            if (zf <= 0.0f || zf >= 1.0f) {
                fprintf(stderr,
                        "zoom_factor must be in (0,1), got %g — using 0.5\n", zf);
                zf = 0.5f;
            }
            args->zoom_factor = zf;
        }

        /* ---- Rotation ---- */
        else if (strcmp(argv[i], "-rot_dangle") == 0)
            args->rot_dangle = (float)atof(argv[++i]);
        else if (strcmp(argv[i], "-rot_axis") == 0) {
            float x = 0.0f, y = 0.0f, z = 0.0f;
            sscanf(argv[++i], "%f,%f,%f", &x, &y, &z);
            float len = sqrtf(x*x + y*y + z*z);
            if (len > 0.0f) { x /= len; y /= len; z /= len; }
            args->rot_axis[0] = x;
            args->rot_axis[1] = y;
            args->rot_axis[2] = z;
        }

        /* ---- Particle type selection ---- */
        else if (strcmp(argv[i], "-all_types") == 0)
            args->ptype_mask = -1;
        else if (strcmp(argv[i], "-gas") == 0) {
            if (args->ptype_mask == (1 << 1)) args->ptype_mask = 0;
            args->ptype_mask |= (1 << 0);
        }
        else if (strcmp(argv[i], "-dark_matter") == 0)
            args->ptype_mask |= (1 << 1);
        else if (strcmp(argv[i], "-stars") == 0) {
            if (args->ptype_mask == (1 << 1)) args->ptype_mask = 0;
            args->ptype_mask |= (1 << 4);
        }
        else if (strcmp(argv[i], "-ptype") == 0) {
            int t = atoi(argv[++i]);
            if (t >= 0 && t < (int)(sizeof(int) * 8)) {
                if (args->ptype_mask == (1 << 1)) args->ptype_mask = 0;
                args->ptype_mask |= (1 << t);
            }
        }

        /* ---- Scene presets ---- */
        else if (strcmp(argv[i], "-scene") == 0)
            apply_scene_preset(args, argv[++i]);

        /* ---- Colour / opacity ---- */
        else if (strcmp(argv[i], "-colormap") == 0)
            render_config_parse_arg(&args->rcfg, "colormap", argv[++i]);
        else if (strcmp(argv[i], "-reverse_colormap") == 0)
            args->rcfg.reverse = 1;
        else if (strcmp(argv[i], "-opacity") == 0)
            render_config_parse_arg(&args->rcfg, "opacity", argv[++i]);
        else if (strcmp(argv[i], "-opacity_func") == 0)
            render_config_parse_arg(&args->rcfg, "opacity_func", argv[++i]);
        else if (strcmp(argv[i], "-opacity_gamma") == 0)
            render_config_parse_arg(&args->rcfg, "opacity_gamma", argv[++i]);
        else if (strcmp(argv[i], "-opacity_threshold") == 0)
            render_config_parse_arg(&args->rcfg, "opacity_threshold", argv[++i]);
        else if (strcmp(argv[i], "-vmin") == 0)
            render_config_parse_arg(&args->rcfg, "vmin", argv[++i]);
        else if (strcmp(argv[i], "-vmax") == 0)
            render_config_parse_arg(&args->rcfg, "vmax", argv[++i]);
        else if (strcmp(argv[i], "-linear_scale") == 0)
            args->rcfg.log_scale = 0;
        else if (strcmp(argv[i], "-colormap_file") == 0)
            render_config_parse_arg(&args->rcfg, "colormap_file", argv[++i]);
        else if (strcmp(argv[i], "-bg_color") == 0)
            render_config_parse_arg(&args->rcfg, "bg_color", argv[++i]);

        /* ---- Auto-levels ---- */
        else if (strcmp(argv[i], "-no_auto_levels") == 0)
            args->rcfg.auto_levels = 0;
        else if (strcmp(argv[i], "-lock_levels") == 0)
            args->lock_levels = 1;
        else if (strcmp(argv[i], "-auto_pct_lo") == 0)
            render_config_parse_arg(&args->rcfg, "auto_pct_lo", argv[++i]);
        else if (strcmp(argv[i], "-auto_pct_hi") == 0)
            render_config_parse_arg(&args->rcfg, "auto_pct_hi", argv[++i]);

        /* ---- Kernel smoothing ---- */
        else if (strcmp(argv[i], "-num_ngb") == 0)
            args->num_ngb = atoi(argv[++i]);
        else if (strcmp(argv[i], "-sph_eta") == 0)
            args->sph_eta = (float)atof(argv[++i]);
        else if (strcmp(argv[i], "-sph_cache") == 0)
            snprintf(args->sph_cache, sizeof(args->sph_cache),
                     "%s", argv[++i]);
        else if (strcmp(argv[i], "-fast_smooth") == 0)
            args->fast_smooth = 1;

        /* ---- Time interpolation ---- */
        else if (strcmp(argv[i], "-interp_frac") == 0)
            args->interp_frac = (float)atof(argv[++i]);
        else if (strcmp(argv[i], "-snap_dt") == 0)
            args->snap_dt = (float)atof(argv[++i]);

        /* ---- 3-D grid ---- */
        else if (strcmp(argv[i], "-ngrid_z") == 0)
            args->ngrid_z = atoi(argv[++i]);

        else
            fprintf(stderr, "Warning: unknown argument '%s' — ignored\n", argv[i]);

        i++;
    }

    return 0;
}
