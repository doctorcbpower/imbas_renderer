#ifndef ARGS_H
#define ARGS_H

/*
 * args.h — CLI argument struct and parser for render_image.
 *
 * Provides cli_args_t (all tuneable parameters in one place),
 * cli_args_default() (fills in the same defaults previously scattered
 * through main()), and parse_args() (the extracted parse loop).
 *
 * Typical call sequence in main():
 *
 *   cli_args_t args;
 *   cli_args_default(&args);
 *   load_yaml_config(args.config_path, &args);  // YAML baseline (config.h)
 *   parse_args(argc, argv, &args);              // CLI wins
 */

#include "colormap.h"   /* render_config_t */

typedef struct {

    /* ----- Config file ----- */
    char config_path[512];      /* -config <file> */

    /* ----- I/O ----- */
    char file_root[256];        /* -input  <root> */
    char image_file_root[256];  /* -output <root> */
    int  isHDF5;                /* -isHDF5        */

    /* ----- View / simulation ----- */
    int    boxunits;            /* -units 0|1     */
    double xcen, ycen, zcen;   /* -xc -yc -zc    */
    double lbox;                /* -lbox          */
    int    itmax;               /* -itmax N       */
    int    zoom;                /* -zoom  N       */
    float  zoom_factor;         /* -zoom_factor f */

    /* ----- Rotation ----- */
    float rot_axis[3];          /* -rot_axis x,y,z */
    float rot_dangle;           /* -rot_dangle deg  */

    /* ----- Particle selection ----- */
    int ptype_mask;             /* -ptype / -gas / -dark_matter / -stars / -all_types */

    /* ----- Colour / opacity ----- */
    render_config_t rcfg;       /* -colormap, -opacity_*, -vmin/vmax, -scene, … */
    int lock_levels;            /* -lock_levels */

    /* ----- Kernel smoothing ----- */
    int   num_ngb;              /* -num_ngb N      */
    float sph_eta;              /* -sph_eta val    */
    char  sph_cache[512];       /* -sph_cache file */
    int   fast_smooth;          /* -fast_smooth    */

    /* ----- Time interpolation ----- */
    float interp_frac;          /* -interp_frac f */
    float snap_dt;              /* -snap_dt val   */

    /* ----- 3-D grid ----- */
    int ngrid_z;                /* -ngrid_z N */

} cli_args_t;

/*
 * cli_args_default() — populate *args with the defaults that were
 * previously scattered through main()'s variable declarations.
 * Always call this before parse_args() or load_yaml_config().
 */
void cli_args_default(cli_args_t *args);

/*
 * Scene preset helper
 */
void apply_scene_preset(cli_args_t *args, const char *scene);

/*
 * parse_args() — walk argv and populate *args.
 *
 * -config is stored in args->config_path but config loading is left
 * to the caller, preserving the explicit precedence order:
 *
 *   cli_args_default()  →  load_yaml_config()  →  parse_args()
 *
 * Unknown flags print a warning to stderr but do not abort.
 * Returns 0 on success.
 */
int parse_args(int argc, char *argv[], cli_args_t *args);

#endif // !ARGS_H
