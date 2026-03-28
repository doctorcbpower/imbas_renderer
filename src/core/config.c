/*
 * config.c — YAML configuration file loader for render_image.
 *
 * Depends on libyaml (linked via the Makefile's YAML_OPTS).
 *
 * Design notes
 * ------------
 * - Walks a flat YAML key→value mapping via the libyaml event API.
 * - All keys match CLI flag names (minus the leading '-'), so the same
 *   mental model applies whether you are writing a config file or a
 *   shell script.
 * - Booleans accept: true/false, yes/no, on/off, 1/0.
 * - 'ptype' accepts a scalar ("1") or a sequence ([0, 1, 4]).
 * - Unknown keys warn and are skipped; the file never aborts the run.
 * - apply_scene() is defined in args.c and declared via args.h so the
 *   scene preset logic lives in exactly one place.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <yaml.h>

#include "config.h"     /* load_yaml_config() declaration  */
#include "args.h"       /* cli_args_t, apply_scene() via colormap.h chain */
#include "colormap.h"   /* render_config_parse_arg(), OPACITY_*, CMAP_*   */

/* ------------------------------------------------------------------ */
/* Internal helpers                                                     */
/* ------------------------------------------------------------------ */

static int yaml_bool(const char *val)
{
    return (strcmp(val, "true")  == 0 ||
            strcmp(val, "yes")   == 0 ||
            strcmp(val, "on")    == 0 ||
            strcmp(val, "1")     == 0);
}

/* Forward-declare the scene helper that lives in args.c */
void apply_scene_preset(cli_args_t *args, const char *scene);

/*
 * apply_kv() — write a single key=value pair into args.
 * Returns 0 on success, -1 if the key is unrecognised.
 */
static int apply_kv(cli_args_t *args, const char *key, const char *val)
{
    /* ---- I/O ---- */
    if      (strcmp(key, "input")  == 0)
        snprintf(args->file_root,       sizeof(args->file_root),       "%s", val);
    else if (strcmp(key, "output") == 0)
        snprintf(args->image_file_root, sizeof(args->image_file_root), "%s", val);

    /* ---- View / simulation ---- */
    else if (strcmp(key, "isHDF5")  == 0) args->isHDF5   = yaml_bool(val);
    else if (strcmp(key, "units")   == 0) args->boxunits = atoi(val);
    else if (strcmp(key, "xc")      == 0) args->xcen     = atof(val);
    else if (strcmp(key, "yc")      == 0) args->ycen     = atof(val);
    else if (strcmp(key, "zc")      == 0) args->zcen     = atof(val);
    else if (strcmp(key, "lbox")    == 0) args->lbox     = atof(val);
    else if (strcmp(key, "itmax")   == 0) args->itmax    = atoi(val);
    else if (strcmp(key, "zoom")    == 0) args->zoom     = atoi(val);
    else if (strcmp(key, "zoom_factor") == 0) {
        float zf = (float)atof(val);
        if (zf > 0.0f && zf < 1.0f)
            args->zoom_factor = zf;
        else
            fprintf(stderr,
                    "config: zoom_factor must be in (0,1), got %s — ignored\n", val);
    }

    /* ---- Rotation ---- */
    else if (strcmp(key, "rot_dangle") == 0)
        args->rot_dangle = (float)atof(val);
    else if (strcmp(key, "rot_axis") == 0) {
        float x = 0.0f, y = 0.0f, z = 0.0f;
        if (sscanf(val, "%f,%f,%f", &x, &y, &z) == 3) {
            float len = sqrtf(x*x + y*y + z*z);
            if (len > 0.0f) {
                args->rot_axis[0] = x / len;
                args->rot_axis[1] = y / len;
                args->rot_axis[2] = z / len;
            }
        } else {
            fprintf(stderr,
                    "config: rot_axis expects 'x,y,z', got '%s' — ignored\n", val);
        }
    }

    /* ---- Particle type ---- */
    else if (strcmp(key, "all_types") == 0) {
        if (yaml_bool(val)) args->ptype_mask = -1;
    }
    else if (strcmp(key, "gas") == 0) {
        if (yaml_bool(val)) {
            if (args->ptype_mask == (1 << 1)) args->ptype_mask = 0;
            args->ptype_mask |= (1 << 0);
        }
    }
    else if (strcmp(key, "dark_matter") == 0) {
        if (yaml_bool(val)) args->ptype_mask |= (1 << 1);
    }
    else if (strcmp(key, "stars") == 0) {
        if (yaml_bool(val)) {
            if (args->ptype_mask == (1 << 1)) args->ptype_mask = 0;
            args->ptype_mask |= (1 << 4);
        }
    }
    else if (strcmp(key, "ptype") == 0) {
        /* Scalar form; sequence items are handled in the event loop */
        int t = atoi(val);
        if (t >= 0 && t < (int)(sizeof(int) * 8)) {
            if (args->ptype_mask == (1 << 1)) args->ptype_mask = 0;
            args->ptype_mask |= (1 << t);
        }
    }

    /* ---- Scene preset ---- */
    else if (strcmp(key, "scene") == 0)
        apply_scene_preset(args, val);   /* defined in args.c */

    /* ---- Colour / opacity ---- */
    else if (strcmp(key, "colormap")          == 0) render_config_parse_arg(&args->rcfg, "colormap",          val);
    else if (strcmp(key, "reverse_colormap")  == 0) { if (yaml_bool(val)) args->rcfg.reverse = 1; }
    else if (strcmp(key, "opacity")           == 0) render_config_parse_arg(&args->rcfg, "opacity",           val);
    else if (strcmp(key, "opacity_func")      == 0) render_config_parse_arg(&args->rcfg, "opacity_func",      val);
    else if (strcmp(key, "opacity_gamma")     == 0) render_config_parse_arg(&args->rcfg, "opacity_gamma",     val);
    else if (strcmp(key, "opacity_threshold") == 0) render_config_parse_arg(&args->rcfg, "opacity_threshold", val);
    else if (strcmp(key, "vmin")              == 0) render_config_parse_arg(&args->rcfg, "vmin",              val);
    else if (strcmp(key, "vmax")              == 0) render_config_parse_arg(&args->rcfg, "vmax",              val);
    else if (strcmp(key, "linear_scale")      == 0) { if (yaml_bool(val)) args->rcfg.log_scale = 0; }
    else if (strcmp(key, "colormap_file")     == 0) render_config_parse_arg(&args->rcfg, "colormap_file",     val);
    else if (strcmp(key, "bg_color")          == 0) render_config_parse_arg(&args->rcfg, "bg_color",          val);

    /* ---- Auto-levels ---- */
    else if (strcmp(key, "no_auto_levels") == 0) { if (yaml_bool(val)) args->rcfg.auto_levels = 0; }
    else if (strcmp(key, "lock_levels")    == 0) args->lock_levels = yaml_bool(val);
    else if (strcmp(key, "auto_pct_lo")    == 0) render_config_parse_arg(&args->rcfg, "auto_pct_lo", val);
    else if (strcmp(key, "auto_pct_hi")    == 0) render_config_parse_arg(&args->rcfg, "auto_pct_hi", val);

    /* ---- Kernel smoothing ---- */
    else if (strcmp(key, "num_ngb")     == 0) args->num_ngb    = atoi(val);
    else if (strcmp(key, "sph_eta")     == 0) args->sph_eta    = (float)atof(val);
    else if (strcmp(key, "sph_cache")   == 0)
        snprintf(args->sph_cache, sizeof(args->sph_cache), "%s", val);
    else if (strcmp(key, "fast_smooth") == 0) args->fast_smooth = yaml_bool(val);

    /* ---- Time interpolation ---- */
    else if (strcmp(key, "interp_frac") == 0) args->interp_frac = (float)atof(val);
    else if (strcmp(key, "snap_dt")     == 0) args->snap_dt     = (float)atof(val);

    /* ---- 3-D grid ---- */
    else if (strcmp(key, "ngrid_z") == 0) args->ngrid_z = atoi(val);

    else {
        fprintf(stderr, "config: unknown key '%s' — ignored\n", key);
        return -1;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public API                                                           */
/* ------------------------------------------------------------------ */

int load_yaml_config(const char *path, cli_args_t *args)
{
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "config: cannot open '%s': ", path);
        perror(NULL);
        return -1;
    }

    yaml_parser_t parser;
    if (!yaml_parser_initialize(&parser)) {
        fprintf(stderr, "config: failed to initialise YAML parser\n");
        fclose(fp);
        return -1;
    }
    yaml_parser_set_input_file(&parser, fp);

    /*
     * State machine for a flat YAML mapping:
     *   STATE_KEY   — next scalar is a key
     *   STATE_VALUE — next scalar is the value for pending_key
     *   STATE_SEQ   — inside a sequence for pending_key (e.g. ptype: [0,1,4])
     */
    typedef enum { STATE_KEY, STATE_VALUE, STATE_SEQ } parse_state_t;
    parse_state_t state = STATE_KEY;
    char pending_key[128] = "";

    yaml_event_t ev;
    int ok = 1;

    while (ok) {
        if (!yaml_parser_parse(&parser, &ev)) {
            fprintf(stderr,
                    "config: YAML parse error in '%s' at line %zu\n",
                    path, parser.problem_mark.line + 1);
            ok = 0;
            break;
        }

        switch (ev.type) {

        case YAML_STREAM_END_EVENT:
        case YAML_DOCUMENT_END_EVENT:
            yaml_event_delete(&ev);
            goto done;

        case YAML_SCALAR_EVENT: {
            const char *val = (const char *)ev.data.scalar.value;
            if (state == STATE_KEY) {
                snprintf(pending_key, sizeof(pending_key), "%s", val);
                state = STATE_VALUE;
            } else if (state == STATE_VALUE) {
                apply_kv(args, pending_key, val);
                state = STATE_KEY;
            } else if (state == STATE_SEQ) {
                apply_kv(args, pending_key, val);   /* one item per call */
            }
            break;
        }

        case YAML_SEQUENCE_START_EVENT:
            if (state == STATE_VALUE) state = STATE_SEQ;
            break;

        case YAML_SEQUENCE_END_EVENT:
            state = STATE_KEY;
            break;

        case YAML_MAPPING_START_EVENT:
        case YAML_MAPPING_END_EVENT:
        case YAML_DOCUMENT_START_EVENT:
        case YAML_STREAM_START_EVENT:
        case YAML_NO_EVENT:
        default:
            break;
        }

        yaml_event_delete(&ev);
    }

done:
    yaml_parser_delete(&parser);
    fclose(fp);
    return ok ? 0 : -1;
}
