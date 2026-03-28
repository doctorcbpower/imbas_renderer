#ifndef CONFIG_H
#define CONFIG_H

/*
 * config.h — YAML configuration file support for render_image.
 *
 * load_yaml_config() reads a YAML file and merges its values into an
 * already-initialised cli_args_t.  It is designed to sit between
 * cli_args_default() and parse_args() so the precedence is:
 *
 *   compiled defaults  →  YAML file  →  command-line flags
 *
 * Every YAML key matches the equivalent CLI flag name (without '-').
 * Keys absent from the file are left at their current value.
 * Unknown keys print a warning and are ignored.
 *
 * Example YAML (see render.yaml for a fully-annotated template):
 *
 *   input:       snapshots/snap_100
 *   output:      frames/frame
 *   isHDF5:      true
 *   xc:          25.0
 *   yc:          25.0
 *   zc:          25.0
 *   lbox:        10.0
 *   itmax:       36
 *   rot_dangle:  10.0
 *   scene:       cluster
 *   colormap:    magma
 *   ptype:       [0, 1]
 */

#include "args.h"   /* cli_args_t — single definition, owned by args.h */

/*
 * load_yaml_config() — parse `path` and populate `args`.
 *
 * Returns  0 on success.
 * Returns -1 if the file cannot be opened or the YAML is malformed;
 *            a message is printed to stderr and args is left unchanged.
 */
int load_yaml_config(const char *path, cli_args_t *args);

#endif /* CONFIG_H */
