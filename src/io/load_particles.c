/*
 * load_particles.c — snapshot reading and position interpolation.
 *
 * Contains:
 *   load_snapshot_header()      — filename resolution, header read, boxunits scaling
 *   load_particles()            — single-snapshot read + linear extrapolation (legacy)
 *   load_particles_hermite()    — 2-snapshot cubic Hermite, index-matched,
 *                                 finite-difference tangents, periodic wrapping
 *   load_particles_hermite3()   — 3 or 4-snapshot Catmull-Rom Hermite, ID-matched
 *                                 via radix sort, periodic wrapping
 *
 * Hermite scheme (2-snapshot)
 * ---------------------------
 * Tangent vectors are derived from the finite difference of the two snapshot
 * positions rather than from the stored velocity field.  Velocities in SWIFT
 * outputs are peculiar velocities in km/s; converting them to comoving position
 * units requires knowing H(a) and the internal time unit, and the Hubble drag
 * term anti-correlates with the peculiar displacement making a naive conversion
 * produce negative effective dt.  The finite-difference approach is unit-agnostic.
 *
 *   tangent[k] = wrap(xB[k] - xA[k])   — periodically wrapped displacement
 *
 * With dt = 1 the Hermite formula is:
 *   p(t) = h00*xA + h10*tangent + h01*xB + h11*tangent
 *
 * Catmull-Rom scheme (3 or 4-snapshot)
 * -------------------------------------
 * Uses the Catmull-Rom centred-difference tangents for smoother motion:
 *
 *   m0 = wrap(xA - xPrev) / 2          — tangent at left endpoint
 *   m1 = wrap(xNext - xA) / 2          — tangent at right endpoint (4-snap)
 *      = wrap(xB   - xA)  / 2          — one-sided backward diff (3-snap only)
 *
 * Particles are matched by 64-bit ParticleID using a radix sort, making this
 * robust to the domain-decomposition reordering that SWIFT applies between
 * snapshots.  Only DM (type 1) particles are fully conserved; gas/star counts
 * vary.  The ptype_mask in cfg selects which types are rendered.
 *
 * Periodic boundary handling
 * --------------------------
 * Raw displacements xB-xA can be ~BoxSize for particles that crossed a
 * periodic boundary.  All displacements are wrapped into (-L/2, +L/2)
 * before use as tangents.  Interpolated positions are wrapped back into
 * [0, BoxSize) after the Hermite evaluation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "load_particles.h"
#include "args.h"
#include "io.h"
#include "globals.h"   /* ThisTask */

/* ================================================================== */
/* Internal helpers                                                     */
/* ================================================================== */

static int alloc_particle_arrays(long long n,
                                 float **x, float **y, float **z, int **ptype)
{
    *x     = (float *)malloc(sizeof(float) * n);
    *y     = (float *)malloc(sizeof(float) * n);
    *z     = (float *)malloc(sizeof(float) * n);
    *ptype = (int   *)malloc(sizeof(int)   * n);
    if (!*x || !*y || !*z || !*ptype) {
        fprintf(stderr, "load_particles: malloc failed for position arrays\n");
        free(*x); free(*y); free(*z); free(*ptype);
        *x = *y = *z = NULL; *ptype = NULL;
        return -1;
    }
    return 0;
}

static int alloc_velocity_arrays(long long n,
                                 float **vx, float **vy, float **vz)
{
    *vx = (float *)malloc(sizeof(float) * n);
    *vy = (float *)malloc(sizeof(float) * n);
    *vz = (float *)malloc(sizeof(float) * n);
    if (!*vx || !*vy || !*vz) {
        fprintf(stderr, "load_particles: malloc failed for velocity arrays\n");
        free(*vx); free(*vy); free(*vz);
        *vx = *vy = *vz = NULL;
        return -1;
    }
    return 0;
}

static long long count_header_particles(const struct sim_info *header)
{
    long long n = 0;
    for (int t = 0; t < header->num_types; t++)
        n += header->nall[t];
    return n;
}

/* Wrap position into [0, L) */
static inline float wrap_pos(float p, float L)
{
    if (p <  0.0f) p += L;
    if (p >= L)    p -= L;
    return p;
}

/* Wrap displacement into (-L/2, +L/2] */
static inline float wrap_disp(float d, float half, float L)
{
    if (d >  half) d -= L;
    if (d < -half) d += L;
    return d;
}

/* ================================================================== */
/* Radix sort helpers for 64-bit ParticleID matching                   */
/* ================================================================== */

/*
 * id_index_t — pairs a ParticleID with its original array index.
 * We sort an array of these by ID, then use the stored index to
 * reorder the position arrays — avoids touching the (large) float
 * arrays during the sort itself.
 */
typedef struct {
    uint64_t id;
    uint32_t idx;   /* original index; safe for up to 4G particles */
} id_index_t;

/*
 * radix_sort_id_index() — stable 8-pass radix sort on the 64-bit id field.
 *
 * Sorts arr[0..n-1] in ascending order of arr[i].id.
 * Uses a temporary buffer of the same size (allocated internally).
 * Returns 0 on success, -1 on malloc failure.
 *
 * Cost: O(8*n) passes, ~8n bytes working space.
 * For n=16M: ~128M iterations, ~128 MB working buffer, ~0.5 s.
 */
static int radix_sort_id_index(id_index_t *arr, long long n)
{
    id_index_t *tmp = (id_index_t *)malloc(sizeof(id_index_t) * n);
    if (!tmp) {
        fprintf(stderr, "radix_sort: malloc failed for %lld elements\n", n);
        return -1;
    }

    /* 8 passes of 8-bit radix sort (8 * 8 = 64 bits) */
    for (int pass = 0; pass < 8; pass++) {
        int shift = pass * 8;

        /* Count occurrences of each byte value */
        long long count[256] = {0};
        for (long long i = 0; i < n; i++)
            count[(arr[i].id >> shift) & 0xFF]++;

        /* Prefix sum → starting positions */
        long long pos[256];
        pos[0] = 0;
        for (int b = 1; b < 256; b++)
            pos[b] = pos[b-1] + count[b-1];

        /* Scatter into tmp */
        for (long long i = 0; i < n; i++) {
            int bucket = (arr[i].id >> shift) & 0xFF;
            tmp[pos[bucket]++] = arr[i];
        }

        /* Swap arr and tmp for next pass */
        id_index_t *swap = arr;
        arr = tmp;
        tmp = swap;

        /*
         * After an even number of swaps the result is back in the
         * original arr pointer.  After odd passes it is in tmp
         * (which now points to what was arr).  We track this with
         * the pass counter and do a final copy if needed.
         */
    }

    /* After 8 (even) passes the result is in the original arr buffer.
     * tmp now points to the working buffer — free it. */
    free(tmp);
    return 0;
}

/*
 * reorder_by_id() — read positions for a snapshot, then reorder them
 * to match the ID ordering of the reference snapshot (snap_a / middle).
 *
 * Parameters:
 *   snap_root    — HDF5 file root path
 *   num_files    — number of distributed HDF5 files
 *   ref_ids      — sorted id_index_t array from the reference snapshot
 *   n_ref        — number of particles in reference
 *   BoxSize      — for periodic position wrapping of input positions
 *   x_out etc.   — pre-allocated output arrays, length n_ref
 *   n_matched    — set to the number of IDs found in both snapshots
 *
 * Particles present in the reference but absent in this snapshot
 * (e.g. converted gas particles) are left at position (0,0,0) and
 * ptype -1.  The caller should filter by ptype_mask.
 *
 * Returns 0 on success, -1 on allocation failure.
 */
static int reorder_by_id(const char    *snap_root,
                         int            num_files,
                         const id_index_t *ref_sorted,  /* sorted by ID */
                         long long      n_ref,
                         float         *x_out,
                         float         *y_out,
                         float         *z_out,
                         int           *pt_out,
                         long long     *n_matched)
{
    /* ---- Read this snapshot's positions and IDs ---- */
    long long nRead = 0;

    /* Allocate temporary position + ID arrays for this snapshot */
    float    *xR = (float    *)malloc(sizeof(float)    * n_ref);
    float    *yR = (float    *)malloc(sizeof(float)    * n_ref);
    float    *zR = (float    *)malloc(sizeof(float)    * n_ref);
    int      *ptR= (int      *)malloc(sizeof(int)      * n_ref);
    uint64_t *idR= (uint64_t *)malloc(sizeof(uint64_t) * n_ref);

    if (!xR || !yR || !zR || !ptR || !idR) {
        fprintf(stderr, "reorder_by_id: malloc failed\n");
        free(xR); free(yR); free(zR); free(ptR); free(idR);
        return -1;
    }

    /*
     * read_particles_from_hdf5_ids_only() reads Coordinates + ParticleIDs
     * only — no velocities — which is all reorder_by_id needs.
     */
    read_particles_from_hdf5(snap_root,
                             xR, yR, zR,
                             (unsigned long long *)idR,
                             ptR,
                             -1,      /* ptype_mask */
                             &nRead);

    /* ---- Sort this snapshot by ID ---- */
    id_index_t *sorted = (id_index_t *)malloc(sizeof(id_index_t) * nRead);
    if (!sorted) {
        fprintf(stderr, "reorder_by_id: malloc failed for sort buffer\n");
        free(xR); free(yR); free(zR); free(ptR); free(idR);
        return -1;
    }
    for (long long i = 0; i < nRead; i++) {
        sorted[i].id  = idR[i];
        sorted[i].idx = (uint32_t)i;
    }
    free(idR);   /* no longer needed after packing into sorted[] */

    if (radix_sort_id_index(sorted, nRead) != 0) {
        free(xR); free(yR); free(zR); free(ptR); free(sorted);
        return -1;
    }

    /* ---- Match against reference IDs and scatter to output ---- */
    /*
     * Both ref_sorted and sorted are now sorted by ID.
     * Walk them with two pointers — O(n_ref + nRead).
     */
    *n_matched = 0;

    /* Initialise output to sentinel values for unmatched particles */
    memset(x_out,  0, sizeof(float) * n_ref);
    memset(y_out,  0, sizeof(float) * n_ref);
    memset(z_out,  0, sizeof(float) * n_ref);
    for (long long i = 0; i < n_ref; i++) pt_out[i] = -1;

    long long j = 0;   /* pointer into sorted (this snapshot) */
    for (long long i = 0; i < n_ref && j < nRead; i++) {
        uint64_t ref_id = ref_sorted[i].id;

        /* Advance j until sorted[j].id >= ref_id */
        while (j < nRead && sorted[j].id < ref_id) j++;

        if (j < nRead && sorted[j].id == ref_id) {
            uint32_t src = sorted[j].idx;
            /* ref_sorted[i].idx is the output position in the reference ordering */
            uint32_t dst = ref_sorted[i].idx;
            x_out[dst]  = xR[src];
            y_out[dst]  = yR[src];
            z_out[dst]  = zR[src];
            pt_out[dst] = ptR[src];
            (*n_matched)++;
            j++;
        }
        /* else: particle i exists in ref but not in this snap — leave as 0/-1 */
    }

    free(xR); free(yR); free(zR); free(ptR); free(sorted);
    return 0;
}

/* ================================================================== */
/* load_snapshot_header                                                 */
/* ================================================================== */

int load_snapshot_header(cli_args_t      *cfg,
                         struct sim_info *header,
                         long long       *NumPart,
                         char            *filename_out,
                         int             *isDistributed_out)
{
    const char *root = hermite_mode(cfg) ? cfg->snap_a : cfg->file_root;
    check_input_filenames(filename_out, root, cfg->isHDF5, isDistributed_out);

    if (ThisTask == 0) {
        if (cfg->xcen > 0 || cfg->lbox > 0) {
            fprintf(stdout, "Centre: (%g|%g|%g)\n",
                    cfg->xcen, cfg->ycen, cfg->zcen);
            fprintf(stdout, "Box Length: %g\n", cfg->lbox);
        }
        fprintf(stdout, "Render configuration:\n");
        render_config_print(&cfg->rcfg);
        if (hermite3_mode(cfg))
            fprintf(stdout,
                    "Interpolation: Catmull-Rom Hermite\n"
                    "  snap_prev = %s\n  snap_a    = %s\n"
                    "  snap_b    = %s\n  snap_next = %s\n",
                    cfg->snap_prev,
                    cfg->snap_a, cfg->snap_b,
                    cfg->snap_next[0] ? cfg->snap_next : "(none — one-sided)");
        else if (hermite_mode(cfg))
            fprintf(stdout, "Interpolation: Hermite (snap_a=%s  snap_b=%s)\n",
                    cfg->snap_a, cfg->snap_b);
        fprintf(stdout, "Reading header...\n");
        fflush(stdout);
    }

    if (cfg->isHDF5)
        read_hdf5_header(filename_out, header, NumPart);
    else
        read_gadget_binary_header(filename_out, header, NumPart);

    if (ThisTask == 0) {
        fprintf(stdout, "Number of particle types in file: %d\n",
                header->num_types);
        for (int i = 0; i < (int)(sizeof(int) * 8); i++) {
            if (cfg->ptype_mask != -1 && !(cfg->ptype_mask & (1 << i)))
                continue;
            if (i >= header->num_types || header->nall[i] == 0)
                fprintf(stdout,
                        "Warning: type %d requested but not present\n", i);
        }
        fflush(stdout);
    }

    if (header->BoxSize == 0)
        header->BoxSize = 1.e6;
    else if (cfg->boxunits == 1) {
        cfg->xcen *= header->BoxSize;
        cfg->ycen *= header->BoxSize;
        cfg->zcen *= header->BoxSize;
        cfg->lbox *= header->BoxSize;
    }

    if (ThisTask == 0) {
        fprintf(stdout, "Number of files: %d\n",       header->NumFiles);
        fprintf(stdout, "Number of particles: %lld\n", *NumPart);
        fprintf(stdout, "Time/Expansion Factor: %g\n", header->time);
        fflush(stdout);
    }
    return 0;
}

/* ================================================================== */
/* load_particles — single-snapshot, linear extrapolation (legacy)     */
/* ================================================================== */

int load_particles(const cli_args_t      *cfg,
                   const struct sim_info  *header,
                   const char             *filename,
                   float                **x,
                   float                **y,
                   float                **z,
                   int                  **ptype,
                   long long             *NumPartRead_out)
{
    long long NumPart = count_header_particles(header);

    if (alloc_particle_arrays(NumPart, x, y, z, ptype) != 0)
        return -1;

    int do_interp = (cfg->interp_frac != 0.0f && cfg->isHDF5);
    float *vx = NULL, *vy = NULL, *vz = NULL;
    if (do_interp) {
        if (alloc_velocity_arrays(NumPart, &vx, &vy, &vz) != 0) {
            free(*x); free(*y); free(*z); free(*ptype);
            return -1;
        }
    }

    if (ThisTask == 0) { fprintf(stdout, "Reading particles...\n"); fflush(stdout); }

    if (cfg->isHDF5)
        read_particles_from_hdf5(cfg->file_root, *x, *y, *z,
                                  NULL,           /* partid not needed here */
                                  *ptype,
                                  header->NumFiles, NumPartRead_out);
    else
        read_particles_from_gadget_binary(cfg->file_root, *x, *y, *z, *ptype,
                                           header->NumFiles, NumPartRead_out);

    fprintf(stdout, "NumPart: %llu\tNumPartRead: %llu\n",
            NumPart, *NumPartRead_out);
    fflush(stdout);

    if (do_interp) {
        float scale = (cfg->snap_dt > 0.0f)
                    ? cfg->interp_frac * cfg->snap_dt
                    : cfg->interp_frac;
        fprintf(stdout, "Linear interpolation: frac=%.4f  dt=%.4g  scale=%.4g\n",
                cfg->interp_frac, cfg->snap_dt, scale);
        fflush(stdout);
        for (long long k = 0; k < *NumPartRead_out; k++) {
            (*x)[k] += scale * vx[k];
            (*y)[k] += scale * vy[k];
            (*z)[k] += scale * vz[k];
        }
        free(vx); free(vy); free(vz);
    }
    return 0;
}

/* ================================================================== */
/* load_particles_hermite — 2-snapshot, index-matched                  */
/* ================================================================== */

int load_particles_hermite(const cli_args_t      *cfg,
                           const struct sim_info  *header,
                           float                **x,
                           float                **y,
                           float                **z,
                           int                  **ptype,
                           long long             *NumPartRead_out)
{
    long long NumPart = count_header_particles(header);

    float *xA = NULL, *yA = NULL, *zA = NULL;
    float *xB = NULL, *yB = NULL, *zB = NULL;
    float *vxA = NULL, *vyA = NULL, *vzA = NULL;
    float *vxB = NULL, *vyB = NULL, *vzB = NULL;
    int   *ptA = NULL, *ptB = NULL;

    if (alloc_particle_arrays(NumPart, &xA, &yA, &zA, &ptA) != 0) return -1;
    if (alloc_velocity_arrays(NumPart, &vxA, &vyA, &vzA) != 0) {
        free(xA); free(yA); free(zA); free(ptA); return -1;
    }
    if (alloc_particle_arrays(NumPart, &xB, &yB, &zB, &ptB) != 0) {
        free(xA); free(yA); free(zA); free(ptA);
        free(vxA); free(vyA); free(vzA); return -1;
    }
    if (alloc_velocity_arrays(NumPart, &vxB, &vyB, &vzB) != 0) {
        free(xA); free(yA); free(zA); free(ptA);
        free(vxA); free(vyA); free(vzA);
        free(xB); free(yB); free(zB); free(ptB); return -1;
    }
    if (alloc_particle_arrays(NumPart, x, y, z, ptype) != 0) {
        free(xA); free(yA); free(zA); free(ptA);
        free(vxA); free(vyA); free(vzA);
        free(xB); free(yB); free(zB); free(ptB); return -1;
    }
    // First read snapshot A
    long long nReadA = 0, nReadB = 0;

    uint64_t *idA = (uint64_t *)malloc(sizeof(uint64_t) * NumPart);
    if (!idA) {
        fprintf(stderr, "hermite3: malloc failed for idA\n");
        return -1;
    }
    
    if (ThisTask == 0) {
        fprintf(stdout, "Reading snapshot A: %s\n", cfg->snap_a); fflush(stdout);
    }
    read_particles_from_hdf5(cfg->snap_a, xA, yA, zA,
                              (unsigned long long *)idA,
                              ptA,
                              header->NumFiles, &nReadA);

    /* Build and sort the reference id_index array */
    id_index_t *ref_sorted = (id_index_t *)malloc(sizeof(id_index_t) * nReadA);
    if (!ref_sorted) {
        fprintf(stderr, "hermite3: malloc failed for ref_sorted\n");
        free(idA); return -1;
    }
    for (long long i = 0; i < nReadA; i++) {
        ref_sorted[i].id  = idA[i];
        ref_sorted[i].idx = (uint32_t)i;
    }
    free(idA);
    
    if (ThisTask == 0) {
        fprintf(stdout, "Sorting %lld particles by ID...\n", nReadA);
        fflush(stdout);
    }
    if (radix_sort_id_index(ref_sorted, nReadA) != 0) {
        free(ref_sorted); return -1;
    }

    // Now do snapshot B
    struct sim_info headerB;
    memset(&headerB, 0, sizeof(headerB));
    long long NumPartB = 0;
    char filenameB[256] = "";
    int isDistB = 0;
    check_input_filenames(filenameB, cfg->snap_b, cfg->isHDF5, &isDistB);
    read_hdf5_header(filenameB, &headerB, &NumPartB);

    if (NumPartB != NumPart)
        fprintf(stderr,
                "load_particles_hermite: particle count mismatch — "
                "snap_a has %lld, snap_b has %lld\n", NumPart, NumPartB);

    if (ThisTask == 0) {
        fprintf(stdout, "Reading snapshot B: %s\n", cfg->snap_b); fflush(stdout);
    }
    
    long long nMatchB = 0;
    if (reorder_by_id(cfg->snap_b, headerB.NumFiles,
                      ref_sorted, nReadA,
                      xB, yB, zB, ptB, &nMatchB) != 0) {
        free(ref_sorted); return -1;
    }
    
    //*NumPartRead_out = (nReadA < nReadB) ? nReadA : nReadB;

    *NumPartRead_out = nReadA;
    
    if (ThisTask == 0) {
        fprintf(stdout, "  matched %lld / %lld particles from B\n",
                nMatchB, nReadA);
        fflush(stdout);
    }
    
    if (ThisTask == 0) {
        fprintf(stdout, "snap_a read: %lld  snap_b read: %lld\n", nReadA, nMatchB);
        fflush(stdout);
    }

    /* Finite-difference tangents with periodic wrapping */
    float BoxSize = (float)header->BoxSize;
    float half    = 0.5f * BoxSize;

    /* Hermite basis polynomials */
    float t  = cfg->interp_frac;
    float t2 = t * t, t3 = t2 * t;
    float h00 =  2.0f*t3 - 3.0f*t2 + 1.0f;
    float h10 =       t3 - 2.0f*t2 + t;
    float h01 = -2.0f*t3 + 3.0f*t2;
    float h11 =       t3 -       t2;
    
    if (ThisTask == 0) {
        fprintf(stdout,
                "Hermite: t=%.4f  h00=%.4f h10=%.4f h01=%.4f h11=%.4f\n",
                t, h00, h10, h01, h11);
        fflush(stdout);
    }
    
    for (long long k = 0; k < *NumPartRead_out; k++) {
        if (ptB[k] == -1) {
            (*x)[k] = xA[k];
            (*y)[k] = yA[k];
            (*z)[k] = zA[k];
            (*ptype)[k] = ptA[k];
            continue;
        }
        
        float dx = wrap_disp(xB[k] - xA[k], half, BoxSize);
        float xB_local = xA[k] + dx;
        float px = h00*xA[k] + h10*dx + h01*xB_local + h11*dx;
        (*x)[k] = wrap_pos(px, BoxSize);
        
        float dy = wrap_disp(yB[k] - yA[k], half, BoxSize);
        float yB_local = yA[k] + dy;
        float py = h00*yA[k] + h10*dy + h01*yB_local + h11*dy;
        (*y)[k] = wrap_pos(py, BoxSize);
        
        float dz = wrap_disp(zB[k] - zA[k], half, BoxSize);
        float zB_local = zA[k] + dz;
        float pz = h00*zA[k] + h10*dz + h01*zB_local + h11*dz;
        (*z)[k] = wrap_pos(pz, BoxSize);
        
        (*ptype)[k] = ptA[k];
    }

    free(xA);  free(yA);  free(zA);  free(ptA);
    free(vxA); free(vyA); free(vzA);
    free(xB);  free(yB);  free(zB);  free(ptB);
    free(vxB); free(vyB); free(vzB);
    sim_info_free(&headerB);
    return 0;
}

/* ================================================================== */
/* load_particles_hermite3 — 3 or 4-snapshot Catmull-Rom, ID-matched  */
/* ================================================================== */

/*
 * read_snap_header_only() — read header for an arbitrary snapshot path.
 * Used to get NumFiles for snap_prev / snap_next.
 */
static int read_snap_header_only(const char      *root,
                                 int              isHDF5,
                                 struct sim_info *hdr,
                                 long long       *npart)
{
    char fname[512] = "";
    int  isDist = 0;
    check_input_filenames(fname, root, isHDF5, &isDist);
    memset(hdr, 0, sizeof(*hdr));
    *npart = 0;
    if (isHDF5)
        read_hdf5_header(fname, hdr, npart);
    else
        read_gadget_binary_header(fname, hdr, npart);
    return 0;
}

int load_particles_hermite3(const cli_args_t      *cfg,
                            const struct sim_info  *header,
                            float                **x,
                            float                **y,
                            float                **z,
                            int                  **ptype,
                            long long             *NumPartRead_out)
{
    /*
     * Snapshots:
     *   P = snap_prev  (always required)
     *   A = snap_a     (left  endpoint of interpolation interval)
     *   B = snap_b     (right endpoint of interpolation interval)
     *   N = snap_next  (optional; improves tangent at B)
     *
     * Tangents (Catmull-Rom):
     *   m0 = wrap(xA - xP) / 2         at t=0 (left)
     *   m1 = wrap(xN - xA) / 2         at t=1 (right, 4-snap)
     *      = wrap(xB - xA) / 2         at t=1 (right, 3-snap fallback)
     *
     * Memory strategy:
     *   1. Read A positions + IDs → sort A by ID → this is the reference order.
     *   2. Read P positions + IDs → match to A order → free P IDs.
     *   3. Read B positions + IDs → match to A order → free B IDs.
     *   4. If snap_next: read N positions + IDs → match → free N IDs.
     *   5. Compute tangents (in-place into reused arrays).
     *   6. Hermite loop → output.
     *
     * Peak memory: 4 position arrays (xA,yA,zA, xP,yP,zP, xB,yB,zB, [xN,yN,zN])
     *   + 1 ID array for the snapshot being read
     *   + sort buffer (same size as ID array)
     *   + ref_sorted array (id_index_t, 12 bytes/particle)
     * For 16M particles: 4×3×16M×4 + 16M×8 + 16M×12 ≈ 896 MB.
     * With snap_next: add another 3×16M×4 = 192 MB → ~1.1 GB total.
     * On memory-constrained systems omit snap_next (3-snap mode).
     */

    long long NumPart = count_header_particles(header);
    float     BoxSize = (float)header->BoxSize;
    float     half    = 0.5f * BoxSize;
    int       has_next = (cfg->snap_next[0] != '\0');

    if (ThisTask == 0) {
        fprintf(stdout,
                "Catmull-Rom Hermite3: %s mode\n"
                "  P=%s\n  A=%s\n  B=%s\n  N=%s\n",
                has_next ? "4-snapshot" : "3-snapshot",
                cfg->snap_prev, cfg->snap_a, cfg->snap_b,
                has_next ? cfg->snap_next : "(none)");
        fflush(stdout);
    }

    /* ---- Allocate position arrays for all snapshots ---- */
    float *xP = NULL, *yP = NULL, *zP = NULL;  int *ptP = NULL;
    float *xA = NULL, *yA = NULL, *zA = NULL;  int *ptA = NULL;
    float *xB = NULL, *yB = NULL, *zB = NULL;  int *ptB = NULL;
    float *xN = NULL, *yN = NULL, *zN = NULL;  int *ptN = NULL;

    /* Macro to free everything on error */
#define FREE_ALL() do { \
    free(xP); free(yP); free(zP); free(ptP); \
    free(xA); free(yA); free(zA); free(ptA); \
    free(xB); free(yB); free(zB); free(ptB); \
    free(xN); free(yN); free(zN); free(ptN); \
    } while(0)

    if (alloc_particle_arrays(NumPart, &xA, &yA, &zA, &ptA) != 0) return -1;
    if (alloc_particle_arrays(NumPart, &xP, &yP, &zP, &ptP) != 0) { FREE_ALL(); return -1; }
    if (alloc_particle_arrays(NumPart, &xB, &yB, &zB, &ptB) != 0) { FREE_ALL(); return -1; }
    if (has_next) {
        if (alloc_particle_arrays(NumPart, &xN, &yN, &zN, &ptN) != 0) { FREE_ALL(); return -1; }
    }

    /* Output arrays */
    if (alloc_particle_arrays(NumPart, x, y, z, ptype) != 0) { FREE_ALL(); return -1; }

    /* ----------------------------------------------------------------
     * Step 1: Read snapshot A + IDs, sort by ID → reference ordering
     * ---------------------------------------------------------------- */
    if (ThisTask == 0) {
        fprintf(stdout, "Reading snapshot A (reference): %s\n", cfg->snap_a);
        fflush(stdout);
    }

    uint64_t *idA = (uint64_t *)malloc(sizeof(uint64_t) * NumPart);
    if (!idA) {
        fprintf(stderr, "hermite3: malloc failed for idA\n");
        FREE_ALL(); return -1;
    }

    long long nReadA = 0;
    
    read_particles_from_hdf5(cfg->snap_a,
                             xA, yA, zA,
                             (unsigned long long *)idA,  /* Add this cast */
                             ptA,
                             header->NumFiles, &nReadA);
//    read_particles_from_hdf5(cfg->snap_a,
//                             xA, yA, zA,
//                             idA, ptA,
//                             header->NumFiles, &nReadA);

    /* Build and sort the reference id_index array */
    id_index_t *ref_sorted = (id_index_t *)malloc(sizeof(id_index_t) * nReadA);
    if (!ref_sorted) {
        fprintf(stderr, "hermite3: malloc failed for ref_sorted\n");
        free(idA); FREE_ALL(); return -1;
    }
    for (long long i = 0; i < nReadA; i++) {
        ref_sorted[i].id  = idA[i];
        ref_sorted[i].idx = (uint32_t)i;
    }
    free(idA);

    if (ThisTask == 0) {
        fprintf(stdout, "Sorting %lld particles by ID...\n", nReadA);
        fflush(stdout);
    }
    if (radix_sort_id_index(ref_sorted, nReadA) != 0) {
        free(ref_sorted); FREE_ALL(); return -1;
    }

    /*
     * ref_sorted now maps sorted-ID position → original array index in A.
     * We need the inverse: for each position i in A's array, what is its
     * rank in the sorted order?
     *
     * Actually we don't need the inverse.  reorder_by_id() uses ref_sorted
     * to scatter other snapshots into A's original order (ref_sorted[i].idx
     * is the destination index in A's array space).  xA itself is already in
     * A's original order — we just need to verify it after the sort by
     * accessing xA[ref_sorted[i].idx] for the i-th sorted particle.
     * The Hermite loop runs over A's original index k, which is correct.
     */

    *NumPartRead_out = nReadA;

    /* ----------------------------------------------------------------
     * Step 2: Read P, match to A's order
     * ---------------------------------------------------------------- */
    if (ThisTask == 0) {
        fprintf(stdout, "Reading snapshot P (prev): %s\n", cfg->snap_prev);
        fflush(stdout);
    }

    struct sim_info hdrP;
    long long npartP = 0;
    read_snap_header_only(cfg->snap_prev, cfg->isHDF5, &hdrP, &npartP);

    long long nMatchP = 0;
    if (reorder_by_id(cfg->snap_prev, hdrP.NumFiles,
                      ref_sorted, nReadA,
                      xP, yP, zP, ptP, &nMatchP) != 0) {
        free(ref_sorted); FREE_ALL(); return -1;
    }
    sim_info_free(&hdrP);

    if (ThisTask == 0) {
        fprintf(stdout, "  matched %lld / %lld particles from P\n",
                nMatchP, nReadA);
        fflush(stdout);
    }

    /* ----------------------------------------------------------------
     * Step 3: Read B, match to A's order
     * ---------------------------------------------------------------- */
    if (ThisTask == 0) {
        fprintf(stdout, "Reading snapshot B: %s\n", cfg->snap_b);
        fflush(stdout);
    }

    struct sim_info hdrB;
    long long npartB = 0;
    read_snap_header_only(cfg->snap_b, cfg->isHDF5, &hdrB, &npartB);

    long long nMatchB = 0;
    if (reorder_by_id(cfg->snap_b, hdrB.NumFiles,
                      ref_sorted, nReadA,
                      xB, yB, zB, ptB, &nMatchB) != 0) {
        free(ref_sorted); FREE_ALL(); return -1;
    }
    sim_info_free(&hdrB);

    if (ThisTask == 0) {
        fprintf(stdout, "  matched %lld / %lld particles from B\n",
                nMatchB, nReadA);
        fflush(stdout);
    }

    /* ----------------------------------------------------------------
     * Step 4 (optional): Read N, match to A's order
     * ---------------------------------------------------------------- */
    if (has_next) {
        if (ThisTask == 0) {
            fprintf(stdout, "Reading snapshot N (next): %s\n", cfg->snap_next);
            fflush(stdout);
        }

        struct sim_info hdrN;
        long long npartN = 0;
        read_snap_header_only(cfg->snap_next, cfg->isHDF5, &hdrN, &npartN);

        long long nMatchN = 0;
        if (reorder_by_id(cfg->snap_next, hdrN.NumFiles,
                          ref_sorted, nReadA,
                          xN, yN, zN, ptN, &nMatchN) != 0) {
            free(ref_sorted); FREE_ALL(); return -1;
        }
        sim_info_free(&hdrN);

        if (ThisTask == 0) {
            fprintf(stdout, "  matched %lld / %lld particles from N\n",
                    nMatchN, nReadA);
            fflush(stdout);
        }
    }

    /* ref_sorted no longer needed — free before the tangent/interp loop */
    free(ref_sorted);

    /* ----------------------------------------------------------------
     * Step 5: Compute Catmull-Rom tangents
     *
     * We reuse the vxA/vyA/vzA and vxB/vyB/vzB names conceptually;
     * to avoid extra allocations we compute tangents into xP/yP/zP
     * (m0, tangent at A) and xN/yN/zN (m1, tangent at B) in-place,
     * overwriting the raw positions which are no longer needed after
     * this point.
     *
     * m0[k] = wrap(xA[k] - xP[k]) / 2
     * m1[k] = wrap(xN[k] - xA[k]) / 2   (4-snap)
     *       = wrap(xB[k] - xA[k]) / 2   (3-snap fallback)
     *
     * After this step:
     *   xP/yP/zP  hold  m0  (tangent at left  endpoint A)
     *   xN/yN/zN  hold  m1  (tangent at right endpoint B)
     *   xA/yA/zA  hold positions of A  (unchanged)
     *   xB/yB/zB  hold positions of B  (unchanged)
     * ---------------------------------------------------------------- */
    if (ThisTask == 0) {
        fprintf(stdout, "Computing Catmull-Rom tangents...\n");
        fflush(stdout);
    }

    long long nwx0=0, nwy0=0, nwz0=0;
    long long nwx1=0, nwy1=0, nwz1=0;

    for (long long k = 0; k < nReadA; k++) {
        /* Tangent at A (m0): centred difference over [P, A] */
        float dx0 = xA[k] - xP[k];
        float dy0 = yA[k] - yP[k];
        float dz0 = zA[k] - zP[k];
        if (dx0 >  half) { dx0 -= BoxSize; nwx0++; }
        else if (dx0 < -half) { dx0 += BoxSize; nwx0++; }
        if (dy0 >  half) { dy0 -= BoxSize; nwy0++; }
        else if (dy0 < -half) { dy0 += BoxSize; nwy0++; }
        if (dz0 >  half) { dz0 -= BoxSize; nwz0++; }
        else if (dz0 < -half) { dz0 += BoxSize; nwz0++; }
        /* Store m0 in xP (overwrite — raw P positions no longer needed) */
        xP[k] = dx0 * 0.5f;
        yP[k] = dy0 * 0.5f;
        zP[k] = dz0 * 0.5f;

        /* Tangent at B (m1) */
        float dx1, dy1, dz1;
        if (has_next) {
            /* 4-snap: centred difference over [A, N] */
            dx1 = xN[k] - xA[k];
            dy1 = yN[k] - yA[k];
            dz1 = zN[k] - zA[k];
        } else {
            /* 3-snap: one-sided backward difference over [A, B] */
            dx1 = xB[k] - xA[k];
            dy1 = yB[k] - yA[k];
            dz1 = zB[k] - zA[k];
        }
        if (dx1 >  half) { dx1 -= BoxSize; nwx1++; }
        else if (dx1 < -half) { dx1 += BoxSize; nwx1++; }
        if (dy1 >  half) { dy1 -= BoxSize; nwy1++; }
        else if (dy1 < -half) { dy1 += BoxSize; nwy1++; }
        if (dz1 >  half) { dz1 -= BoxSize; nwz1++; }
        else if (dz1 < -half) { dz1 += BoxSize; nwz1++; }

        if (has_next) {
            /* Store m1 in xN (overwrite — raw N positions no longer needed) */
            xN[k] = dx1 * 0.5f;
            yN[k] = dy1 * 0.5f;
            zN[k] = dz1 * 0.5f;
        } else {
            /*
             * 3-snap: allocate xN/yN/zN now if not already done,
             * or just reuse xB temporarily — but xB is still needed
             * in the Hermite loop.  Safest: store m1 in a dedicated
             * set.  We allocate it here on first use.
             */
            if (!xN) {
                xN = (float *)malloc(sizeof(float) * nReadA);
                yN = (float *)malloc(sizeof(float) * nReadA);
                zN = (float *)malloc(sizeof(float) * nReadA);
                if (!xN || !yN || !zN) {
                    fprintf(stderr, "hermite3: malloc failed for m1 arrays\n");
                    FREE_ALL(); return -1;
                }
            }
            xN[k] = dx1 * 0.5f;
            yN[k] = dy1 * 0.5f;
            zN[k] = dz1 * 0.5f;
        }
    }

    /* ptP/ptN no longer needed */
    free(ptP); ptP = NULL;
    if (ptN) { free(ptN); ptN = NULL; }

    if (ThisTask == 0) {
        fprintf(stdout,
                "Wraps m0: x=%lld y=%lld z=%lld\n"
                "Wraps m1: x=%lld y=%lld z=%lld\n",
                nwx0, nwy0, nwz0, nwx1, nwy1, nwz1);
        fflush(stdout);
    }

    /* ----------------------------------------------------------------
     * Step 6: Hermite interpolation
     *
     * Arrays at this point:
     *   xA/yA/zA  — positions of snapshot A
     *   xB/yB/zB  — positions of snapshot B
     *   xP/yP/zP  — tangent m0 at A  (= wrap(xA-xPrev)/2)
     *   xN/yN/zN  — tangent m1 at B  (= wrap(xNext-xA)/2 or wrap(xB-xA)/2)
     * ---------------------------------------------------------------- */
    float t  = cfg->interp_frac;
    float t2 = t * t, t3 = t2 * t;
    float h00 =  2.0f*t3 - 3.0f*t2 + 1.0f;
    float h10 =       t3 - 2.0f*t2 + t;
    float h01 = -2.0f*t3 + 3.0f*t2;
    float h11 =       t3 -       t2;

    if (ThisTask == 0) {
        fprintf(stdout,
                "Catmull-Rom Hermite: t=%.4f\n"
                "  h00=%.4f h10=%.4f h01=%.4f h11=%.4f\n",
                t, h00, h10, h01, h11);
        fflush(stdout);
    }

    for (long long k = 0; k < nReadA; k++) {
        float px = h00*xA[k] + h10*xP[k] + h01*xB[k] + h11*xN[k];
        float py = h00*yA[k] + h10*yP[k] + h01*yB[k] + h11*yN[k];
        float pz = h00*zA[k] + h10*zP[k] + h01*zB[k] + h11*zN[k];
        (*x)[k]     = wrap_pos(px, BoxSize);
        (*y)[k]     = wrap_pos(py, BoxSize);
        (*z)[k]     = wrap_pos(pz, BoxSize);
        (*ptype)[k] = ptA[k];
    }

    /* Free everything */
    free(xA); free(yA); free(zA); free(ptA);
    free(xB); free(yB); free(zB); free(ptB);
    free(xP); free(yP); free(zP);   /* m0 */
    free(xN); free(yN); free(zN);   /* m1 */

#undef FREE_ALL
    return 0;
}
