"""
fill_grid_holes.py
===================

Standalone post-processing step: fills missing nodes in an already-
COLLECTED spectra grid (the .npz written by turbogrid.collect_grid())
via Radial Basis Function (RBF) interpolation from neighboring spectra
in parameter space.

This is the turbogrid/Turbospectrum equivalent of Synple's
fill_synth()/rbf_get()/rbf_apply() (Synple User's Guide, Sect. 5), which
does the same thing for a Synspec-built FERRE grid -- and, like that
one, operates at the SPECTRA stage: it fills gaps in the already-
synthesized flux/cont grid, whatever caused them (a point that failed
synthesis, one still queued/running, or one dropped upstream for lacking
a matching atmosphere). It is NOT the same as atmos.fill_holes(), which
interpolates missing ATMOSPHERE models, upstream, before synthesis ever
runs.

Deliberately has NO dependency on turbogrid.py (and therefore none on
synthesis.py/atmos.py/utils.py either) -- only numpy and scipy -- so it
can run anywhere the collected .npz lives, on any node, without the
Turbospectrum synthesis toolchain being set up at all. It never
re-synthesizes anything; every filled node is an approximate,
interpolated stand-in, clearly flagged as such (see `interpolated` in
the output), not a real Turbospectrum spectrum.

FERRE itself cannot do this: the SYNTH grid format it reads requires
exactly N_P(1)*N_P(2)*...*N_P(ndim) lines with no missing-node concept,
and its own interpolation (INTER keyword) operates at query time between
EXISTING nodes for an observation being fit -- it does not repair gaps
in the grid's own structure. The grid must already be complete before
FERRE (or turbogrid's own write_ferre_grid()) will touch it.

Parallelism (single-job mode): the RBF fit/query is split across `nproc`
worker processes via a fork-based process pool. flux/cont are set as
module globals BEFORE the pool is created, so on Linux each forked
worker shares that memory copy-on-write instead of it being pickled and
duplicated per worker. HOWEVER: scipy's RBFInterpolator unconditionally
upcasts its `d` (data value) argument to float64 on construction
(`d = np.asarray(d, dtype=float, ...)`), regardless of the `neighbors`
argument or the input dtype -- so each worker independently re-inflates
its own full float64 copy of flux and of cont the moment it builds its
interpolator, and that copy is NOT shared across workers (it's a fresh
allocation, not the forked COW buffer). Passing --dtype float32 still
halves the ONE shared parent-process copy and the final written grid,
but does nothing for this per-worker cost. Rough per-worker memory:

    2 * n_real_points * npix * 8 bytes    (flux + cont, upcast to float64)

e.g. ~857k real points x 2401 pixels is ~30.7 GB PER WORKER, on top of
the shared parent-process copy. Size nproc accordingly -- this is the
dominant memory cost, not thread/process overhead.

Multi-node: ProcessPoolExecutor's fork pool only shares memory within a
single node, so it cannot itself span nodes. To spread the work across
several nodes, use fill_grid_holes_chunk() to process a slice of the
missing nodes per SLURM array task/node (each still loads the full real
-point set, since RBF neighbor search needs it, but only needs a small
nproc and modest --mem), then merge_filled_chunks() to stitch the
per-chunk outputs back into one completed grid.

Run this as ITS OWN job, separate from the synthesis array -- it wants
many cores on one node (not many nodes, per job), and needs the entire
collected flux/cont array resident in memory at once.

Usage
-----
As a library:

    from fill_grid_holes import fill_grid_holes
    fill_grid_holes('grid1/grid.npz',
                     labels=('teff', 'logg', 'mh', 'cm', 'nm', 'am'),
                     nproc=16)

From the command line, single job:

    python fill_grid_holes.py grid1/grid.npz \\
        --labels cm nm am mh logg teff \\
        --kernel thin_plate_spline --neighbors 100 --degree 0 --nproc 16 \\
        --outfile grid1/grid_filled.npz

Split across N nodes (e.g. a SLURM array of size N), then merge:

    # one of these per array task, $SLURM_ARRAY_TASK_ID = 0..N-1
    python fill_grid_holes.py grid1/grid.npz --labels cm nm am mh logg teff \\
        --dtype float32 --nproc 4 \\
        --n-chunks 8 --chunk-index $SLURM_ARRAY_TASK_ID \\
        --outfile grid1/chunk${SLURM_ARRAY_TASK_ID}.npz

    # once all chunks are done
    python fill_grid_holes.py grid1/grid.npz --labels cm nm am mh logg teff \\
        --merge grid1/chunk0.npz grid1/chunk1.npz ... grid1/chunk7.npz \\
        --outfile grid1/grid_filled.npz
"""

import os
import sys
import argparse
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import numpy as np


def _default_nproc():
    """Cores actually available to this job (respects Slurm/cgroup cpu
    binding), falling back to os.cpu_count() where that's not available
    (e.g. non-Linux). NOTE: on some Slurm configurations this does NOT
    respect --cpus-per-task -- if you see nproc reported higher than
    what you requested, always pass --nproc explicitly rather than
    relying on this default."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


# Set (as plain module globals) right before the worker pool is created,
# so forked workers on Linux inherit them via copy-on-write instead of
# having them pickled and duplicated per worker -- see module docstring.
# NOT meaningful outside of a fill_grid_holes()/_run_rbf() call in
# progress.
_have_coords = None
_flux = None
_cont = None
_kernel = None
_n_neighbors = None
_degree = None
# built lazily, once, inside each worker process on its first chunk
_worker_rbf_flux = None
_worker_rbf_cont = None


def _fill_chunk(missing_coords_chunk):
    """Runs inside a worker process. Reads the module globals set by the
    parent just before the pool was created (inherited via fork, not
    pickled) -- builds this worker's own RBFInterpolator once, then
    evaluates it for this chunk of missing coordinates.

    NOTE: RBFInterpolator upcasts flux/cont to float64 internally on
    construction (see module docstring) -- this is where that per-worker
    memory cost is actually paid, not here in the query."""
    global _worker_rbf_flux, _worker_rbf_cont
    from scipy.interpolate import RBFInterpolator
    if _worker_rbf_flux is None:
        _worker_rbf_flux = RBFInterpolator(_have_coords, _flux, kernel=_kernel,
                                            neighbors=_n_neighbors, degree=_degree)
        _worker_rbf_cont = RBFInterpolator(_have_coords, _cont, kernel=_kernel,
                                            neighbors=_n_neighbors, degree=_degree)
    return _worker_rbf_flux(missing_coords_chunk), _worker_rbf_cont(missing_coords_chunk)


def read_grid(gridfile):
    """Load a grid file written by turbogrid.collect_grid()."""
    d = np.load(gridfile, allow_pickle=True)
    return d['params'], d['wave'], d['flux'], d['cont']


def _locate_holes(params, labels, rtol=1e-4):
    """Shared setup: figure out which regular-grid nodes are missing.

    Returns a dict with everything fill_grid_holes(), fill_grid_holes_chunk()
    and merge_filled_chunks() all need, so the "what's missing and where"
    logic exists in exactly one place.

    Coordinates are rescaled per axis to (value - llimits)/steps -- i.e.
    integer grid-index units -- before being handed to the RBF
    interpolator. Without this, an axis like teff (spacing of hundreds of
    K) would completely dominate the distance metric over an axis like
    logg or [M/H] (spacing of a few tenths).

    Each missing node is classified `edge` (touches the min/max index of
    at least one varying axis) or `interior` -- edge fills are
    extrapolations (real data on only one side along that axis) and are
    also the points most likely to hit a degenerate/singular local fit
    (see the `degree` parameter on fill_grid_holes()).
    """
    axes = []
    for lab in labels:
        vals = np.unique(np.round(params[lab], 6))
        if len(vals) > 1:
            axes.append((lab, vals))
    if not axes:
        raise ValueError('none of the requested labels vary in this grid')

    n_p, steps, llimits = [], [], []
    for lab, vals in axes:
        d_v = np.diff(vals)
        if not np.allclose(d_v, d_v[0], rtol=rtol):
            raise ValueError(f"axis '{lab}' is not evenly spaced ({vals}); "
                              'a regular grid is required to locate holes '
                              '(this must match the completeness check '
                              'write_ferre_grid() will later apply)')
        n_p.append(len(vals))
        steps.append(float(d_v[0]))
        llimits.append(float(vals[0]))

    strides = [1] * len(axes)
    for i in range(1, len(axes)):
        strides[i] = strides[i - 1] * n_p[i - 1]
    ntotal = int(np.prod(n_p))

    def axis_index(lab_vals, value):
        return int(np.argmin(np.abs(lab_vals - value)))

    have_flat = np.zeros(ntotal, dtype=bool)
    have_coords = np.empty((len(params), len(axes)))
    for row in range(len(params)):
        flat = 0
        for j, (lab, vals) in enumerate(axes):
            aidx = axis_index(vals, params[lab][row])
            have_coords[row, j] = (params[lab][row] - llimits[j]) / steps[j]
            flat += aidx * strides[j]
        have_flat[flat] = True

    missing_flat = np.where(~have_flat)[0]
    missing_coords_scaled = np.empty((len(missing_flat), len(axes)))
    missing_coords_real = np.empty((len(missing_flat), len(axes)))
    is_edge = np.zeros(len(missing_flat), dtype=bool)
    for k, flat in enumerate(missing_flat):
        rem = flat
        for j in range(len(axes) - 1, -1, -1):
            aidx = rem // strides[j]
            rem = rem % strides[j]
            missing_coords_scaled[k, j] = aidx
            missing_coords_real[k, j] = axes[j][1][aidx]
            if aidx == 0 or aidx == n_p[j] - 1:
                is_edge[k] = True

    return dict(axes=axes, n_p=n_p, steps=steps, llimits=llimits,
                strides=strides, ntotal=ntotal, axis_index=axis_index,
                have_flat=have_flat, have_coords=have_coords,
                missing_flat=missing_flat,
                missing_coords_scaled=missing_coords_scaled,
                missing_coords_real=missing_coords_real, is_edge=is_edge)


def _run_rbf(have_coords, flux, cont, coords_scaled, kernel, neighbors, nproc, degree=0):
    """Runs the (possibly pooled) RBF fit/query for a set of query coords.

    `degree` is the RBF's added-polynomial degree (passed straight to
    RBFInterpolator). Default 0 (constant term only) rather than scipy's
    own default (1, for thin_plate_spline) -- degree 1 needs the local
    neighborhood to vary independently along every axis, which can fail
    with a Singular matrix / rank-deficient error for edge/corner missing
    points where the 100 nearest real neighbors don't span all axes
    independently. Degree 0 is always solvable given >=1 neighbor, at
    the cost of the local linear trend term -- a reasonable trade for
    exactly the points that were already the least-trustworthy fills.
    """
    n_neighbors = min(neighbors, len(have_coords))
    global _have_coords, _flux, _cont, _kernel, _n_neighbors, _degree
    _have_coords, _flux, _cont, _kernel, _n_neighbors, _degree = (
        have_coords, flux, cont, kernel, n_neighbors, degree)

    if nproc <= 1 or len(coords_scaled) < nproc:
        return _fill_chunk(coords_scaled)

    chunks = np.array_split(coords_scaled, nproc)
    ctx = mp.get_context('fork')
    with ProcessPoolExecutor(max_workers=nproc, mp_context=ctx) as ex:
        results = list(ex.map(_fill_chunk, chunks))
    filled_flux = np.vstack([r[0] for r in results])
    filled_cont = np.vstack([r[1] for r in results])
    return filled_flux, filled_cont


def _assemble_real_points(params, flux, cont, wave, loc, npix):
    """Builds all_params/all_flux/all_cont/interpolated/edge_fill arrays,
    populated with the REAL points only (flat-index order). Shared by
    fill_grid_holes() and merge_filled_chunks()."""
    dtype = params.dtype
    ntotal = loc['ntotal']
    all_params = np.zeros(ntotal, dtype=dtype)
    all_flux = np.empty((ntotal, npix), dtype=flux.dtype)
    all_cont = np.empty((ntotal, npix), dtype=cont.dtype)
    interpolated = np.zeros(ntotal, dtype=bool)
    edge_fill = np.zeros(ntotal, dtype=bool)

    axis_labels = [lab for lab, _ in loc['axes']]
    axis_index = loc['axis_index']
    strides = loc['strides']

    for row in range(len(params)):
        flat = 0
        for j, (lab, vals) in enumerate(loc['axes']):
            flat += axis_index(vals, params[lab][row]) * strides[j]
        all_params[flat] = params[row]
        all_flux[flat] = flux[row]
        all_cont[flat] = cont[row]

    have_rows_coords_real = np.array(
        [[params[lab][row] for lab in axis_labels] for row in range(len(params))])

    return dict(all_params=all_params, all_flux=all_flux, all_cont=all_cont,
                interpolated=interpolated, edge_fill=edge_fill,
                axis_labels=axis_labels, have_rows_coords_real=have_rows_coords_real,
                dtype=dtype)


def _place_filled(assembled, flat_idx, coords_real, filled_flux, filled_cont,
                   edge_flags, params, seen=None):
    """Writes one batch of interpolated points into the assembled arrays
    (in place). Raises if a flat index was already placed (only relevant
    when called repeatedly across chunk files, see merge_filled_chunks).
    Shared by fill_grid_holes() and merge_filled_chunks()."""
    all_params = assembled['all_params']
    all_flux = assembled['all_flux']
    all_cont = assembled['all_cont']
    interpolated = assembled['interpolated']
    edge_fill = assembled['edge_fill']
    dtype = assembled['dtype']
    axis_labels = assembled['axis_labels']
    have_rows_coords_real = assembled['have_rows_coords_real']

    for k, flat in enumerate(flat_idx):
        if seen is not None:
            if seen[flat]:
                raise ValueError(f'flat index {flat} appears in more than one chunk '
                                  f'file -- chunks must be non-overlapping (from the '
                                  f'same --n-chunks run)')
            seen[flat] = True
        nearest_row = int(np.argmin(
            ((have_rows_coords_real - coords_real[k]) ** 2).sum(axis=1)))
        rec = list(params[nearest_row])
        for j, lab in enumerate(axis_labels):
            rec[dtype.names.index(lab)] = coords_real[k, j]
        all_params[flat] = tuple(rec)
        all_flux[flat] = filled_flux[k]
        all_cont[flat] = filled_cont[k]
        interpolated[flat] = True
        edge_fill[flat] = edge_flags[k]


def fill_grid_holes(gridfile, outfile=None,
                     labels=('teff', 'logg', 'mh', 'cm', 'am', 'nm'),
                     kernel='thin_plate_spline', neighbors=100, degree=0, rtol=1e-4,
                     nproc=None, dtype=None):
    """Single-node, all-in-one fill. See module docstring for the memory
    caveats (RBFInterpolator's float64 upcast) and the `degree` default.

    Parameters mostly as before; see module docstring's Usage section.
    """
    d = np.load(gridfile, allow_pickle=True)
    params, wave, flux, cont = d['params'], d['wave'], d['flux'], d['cont']
    if dtype is not None:
        flux = flux.astype(dtype)
        cont = cont.astype(dtype)
    npix = len(wave)

    loc = _locate_holes(params, labels, rtol=rtol)
    missing_flat = loc['missing_flat']
    if len(missing_flat) == 0:
        print('no holes to fill -- grid is already complete')
        dest = outfile or gridfile
        if dest != gridfile:
            np.savez(dest, params=params, wave=wave, flux=flux, cont=cont,
                     interpolated=np.zeros(len(params), dtype=bool),
                     edge_fill=np.zeros(len(params), dtype=bool))
        return dest

    n_edge = int(loc['is_edge'].sum())
    n_interior = len(missing_flat) - n_edge
    nproc = nproc or _default_nproc()
    print(f'{len(missing_flat)}/{loc["ntotal"]} grid nodes are missing '
          f'({n_interior} interior, {n_edge} on a grid edge) -- '
          f'filling via RBF ({kernel}, neighbors={neighbors}, degree={degree}, '
          f'nproc={nproc})')
    if n_edge:
        print(f'  note: edge fills are extrapolated -- treat those {n_edge} '
              f'points with extra caution; see the edge_fill flag in the output')

    filled_flux, filled_cont = _run_rbf(
        loc['have_coords'], flux, cont, loc['missing_coords_scaled'],
        kernel, neighbors, nproc, degree=degree)

    assembled = _assemble_real_points(params, flux, cont, wave, loc, npix)
    _place_filled(assembled, missing_flat, loc['missing_coords_real'],
                  filled_flux, filled_cont, loc['is_edge'], params)

    outfile = outfile or (os.path.splitext(gridfile)[0] + '_filled.npz')
    np.savez(outfile, params=assembled['all_params'], wave=wave,
              flux=assembled['all_flux'], cont=assembled['all_cont'],
              interpolated=assembled['interpolated'], edge_fill=assembled['edge_fill'])
    n_interp = int(assembled['interpolated'].sum())
    print(f'wrote completed grid ({loc["ntotal"]} points, {n_interp} interpolated '
          f'[{n_edge} on an edge] / {loc["ntotal"] - n_interp} real) to {outfile}')
    return outfile


def fill_grid_holes_chunk(gridfile, chunk_index, n_chunks, outfile=None,
                           labels=('teff', 'logg', 'mh', 'cm', 'am', 'nm'),
                           kernel='thin_plate_spline', neighbors=100, degree=0,
                           rtol=1e-4, nproc=None, dtype=None):
    """Multi-node step 1: process only this job's slice of the missing
    points. Run one of these per SLURM array task / per node, then feed
    all the outputs to merge_filled_chunks(). Each task still loads the
    full real-point grid (needed for the neighbor search), but only needs
    `nproc` workers and modest --mem for its own slice -- see the module
    docstring for how to size nproc against --mem (RBFInterpolator's
    float64 upcast, ~30GB/worker at 857k real points x 2401 pixels)."""
    if not (0 <= chunk_index < n_chunks):
        raise ValueError('chunk_index must be in [0, n_chunks)')

    d = np.load(gridfile, allow_pickle=True)
    params, wave, flux, cont = d['params'], d['wave'], d['flux'], d['cont']
    if dtype is not None:
        flux = flux.astype(dtype)
        cont = cont.astype(dtype)

    loc = _locate_holes(params, labels, rtol=rtol)
    missing_flat = loc['missing_flat']
    if len(missing_flat) == 0:
        print('no holes to fill -- grid is already complete; nothing for this chunk to do')
        outfile = outfile or (os.path.splitext(gridfile)[0] + f'_chunk{chunk_index}.npz')
        np.savez(outfile, missing_flat=np.array([], dtype=int),
                 filled_flux=np.empty((0, len(wave))), filled_cont=np.empty((0, len(wave))),
                 is_edge=np.array([], dtype=bool),
                 missing_coords_real=np.empty((0, 0)), labels=np.array(labels))
        return outfile

    idx_slices = np.array_split(np.arange(len(missing_flat)), n_chunks)
    sl = idx_slices[chunk_index]
    nproc = nproc or _default_nproc()
    print(f'chunk {chunk_index}/{n_chunks}: filling {len(sl)}/{len(missing_flat)} '
          f'missing nodes, nproc={nproc}')

    filled_flux, filled_cont = _run_rbf(
        loc['have_coords'], flux, cont, loc['missing_coords_scaled'][sl],
        kernel, neighbors, nproc, degree=degree)

    outfile = outfile or (os.path.splitext(gridfile)[0] + f'_chunk{chunk_index}.npz')
    np.savez(outfile,
              missing_flat=missing_flat[sl],
              filled_flux=filled_flux,
              filled_cont=filled_cont,
              is_edge=loc['is_edge'][sl],
              missing_coords_real=loc['missing_coords_real'][sl],
              labels=np.array(labels))
    print(f'wrote chunk {chunk_index} ({len(sl)} points) to {outfile}')
    return outfile


def merge_filled_chunks(gridfile, chunk_files, outfile=None,
                         labels=('teff', 'logg', 'mh', 'cm', 'am', 'nm'),
                         rtol=1e-4):
    """Multi-node step 2: stitch chunk outputs from fill_grid_holes_chunk()
    back into one completed grid, identical in format/content to what
    fill_grid_holes() would have produced in one shot."""
    d = np.load(gridfile, allow_pickle=True)
    params, wave, flux, cont = d['params'], d['wave'], d['flux'], d['cont']
    npix = len(wave)

    loc = _locate_holes(params, labels, rtol=rtol)
    assembled = _assemble_real_points(params, flux, cont, wave, loc, npix)

    seen = np.zeros(loc['ntotal'], dtype=bool)
    total_placed = 0
    for cf in chunk_files:
        c = np.load(cf, allow_pickle=True)
        flat_idx = c['missing_flat']
        if len(flat_idx) == 0:
            continue
        _place_filled(assembled, flat_idx, c['missing_coords_real'],
                      c['filled_flux'], c['filled_cont'], c['is_edge'],
                      params, seen=seen)
        total_placed += len(flat_idx)

    n_missing = len(loc['missing_flat'])
    if total_placed != n_missing:
        raise ValueError(f'chunk files covered {total_placed} missing nodes, '
                          f'expected {n_missing} -- did every chunk run, from the '
                          f'same --n-chunks value, get included?')

    outfile = outfile or (os.path.splitext(gridfile)[0] + '_filled.npz')
    np.savez(outfile, params=assembled['all_params'], wave=wave,
              flux=assembled['all_flux'], cont=assembled['all_cont'],
              interpolated=assembled['interpolated'], edge_fill=assembled['edge_fill'])
    n_edge = int(assembled['edge_fill'].sum())
    print(f'merged {len(chunk_files)} chunks -> {loc["ntotal"]} points, '
          f'{total_placed} interpolated [{n_edge} on an edge] / '
          f'{loc["ntotal"] - total_placed} real, wrote {outfile}')
    return outfile


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('gridfile', help='.npz from turbogrid.collect_grid()')
    p.add_argument('--outfile', default=None,
                    help="default: '<gridfile stem>_filled.npz' (or "
                         "'_chunkN.npz' in chunk mode)")
    p.add_argument('--kernel', default='thin_plate_spline',
                    help='RBF kernel (default: thin_plate_spline)')
    p.add_argument('--neighbors', type=int, default=100,
                    help='local neighborhood size per RBF fit (default: 100)')
    p.add_argument('--degree', type=int, default=0,
                    help='RBF added-polynomial degree (default: 0, i.e. '
                         'constant term only -- always solvable; scipy\'s own '
                         'default of 1 for thin_plate_spline can raise a '
                         'Singular matrix error for edge/corner missing '
                         'points whose local neighborhood does not vary '
                         'independently along every axis)')
    p.add_argument('--nproc', type=int, default=None,
                    help='worker processes for the RBF fit/query. IMPORTANT: '
                         'each worker independently upcasts flux+cont to '
                         'float64 on construction (~30GB/worker at 857k real '
                         'points x 2401 pixels) -- size this against --mem, '
                         'not against core count. Always pass explicitly on '
                         'clusters where cgroup CPU affinity may not reflect '
                         '--cpus-per-task.')
    p.add_argument('--dtype', default=None, choices=['float32', 'float64'],
                    help='downcast flux/cont before fitting -- halves the '
                         'shared parent-process copy and the final written '
                         'grid, but NOT the per-worker RBFInterpolator cost '
                         '(scipy upcasts to float64 internally regardless)')
    p.add_argument('--n-chunks', type=int, default=1,
                    help='split the missing nodes into this many independent '
                         'jobs (e.g. one per SLURM array task / node)')
    p.add_argument('--chunk-index', type=int, default=0,
                    help='which chunk this job processes, 0-indexed '
                         '(e.g. $SLURM_ARRAY_TASK_ID)')
    p.add_argument('--merge', nargs='+', default=None, metavar='CHUNK_NPZ',
                    help='merge these fill_grid_holes_chunk() outputs into '
                         'one completed grid, instead of filling anything')
    p.add_argument('--labels', nargs='+',
                default=['teff', 'logg', 'mh', 'cm', 'am', 'nm'],
                help='grid axes to treat as the regular grid (default: '
                     'teff logg mh cm am nm)')
    args = p.parse_args()
    dtype = {'float32': np.float32, 'float64': np.float64}.get(args.dtype)

    if args.merge:
        merge_filled_chunks(args.gridfile, args.merge, outfile=args.outfile,
                             labels=tuple(args.labels))
    elif args.n_chunks > 1:
        fill_grid_holes_chunk(args.gridfile, args.chunk_index, args.n_chunks,
                               outfile=args.outfile, labels=tuple(args.labels),
                               kernel=args.kernel, neighbors=args.neighbors,
                               degree=args.degree, nproc=args.nproc, dtype=dtype)
    else:
        fill_grid_holes(args.gridfile, outfile=args.outfile, labels=tuple(args.labels),
                         kernel=args.kernel, neighbors=args.neighbors,
                         degree=args.degree, nproc=args.nproc, dtype=dtype)