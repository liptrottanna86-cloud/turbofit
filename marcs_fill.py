"""
marcs_fill.py
=============

Fill holes in a MARCS atmosphere grid with APOGEE-style file names
(s5000_g+2.0_m1.0_t02_x3_z+0.00_a+0.00_c+0.00_..._s+0.00.mod), writing the
filled models to a separate directory so the original grid is untouched.

For each missing node:

1. **Interpolation** between the two closest existing models that bracket the
   node along one axis ([M/H], [a/M], [C/M], log g or Teff; at most two grid
   steps either side, same geometry, mass and microturbulence). Every
   number in the model is interpolated layer by layer on the Rosseland
   optical-depth scale: logarithmically for positive quantities spanning
   more than a factor 10 (pressures, opacities, densities, radius), linearly
   otherwise (T, log quantities, abundances). The two models must share
   the same depth scale; otherwise this step is skipped.
2. **Nearest neighbour** otherwise: a symlink to the closest existing model,
   preferring changes in abundance over log g, and log g over Teff.

Holes are filled only from original models, never from other fills.

Batch use (e.g. as a SLURM job):

    python marcs_fill.py --grid /mnt/.../marcs_apogee --out /mnt/.../marcs_apogee_filled \\
        --teff=2800:4000:100,4250:8000:250 --logg=-0.5:5:0.5 --mh=-2.5:1:0.25 \\
        --am=-1:1:0.25 --cm=-1.5:1:0.25 --nproc 16

On demand (what apogee_gc_test.py does):

    import marcs_fill
    filler = marcs_fill.Filler(griddir, filldir, teff_nodes)
    path, method, sources = filler.fill(teff, logg, mh, am, cm)

A report of every fill (node, method, source models, weight) is written to
<out>/fill_report.csv in batch mode.
"""

import os
import re
import sys
import argparse
import itertools
import numpy as np

import atmos

AXES = ('mh', 'am', 'cm', 'logg', 'teff')          # interpolation preference order
LINEAR = ('Depth', 'Vconv', 'Fconv/F', 'Mu', 'T')   # always interpolated linearly
NUM = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$')


def _axis(spec):
    """'a:b:step[,a:b:step...]' -> sorted array of node values."""
    vals = []
    for part in spec.split(','):
        a, b, d = (float(x) for x in part.split(':'))
        vals.extend(np.arange(a, b + 0.5 * d, d))
    return np.unique(np.round(vals, 4))


def model_name(teff, logg, mh, am, cm, vers='x3'):
    return os.path.basename(atmos.filename(teff, logg, mh, cm, am, model='marcs', vers=vers))


def _kind(name):
    """geometry / mass / microturbulence part of a file name, e.g. ('s', 'm1.0', 't02')."""
    p = name.split('_')
    return p[0][0], p[2], p[3]


# ---------------------------------------------------------------------
# reading / interpolating MARCS files
# ---------------------------------------------------------------------

def _tokens(line):
    return [(m.start(), m.end(), m.group()) for m in re.finditer(r'\S+', line)]


def _fmt_like(tok, val):
    """Format val like tok (same E/fixed style and number of decimals)."""
    if 'e' in tok.lower():
        mant = re.split('[eE]', tok)[0]
        d = len(mant.split('.')[1]) if '.' in mant else 0
        return f'{val:.{d}E}'
    d = len(tok.split('.')[1]) if '.' in tok else 0
    s = f'{val:.{d}f}'
    if tok.endswith('.') and d == 0:
        s += '.'
    if tok.startswith('+') and val >= 0:
        s = '+' + s
    return s


def _combine(a, b, f, logscale):
    """Interpolate two numbers, in log10 if logscale and both are positive."""
    if a == b:
        return a
    if logscale and a > 0 and b > 0:
        return 10 ** ((1 - f) * np.log10(a) + f * np.log10(b))
    return (1 - f) * a + f * b


def interpolate_models(file_a, file_b, f, name):
    """
    Model at weight f between file_a (f=0) and file_b (f=1), as text.
    Returns None if the two files do not have the same layout or depth scale.
    """
    A = open(file_a).read().splitlines()
    B = open(file_b).read().splitlines()
    if len(A) != len(B):
        return None
    out = [name]
    cols = None                                     # column names of the current table
    for la, lb in zip(A[1:], B[1:]):
        ta, tb = _tokens(la), _tokens(lb)
        if len(ta) != len(tb):
            return None
        if ta and ta[0][2] == 'k':
            cols = [t[2] for t in ta]
        line, pos = '', 0
        for ic, ((s0, e0, xa), (_, _, xb)) in enumerate(zip(ta, tb)):
            if NUM.match(xa) and NUM.match(xb):
                va, vb = float(xa), float(xb)
                is_int = re.match(r'^[+-]?\d+$', xa) is not None
                if is_int:
                    if va != vb:
                        return None                 # k index / layer count differ
                    new = xa
                else:
                    # E-format numbers (pressures, opacities, densities, radius,
                    # flux, gravity) in log; fixed-format (T, log tau, logs) linearly
                    lin = cols is not None and ic < len(cols) and cols[ic] in LINEAR
                    new = _fmt_like(xa, _combine(va, vb, f, 'e' in xa.lower() and not lin))
            else:
                new = xa                            # text (e.g. dates): keep A's
            pad = max(e0 - len(line) - len(new), 1 if line else 0)
            line += ' ' * pad + new
        out.append(line)
    # same Rosseland depth scale (first structure table, column lgTauR)
    try:
        i = next(k for k, l in enumerate(A) if l.strip() == 'Model structure')
        n = int(A[i - 1].split()[0])
        ra = np.array([float(A[i + 2 + k].split()[1]) for k in range(n)])
        rb = np.array([float(B[i + 2 + k].split()[1]) for k in range(n)])
        if np.max(np.abs(ra - rb)) > 0.01:
            return None
    except (StopIteration, ValueError, IndexError):
        return None
    return '\n'.join(out) + '\n'


# ---------------------------------------------------------------------
# hole filling
# ---------------------------------------------------------------------

class Filler:
    """
    griddir : original grid (flat directory of *.mod files).
    filldir : where filled models go (created if needed).
    teff_nodes : Teff node values (list, or (n, start, step)).
    steps : grid steps of the other axes.
    """

    def __init__(self, griddir, filldir, teff_nodes, steps=None, vers='x3', max_abund=2):
        self.griddir, self.filldir, self.vers = griddir, filldir, vers
        os.makedirs(filldir, exist_ok=True)
        if isinstance(teff_nodes, tuple) and len(teff_nodes) == 3:
            n, s, d = teff_nodes
            teff_nodes = s + d * np.arange(int(n))
        self.tn = np.sort(np.asarray(teff_nodes, float))
        self.steps = dict(mh=0.25, am=0.25, cm=0.25, logg=0.5)
        self.steps.update(steps or {})
        self.max_abund = max_abund
        self.have = set(f for f in os.listdir(griddir) if f.endswith('.mod'))

    def exists(self, t, g, m, a, c):
        return model_name(t, g, m, a, c, self.vers) in self.have

    def _shift(self, node, ax, n):
        """node moved n grid steps along ax (None if off the Teff list)."""
        v = dict(node)
        if ax == 'teff':
            k = int(np.argmin(np.abs(self.tn - v['teff']))) + n
            if not 0 <= k < len(self.tn):
                return None
            v['teff'] = float(self.tn[k])
        else:
            v[ax] = round(v[ax] + n * self.steps[ax], 3)
        return v

    def _tuple(self, v):
        return (v['teff'], v['logg'], v['mh'], v['am'], v['cm'])

    def _bracket(self, node, kind):
        """Closest bracketing pair along one axis: (span, axis, lo, hi, f)."""
        best = None
        for ia, ax in enumerate(AXES):
            for i, j in sorted(itertools.product((1, 2), (1, 2)), key=lambda x: x[0] + x[1]):
                lo, hi = self._shift(node, ax, -i), self._shift(node, ax, j)
                if lo is None or hi is None:
                    continue
                if not (self.exists(*self._tuple(lo)) and self.exists(*self._tuple(hi))):
                    continue
                nlo, nhi = model_name(*self._tuple(lo), self.vers), model_name(*self._tuple(hi), self.vers)
                if _kind(nlo) != kind or _kind(nhi) != kind:
                    continue
                x0, x1 = lo[ax], hi[ax]
                f = (node[ax] - x0) / (x1 - x0)
                cand = (i + j, ia, ax, lo, hi, f)
                if best is None or cand[:2] < best[:2]:
                    best = cand
                break                                  # smallest span on this axis found
        return best

    def _nearest(self, node):
        k = int(np.argmin(np.abs(self.tn - node['teff'])))
        r = range(-self.max_abund, self.max_abund + 1)
        cands = []
        for dm, da, dc, dg, dt in itertools.product(r, r, r, (0, -1, 1), (0, -1, 1)):
            if not 0 <= k + dt < len(self.tn):
                continue
            cost = abs(dm) + abs(da) + abs(dc) + 2.5 * abs(dg) + 5.0 * abs(dt)
            if cost > 0:
                cands.append((cost, dm, da, dc, dg, dt))
        cands.sort()
        for cost, dm, da, dc, dg, dt in cands:
            v = (float(self.tn[k + dt]), round(node['logg'] + dg * self.steps['logg'], 3),
                 round(node['mh'] + dm * self.steps['mh'], 3),
                 round(node['am'] + da * self.steps['am'], 3),
                 round(node['cm'] + dc * self.steps['cm'], 3))
            if self.exists(*v):
                return v
        return None

    def fill(self, teff, logg, mh, am=0.0, cm=0.0):
        """
        Path of the model for this node: the original if it exists, else a
        filled one in filldir (made now if needed). Returns (path, method,
        sources) with method 'original', 'interp', 'nearest' or 'none'.
        """
        node = dict(teff=float(teff), logg=float(logg), mh=float(mh), am=float(am), cm=float(cm))
        name = model_name(teff, logg, mh, am, cm, self.vers)
        if name in self.have:
            return os.path.join(self.griddir, name), 'original', ''
        out = os.path.join(self.filldir, name)
        note = out + '.how'
        if os.path.lexists(out) and os.path.exists(note):
            how = open(note).read().split('\t')
            return out, how[0], how[1].strip() if len(how) > 1 else ''

        b = self._bracket(node, _kind(name))
        if b is not None:
            _, _, ax, lo, hi, f = b
            fa = os.path.join(self.griddir, model_name(*self._tuple(lo), self.vers))
            fb = os.path.join(self.griddir, model_name(*self._tuple(hi), self.vers))
            txt = interpolate_models(fa, fb, f, name[:-4])
            if txt is not None:
                tmp = f'{out}.tmp{os.getpid()}'
                open(tmp, 'w').write(txt)
                os.replace(tmp, out)
                src = f'{os.path.basename(fa)} ({1 - f:.3f}) + {os.path.basename(fb)} ({f:.3f}) along {ax}'
                open(note, 'w').write(f'interp\t{src}\n')
                return out, 'interp', src
        v = self._nearest(node)
        if v is None:
            return None, 'none', ''
        src_name = model_name(*v, self.vers)
        tmp = f'{out}.tmp{os.getpid()}'
        if os.path.lexists(tmp):
            os.remove(tmp)
        os.symlink(os.path.abspath(os.path.join(self.griddir, src_name)), tmp)
        os.replace(tmp, out)
        open(note, 'w').write(f'nearest\t{src_name}\n')
        return out, 'nearest', src_name

    def lookup(self, base_orig, base_fill, verbose=True):
        """
        An atmod_lookup(teff, logg, mh, am, cm) for turbofit: originals through
        base_orig (turbogrid.make_atmod_lookup(griddir, flat=True)), fills
        through base_fill (make_atmod_lookup(filldir, flat=True)).
        """
        seen = {}

        def lookup(teff, logg, mh, am=0.0, cm=0.0):
            key = (round(teff, 1), round(logg, 2), round(mh, 2), round(am, 2), round(cm, 2))
            if key in seen:
                return seen[key]
            path, how, src = self.fill(teff, logg, mh, am, cm)
            if how == 'original':
                seen[key] = base_orig(teff, logg, mh, am, cm)
            elif how == 'none':
                if verbose:
                    print(f'hole ({teff:.0f}, {logg:+.1f}, {mh:+.2f}, {am:+.2f}, {cm:+.2f}): nothing nearby')
                seen[key] = None
            else:
                if verbose:
                    print(f'hole ({teff:.0f}, {logg:+.1f}, {mh:+.2f}, {am:+.2f}, {cm:+.2f}) {how}: {src}')
                seen[key] = base_fill(teff, logg, mh, am, cm)
            return seen[key]
        return lookup


# ---------------------------------------------------------------------
# batch
# ---------------------------------------------------------------------

_F = None


def _init(args):
    global _F
    _F = Filler(*args)


def _one(node):
    path, how, src = _F.fill(*node)
    return node, how, src


def main(argv=None):
    ap = argparse.ArgumentParser(description='Fill holes in a MARCS grid.')
    ap.add_argument('--grid', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--teff', required=True, help="e.g. --teff=2800:4000:100,4250:8000:250 (use = for negative starts)")
    ap.add_argument('--logg', required=True)
    ap.add_argument('--mh', required=True)
    ap.add_argument('--am', required=True)
    ap.add_argument('--cm', required=True)
    ap.add_argument('--vers', default='x3')
    ap.add_argument('--nproc', type=int, default=1)
    a = ap.parse_args(argv)

    ax = {k: _axis(getattr(a, k)) for k in ('teff', 'logg', 'mh', 'am', 'cm')}
    steps = {k: float(np.median(np.diff(ax[k]))) for k in ('logg', 'mh', 'am', 'cm') if len(ax[k]) > 1}
    F = Filler(a.grid, a.out, ax['teff'], steps=steps, vers=a.vers)
    nodes = [n for n in itertools.product(ax['teff'], ax['logg'], ax['mh'], ax['am'], ax['cm'])
             if not F.exists(*n)]
    total = np.prod([len(v) for v in ax.values()])
    print(f'{total} nodes, {total - len(nodes)} in the grid, {len(nodes)} to fill', flush=True)

    rows = []
    if a.nproc > 1:
        from multiprocessing import Pool
        with Pool(a.nproc, _init, ((a.grid, a.out, ax['teff'], steps, a.vers),)) as pool:
            for i, r in enumerate(pool.imap_unordered(_one, nodes, chunksize=64)):
                rows.append(r)
                if (i + 1) % 5000 == 0:
                    print(f'{i + 1}/{len(nodes)}', flush=True)
    else:
        _init((a.grid, a.out, ax['teff'], steps, a.vers))
        rows = [_one(n) for n in nodes]

    rows.sort()
    with open(os.path.join(a.out, 'fill_report.csv'), 'w') as f:
        f.write('teff,logg,mh,am,cm,method,source\n')
        for (t, g, m, al, c), how, src in rows:
            f.write(f'{t:.0f},{g:+.2f},{m:+.2f},{al:+.2f},{c:+.2f},{how},"{src}"\n')
    n = {k: sum(1 for r in rows if r[1] == k) for k in ('interp', 'nearest', 'none')}
    print(f"filled: {n['interp']} interpolated, {n['nearest']} nearest neighbour, "
          f"{n['none']} with nothing nearby; report in {a.out}/fill_report.csv")


if __name__ == '__main__':
    main()
