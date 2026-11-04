"""
nlte.py
=======

NLTE support for synthesis.synthesize(..., nlte=cfg, nlte_elems=[...]).

In NLTE mode the model atmosphere is NOT taken from the APOGEE C/alpha/M
grid: Turbospectrum's own Fortran interpolator (interpol_modeles_nlte,
T. Masseron's MARCS interpolator extended by Gerber, Magg & Bergemann)
interpolates the 8 bracketing STANDARD MARCS models and, with the same
weights, the grid of NLTE departure coefficients for each requested element.
Those departure-coefficient grids were computed on standard MARCS models, so
this is the only self-consistent combination. Elements with no departure
grid fall back to LTE with a warning.

The interpolator input, the species-file format and the NaN clean-up of the
departure files follow TSFitPy (Gerber et al. 2023; Storm & Bergemann 2023,
https://github.com/TSFitPy-developers/TSFitPy, GPL-3).

Setup (once)
------------
1. Compile the NLTE interpolator that ships with Turbospectrum_NLTE:
       cd $TURBOSPEC_DIR/interpolator
       gfortran -o interpol_modeles_nlte interpol_modeles_nlte_gfort.f
   (use interpol_modeles_nlte.f with ifort/ifx).
2. Download from https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/ the
   MARCS-type departure grids (NLTEgrid*MARCS*.bin + auxData*MARCS*.dat) and
   model atoms (atom.*) for the elements you want, laid out as
       nlte_dir/<El>/<bin and aux files>      (paths as in nlte_filenames.cfg)
       nlte_dir/model_atoms/atom.*
3. Have the standard MARCS grid (ASCII web format, *.mod) in one directory.

    import nlte
    cfg = nlte.NLTEConfig(marcs_dir='/data/marcs_standard',
                          nlte_dir='/data/nlte_data')
    print(cfg.available_elements())
"""

import os
import re
import glob
import subprocess
import warnings
import configparser
import numpy as np

import atomic

_HERE = os.path.dirname(os.path.abspath(__file__))

# MARCS file name, e.g. p5000_g+4.0_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod
_MARCS_RE = re.compile(r'^([ps])(\d{4,5})_g([+-]\d+\.\d+)_m(\d+\.\d+)_t(\d+)_(\w\w)_'
                       r'z([+-]\d+\.\d+)_a([+-]\d+\.\d+)_c([+-]\d+\.\d+)_n([+-]\d+\.\d+)_'
                       r'o([+-]\d+\.\d+)_r([+-]\d+\.\d+)_s([+-]\d+\.\d+)\.mod$')


class NLTEConfig:
    """
    Everything NLTE mode needs to know, set up once and passed to
    synthesis.synthesize(..., nlte=cfg, nlte_elems=[...]).

    Parameters
    ----------
    marcs_dir : str
        Directory holding the standard MARCS grid (ASCII *.mod files, flat).
    nlte_dir : str
        Root of the departure-coefficient data; bin/aux paths in the config
        file are relative to this.
    config_file : str, optional
        TSFitPy-style nlte_filenames.cfg. Default: nlte_dir/nlte_filenames.cfg
        if present, else the copy shipped with turbofit.
    model_atom_dir : str, optional
        Where the atom.* files live. Default nlte_dir/model_atoms.
    interpol_dir : str, optional
        Where the compiled interpol_modeles_nlte lives. Default
        $TURBOSPEC_DIR/interpolator.
    comp : str
        MARCS chemical-composition class to use ('st' = standard). Default 'st'.
    """

    def __init__(self, marcs_dir, nlte_dir, config_file=None, model_atom_dir=None,
                 interpol_dir=None, comp='st'):
        self.marcs_dir = os.path.abspath(marcs_dir)
        self.nlte_dir = os.path.abspath(nlte_dir)
        if config_file is None:
            config_file = os.path.join(self.nlte_dir, 'nlte_filenames.cfg')
            if not os.path.exists(config_file):
                config_file = os.path.join(_HERE, 'nlte_filenames.cfg')
        self.config_file = config_file
        self.model_atom_dir = os.path.abspath(model_atom_dir or
                                              os.path.join(self.nlte_dir, 'model_atoms'))
        if interpol_dir is None:
            interpol_dir = os.path.join(os.environ.get('TURBOSPEC_DIR', ''), 'interpolator')
        self.interpol_exe = os.path.join(os.path.abspath(interpol_dir), 'interpol_modeles_nlte')
        self.comp = comp

        cp = configparser.ConfigParser()
        cp.read(config_file)
        self.files = {}
        for el in cp.sections():
            if '1d_bin' in cp[el] and '1d_aux' in cp[el] and 'atom_file' in cp[el]:
                self.files[el] = dict(bin=os.path.join(self.nlte_dir, cp[el]['1d_bin']),
                                      aux=os.path.join(self.nlte_dir, cp[el]['1d_aux']),
                                      atom=cp[el]['atom_file'])
        self._aux_len = {}
        self._grid = None

    # -- data checks ----------------------------------------------------
    def available_elements(self):
        """Elements whose bin, aux and model-atom files are all on disk."""
        ok = []
        for el, f in self.files.items():
            if (os.path.exists(f['bin']) and os.path.exists(f['aux']) and
                    os.path.exists(os.path.join(self.model_atom_dir, f['atom']))):
                ok.append(el)
        return ok

    def aux_length(self, el):
        if el not in self._aux_len:
            self._aux_len[el] = len(np.loadtxt(self.files[el]['aux'], dtype=str, ndmin=2))
        return self._aux_len[el]

    # -- MARCS grid ------------------------------------------------------
    @property
    def grid(self):
        """{(geom, mass, vt): {(teff, logg, z): filename}} for comp==self.comp."""
        if self._grid is None:
            grid = {}
            for path in glob.glob(os.path.join(self.marcs_dir, '*.mod')):
                m = _MARCS_RE.match(os.path.basename(path))
                if not m or m.group(6) != self.comp:
                    continue
                geom, teff, logg, mass, vt = (m.group(1), int(m.group(2)), float(m.group(3)),
                                              float(m.group(4)), int(m.group(5)))
                z = float(m.group(7))
                grid.setdefault((geom, mass, vt), {})[(teff, logg, z)] = path
            if not grid:
                raise FileNotFoundError(f'no standard-composition ({self.comp}) MARCS *.mod '
                                        f'files found in {self.marcs_dir}')
            self._grid = grid
        return self._grid


def _family(cfg, logg):
    """Spherical 1 Msun vt=2 models for logg<3, plane-parallel vt=1 (else 2) above,
    following the MARCS grid conventions (and TSFitPy)."""
    if logg < 3.0:
        prefs = [('s', 1.0, 2)]
    else:
        prefs = [('p', 0.0, 1), ('p', 0.0, 2)]
    for key in prefs:
        if key in cfg.grid:
            return key
    raise FileNotFoundError(f'no MARCS models of type {prefs} in {cfg.marcs_dir}')


def bracket_models(cfg, teff, logg, mh):
    """
    The 8 MARCS models bracketing (teff, logg, [M/H]), in the order the
    interpolator wants (bit0=[M/H], bit1=logg, bit2=Teff; low before up).
    If a corner model is missing, the face with most missing corners is moved
    out by one grid step and the search repeated (as in TSFitPy).
    """
    fam = _family(cfg, logg)
    models = cfg.grid[fam]
    axes = [sorted({k[i] for k in models}) for i in range(3)]   # teff, logg, z
    target = [teff, logg, mh]
    lohi = []
    for vals, x, name in zip(axes, target, ('Teff', 'logg', '[M/H]')):
        if x < vals[0] or x > vals[-1]:
            raise ValueError(f'{name}={x} outside the MARCS grid range {vals[0]}..{vals[-1]}')
        i = int(np.searchsorted(vals, x, side='right')) - 1
        i = min(i, len(vals) - 1)
        if np.isclose(vals[i], x):
            lohi.append([i, i])
        else:
            lohi.append([i, i + 1])

    while True:
        corners, missing = [], []
        for v in range(8):
            it = lohi[0][(v >> 2) & 1]
            ig = lohi[1][(v >> 1) & 1]
            iz = lohi[2][v & 1]
            key = (axes[0][it], axes[1][ig], axes[2][iz])
            corners.append(models.get(key))
            if key not in models:
                missing.append(v)
        if not missing:
            return corners, fam
        # count missing corners per face: (axis, side)
        bits = {0: 2, 1: 1, 2: 0}     # axis -> bit position
        counts = []
        for ax in range(3):
            for side in (0, 1):
                n = sum(1 for v in missing if ((v >> bits[ax]) & 1) == side)
                counts.append((n, ax, side))
        n, ax, side = max(counts)
        if side == 0:
            lohi[ax][0] -= 1
            if lohi[ax][0] < 0:
                raise ValueError(f'cannot bracket ({teff},{logg},{mh}): MARCS grid holes')
        else:
            lohi[ax][1] += 1
            if lohi[ax][1] >= len(axes[ax]):
                raise ValueError(f'cannot bracket ({teff},{logg},{mh}): MARCS grid holes')


def _fix_nan_departures(path):
    """Replace NaN departure coefficients (port of TSFitPy's
    check_nan_in_coefficients): NaNs in the deepest layer become 1 (LTE),
    others take the value of the layer below."""
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while lines[i].startswith('#'):
        i += 1
    i += 1
    ndep = int(lines[i])
    i += 2 + ndep
    block = lines[i:i + ndep][::-1]
    changed = False
    for k, line in enumerate(block):
        toks = line.split()
        if not any(t.lower() == 'nan' for t in toks):
            continue
        changed = True
        below = block[k - 1].split() if k > 0 else None
        toks = [('1.00000E+00' if below is None else below[j]) if t.lower() == 'nan' else t
                for j, t in enumerate(toks)]
        block[k] = ' '.join(f'{t:>12s}' for t in toks)
    if changed:
        lines = lines[:i] + block[::-1] + lines[i + ndep:]
        with open(path, 'w') as f:
            f.write('\n'.join(lines) + '\n')


def prepare(cfg, teff, logg, mh, abundances, nlte_elems, workdir='.', verbose=False):
    """
    Interpolate the MARCS atmosphere and departure coefficients for each NLTE
    element and write the Turbospectrum species file.

    Parameters
    ----------
    abundances : array
        Absolute abundances A(X), index Z-1 (as built in synthesis.synthesize);
        the departure grids are interpolated to A(X) of each element too.
    nlte_elems : list of str

    Returns
    -------
    (atmod, nlte_info_file) : absolute paths; atmod is TS 'interpol' format
    (MARCS-FILE .false.).
    """
    workdir = os.path.abspath(workdir)
    if not os.path.exists(cfg.interpol_exe):
        raise FileNotFoundError(f'{cfg.interpol_exe} not found -- compile it, see nlte.py docstring')
    avail = cfg.available_elements()
    elems = []
    for el in nlte_elems:
        if el in avail:
            elems.append(el)
        else:
            warnings.warn(f'no NLTE data for {el} (need bin/aux/model atom); treating it in LTE')

    models, fam = bracket_models(cfg, teff, logg, mh)
    root = os.path.join(workdir, 'nlte_t{:.0f}g{:.2f}m{:.2f}'.format(teff, logg, mh))
    atmod = root + '.interpol'

    if not elems:
        raise RuntimeError('none of the requested nlte_elems have NLTE data available: '
                           f'{nlte_elems}; available: {avail}')

    for el in elems:
        z = int(np.atleast_1d(atomic.periodic(el))[0])
        # rounded as synthesize() writes it (8.3f): bsyn checks that the departure
        # coefficients were computed for exactly the abundance it uses
        a_x = round(float(abundances[z - 1]), 3) if el.upper() not in ('H',) else 12.0
        cfgtxt = ''.join(f"'{m}'\n" for m in models)
        cfgtxt += f"'{atmod}'\n'{root}.alt'\n'{root}_{el}_coef.dat'\n"
        cfgtxt += f"'{cfg.files[el]['bin']}'\n'{cfg.files[el]['aux']}'\n"
        cfgtxt += f"{cfg.aux_length(el)}\n{teff}\n{logg}\n{mh:.6f}\n{a_x:.6f}\n"
        cfgtxt += ".false.\n.false.\n'{}.test'\n".format(root)
        res = subprocess.run([cfg.interpol_exe], input=cfgtxt.encode(),
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        coef = f'{root}_{el}_coef.dat'
        if res.returncode != 0 or not os.path.exists(coef) or not os.path.exists(atmod):
            raise RuntimeError(f'NLTE interpolation failed for {el}:\n{res.stdout.decode()[-2000:]}')
        if verbose:
            print(res.stdout.decode())
        _fix_nan_departures(coef)

    info = root + '_SPECIES_LTE_NLTE.dat'
    with open(info, 'w') as f:
        f.write('# written by turbofit nlte.py\n#\n')
        f.write("# path for model atom files     ! don't change this line !\n")
        f.write(cfg.model_atom_dir.rstrip('/') + '/\n#\n')
        f.write("# path for departure files      ! don't change this line !\n")
        f.write(workdir.rstrip('/') + '/\n#\n')
        f.write('# atomic (N)LTE setup\n')
        for el in elems:
            z = int(np.atleast_1d(atomic.periodic(el))[0])
            f.write(f"{z}  '{el}'  'nlte'  '{cfg.files[el]['atom']}'  "
                    f"'{os.path.basename(root)}_{el}_coef.dat'  'ascii'\n")
    return atmod, info
