"""Build a CNO grid for FERRE (Bob's build_turbogrid.py, with the LSF step).
See turbogrid_user_guide.pdf for the whole workflow (collect, fill holes, write)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import turbogrid as tg

griddir = '/path/to/marcs_apogee'
linelists = ['/path/to/molecular_vac.list', '/path/to/atomic_vac.list', '/path/to/Hlinedata']

lookup = tg.make_atmod_lookup(griddir=griddir, model='marcs', flat=True)
params = tg.make_grid_params(teff_grid=(15, 3000., 250.), logg_grid=(12, -0.5, 0.5),
                             mh_grid=(14, -2.5, 0.25), cm=(11, -1.5, 0.25),
                             am=(7, -1.0, 0.25), nm=(6, -0.5, 0.5),
                             vmicro='holtzman', atmod_lookup=lookup)
print(len(params), 'grid points')

# broaden on the fly with broaden.convolve: constant R on a log-lambda grid ...
bk = dict(R=22500., dlogw_out=6e-6)
# ... or a wavelength-dependent R:      bk = dict(R=(wave_R, R_of_wave), dlogw_out=6e-6)
# ... or an instrument LSF per pixel:    bk = dict(lsf=lsf_array, lsf_offsets=np.arange(-7, 8),
#                                                  wave_out=instrument_wave)
tg.write_grid_table(params, outdir='gridm31', wrange=(15050., 17050.), dw=0.01,
                    linelists=linelists, atmos_type='marcs',
                    outputs_dir='/path/to/big/storage/m31grid', shard_size=1000,
                    broaden_kw=bk)
tg.write_slurm_array('gridm31', jobname='m31grid', time='24:00:00', mem='8G',
                     partition='nodes', array_step=1011, max_concurrent=30,
                     max_array_size=999)

# --- an NLTE grid instead (standard MARCS; [a/M],[C/M] then only change the synthesis):
# import nlte
# cfg = nlte.NLTEConfig('/path/to/marcs_standard', '/path/to/nlte_data')
# params = tg.make_grid_params(..., atmod_lookup=tg.nlte_atmod_lookup)
# tg.write_grid_table(params, ..., nlte=cfg, nlte_elems=['Mg', 'Ca', 'Fe'])
