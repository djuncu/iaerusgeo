#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import zip
from AerusL3_BasicImports import *
from AerusL3_Reading import GetProductCfg

SDS_NAME = 'AOD_VIS06'
SDS_NAMES = [ SDS_NAME, 'CM_VIS06', 'Q-Flag' ]
IN_AOD, MAX_AOD = 0.0, 30000.0


def Smoothing(cfg, fpath):
    if cfg['VERBOSE']: print('--> Smoothing of %s' % SDS_NAME)
    gnf = path.splitext(fpath)[0]
    h5f = h5py.File(fpath, 'r+')
    sat = h5f.attrs['Sensors'].split('/')[0]
    dsets = dict(list(zip(SDS_NAMES, [ h5f[key] for key in SDS_NAMES ])))
    dset = h5f[SDS_NAME]
    nx, ny = dset.shape

    if sat[:3] == "MSG": win_size1, win_size2 = 3, 3
    else               : win_size1, win_size2 = 5, 5

    fdsets = dict(list(zip(SDS_NAMES, [ ravel(dsets[key][:]) for key in SDS_NAMES ])))
    fcst   = ravel(cfg['cst_pix'])

    status = smoothing(fdsets['AOD_VIS06'], fdsets['CM_VIS06'], fdsets['Q-Flag'], fcst,
                       nx, ny, win_size1, win_size2, 1., MIN_AOD, MAX_AOD)

    newdset = h5f.create_dataset(SDS_NAME + '_smooth', (nx, ny), dtype=fdsets[SDS_NAME].dtype,
                                 compression=dset.compression, compression_opts=dset.compression_opts)
    for key, val in list(dset.attrs.items()): newdset.attrs[key] = val
    newdset[:] = fdsets[SDS_NAME].reshape((nx, ny))

    h5f.close()
    try   : del h5f, dsets
    except: pass
    gc.collect()
