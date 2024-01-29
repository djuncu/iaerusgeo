#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
from os import path
from glob import glob

from ConfigFile import *
from IcareUtil import *

inputs, geofil = ("inputs", "hdf_sat_reg.cfg")
slot_freq = 15

ecmwf_parm = [[ 6, 0 ], [ 6, 3 ], [ 3, 0 ]] # Analysis at 0, 6, 12, 18 - Forecast at 3, 9, 15, 21

STEPS = [ 24, 3, 24]

def PCF_SetFilesListOld(key, comment, dates, basedir, proddir, pattern, geocorrect=False, datedir=True, yeardir=True):
    if not isinstance(dates, list): dates = [dates]
    comment_lines = comment.split('\n')
    print()
    for l in comment_lines: print('# ' + l)
    print(key + ' = ', end=' ')
    flist = []
    for dt in sorted(dates):
        rep = path.join(path.abspath(basedir), proddir)
        if yeardir: rep = path.join(rep, dt.strftime("%Y"))
        if datedir: rep = path.join(rep, dt.strftime("%Y_%m_%d"))
        for f in sorted(glob(path.join(rep, pattern.replace('<SLOT>', dt.strftime("%Y%m%d%H%M"))))): flist.append(f)
    for f in flist: print("\\\n  " + f, end=' ')    
    print()

                    
def PCF_SetFilesList(key, comment, dt_from, dt_to, basedir, proddir, pattern, freq, geocorrect=False, datedir=True, yeardir=True):
    comment_lines = comment.split('\n')
    print()
    for l in comment_lines: print('# ' + l)
    print(key + ' = ', end=' ')
    flist = []
    dt = dt_from
    while dt <= dt_to:
        rep = path.join(path.abspath(basedir), proddir)
        if yeardir: rep = path.join(rep, dt.strftime("%Y"))
        if datedir: rep = path.join(rep, dt.strftime("%Y_%m_%d"))
        for f in glob(path.join(rep, pattern.replace('<SLOT>', IcareHdfFileNameUTC(dt)))): flist.append(f)
#        dt += timedelta(minutes = slot_freq)
        dt += timedelta(minutes = freq)
    for f in sorted(set(flist)): print("\\\n  " + f, end=' ')    
    print()

                    
def main():
    # Parameters definition
    req_cfg = {
        'ICARE_ID'  : (str, 'i:', 'RFU'),
        'PROD_VER'  : (str, 'v:', 'DEFAULT'),
        'OUTDIR'    : (str, 'O:',  '.'),
        'TMPDIR'    : (str, 'T:',  '.'),
        'MAINDIR'   : (str, 'M:',  '.'),
        'DATE'      : (str, 'd:', None),
        'GEOREG'    : (str, 'g:', None),
        'TRIHOR'    : (int, 't:',    0),
#        'HIST_DAY': (int, 'p:',    0),
        'HIST_DAY': (int, 'p:',    1),
        'MAXHIST'   : (int, 'm:',    1),
        'START'     : (int, 's:',    0),
        'GEOBASEDIR': (str, 'G:',   ''),
        'ANGLEDIR'  : (str, 'a:',   ''),
        'AE1BASEDIR': (str, 'A:',   ''),
        'AE2BASEDIR': (str, 'B:',   ''),
        'AE1DIR'    : (str, '1:',   ''),
        'AE2DIR'    : (str, '2:',   ''),
        'YEARDIR'   : (int, 'Y:',    3),        
        'DATEDIR'   : (int, 'D:',    3),        
        'NRT'       : (int, 'n:', None),        
        'SUBSET'    : (str, 'S:',   ''),
        'EXTENSION' : (str, 'e:',   ''),
        'INST_OUT'  : (int, 'I:',    0),
        'NSLOTS'    : (int, 'N:',    0),
        'ACCELERATE': (int, 'x:',    0),
        'FLOTSAM_FLAG': (int, 'F:',    0),
        'VERBOSE'   : (int, 'V:',    1),        
        }

    # Settings
    cfg = Config(None, [""]).get(req_cfg)

    # Make necessary directories
    if not MakeNewDir(cfg['TMPDIR']): exit(1)

    if cfg['AE2DIR'] is None: cfg['AE2DIR'] = ''
    if cfg['START'] == 0: cfg['HIST_DAY'] = 1
    if cfg['HIST_DAY']: cfg['HIST_DAY'] = 1

    instant = False
    try:
        if len(cfg['DATE']) < 12: raise ValueError()
        dt = datetime.strptime(cfg['DATE'],'%Y%m%d%H%M')
#        if not cfg['TRIHOR']: exit("Daily run: date must be YYYYMMDD")
        if cfg['TRIHOR']: 
            step, priori_step = STEPS[cfg['TRIHOR']], STEPS[1 - cfg['HIST_DAY']]
        else:
            instant = True
            step, priori_step = -1,24
        title_dt = dt
        daily = 0
    except ValueError:
        dt = datetime.strptime(cfg['DATE'],'%Y%m%d')
        if cfg['TRIHOR']: exit("3 hours run: date must be YYYYMMDDhhmm")
        title_dt, cfg['HIST_DAY'] = dt.date(), 1
        daily = 1
        step, priori_step = STEPS[cfg['TRIHOR']], STEPS[1 - cfg['HIST_DAY']]
#    step, priori_step = STEPS[cfg['TRIHOR']], STEPS[1 - cfg['HIST_DAY']]

    print('#####################################################')
    print('# FICHIER DE COMMANDE POUR LE FRAMEWORK FMK_AERUSL3 #')
    print('#####################################################')
    print('# EXECUTION : python FMK_AerusL3.py <this file>     #')
    print('#####################################################')
    print('# - Date : %s' % IcareHdfFileNameUTC(title_dt))
    print('#####################################################')

    print()
    print('# Identificateur interne CGTD. Requis.')
    print('ICARE_ID = ' + cfg['ICARE_ID'])

    print()
    print('# Version produit. Requis.')
    print('PROD_VER = ' + cfg['PROD_VER'])

    print()
    print('# Date de traitement. Requis.')
    print('DATE = ' + cfg['DATE'])

    print()
    print('# Niveau des affichage.')
    print('VERBOSE = %d' % cfg['VERBOSE'])

    print()
    print('# Chemin du répertoire principal de la chaîne')
    print('MAINDIR = ' + cfg['MAINDIR'])

    print()
    print('# Chemin du répertoire de sortie')
    print('OUTDIR = ' + cfg['OUTDIR'])

    print()
    print('# Répertoire de travail temporaire. Facultatif. Valeur par défaut: "."')
    print('TMPDIR = ' + path.abspath(cfg['TMPDIR']))

    print()
    print('# Exécution en mode quasi temps réel ou non (0: non, 1: oui). Requis.')
    print('NRT = ', cfg['NRT'])

    print()
    print('# Région géostationnaire à traiter (MSG+0000 ou MSG+0415). Requis.')
    print('GEOREG = ', cfg['GEOREG'])

    print()
    print('# Nom du fichier modèle de subsetting (recherché dans ancillary).')
    print('# Facultatif. Valeur par défaut: "" (pas de subsetting)')
    if cfg['SUBSET'] is None: cfg['SUBSET'] = ''
    if cfg['SUBSET'] is not None: print('SUBSET = ', cfg['SUBSET'])

    print()
    print('# Extension à ajouter en fin de nom de fichier (utile en cas de subsetting).')
    print("# Facultatif. Valeur par défaut: '' (pas d'extension)")
    if cfg['EXTENSION'] is None: cfg['EXTENSION'] = ''
    print('EXTENSION = ', cfg['EXTENSION'])

    print()
    print('# Start?. Facultatif. Valeur par défaut: 0')
    print('START = %d' % cfg['START'])

    print()
    print('# Mode de traitement (daily ou instantaneous). Facultatif. Valeur par défaut: 1 (daily).')
    print('DAILY_MODE = %d' % daily)

    print()
    print('# Run tri-horaire (1) ou daily (0)?. Facultatif. Valeur par défaut: 0')
    print('TRIHOR = %d' % cfg['TRIHOR'])

    print()
    print('# Information à priori dans les daily plutôt que dans les 3h. Facultatif. Valeur par défaut: 0')
    print('HIST_DAY = %d' % cfg['HIST_DAY'])

    print()
    print("# Délai en jours admissible pour l'information à priori. Facultatif. Valeur par défaut: 1")
    print('MAXHIST = %d' % cfg['MAXHIST'])

    print()
    print("# Sorties instantanées. Facultatif. Valeur par défaut: 0")
    print('INST_OUT = %d' % cfg['INST_OUT'])

    print()
    print("# Mode accéléré. Facultatif. Valeur par défaut: 0")
    print('ACCELERATE = %d' % cfg['ACCELERATE'])

    print()
    print('# Use FLOTSAM as RTM 1 - yes 0 - MSA.')
    print('FLOTSAM_FLAG = ', cfg['FLOTSAM_FLAG'])

    datedir1, yeardir1 = cfg['DATEDIR'] not in (0,2), cfg['YEARDIR'] not in (0,2)
    datedir2, yeardir2 = cfg['DATEDIR'] >= 2, cfg['YEARDIR'] >= 2

    IcareDate = IcareHdfFileNameUTC(cfg['DATE'])

#    dt_from = dt
#    dt_to  = dt + timedelta(hours=step) - timedelta(minutes=1)
#    if instant: tol = slot_freq
#    if instant: tol = 0
#    else: tol = 1
#    dt_from, dt_to = sorted([dt, dt + timedelta(hours=step) - timedelta(minutes=tol)])
    if instant:
        dt_from, dt_to = sorted([dt, dt - timedelta(minutes=cfg['NSLOTS']*slot_freq)])
    else:
        dt_from, dt_to = sorted([dt, dt + timedelta(hours=step) - timedelta(minutes=1)])

    if not cfg['START']:
        dt_ae2 = dt     - timedelta(hours=priori_step)
        dt_ae3 = dt_ae2 - timedelta(hours=priori_step)
        dt_ae1 = dt_ae2 - timedelta(hours=priori_step*(cfg['MAXHIST']-1))
        if cfg['HIST_DAY'] : dt_ae1, dt_ae2, dt_ae3 = dt_ae1.date(), dt_ae2.date(), dt_ae3.date()

        PCF_SetFilesList('CK_K012_IN', 'Fichier "BRDF Model et covariance matrix" du run précédent. Requis sauf si start.',
                         dt_ae2, dt_ae2, cfg['AE2BASEDIR'], 
                         cfg['AE2DIR'], 'SEV_AERUS-INTERNAL*_<SLOT>_*.h5', priori_step*60, datedir=datedir2, yeardir=yeardir2)
        PCF_SetFilesList('ALT_CK_K012_IN', 'Fichiers "BRDF Model et covariance matrix" des runs précédents (-->J-2). Facultatif.',
                         dt_ae1, dt_ae3, cfg['AE2BASEDIR'], 
                         cfg['AE2DIR'], 'SEV_AERUS-INTERNAL*_<SLOT>_*.h5', priori_step*60, datedir=datedir2, yeardir=yeardir2)

    if cfg['ACCELERATE']: channels = ('VIS06', )
    else                : channels = ('VIS06', 'VIS08', 'IR016')

    for chn in channels:
        PCF_SetFilesList(chn, 'Fichiers AE1 de la période à traiter pour la bande %s. Requis.' % chn, dt_from, dt_to, cfg['AE1BASEDIR'], 
                         cfg['AE1DIR'  ], 'GEO_AERUS*L1*_<SLOT>_%s*.h5' % chn, slot_freq, datedir=datedir1, yeardir=yeardir1)
    PCF_SetFilesList('ECMWF', 'Fichiers climatos ECMWF de la période à traiter. Requis.', dt_from, dt_to, cfg['AE1BASEDIR'], 
                     cfg['AE1DIR'  ], 'GEO_AERUS*L1*_<SLOT>_%s*.h5' % 'ECMWF', slot_freq, datedir=datedir1, yeardir=yeardir1)
    PCF_SetFilesList('ANGLES', 'Fichiers des angles solaires et de visée de la période à traiter. Requis.', dt_from, dt_to, cfg['GEOBASEDIR'],
#                     cfg['ANGLEDIR'], 'GEO_L1B-ANGLES*<SLOT>*.hdf5', slot_freq, datedir=datedir1, yeardir=yeardir1)
                     cfg['ANGLEDIR'], 'GEO_L1B-ANGLES*<SLOT>*.hdf*', slot_freq, datedir=datedir1, yeardir=yeardir1)


if __name__== "__main__":
    main()
