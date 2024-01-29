#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from os import path

from AerusL3_BasicImports import *
from AerusL3_Reading import *
from AerusL3_Writing import *

ReadLatLon = ReadLatLonRaw

_OUT, _INP, _IN, _INS, _NOPE = range(5)

#FIXED_PRV_FILES = [ 'LAT', 'LON' ] + [ 'AOD_CLM' ]*6 + [ 'CST_PIX' ]
FIXED_PRV_FILES = [ 'LAT', 'LON' ] + [ ('AOD_CLM' + str(i)) for i in range(6) ] + [ 'CST_PIX' ]
#FIXED_PRV_FILES = [ 'LAT', 'LON' ] + [ ( 'AOD_CLM', )*6 ] + [ 'CST_PIX' ]

# This function Checks the impact of dated infos on a slot
def CheckImpact(bad_dates, dt, duration, counter, pos):
    bad = filter(lambda x: x[1] > dt and x[0] < dt + timedelta(minutes=duration), bad_dates)
    chk = len(bad) > 0
    if chk: counter[pos] += 1
    return chk


# This function prepares the lists of input, previous and output files
def SetFilesConfig(cfg):
    global dtchk          # global for use in map() function
    dtchk = 0

    # First check input AERUS-L1 and angles files to determine the available scenes.
    # This checking searches slot date and time in file names: not the best way because it
    # requires a standard naming of these files (what is done in Icare archive) but very more
    # efficient than opening and closing up to 672 HDF files just to read an attribute
    dt0, keys = cfg['DATETIME'], cfg['VIS_CHN'] + [ 'ANGLES', 'ECMWF' ]

    if cfg['VERBOSE']: print "Check for available scenes. "

    # Abort if no available scene
    if None in [cfg[key] for key in keys]: exit("No L1 scene available. Abort")

    files_ok = dict.fromkeys(keys)
    for key in keys: files_ok[key] = []

    sc_times, bad_slts = [], [0, 0]
    for hh in range(cfg['NHOURS']):
#        for mm in map(lambda x:x*SLOT_FREQ, range(60/SLOT_FREQ)):
        for mm in map(lambda x:x*SLOT_FREQ, set(range(60/SLOT_FREQ) + [(1-cfg['SENS'])*2])):
            # For INST, we do not want to repeat XX-00 if temporal window is longer than one hour
            if not cfg['TRIHOR'] and not cfg['DAILY_MODE'] and hh > 0 and mm == 0:
                continue
            else:
                # Build datetime string to search in file names
                #scdt = dt0 + timedelta(hours=hh, minutes=mm)
                scdt = dt0 + cfg['SENS']*timedelta(hours=hh, minutes=mm)
                dtchk = IcareHdfFileNameUTC(scdt)

                # Check impact of solar eclipses and corrupted slots
                if CheckImpact(cfg['ECLIPSES' ], scdt, SLOT_FREQ, bad_slts, 0): continue
                if CheckImpact(cfg['CORRUPTED'], scdt, SLOT_FREQ, bad_slts, 1): continue

                nok, tmpfiles = 0, dict.fromkeys(keys)
                for key in keys:
                    check = map(lambda x:dtchk in x, cfg[key])
                    if any(check): # Datetime string found in at least one file of the list
                        nok += 1
                        tmpfiles[key] = cfg[key][check.index(True)]
                if nok == len(keys): # Datetime string found in all input files list
                    for key in keys: files_ok[key].append(tmpfiles[key]) # Remember matching files
                    sc_times.append(scdt)

    # Dans le cas instantané, il faut que les données du slot traité soitent disponible
    if cfg['DAILY_MODE'] == 0 and cfg['DATETIME'] not in sc_times:
        exit("No L1 scene for instantaneous mode processed slot %s. Abort" % cfg['DATE'])

    for key in keys: cfg[key] = files_ok[key] # Set the actual lists of files to process

    if cfg['VERBOSE']:
        if bad_slts[0] > 0: print "Warning: %d possibly eliminated slots due to solar eclipse"    % bad_slts[0]
        if bad_slts[1] > 0: print "Warning: %d possibly eliminated slots due to known corruption" % bad_slts[1]

    # Now set the generic names of the temporary files to be used by fortran code
    for key in cfg['VIS_CHN']: cfg['tmp' + key] = map(lambda x: path.splitext(x)[0], files_ok[key])
    for key in ANGLES_KEYS   : cfg['tmp' + key] = map(lambda x: path.splitext(x)[0] + '_' + key, cfg['ANGLES'])
    for key in ECMWF_KEYS    : cfg['tmp' + key] = map(lambda x: path.splitext(x)[0] + '_' + key, cfg['ECMWF'])


    try   : nhistfiles = len(cfg['CK_K012_IN'])
    except: nhistfiles = 0
    if cfg['K012_IN'] is None: cfg['K012_IN'] = [cfg['CK_K012_IN']]*(py_ifc.n_channels + 0)
    if cfg['CK_IN'  ] is None: cfg['CK_IN'  ] = [cfg['CK_K012_IN']]*(py_ifc.n_channels + 0)

    for ichn in range(len(cfg['VIS_CHN'])):
        for prd in ['K012', 'CK']:
            for hstf in range(nhistfiles):
                key = prd + '_' + cfg['VIS_CHN'][ichn] + '_' + ('%02d' % hstf)
                if cfg[prd + '_IN'][ichn] is None: cfg['tmp' + key] = None
#                else: cfg['tmp' + key] = path.splitext(cfg[prd + '_IN'][ichn])[0] + '_' + key
                else: cfg['tmp' + key] = path.splitext(cfg[prd + '_IN'][ichn][hstf])[0] + '_' + key

    if not cfg['DAILY_MODE'] : lvl, fdt = 'L3', cfg['DATETIME']
    else                     : lvl, fdt = 'D3', cfg['DATETIME'].date()
    lvl = NRT_STRS[cfg['NRT']] + lvl

    gen_out = IcareHdfFileName('SEV', 'AERUS', '%s' , fdt, None, cfg['PROD_VER'], None)
    alt_out = IcareHdfFileName('SEV', 'AERUS', '%s' , fdt, '%s', cfg['PROD_VER'], None)
    aer_out, unf_out = gen_out % ('AEROSOL-' + lvl), gen_out % ('INTERNAL-' + lvl)

    gen_ins = IcareHdfFileName('SEV', 'AERUS', 'INST-' + lvl , '%s', None, cfg['PROD_VER'], None)

    # Output Albedos as a intrinsic product ?
    if not cfg['PROD_MODE']: cfg['ALB_OFF_PRD'] = True
    if   cfg['ALB_OFF_PRD']: alb_out = gen_out % ('ALBEDO-' + lvl)
    else                   : alb_out = unf_out

    for chn in cfg['VIS_CHN']:
        for prd in [ 'ALK012', 'ALCK', 'AL' ]:
            prd_out = unf_out
            if prd == 'AL': prd_out = alb_out
            key = prd + '_' + chn
            cfg['out' + key] = path.splitext(prd_out)[0] + '_' + key
    cfg['outALBEDO' ] =  path.splitext(unf_out)[0] + '_' + 'ALBEDO'

    for outkey in [ 'AL-K012', 'AL-CK', 'AL' ]:
        cfg[outkey] = []
        prd_out = unf_out
        if outkey == 'AL': prd_out = alb_out
        for key in cfg['VIS_CHN']:
            if cfg['PROD_MODE']: cfg[outkey].append(prd_out)
            else: cfg[outkey].append(alt_out % (outkey + '-L2', key))

    cfg['INSTANT'], cfg['INST_TIMES'] = [], sc_times
    for scn in range(len(sc_times)):
        cfg['INSTANT'].append(gen_ins % IcareHdfFileNameUTC(sc_times[scn]))
        
    cfg['AEROSOL'] = aer_out
    cfg['ALBEDO']  = alb_out

    # Return the number of daily available scenes
    return len(cfg[keys[0]])


def FilesPrvReduce(cfg):
    files = []
    for l in range(py_ifc.nhistfiles):
        for k in cfg['VIS_CHN']:
            for key in ('K012', 'CK'):
                f = 'tmp' + key + '_' + k + ('_%02d' % l)
                files.append(f)
    return files


def ReordonInFiles(files, offset):
    newfiles = files[:offset]
    oldfiles = files[offset:]
    for k in range(py_ifc.nhistfiles):
        for i in range(py_ifc.n_channels):
            for j in range(2):
                kkk = k + py_ifc.nhistfiles*(2*i +j)
                newfiles.append(oldfiles[kkk])
    return newfiles


# This function sends algorithm and product configuration to the Fortran 90 code part
def interface_algo(alg_cfg, cfg):

    # Send algorithm configuation
    for key, val in alg_cfg.items():
        if isinstance(val, str):
#            exec("py_ifc.%s[:%d] = '%s'" % (key.lower(), len(val), val))
#            exec("py_ifc.%s[%d:] = ' '"  % (key.lower(), len(val)))
            exec("strlth=py_ifc.%s.dtype.itemsize" % key.lower())
            exec("py_ifc.%s = '%s%s'" % (key.lower(), val, ' '*(strlth - len(val))))
        elif isinstance(val, list):
            exec("py_ifc.%s[:%d] = %s" % (key.lower(), len(val), str(val)))
        else:
            exec("py_ifc.%s = %s" % (key.lower(), str(val)))

    # Send product configuration including temporary file names arrays
    py_ifc.mode          = cfg['MODE']
    py_ifc.ocean_ok      = cfg['OCEAN_OK']
    py_ifc.land_ok       = cfg['LAND_OK']
    py_ifc.coast_ok      = cfg['COAST_OK']
    py_ifc.instoutput    = cfg['INST_OUT'] != 0
    py_ifc.instantaneous = cfg['DAILY_MODE'] == 0
    try   : py_ifc.nhistfiles = len(cfg['CK_K012_IN'])
    except: py_ifc.nhistfiles = 0

    files_prv_ok  = SetFileArray( _INP, cfg, FIXED_PRV_FILES + FilesPrvReduce(cfg))
    files_prv_in  = SetFileArray(_NOPE, cfg, FIXED_PRV_FILES + [ ('K012_IN', 'CK_IN') ])
    files_prv_in  = ReordonInFiles(files_prv_in, len(FIXED_PRV_FILES))

    files_inp_ok  = SetFileArray(  _IN, cfg, [ tuple( map(lambda x: 'tmp' + x, cfg['VIS_CHN'] + ANGLES_KEYS + ECMWF_KEYS) ) ])
    files_inp_in  = SetFileArray(_NOPE, cfg, [ tuple( cfg['VIS_CHN'] + ['ANGLES']*len(ANGLES_KEYS) + ['ECMWF']*len(ECMWF_KEYS) ) ])

    files_out_ok  = SetFileArray( _OUT, cfg, reduce(lambda x,y: x+y, [['outALK012_' + k, 'outALCK_' + k, 'outAL_' + k]
                                                                      for k in cfg['VIS_CHN'] ]) + [ 'outALBEDO', 'AEROSOL' ])
    files_out_out = SetFileArray(_NOPE, cfg, [ ('AL-K012', 'AL-CK', 'AL'), 'ALBEDO', 'AEROSOL' ])

    files_ins_ok  = SetFileArray( _INS, cfg, [ ('INSTANT',) ])

    # Return actual and temporary file names arrays for previous, input and output files
    return (files_prv_ok, files_prv_in, files_inp_ok, files_inp_in, files_out_ok, files_out_out, files_ins_ok)


# This function sets temporary files names to be used by the Fortran 90 code part
def SetFileArray(inputs, cfg, filseq):
    lstfiles = []
    for ils in range(len(filseq)):
        lsf = filseq[ils]
        if not isinstance(lsf, tuple): lstfiles.append(cfg[lsf])
        else:
            for nseries in range(len(cfg[lsf[0]])):
                for ifil in range(len(lsf)):
#                    lstfiles.append(cfg[lsf[ifil]][nseries])
                    fff = cfg[lsf[ifil]][nseries]
                    if not isinstance(fff, list): fff = [fff]
                    lstfiles += fff
    if inputs != _NOPE:
        for ils in range(len(lstfiles)): py_ifc.setfilearray(ils, inputs, lstfiles[ils])
    return lstfiles


# This function sets some essential variables for the Fortran 90 code part
def interface_code(cfg, gattr):
    if cfg['SUBSET'] is not None :
        ok = cfg['SUBSET_GRID'][0]
        py_ifc.msgpixy, py_ifc.msgpixx = ok.shape
    else:
        py_ifc.msgpixx, py_ifc.msgpixy = gattr['NC'], gattr['NL']
    start_mod.n_pix       = py_ifc.msgpixx *  py_ifc.msgpixy
    start_mod.n_blocks_f  = py_ifc.msgpixy // py_ifc.linesblock
    start_mod.linesrest   = py_ifc.msgpixy %  py_ifc.linesblock
    start_mod.n_blocks    = start_mod.n_blocks_f
    if start_mod.linesrest > 0: start_mod.n_blocks += 1


# This function sets some global output attributes from selected input files global attributes
def interface_attrs(gattl, gattr, gatta, gatte):
    for key in ('Sensors', 'L3_Input_Files', 'L1_Input_Files', 'Angles_Input_Files', 'ECMWF_Input_Files'):
        if len(g_attr_r[key]) == 0 : g_attr_r[key] = 'None'
        else                       : g_attr_r[key] = ','.join(sorted(g_attr_r[key]))
    g_attr_r.update(gattl)
    for key in set(['Sensors', 'Ancillary_Files', 'Beginning_Acquisition_Date', 'End_Acquisition_Date']) & set(gattr.keys()): del gattr[key]
    for key in set(['Sensors', 'Ancillary_Files', 'Beginning_Acquisition_Date', 'End_Acquisition_Date']) & set(gatta.keys()): del gatta[key]
    g_attr_r.update(gattr)
    for key in set(['Date', 'Ancillary_Files', 'Beginning_Acquisition_Date', 'End_Acquisition_Date', 'Contact' ]) & set(gatta.keys()): del gatta[key]
    for key in set(['Date', 'Ancillary_Files', 'Beginning_Acquisition_Date', 'End_Acquisition_Date', 'Contact' ]) & set(gatta.keys()): del gatta[key]
    g_attr_r.update(gatta)


def AerusL3(cfg):

    # INITIALIZE
    if cfg['VERBOSE']: print "Initialisations"

    # Try to build temporary directory if none
    if not MakeNewDir(cfg['TMPDIR']): exit(1)

    # Determine algorithm configuration mode
    alg_mode = cfg['RECURSION'] + 3*cfg['COMPO']
    if alg_mode == 0 or alg_mode >= len(ALG_CFG_NAMES): exit("Error: invalid algorithm configuration mode: %d" % alg_mode)
    if cfg['VERBOSE']: print "Algorithm configuration: " + ALG_CFG_MODES[alg_mode]

    # Check parameters and algorithm mode consistency
    if cfg['RECURSION'] > 1:
        if cfg['CK_K012_IN'] is None:
            for key in [ 'CK_IN', 'K012_IN' ]:
                if cfg[key] is None: exit ('Error: CK_K012_IN or %s parameter must exist' % key)
                cfg[key] = cfg[key].split()
                cfg[key].sort(compare_inputs)
        else:
            for key in [ 'CK_IN', 'K012_IN' ]:
                if cfg[key] is not None: exit ('Error: only one of CK_K012_IN and %s parameter can exist' % key)

    if cfg['VERBOSE']: print "Read configuration files"

    # Read algo technical parameters
    alg_cfg_pth = path.join(cfg['MAINDIR'], CFG_DIR_NAM, ALG_CFG_NAMES[alg_mode])
    algo_cfg = GetAlgoCfg(cfg, alg_cfg_pth)
    if algo_cfg is None: exit(1)
    algo_cfg['path_tmp'] = cfg['TMPDIR'] + '/' # Pour les fichiers temporaires Fortran

    # Get region name from first band
#    georeg = GetRegionFromBandFile(cfg['VIS06'][0])

    # Read product configuration
#    product_cfg = GetProductCfg(cfg, path.join(cfg['MAINDIR'], CFG_DIR_NAM, PRD_CFG_NAM % georeg))
    product_cfg = GetProductCfg(cfg, path.join(cfg['MAINDIR'], CFG_DIR_NAM, PRD_CFG_NAM % cfg['GEOREG']))
    if product_cfg is None: exit(1)
    cfg.update(product_cfg)
    if   cfg['TRIHOR']    : sens, nh, dt =  1,  3, datetime.strptime(cfg['DATE'], "%Y%m%d%H%M")
    elif cfg['DAILY_MODE']: sens, nh, dt =  1, 24, datetime.strptime(cfg['DATE'], "%Y%m%d")
    # temporal window set to 2 hours: TO BE IMPROVED as it should depend on NSLOTS
    else                  : sens, nh, dt = -1,  2, datetime.strptime(cfg['DATE'], "%Y%m%d%H%M")
    cfg.update(dict(Input_Files=[], Ancillary_Files=[], DATETIME=dt, NHOURS=nh, SENS=sens))

    # Complete ancillary files pathes
    for key in ['LAT', 'LON', 'CST_PIX']: cfg[key] = path.join(cfg['ANC_DIR'], cfg[key])
    cfg['AOD_CLM'] = path.join(cfg['ANC_DIR'], AOD_CLM_NAM % cfg['GEOREG'])
    for i in range(6): cfg['AOD_CLM' + str(i)] = path.join(cfg['ANC_DIR'], (AOD_CLM_NAM % cfg['GEOREG']) + '.' + str(i))

    # Set input and output files configuration
    py_ifc.n_channels = cfg['NCHANNELS']
    py_ifc.n_files_in_scene = py_ifc.n_channels + 4 + 2
    py_ifc.n_scenes = SetFilesConfig(cfg)
    if py_ifc.n_scenes == 0: exit("No complete input scenes available! No product generated!")
    if cfg['VERBOSE']: print 'Number of complete input scenes: %d' % py_ifc.n_scenes

    # Read or initialize subsetting mask
    cfg['SUBSET_MSK'] = ReadSubsetMask(cfg)
    msk_nlin = cfg['SUBSET_GRID'][0].shape[0]
    if prod(cfg['SUBSET_GRID'][0].shape) < algo_cfg['LinesBlock']*NC : algo_cfg['LinesBlock'] = msk_nlin
    elif msk_nlin < algo_cfg['LinesBlock']: algo_cfg['LinesBlock'] = msk_nlin

    # Initialize algo configuration Fortran90 interface
    if cfg['VERBOSE']: print 'Initialize algo configuration Fortran90 interface'
    files_prv_ok, files_prv_in, files_inp_ok, files_inp_in, files_out_ok, files_out_out, files_inst = interface_algo(algo_cfg, cfg)

    if cfg['VERBOSE']: print "Read ancillary files"

    # Read latitudes and longitudes files
    gattlat, lats, lons = ReadLatLon(cfg)
    if lats is None: exit(1)
    lats[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[0])[1] + '.gb2'))
    lons[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[1])[1] + '.gb2'))

    # Read AOD climatology file
#    aer_mod, tot_aod = ReadAodClimato(cfg)
#    aer_mod[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[2])[1] + '.gb1'))
#    tot_aod[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[3])[1] + '.gb2'))
    aer_mod1, aer_mod2, aer_mod3, tot_aod, bdw_mod2, bdw_mod3 = ReadAodClimato(cfg)
    aer_mod1[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[2])[1] + '.gb1'))
    aer_mod2[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[3])[1] + '.gb1'))
    aer_mod3[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[4])[1] + '.gb1'))
    tot_aod [cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[5])[1] + '.gb2'))
    bdw_mod2[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[6])[1] + '.gb2'))
    bdw_mod3[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[7])[1] + '.gb2'))

    # Read coast pixels file
    cst_pix = ReadCoastPixels(cfg, algo_cfg['coast_size'])
#    cst_pix[cfg['SUBSET_GRID']].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[4])[1] + '.gb1'))
    cfg['cst_pix'] = cst_pix[cfg['SUBSET_GRID']]
    cfg['cst_pix'].tofile(path.join(cfg['TMPDIR'], path.split(files_prv_in[8])[1] + '.gb1'))

    prv_off = len(FIXED_PRV_FILES)
    py_ifc.prv_off = prv_off + 1
    g_attr_r['Ancillary_Files'] = ','.join(sorted(set(map(IcarePath, files_prv_in[:prv_off]))))

    # Initialize Fortran90 code interface
    interface_code(cfg, gattlat)
    data_p = [brf_p]*(py_ifc.n_channels + 0) + angles_p + ecmwf_p
    data_d = [brf_d]*(py_ifc.n_channels + 0) + angles_d + ecmwf_d
    data_t = [brf_t]*(py_ifc.n_channels + 0) + angles_t + ecmwf_t

    # Read previous estimate for model parameters and covariance matrix if any
    proc_mode = PROCESSING_MODES[cfg['NRT']]
    current_prodver = "%d.%d" % GetIcareVersionFromStr(cfg['PROD_VER'])
    if py_ifc.recursion and not py_ifc.startseries:

        if cfg['VERBOSE']: print "Read previous memory product files"

        for k in range(py_ifc.nhistfiles):

            for i in range(py_ifc.n_channels):

                for j in range(py_ifc.mmp1):

                    # Read and write parameter
                    ret = ReadAndWriteBin(cfg, files_prv_in, files_prv_ok, 2*(k*py_ifc.n_channels + i) + prv_off,
                                          par_d[j], par_p[i], par_t, 'K%cin' % chr(48+j))
                    if ret is None: exit(1)
                    py_ifc.scale_par_in[i,j], py_ifc.offset_par_in[i,j], py_ifc.missing_par_in[i,j], gatt = ret

                    if j != py_ifc.mmp1-1:
                       ret = ReadAndWriteBin(cfg, files_prv_in, files_prv_ok, 2*(k*py_ifc.n_channels + i) + prv_off,
                                             par_d_[j], par_p[i], par_t, 'R%cin' % chr(48+j))
                       if ret is None: exit(1)

                    # Check processing mode consistency
                    if proc_mode is None: proc_mode = gatt['Processing_Mode']
                    elif gatt['Processing_Mode'] != proc_mode:
                        exit("Error: Bad processing mode %s - %s expected: %s" % (gatt['Processing_Mode'], proc_mode,
                                                                                  files_prv_in[2*(k*py_ifc.n_channels + i) + prv_off]))

                    # Check product version consistency
                    if gatt['Product_Version'] != current_prodver:
                        exit("Error: Bad product version %s - %s expected: %s" % (gatt['Product_Version'], current_prodver,
                                                                                  files_prv_in[2*(k*py_ifc.n_channels + i) + prv_off]))

                # Read and write processing flag
                ret = ReadAndWriteBin(cfg, files_prv_in, files_prv_ok, 2*(k*py_ifc.n_channels + i) + prv_off,
                                      qua_d, par_p[i], par_t, 'K012.qua_in')
                if ret is None: exit(1)
                py_ifc.missing_qua = ret[2]

                # Read and write age of information
                ret = ReadAndWriteBin(cfg, files_prv_in, files_prv_ok, 2*(k*py_ifc.n_channels + i) + prv_off,
                                      age_d, par_p[i], par_t, 'K012.age_in')
                if ret is None: exit(1)
                py_ifc.missing_age = ret[2]

                g_attr_r['L3_Input_Files'].add(IcarePath(files_prv_in[2*(k*py_ifc.n_channels + i) + prv_off]))

                j = 0
                for n in range(py_ifc.mmp1):
                    for nn in range(n, py_ifc.mmp1):
                        # Read and write covariance matrix element
                        ret = ReadAndWriteBin(cfg, files_prv_in, files_prv_ok, 2*(k*py_ifc.n_channels + i) + prv_off + 1,
                                              cov_d[j], cov_p[i], cov_t, 'E%cin' % chr(48+j+1))
                        if ret is None: exit(1)
                        if k == 0: py_ifc.scale_cov_in[i,j], py_ifc.offset_cov_in[i,j], py_ifc.missing_cov_in[i,j], gatt = ret
                        g_attr_r['L3_Input_Files'].add(IcarePath(files_prv_in[2*(k*py_ifc.n_channels + i) + prv_off + 1]))
                        j = j + 1

            if k == 0:
                # Extract the date of the previous result from the image aquisition time attribute
                try:
                    if cfg['TRIHOR'] and not cfg['HIST_DAY']:
                        prev_date = FromIcareHdfUTC(gatt['Date'])
                        if prev_date < MINDATE or prev_date > MAXDATE: raise
                    else:
                        prev_date = FromIcareHdfDate(gatt['Date'])
                        if prev_date < MINDATE.date() or prev_date > MAXDATE.date(): raise
                except:
                    print "Warning: invalid date of previous result: " + gatt['Date']
                    set_fortran_string(py_ifc.o_qua_flag, "NOK")

                start_mod.year_in, start_mod.month_in, start_mod.day_of_month_in = prev_date.year, prev_date.month, prev_date.day

    # Send scale, offset and missing_value parameters through the py_ifc interface
    py_ifc.scale_ref   = zeros((py_ifc.n_scenes,), dtype=float32)
    py_ifc.offset_ref  = zeros((py_ifc.n_scenes,), dtype=float32)
    py_ifc.missing_ref = zeros((py_ifc.n_scenes,), dtype=int32)
    for i in range(py_ifc.n_channels, py_ifc.n_channels + 6):
        exec ("py_ifc.scale_%s   = zeros((py_ifc.n_scenes,), dtype=float32)" % (data_p[i].lower()))
        exec ("py_ifc.offset_%s  = zeros((py_ifc.n_scenes,), dtype=float32)" % (data_p[i].lower()))
        exec ("py_ifc.missing_%s = zeros((py_ifc.n_scenes,), dtype=int32  )" % (data_p[i].lower()))

    if cfg['VERBOSE']: print "Read complete scenes input data files"

    # Read input files and write temporary files used by F90 code
    for n in range(py_ifc.n_scenes):

        for i in range(py_ifc.n_channels):

            # Read and write reflectances
            ret = ReadAndWriteBin(cfg, files_inp_in, files_inp_ok, n*py_ifc.n_files_in_scene + i,
                                  data_d[i], data_p[i], data_t[i], 'c%d' % (i+1))
            if ret is None: exit(1)
            py_ifc.scale_ref[n], py_ifc.offset_ref[n], py_ifc.missing_ref[n], gattref = \
                (ret[0] * py_ifc.scale_conv_ref,) + ret[1:]

            # Check processing mode consistency
            if proc_mode is None: proc_mode = gattref['Processing_Mode']
            elif gattref['Processing_Mode'] != proc_mode:
                exit("Error: Bad processing mode %s - %s expected: %s" %
                     (gattref['Processing_Mode'], proc_mode, files_inp_in[n*py_ifc.n_files_in_scene + i]))

            # Read and write quality flags
            if ReadAndWriteBin(cfg, files_inp_in, files_inp_ok, n*py_ifc.n_files_in_scene + i,
                               brq_d, None, None, 'c%d.qua' % (i+1)) is None: exit(1)

            g_attr_r['L1_Input_Files'].add(IcarePath(files_inp_in[n*py_ifc.n_files_in_scene + i]))
            if i== 0:
                if n == 0: g_attr_r['Beginning_Acquisition_Date'] = gattref['Beginning_Acquisition_Date']
                elif n + 1 == py_ifc.n_scenes: g_attr_r['End_Acquisition_Date'] = gattref['End_Acquisition_Date']
                g_attr_r['Sensors'].add(gattref['Sensors'])

        # Read and write angles
        for i in range(py_ifc.n_channels, py_ifc.n_channels + 4):
            ret = ReadAndWriteBin(cfg, files_inp_in, files_inp_ok, n*py_ifc.n_files_in_scene + i, data_d[i], data_p[i], data_t[i], '')
            if ret is None: exit(1)
            exec ("py_ifc.scale_%s[n], py_ifc.offset_%s[n], py_ifc.missing_%s[n] = ret[:3]" % 
                  tuple(map(lambda x:x.lower(), (data_p[i], data_p[i], data_p[i]))))
            gattang = ret[3]
            g_attr_r['Angles_Input_Files'].add(IcarePath(files_inp_in[n*py_ifc.n_files_in_scene + i]))

        # Read and write wind data
        for i in range(py_ifc.n_channels + 4, py_ifc.n_channels + 6):
            ret = ReadAndWriteBin(cfg, files_inp_in, files_inp_ok, n*py_ifc.n_files_in_scene + i, data_d[i], data_p[i], data_t[i], '')
            if ret is None: exit(1)
            exec ("py_ifc.scale_%s[n], py_ifc.offset_%s[n], py_ifc.missing_%s[n] = ret[:3]" % 
                  tuple(map(lambda x:x.lower(), (data_p[i], data_p[i], data_p[i]))))
            gattecm = ret[3]
            g_attr_r['ECMWF_Input_Files'].add(IcarePath(files_inp_in[n*py_ifc.n_files_in_scene + i]))

    # set various attributes according to the last scene
    interface_attrs(gattlat, gattref, gattang, gattecm)
    ladate = datetime.strptime(cfg['DATE'][:8], "%Y%m%d")
    start_mod.year, start_mod.month, start_mod.day_of_month = ladate.year, ladate.month, ladate.day

    # GO! Lauch F90 scientific software
    if cfg['VERBOSE']: print "Lauch Aerosol and Albedo retrieval (blocks of %d rows)" % py_ifc.linesblock
    if cfg['M_A_P'] < 2:
        status = start_exe()
        if cfg['VERBOSE']: print '\n  Retrieval completed'

    # Write outputs
    if py_ifc.recursion :
        stat_type = "recursive, timescale: %ddays" % py_ifc.timescale

        if cfg['TRIHOR'] or cfg['DAILY_MODE']:
            if cfg['VERBOSE']: print "Write official output products"
            WriteOfficialProducts  (cfg, stat_type, files_out_out, files_out_ok)
            if cfg['VERBOSE']: print "Write unofficial and/or memory output products"
            WriteUnofficialProducts(cfg, stat_type, files_out_out, files_out_ok, compo=False)

        if py_ifc.instoutput:
            if cfg['VERBOSE']: print "Write instantaneous estimation products"
            files_inst = sorted(files_inst)
            fidx = range(len(files_inst))
            if not cfg['TRIHOR'] and not cfg['DAILY_MODE']: fidx = fidx[-1:]
            for k in fidx:
                WriteInstantaneous(cfg, stat_type, files_inst[k], files_inst[k], cfg['INST_TIMES'][k], do_close=True)

# ############################ R F U ############################
#    if py_ifc.composition :
#        stat_type = "composition period: 1day"
#        # Write BRDF-model output files
#        WriteBRDFModel(cfg, stat_type, files_out, py_ifc.n_files_out_channel_compo, py_ifc.n_fileout_rec, True)
#        # Write covariance matrix output files
#        WriteCovarianceMatrix(cfg, stat_type, files_out, py_ifc.n_files_out_channel_compo, py_ifc.n_fileout_rec)

    if cfg['VERBOSE']: print "Execution OK"


def main():
    # Parameters definition
    req_cfg = {
        'ICARE_ID'      : (str, 'i:',  "RFU"),
        'PROD_VER'      : (str, 'v:',   None),
        'MAINDIR'       : (str, 'M:',    '.'),
        'OUTDIR'        : (str, 'O:',    '.'),
        'TMPDIR'        : (str, 'T:', '/tmp'),
        'VERBOSE'       : (int, 'V:',      1),
        'DATE'          : (str, 'd:',   None),
        'GEOREG'        : (str, 'g:',   None),
        'NRT'           : (int, 'n:',   None),
        'START'         : (int, 's:',      0),
        'DAILY_MODE'    : (int, 'D:',      1),
        'TRIHOR'        : (int, 't:',      0),
        'CK_K012_IN'    : (str, 'c:',     ''),
        'ALT_CK_K012_IN': (str, 'C:',     ''),
        'MAXHIST'       : (int, 'm:',      1),
        'HIST_DAY'      : (int, 'h:',      1),
        'SUBSET'        : (str, 'S:',     ''),
        'EXTENSION'     : (str, 'E:',     ''),
        'INST_OUT'      : (int, 'I:',      0),
        'SUBSUBSET'     : (str, 'Y:',     ''),
        'SUBTRANSP'     : (bool,'X:',  False),
        'ANGLES'        : (str, 'a:',   None),
        'VIS06'         : (str, '6:',     ''),
        'VIS08'         : (str, '8:',     ''),
        'IR016'         : (str, '1:',     ''),
        'ECMWF'         : (str, 'e:',     ''),
        'SMOOTH'        : (bool,'l:',   True),
        'ACCELERATE'    : (bool,'x:',  False),
#        'RECURSION'     : (int, 'R:',      0),
#        'COMPO'         : (int, 'Z:',      0),
#        'M_A_P'         : (int, 'm:',      0),
        }

    # Settings
    config = Config(None, [''])
    cfg = config.get(req_cfg)

   # In this version:
    # - activate 3-hourly mode
    # - activate instantaneous mode
    # - activate multi-products mode for earliest information
    # - deactivate debug mode for production
    # - deactivate composition mode 

    cfg['M_A_P' ], cfg['COMPO'   ] = 0, 0

    # Complete pathes
    cfg['ANC_DIR'] = path.join(cfg['MAINDIR'], ANC_DIR_NAM)

    # Set product level label, product and product name
    if cfg['DAILY_MODE']: cfg['PRODUCT_LEVEL'] = 'D' + PRODUCT_LEVEL
    else                : cfg['PRODUCT_LEVEL'] = 'L' + PRODUCT_LEVEL

    if not cfg['TRIHOR'] and not cfg['DAILY_MODE']: cfg['INST_OUT'] = True
    elif cfg['START']                             : cfg['INST_OUT'] = False

    if cfg['NRT']:
        cfg['NRT'] = 1
        cfg['PRODUCT'] = '-'.join([PRODUCT_TYPE, 'NRT', cfg['PRODUCT_LEVEL']])
    else         :
        cfg['NRT'] = 0
        cfg['PRODUCT'] = '-'.join([PRODUCT_TYPE, cfg['PRODUCT_LEVEL']])
    cfg['PRODUCT_NAME'] = '_'.join([MISSION, cfg['PRODUCT']])

    # Set channels to process (3 or only 1 if ACCELERATE is True)
    init_vis_chn(cfg)

    if cfg['SUBSUBSET'] is not None:
        cfg['SUBSUBSET'] = map(int, cfg['SUBSUBSET'].split(','))

    # Set product version if not specified
    if cfg['PROD_VER'] == ICA_DFT_VER_KEY: cfg['PROD_VER'] = AERUS_L3_DFT_VER
    if cfg['VERBOSE']: print "\n%s v%d.%d.%d\n" % ((cfg['PRODUCT'],) + AERUS_L3_VER)

    # Check if processing true start
    if cfg['START'] > 0 :
        cfg['RECURSION'] = 1
    else:
        # Check if processing fake start
        if cfg['CK_K012_IN'] is None:
            version  = 'V' + GetIcareFileVersion(GetIcareVersionFromStr(cfg['PROD_VER']))
            prevdate = datetime.strptime(cfg['DATE'], "%Y%m%d") - timedelta(days=1)
            cfg['CK_K012_IN'] = path.join(cfg['ANC_DIR'], FAKE_START_NAME.replace(NRTSTR, NRT_STRS[cfg['NRT']]). \
                                              replace(VERSTR, version).replace(DTSTR, IcareHdfFileNameUTC(prevdate.date())))
            if not path.exists(cfg['CK_K012_IN']):
                exit("Cannot process with START=0 and without CK_K012_IN parameter nor fakestart file for this date")
        else:
            # Find most recent valid supplied history file if possible and check his date
            cfg['CK_K012_IN'] = cfg['CK_K012_IN'].split()
            if cfg['ALT_CK_K012_IN'] is not None: cfg['CK_K012_IN'] += sorted(cfg['ALT_CK_K012_IN'].split())[-1::-1]
            eve_ok, hist_file_ok, delay = False, None, -1
            try   : ladate  = datetime.strptime(cfg['DATE'], "%Y%m%d")
            except: ladate  = datetime.strptime(cfg['DATE'], "%Y%m%d%H%M")
            prevdate = ladate - timedelta(days=cfg['MAXHIST'])
            hist_file_ok, age_hist_file = None, []
            for hist_file in cfg['CK_K012_IN']:
                hist_date, pm_flag = GetHistoryDate(hist_file)
                age_hist_file.append((ladate.date() - hist_date).days)
#                if hist_date == ladate.date() - timedelta(days=1): eve_ok = True
                if age_hist_file[-1] == 1: eve_ok = True
                if pm_flag == 'ERR': continue
                if hist_file_ok is None: 
                    if hist_date < prevdate.date(): exit("Most recent history file is too old")
                    hist_file_ok, delay = hist_file, age_hist_file[-1]
#                break
            if not eve_ok: exit("No history file supplied for the eve (even a fake one)")
            if hist_file_ok is None: exit("No valid history file supplied")
            if delay > 1: print "Warning: history file older than 1 day %s" % hist_file_ok
#            cfg['CK_K012_IN'] = hist_file
            if cfg['VERBOSE']: print "History file: %s\nAge: %d\n" % (hist_file, delay)
            py_ifc.age_hist_file[:len(age_hist_file)] = array(age_hist_file)
            for kkk in range(py_ifc.age_hist_file.size): print kkk, py_ifc.age_hist_file[kkk]

        cfg['RECURSION'] = 2
    cfg['K012_IN'], cfg['CK_IN'] = None, None

    # Make input files lists from parameter strings
    for key in ('ANGLES', 'ECMWF', 'VIS06', 'VIS08', 'IR016') :
        if cfg[key] is not None: cfg[key] = sorted(cfg[key].split())

    # Get solar eclipses and known corrupted slots data
    cfg['ECLIPSES'] = GetSolarEclipseCfg(cfg, path.join(cfg['ANC_DIR'], ECL_CFG_NAM))
    cfg['CORRUPTED'] = GetCorruptedSlots(cfg, path.join(cfg['ANC_DIR'], BAD_CFG_NAM))
    if None in (cfg['ECLIPSES'], cfg['CORRUPTED']): exit("Abort.")

    # Run main code. From now on, any handled failure will produce a fake history file.
    # This prevents from operational production failure.
    try: 
        AerusL3(cfg)

    except SystemExit as err:
        print err
        if cfg['DAILY_MODE']:
            print "Writing fake history file"
            WriteFakeHistoryFile(cfg)


if __name__== '__main__':
    main()
