#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import gc
from datetime import *
from os import path
from numpy import *
from numpy.core.records import fromrecords

try   : import h5py
except:
    print "WARNING - AERUSL3_BasicImports.py: cannot load h5py"
    pass

try:
    from pyhdf.SD import *
    print "Loading pyhdf.SD OK"
except:
    print "WARNING - AERUSL3_BasicImports.py: Cannot load pyhdf.SD"
    pass

from IcareUtil import *
from ConfigFile import *

from AerusL3 import *

AERUS_L3_VER = (1, 1, 8)
AERUS_L3_COLLECTION = "GEO-1.0"

AERUS_L3_DFT_VER = GetIcareFileVersion(AERUS_L3_VER)
MISSION = 'SEV'
PRODUCT_TYPE = 'AERUS'
PRODUCT_LEVEL = '3'

PROCESSING_MODES = ['REA', 'NRT']

ANC_DIR_NAM, CFG_DIR_NAM, PRD_CFG_NAM, ECL_CFG_NAM, BAD_CFG_NAM, AOD_CLM_NAM = \
    "ancillary", "config", "Product.%s.cfg", "solar_eclipses.dat", "invalid_slots.dat", "%s.3km.aer.h5"
ALG_CFG_MODES = [ "", "only_recursion_start", "only_recursion", "only_composition", "recursion_and_composition_start", "recursion_and_composition" ]
ALG_CFG_NAMES = map(lambda x:'Algo_' + x + '.cfg', ALG_CFG_MODES)

FAKE_START_NAME = 'SEV_AERUS-INTERNAL-@NRT@D3_@DATE@_@VERSION@.fake2start.h5'
DTSTR, NRTSTR, VERSTR = '@DATE@', '@NRT@', '@VERSION@'
NRT_STRS = [ '', 'NRT-' ]

# Date and time validity interval for processing
MINDATE, MAXDATE = datetime.strptime("20020101","%Y%m%d"), datetime.strptime("99990101","%Y%m%d")

SLOT_FREQ = 15   # Slot frequency in minutes

AGS_AOD_BNDS = [1, 5000]

SCL_LATLON, OFF_LATLON, BAD_LATLON = 100., 0.0, -32768
FILL_ECMWF = -32768

g_attr_r = dict(L3_Input_Files = set([]), Sensors = set([]), L1_Input_Files = set([]),
                Angles_Input_Files = set([]), ECMWF_Input_Files = set([]))

HdfHandle = None
H4TOH5_CONV_SFX = '_GLOSDS'

# !!! WORKS ONLY WITH MSG+0000 !!!
NL, NC, PRJLON = 3712, 3712, 0.0

ALL_VIS_CHN = [ 'VIS06', 'VIS08', 'IR016']
#VIS_CHN     = [ 'VIS06', 'VIS08', 'IR016']
#VIS_CHN     = ALL_VIS_CHN[:py_ifc.n_channels + 0]
ANGLES_KEYS = [ 'SAA', 'SZA', 'VAA', 'VZA' ]
ECMWF_KEYS = [ 'SWND', 'DWND' ]

prd_brf, prd_saa, prd_sza, prd_vaa, prd_vza = ['BRF'] + ANGLES_KEYS

# For each variable to be read or written set when necessary product (_p), dataset name (_d), title (_t) and data type (_dt)

lat_p, lat_d, lat_t = 'LAT' , 'LAT' , 'latitudes'
lon_p, lon_d, lon_t = 'LON' , 'LON' , 'longitudes'

brf_p, brf_d, brf_t = 'GEO_AERUS-L1', 'BRF-TOC', 'reflectance'
brq_p, brq_d = brf_p, 'BRF_Q_Flag'

angles_p = ANGLES_KEYS
angles_d = [ 'Solar_Azimuth', 'Solar_Zenith', 'View_Azimuth', 'View_Zenith' ]
angles_t = [ 'solar azimuth angle', 'solar zenith angle', 'view azimuth angle', 'view zenith angle' ]

ecmwf_p = ECMWF_KEYS
ecmwf_d = [ 'WIND_SPEED', 'WIND_DIRECTION' ]
ecmwf_t = [ 'wind speed', 'wind direction' ]

par_p = map(lambda x: 'AL-' + x + '-K012', ALL_VIS_CHN)
par_d = [ 'K0', 'K1', 'K2', 'K3' ]
par_d_ = [ 'R0', 'R1', 'R2']
par_t = 'model parameter'
par_dt = int16

cov_p = map(lambda x: 'AL-' + x + '-CK', ALL_VIS_CHN)
cov_d = [ 'C00', 'C01', 'C02', 'C03', 'C11', 'C12', 'C13', 'C22', 'C23', 'C33' ]
cov_t = 'model parameter'
cov_dt = int16

qua_p = "Quality Flag"
qua_d = "Q-Flag"
qua_dt = uint8

age_p = "Age of Information"
age_d = "Z_Age"
age_dt = int8

xi_p = "Scattering angle"
xi_d = "xi"
xi_dt = int16

ref_p = "Surface reflectance VIS06"
ref_d = "ref_VIS06"
ref_dt = int16

jacoAOD_p = "Jacobian AOD VIS06"
jacoAOD_d = "jac_VIS06"
jacoAOD_dt = int16

cm_p = "Confidence measure VIS06"
cm_d = "CM_VIS06"
cm_dt = int8

asp_p = map(lambda x: 'AL-' + x, ALL_VIS_CHN)

aer_p = [ "Aerosol Optical Depth - 635nm", "(beta) Aerosol Optical Depth - 810nm", "(beta) Aerosol Optical Depth - 1640nm",
          "AOD Uncertainty - 635nm"      , "AOD Uncertainty - 810nm"             , "AOD Uncertainty - 1640nm"              ]
aer_d = [ "AOD_VIS06", "AOD_VIS08", "AOD_IR016", "CK_VIS06", "CK_VIS08", "CK_IR016" ]

ags_p = "Angstrom coefficient - 635-810nm"
ags_d = "ANGS_06_08"

alb_p = [ "Spectral Directional-Hemispherical Albedo", "Error of Spectral Directional-Hemispherical Albedo",
          "Spectral Bi-Hemispherical Albedo"         , "Error of Spectral Bi-Hemispherical Albedo"          ]
alb_d = [ "AL-SP-DH", "AL-SP-DH-ERR", "AL-SP-BH", "AL-SP-BH-ERR" ]
alb_dt = int16

bal_p = [ "Shortwave Broadband Directional-Hemispherical Albedo"   , "Error of Shortwave Broadband Directional-Hemispherical Albedo"   ,
          "Shortwave Broadband Bi-Hemispherical Albedo"            , "Error of Shortwave Broadband Bi-Hemispherical Albedo"            ,
          "Visual Broadband Directional-Hemispherical Albedo"      , "Error of Visual Broadband Directional-Hemispherical Albedo"      ,
          "NearInfrared Broadband Directional-Hemispherical Albedo", "Error of NearInfrared Broadband Directional-Hemispherical Albedo" ]
bal_d = [ "AL-BB-DH", "AL-BB-DH-ERR", "AL-BB-BH", "AL-BB-BH-ERR", "AL-VI-DH", "AL-VI-DH-ERR",  "AL-NI-DH", "AL-NI-DH-ERR" ]


# This function copies a python string to a fixed length ' '-padded fortran string
def set_fortran_string(str_out, str_in):
    lth = min(len(str_in), len(str_out))
    str_out[:lth] = str_in[:lth]
    str_out[lth:] = ' '


# This function copies an array of fixed length ' '-padded fortran strings to a single python stripped string
def get_fortran_string(str_arr):
#    return ''.join([str_arr[i] for i in range(str_arr.size)]).strip()
    return str(str_arr)


# This function compares input files names for sorting versus visible channels
def compare_inputs(a,b):
    ia, ib = (None,)*2
    for ikey in range(len(ALL_VIS_CHN)):
        if ALL_VIS_CHN[ikey] in a: ia = ikey
        if ALL_VIS_CHN[ikey] in b: ib = ikey
    return ia - ib


# Utility functions to print uniformly formatted I/O infos 
def print_reading_info(obj, fpath): print_info('Reading', 'from', obj, fpath)
def print_writing_info(obj, fpath): print_info('Writing', 'to  ', obj, fpath)
def print_info(action, direction, obj, fpath):
    print '  - %s %23s %s file %s' % (action, obj, direction, fpath)

def init_vis_chn(cfg):
    if cfg['ACCELERATE']: cfg['NCHANNELS'] = 1
    else                : cfg['NCHANNELS'] = len(ALL_VIS_CHN)
    cfg['VIS_CHN'] = ALL_VIS_CHN[:cfg['NCHANNELS']]
