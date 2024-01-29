#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
#from builtins import range
from AerusL3_BasicImports import *

# Product configuration parameters
PRD_REQ_CFG = {
    'MODE'       : (int  , '', None),
    'OCEAN_OK'   : (bool , '', None),
    'LAND_OK'    : (bool , '', None),
    'COAST_OK'   : (bool , '', None),
    'LAT'        : (str  , '', None),
    'LON'        : (str  , '', None),
    'CST_PIX'    : (str  , '', None),
    'FREQINSECS' : (int  , '', None),
    'PROD_MODE'  : (bool , '', None),
    'ALB_OFF_PRD': (bool , '', None),
    'CONTACT'    : (str  , '', None),
    'REFERENCE'  : (str  , '', None),
    'REF_URL'    : (str  , '', None),
}

# Some useful constants
RD_DAT, RD_GATT, RD_DATT = [set([x]) for x in range(3)]
RD_ATT, RD_SET, RD_ALL = RD_GATT | RD_DATT, RD_DAT | RD_DATT, RD_DAT | RD_GATT | RD_DATT
CHK_PRD, CHK_SIZ = [set([x]) for x in range(2)]
CHK_NOT, CHK_ALL = set([]), CHK_PRD | CHK_SIZ


# This function reads and returns
def GetHistoryDate(hist_file):
    f = h5py.File(hist_file, 'r')
    dt, pm = f.attrs['Date'], f.attrs['Processing_Mode']
    f.close()
    gc.collect()
    f = None
#    return FromIcareHdfDate(dt), pm
    # Pour tri-horaire - à revoir
    try   : return FromIcareHdfDate(dt), pm
    except: return FromIcareHdfUTC(dt).date(), pm


# This function reads and returns algorithm configuration from a text input file
def GetAlgoCfg(cfg, algo_cfg_file):
    if cfg['VERBOSE'] > 1: print_reading_info('algorithm configuration', algo_cfg_file)
    try:
        algo_cfg = Config(None, [""], getall=True).fget({}, algo_cfg_file)
    except:
        print("Error while reading algo configuration file")
        algo_cfg = None
    return algo_cfg


# This function reads and returns product configuration from a text input file
def GetProductCfg(cfg, product_cfg_file):
    if cfg['VERBOSE'] > 1: print_reading_info('product   configuration', product_cfg_file)
    try:
        prd_cfg = Config(None, [""]).fget(PRD_REQ_CFG, product_cfg_file)
    except:
        print("Error while reading product configuration file")
        prd_cfg = None
    return prd_cfg


# This function reads and returns dated infos from a text input file
def ReadInfoTextFile(cfg, txt_file, label):
    if cfg['VERBOSE'] > 1: print_reading_info(label, txt_file)
    try:
        f = open(txt_file, 'r')
        lines = [x for x in f.readlines() if x[0] != '#' and len(x.strip()) > 0]
        info_cfg = [x.split() for x in lines]
        for info in info_cfg: info[0], info[1] = [FromIcareHdfFullUTC(x) for x in info[:2]]
        f.close()
    except:
        print("Error while reading " + label +  "file.")
        info_cfg = None
    return info_cfg


# This function reads and returns solar eclipses data from a text input file
def GetSolarEclipseCfg(cfg, solar_eclipse_file):
    ecl_cfg = ReadInfoTextFile(cfg, solar_eclipse_file, 'solar eclipses data')
    if ecl_cfg is not None:
        for ecl in ecl_cfg: ecl[3] = float(ecl[3])
    return ecl_cfg


# This function reads and returns known corrupted slots from a text input file
def GetCorruptedSlots(cfg, invalid_slot_file):
    return ReadInfoTextFile(cfg, invalid_slot_file, 'corrupted slots data')


# This function reads region name from a band file and returns it
def GetRegionFromBandFile(band):
    ret = None
    try:
        h5f = h5py.File(band, 'r')
        ret = h5f.attrs['Region_Name']
        h5f.close()
    except: print("*** Error while reading from file " + band)
    return ret


# This function reads a subsetting mask form a given binary file
def ReadSubsetMask(cfg): 
    cfg['SUBSET_OUTPUT'] = False
    if cfg['SUBSET'] is not None :
        mskfile = path.join(cfg['ANC_DIR'], cfg['SUBSET'])
        try:
            sdid = SD(mskfile)
            msk = sdid.select('FullDiskMask').get()
            cfg['SUBSET_ROWS'] = sdid.select('GeoGridRow').get()
            cfg['SUBSET_COLS'] = sdid.select('GeoGridCol').get()
            sdid.end()
        except:
            if IsH5File(mskfile):
                sdid = h5py.File(mskfile, 'r')
                sdsid = sdid['FullDiskMask']
                msk = sdsid[:]
                cfg['SUBSET_ROWS'], cfg['SUBSET_COLS'] = sdid['GeoGridRow'][:], sdid['GeoGridCol'][:]
                sdid.close()
                # To prevent buggy behaviour of h5py when closing files
                gc.collect()
            else:
                msk = fromfile(path.join(cfg['ANC_DIR'], cfg['SUBSET']), dtype=uint8)
                msk.shape = (NL, NC)
                cfg['SUBSET_ROWS'], cfg['SUBSET_COLS'] = where(msk)
        if cfg['SUBSUBSET'] is not None :
            cfg['SUBSET_ROWS'] = cfg['SUBSET_ROWS'][:,cfg['SUBSUBSET']]
            cfg['SUBSET_COLS'] = cfg['SUBSET_COLS'][:,cfg['SUBSUBSET']]
            msk *= 2
            msk[cfg['SUBSET_ROWS'], cfg['SUBSET_COLS']] = 1
            msk[where(msk == 2)] = 0
        if cfg['SUBTRANSP']: cfg['SUBSET_GRID'] = (transpose(cfg['SUBSET_ROWS']), transpose(cfg['SUBSET_COLS']))
        else               : cfg['SUBSET_GRID'] = (cfg['SUBSET_ROWS'], cfg['SUBSET_COLS'])
    else :
        cfg['SUBSET_GRID'] = tuple(meshgrid(list(range(NC)), list(range(NL)))[-1::-1])
        cfg['SUBSET_ROWS'], cfg['SUBSET_COLS'] = cfg['SUBSET_GRID']
        msk = ones((NL, NC), dtype=uint8)
        cfg['SUBTRANSP'] = False
    return msk


# This function reads, checks and returns an input variable from HDF4 input file
def ReadHdf4File(cfg, fpath, dset_name, prod, ftyp, todo, check, errwrn=False, DO_CLOSE=True):
    global HdfHandle
    g_attr, d_attr, data = (None, )*3

    try:
        msg = "Error while opening " + fpath
        if HdfHandle is None: HdfHandle = SD(fpath)
        if RD_GATT <= todo:
            g_attr = dict()
            msg = "Error while reading " + fpath + "attributes"
            for key, val in list(HdfHandle.attributes().items()): g_attr[key] = val
            msg = 'Error: not a ' + ftyp + ' file: ' + fpath
            if CHK_PRD <= check and g_attr['Product_Name'] != prod: raise
            msg = 'Error: inconsistent image size for file: ' + fpath
            if CHK_SIZ <= check:
                if cfg['SUBSET'] is None and 'Number_of_Rows_in_Region' in g_attr and \
                        (g_attr['Number_of_Columns_in_Region'] != NC or 
                         g_attr['Number_of_Rows_in_Region'   ] != NL): raise
        if len(RD_SET & todo) > 0:
            dset = HdfHandle.select(dset_name)
            if RD_DATT <= todo:
                d_attr = dict()
                msg = "Error while reading " + fpath + " attributes"
                for key, val in list(dset.attributes().items()): d_attr[key] = val
            if RD_DAT <= todo:
                msg = "Error while reading " + fpath + " data"
                tmpdata = dset.get()
#                if prod in angles_p and prod[1] == 'A': tmpdata = AngleConversion(tmpdata, d_attr)
                if prod in angles_p: tmpdata = AngleConversion(tmpdata, d_attr, azi=(prod[1] == 'A'))
                if cfg['SUBSET'] is not None and tmpdata.shape != (NL, NC):
                    if cfg['SUBSUBSET'] is not None and tmpdata.shape != cfg['SUBSET_ROWS'].shape:
                        tmpdata = tmpdata[:,cfg['SUBSUBSET']]
                    data = ones((NL,NC), dtype=tmpdata.dtype)
                    data[:] = dset.getfillvalue()
                    data[cfg['SUBSET_ROWS'], cfg['SUBSET_COLS']] = tmpdata
                    cfg['SUBSET_OUTPUT'] = True
                else:
                    if cfg['SUBSET_OUTPUT']:
                        msg = 'Inconsistent data size for subset inputs'
                        raise
                    data = tmpdata
            dset.endaccess()
    except:
        if (errwrn): msg = "Warning : " + msg
        print(msg)

    # Close this file?
    if DO_CLOSE:
        if HdfHandle is not None: HdfHandle.end()
        # To prevent buggy behaviour when closing files. No effect with HDF4
        gc.collect()
        HdfHandle = None
    return g_attr, d_attr, data


# This function sets gatt_out with global attributes of gatt_in of an HDF5 file
# If HDF4 file has been converted to HDF5 file, strip specific suffix
def SetHdf5GlobalAttributes(gatt_in, gatt_out, sfx2del):
    for key_in, val in list(gatt_in.items()):
        if len(sfx2del) > 0 and key_in.endswith(sfx2del): key_out = key_in[:-len(sfx2del)]
        else: key_out = key_in
        gatt_out[key_out] = val
            

# This function reads, checks and returns an input variable from HDF5 input file
def ReadHdf5File(cfg, fpath, dset_name, prod, ftyp, todo, check, errwrn=False, DO_CLOSE=True, sfx2del=''):
    global HdfHandle
    g_attr, d_attr, data = (None, )*3

    try:
        msg = "Error while opening " + fpath
        if HdfHandle is None: HdfHandle = h5py.File(fpath, 'r')
        if prod in list(HdfHandle.keys()): handle = HdfHandle[prod]
        else: handle = HdfHandle
        if RD_GATT <= todo:
            g_attr = dict()
            msg = "Error while reading " + fpath + " attributes"
            SetHdf5GlobalAttributes(HdfHandle.attrs, g_attr, sfx2del)
            if handle != HdfHandle: SetHdf5GlobalAttributes(handle.attrs, g_attr, sfx2del)
            msg = 'Error: not a ' + ftyp + ' file: ' + fpath
            if CHK_PRD <= check:
                if 'Product_Name' in g_attr and g_attr['Product_Name'] != prod:
                    if 'ByProduct_Name' in g_attr and g_attr['ByProduct_Name'] != prod: raise
            msg = 'Error: inconsistent image size for file: ' + fpath
            if CHK_SIZ <= check:
                if cfg['SUBSET'] is None and 'Number_of_Rows_in_Region' in g_attr and \
                        (g_attr['Number_of_Columns_in_Region'] != NC or 
                         g_attr['Number_of_Rows_in_Region'   ] != NL): raise
        if len(RD_SET & todo) > 0:
            dset = handle[dset_name]
            if RD_DATT <= todo:
                d_attr = dict()
                msg = "Error while reading " + fpath + " attributes"
                for key, val in list(dset.attrs.items()): d_attr[key] = val
            if RD_DAT <= todo:
                msg = "Error while reading " + fpath + " data"
                tmpdata = dset[:]
#                if prod in angles_p and prod[1] == 'A': tmpdata = AngleConversion(tmpdata, d_attr)
                if prod in angles_p:
                    tmpdata = AngleConversion(tmpdata, d_attr, azi=(prod[1] == 'A'))

                if cfg['SUBSET'] is not None and tmpdata.shape != (NL, NC):
                    if cfg['SUBSUBSET'] is not None and tmpdata.shape != cfg['SUBSET_ROWS'].shape:
                        tmpdata = tmpdata[:,cfg['SUBSUBSET']]
                    data = ones((NL,NC), dtype=tmpdata.dtype)
                    data[:] = dset.attrs['_FillValue']
                    data[cfg['SUBSET_ROWS'], cfg['SUBSET_COLS']] = tmpdata
                    cfg['SUBSET_OUTPUT'] = True
                else:
                    if cfg['SUBSET_OUTPUT']:
                        msg = 'Inconsistent data size for subset inputs'
                        raise
                    data = tmpdata
    except Exception as e:
        if (errwrn): msg = "Warning : " + msg
        print(msg)
        print(e)

    # Close this file?
    if DO_CLOSE:
        if HdfHandle is not None: HdfHandle.close()
        # To prevent buggy behaviour of h5py when closing files
        gc.collect()
        HdfHandle, handle = None, None

    # Return file global attributes, SDS attributes and SDS data
    return g_attr, d_attr, data


# This function converts read angles from (0, 360, uint16) to (-180, 180, int16) with the same origin 0
def AngleConversion(ang_in, datt, azi=False):
    bad = where(ang_in == datt['_FillValue'])
    datt['_FillValue'] = -32768
    if azi:
        ct360 = int32(rint(360/datt['scale_factor'] + datt['add_offset']))
        ct180 = int32(rint(180/datt['scale_factor'] + datt['add_offset']))
        ang_out = int16((int32(ang_in) + ct180) % ct360 - ct180)
    else:
        ang_out = int16(ang_in)
    ang_out[bad] = datt['_FillValue']
    return ang_out


# This function reads an input variable from a HDF4 or HDF5 input file and stores it 
# in a temporary binary file
def ReadAndWriteBin(cfg, path_in, path_out, ipath, sdsnam, prod, ftyp, pfx=''):
    ret = "OK"
    try:
        if   prod is None    : rd, chk = RD_DAT, CHK_NOT
        elif prod in angles_p: rd, chk = RD_ALL, CHK_NOT
        elif prod in ecmwf_p : rd, chk = RD_ALL, CHK_NOT
        else                 : rd, chk = RD_ALL, CHK_ALL
        if cfg['VERBOSE'] > 1:
            if   sdsnam  is None: message = ' '
            elif prod is None: message = sdsnam
            else                : message = prod + '/' + sdsnam + '/' + pfx
            print_reading_info(message, path_in[ipath])
        if IsH5File(path_in[ipath]):
            g_attr, d_attr, data = ReadHdf5File(
                cfg, path_in[ipath], sdsnam, prod, ftyp, rd, chk,
                sfx2del=H4TOH5_CONV_SFX
            )

        else:
            g_attr, d_attr, data = ReadHdf4File(cfg, path_in[ipath], sdsnam, prod, ftyp, rd, chk)
        ext = ('.%s.gb%d' % (pfx, data.dtype.itemsize)).replace('..', '.')
        if rd == RD_ALL:
            if   prod in angles_p:
                ret = (1./d_attr['scale_factor'], -d_attr['add_offset'], d_attr['_FillValue'], g_attr)
            elif prod in ecmwf_p:
                data = EcmwfCalibration(data, d_attr)
                ret = (d_attr['Scaling_Factor'], d_attr['Offset'], d_attr['Missing_Output'], g_attr)
            else:
                ret = (d_attr['Scaling_Factor'], d_attr['Offset'], d_attr['Missing_Output'], g_attr)
        path_out = path.join(cfg['TMPDIR'], path.split(path_out[ipath])[1] + ext)
        data[cfg['SUBSET_GRID']].tofile(path_out)
        print_writing_info(message, path_out)
    except:
        ret = None
        raise
    return ret


def EcmwfCalibration(data, d_attr):
    roffset = FILL_ECMWF + 1
    bad = where(data == d_attr['Missing_Output'])
    d_attr['_FillValue'], d_attr['Missing_Output'], d_attr['Offset'] = FILL_ECMWF, FILL_ECMWF, -roffset
    data = int16(int32(data) + roffset)
    data[bad] = FILL_ECMWF
    return data


# This function reads, scales and returns latitudes from binary float32 input file
def ReadLatitudesRaw(cfg):
    if cfg['VERBOSE'] > 1: print_reading_info('latitudes', cfg['LAT'])
    try:
        data = fromfile(cfg['LAT'], dtype=float32)
        data[where(data <= -100)] = -32768
        data[where(data > -100)] *= 100
        gattr = dict(NL=NL, NC=NC)
        lats = int16(rint(data)).reshape(gattr['NL'], gattr['NC'])
        py_ifc.scale_lat, py_ifc.offset_lat, py_ifc.missing_lat = 100., 0.0, -32768
    except:
        print("Error reading latitudes from file %s" % cfg['LAT'])
        lats, gattr = None, None
    return gattr, lats


def ReadLatLonRaw(cfg):
    gattr = dict(NL=NL, NC=NC)
    lats = GetBinaryLatLon(cfg, 'LAT', 'latitudes' , (NL, NC))
    py_ifc.scale_lat, py_ifc.offset_lat, py_ifc.missing_lat = SCL_LATLON, OFF_LATLON, BAD_LATLON
    lons = GetBinaryLatLon(cfg, 'LON', 'longitudes', (NL, NC))
    py_ifc.scale_lon, py_ifc.offset_lon, py_ifc.missing_lon = SCL_LATLON, OFF_LATLON, BAD_LATLON
    if lats is None or lons is None : lats, lons, gattr = [None]*3
    return gattr, lats, lons


def GetBinaryLatLon(cfg, key, name, shp):
    if cfg['VERBOSE'] > 1: print_reading_info(name, cfg[key])
    try:
        data = fromfile(cfg[key], dtype=float32)
        data[where(data <= -500)] = -32768
        data[where(data > -500)] *= 100
        return int16(rint(data)).reshape(shp)
    except:
        print("Error reading %s from file %s" % (name, cfg[key]))
        return None


def ReadAodClimato(cfg):
    name, key = 'Total AOD, Aerosol model, Boundary weight', cfg['AOD_CLM']
    if cfg['VERBOSE'] > 1: print_reading_info(name, key)
    try:
        h5f = h5py.File(cfg['AOD_CLM'], 'r')
        mth = int(cfg['DATE'][4:6]) - 1
        print('  --> Month:', mth + 1)
        aer_mod1, aer_mod2, aer_mod3, tot_aod, bdw_mod2, bdw_mod3 = h5f['Aerosol model 1'][:],  h5f['Aerosol model 2'][:], \
            h5f['Aerosol model 3'][:],  h5f['Total AOD'][:], h5f['Boundary weight 2'][:], h5f['Boundary weight 3'][:]
        h5f.close()
        py_ifc.missing_aem = -1
        py_ifc.scale_tao, py_ifc.offset_tao, py_ifc.missing_tao = 10000., 0.0,      0
        py_ifc.scale_bdw, py_ifc.offset_bdw, py_ifc.missing_bdw = 10000., 0.0, -32768
        # To prevent buggy behaviour of h5py when closing files
        gc.collect()
        h5f = None
        return aer_mod1[..., mth], aer_mod2[..., mth], aer_mod3[..., mth], tot_aod[..., mth], bdw_mod2[..., mth], bdw_mod3[..., mth]
    except:
        print("Error reading %s from file %s" % (name, key))
        return None


def ReadCoastPixels(cfg, size):
    sdid = SD(cfg['CST_PIX'])
    cst_pix = sdid.select('coasts_lw%d' % size).get()
    sdid.end()
    py_ifc.missing_cst = 0
    return cst_pix
