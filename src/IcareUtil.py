#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################################################################################
###################################################################################################
###################################################################################################
##  Author :
##     ICARE/UDEV - SIX Bruno
##
##  License :
##     This file must be used under the terms of the CCeCILL free software license agreement.
##     This source file is licensed as described in the file COPYING, which you should have
##     received as part of this distribution.
##     The terms are also available at http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
###################################################################################################
###################################################################################################
###################################################################################################


from __future__ import print_function
from __future__ import division
#from builtins import map
#from builtins import str
from datetime import *
import os
from tempfile import mkdtemp
from errno import *

try               : from collections import OrderedDict
except ImportError: from ordereddict import OrderedDict

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from numpy import *

try   : from pyhdf import HDF
except:
#    print "WARNING - IcareUtil.py: cannot load HDF from pyhdf"
    pass

try   : import h5py
except:
#    print "WARNING - IcareUtil.py: cannot load h5py"
    pass


PROD_CENTER = 'AERIS / ICARE Data and Services Center'
CONTACT = 'contact@icare.univ-lille1.fr'
ICA_DFT_VER_KEY = 'DEFAULT'

NetCDF_, HDF4_, HDF5_ = (0, 1, 2)
HDF_ = HDF4_
OUTFILE_EXTENSION = ['nc', 'hdf', 'h5']

HDF5_SIGNATURE = b'\x89HDF\r\n\x1a\n'

EXIT_CODE_WRN = 127
EXIT_CODE_ERR =   1
EXIT_CODE_OK_ =   0

# This function reads and returns dated infos from a text input file
def ReadDatedInfoTextFile(cfg, txt_file, label):
    info_cfg = None
    try:
        f = open(txt_file, 'r')
        try:
            lines = [x for x in f.readlines() if x[0] != '#' and len(x.strip()) > 0]
            info_cfg = [x.split() for x in lines]
            for info in info_cfg: info[0], info[1] = [FromIcareHdfFullUTC(x) for x in info[:2]]
        except:
            print("Error while reading " + label +  "file.")
        f.close()
    except:
        pass
    return info_cfg


def ReadConfigFile(fpath, ordered=False, sep='=', err=True, add=False):
    rep, nam = os.path.split(fpath)
    if ordered: d = OrderedDict()
    else      : d = dict()
    try:
        fcfg = open(fpath, 'r')
        l = ''
        for ll in fcfg:
            ll0 = ll.strip()
            if len(ll0) == 0: continue
            if ll0[0] == '#': continue
            if ll0[0] == '@':
                incl_pth = os.path.join(rep, ll0[1:])
                d.update(ReadConfigFile(incl_pth, ordered=ordered, sep=sep, err=err))
                continue
            l += ll0
            if l[-1] == '\\' :
                l = l[:-1]
                continue
            name, value = l.split(sep)
            key = name.strip()
#            d[key] = value.strip().strip('\n')
            if add:
                if key not in d: d[key] = []
                d[key].append(value.strip().strip('\n'))
            else: d[key] = value.strip().strip('\n')
            l = ''
        fcfg.close()
    except:
        if err: print("Error reading configuration file " + fpath)
        d = None
    return d


def WriteConfigFile(path, d):
    status = 0
    try:
        fcfg = open(path, 'w')
        for k in list(d.keys()):
            fcfg.write("%s = %s\n" %(k, str(d[k])))
        fcfg.close()
    except:
        print('Error writing file ' + path)
        status = -1
    return status


def GetRedhatRelease():
    f = open('/etc/redhat-release')
    ll = f.readline().strip('\n').split()
    f.close()
    return float(ll[ll.index('release')+1])


def GetHDF4Version(): return "%d.%dr%d" % HDF.getlibversion()[:3]


def GetHDF5Version(): return h5py.version.hdf5_version


def IsH5File(fpath):
    # changed 'r' to 'rb' during porting to py3
    f = open(fpath, 'rb')
    signature = f.read(len(HDF5_SIGNATURE))
    f.close()
    return signature == HDF5_SIGNATURE


def NumpyToHdf4DataType(numpy_obj):
    try   : return eval('HDF.HC.' + numpy_obj.dtype.name.upper())
    except:
        if   numpy_obj.dtype.kind == 'b'            : return SDC.UCHAR8
        elif numpy_obj.dtype.kind in ('U', 'S', 'a'): return SDC.CHAR8
    return None


def GetIcareVersionFromStr(ver_str):
    return tuple(map(int, ver_str.strip("Vv ").replace("-",".").split('.')))


def GetIcareFileVersion(version):
    return "%d-%02d" % version[:2]


def GetIcareMetadata():
    return dict(
        Production_Center = PROD_CENTER,
        Production_Date   = IcareHdfUTC(datetime.utcnow()),
        )


def GetIcareCgtdMetadata(pversion, cversion, icare_id, contact=CONTACT):
    d = dict(
        Icare_ID = icare_id,
        Product_Version = "%d.%d" % GetIcareVersionFromStr(pversion)[:2],
        Software_Version = "%d.%d.%d" % cversion,
        )
    if contact is not None: d['Contact'] = contact
    d.update(GetIcareMetadata())
    return d

def IcareAttributeComment(begin, header, text, end):
    return begin + header + ': ' + text + end

def IcareHdfCalibrationEquation():
    return 'physical_value = scale_factor*(SDS_count - add_offset)'

def IcareScanOriginAttrComment(begin = '', header = 'Scan_Origin', end = '.'):
    return IcareAttributeComment(begin, header, 'NW, NE, SW or SE', end)

def IcareAltitudeAttrComment(begin = '', header = 'Altitude', end = '.'):
    return IcareAttributeComment(begin, header, 'unit = km from Earth center', end)

def IcarePixelSizeAttrComment(begin = '', header = 'Nadir_Pixel_Size', end = '.'):
    return IcareAttributeComment(begin, header, 'unit = km', end)

def IcareMaxLatLonAttrComment(begin = '', header = 'Northernmost_Latitude, Southernmost_Latitude, ' +
                              'Westernmost_Longitude, Easternmost_Longitude', end = '.'):
    return IcareAttributeComment(begin, header, 'minimum and maximum latitudes and longitudes in degrees ' +
                                 '(latitudes positive North of the equator, ' +
                                 'longitudes positive East of Greenwich meridian)', end)

def IcareGeoZoneAttrComment(begin = '', header = 'Zone', end = '.'):
    return IcareAttributeComment(begin, header, 'scan area = G (Globe), N (North) or S (South)', end)

def IcareInputFileAttrComment(begin = '', header = 'Input_Files, Ancillary_Files', end = '.'):
    return IcareAttributeComment(begin, header, 'when applicable, filenames are preceded by the parent "date" ' +
                                 'directory (YYYY_MM_DD) used in the %s archive tree structure' % PROD_CENTER, end)

def IcareDateTimeAttrComment(microseconds = False, begin = '',
                             header = 'Beginning_Acquisition_Date, End_Acquisition_Date, Production_Date', end = '.'):
    ms_str = ''
    if microseconds: ms_str = '[.9*]'
    text = 'dates and times comply with ISO-8601 standard:\n' + \
           '  - date: YYYY-MM-DD\n' + \
           '  - time: hh:mm:ss%sZ\n' % ms_str + \
           '  - date & time: YYYY-MM-DDThh:mm:ss%sZ' % ms_str
    return IcareAttributeComment(begin, header, text, end)

def IcareGeoCenteredGridComment(prj_lon, nbl, nbc, begin = '', header = 'Centered_Grid', end = '.'):
    inc = nbl - 2*(nbl//2)
    str = [ 'at the center of %dth row and column' % (nbl//2 + 1),
            'halfway between %dth and %dth row and column' % (nbl//2 + inc , nbl//2 + inc + 1) ]
    if inc : str.reverse()
    text = ('indicates whether the projection center (latitude: 0.0, longitude: %.1f) is located at the center of the SDS grid; ' + \
        '1 = centered, i.e. projection center ' + str[1] + '; 0 = not centered, i.e. projection center ' + str[0]) % prj_lon
    return IcareAttributeComment(begin, header, text, end)

def IcareHdfDateIso(dt):
    if   isinstance(dt, date): return dt.strftime("%Y-%m-%d")
    elif isinstance(dt, str ): return datetime.strptime(dt, "%Y%m%d").strftime("%Y-%m-%d")
    return None

def IcareHdfUTC(dt):
    if   isinstance(dt, datetime): return dt.strftime("%Y-%m-%dT%H:%M:%SZ")
    elif isinstance(dt,     time): return dt.strftime("%H:%M:%SZ")
    else                         : return dt.strftime("%Y-%m-%d")

def IcareHdfIsoDuration(td, detailed=False, prec=-1):
    ret = 'P'
    if td.days != 0: ret += "%dD" % td.days
#    seconds = td.seconds + int(rint(td.microseconds*1.e-6))
    seconds = td.seconds + td.microseconds*1.e-6
    if seconds > 0:
        ret += 'T'
        if detailed:
            dd = datetime(1900,1,1,0,0,0) + td
            if dd.hour > 0: ret += '%dH' % dd.hour
            if dd.minute > 0: ret += '%dM' % dd.minute
#            seconds = dd.second + int(rint(dd.microsecond*1.e-6))
            seconds = dd.second + dd.microsecond*1.e-6
        if   prec == 0: ret += '%dS' % int(rint(seconds))
        elif prec < 0 : ret += '%.3fS' % seconds
        else          : ret += ('%%.%dfS' % prec) % seconds
    return ret

def IcareHdfFullUTC(dt):
    return dt.strftime("%Y-%m-%dT%H:%M:%S") + (".%03dZ" % int(rint((dt.microsecond/1000.))))

def IcareHdfFileNameUTC(dt, utc=False):
    if utc: utc_fmt = "%Y-%m-%dT%H-%M-%SZ"
    else  : utc_fmt = "%Y-%m-%dT%H-%M-%S"
    if   isinstance(dt, datetime): return dt.strftime(utc_fmt)
    elif isinstance(dt, time    ): return dt.strftime("%H-%M-%S")
    elif isinstance(dt, date    ): return dt.strftime("%Y-%m-%d")
    elif isinstance(dt, str     ):
        if dt[:2] == '%s' : return dt
        elif dt[0] == '@' : return dt[1:]
        elif len(dt) ==  8: return datetime.strptime(dt, "%Y%m%d").strftime("%Y-%m-%d")
        elif len(dt) == 12: return datetime.strptime(dt, "%Y%m%d%H%M").strftime(utc_fmt)
        elif len(dt) == 14: return datetime.strptime(dt, "%Y%m%d%H%M%S").strftime(utc_fmt)
    return None


def IcareHdfFileName(mission, prod, subprod, dt, spec, version, extension, dt2=None, duration=None,
                     utc=False, subprod_delim='-'):
    fname = ""
    if mission   is not None: fname += "_"  + mission
    if prod      is not None: fname += "_"  + prod
    if subprod   is not None: fname += subprod_delim + subprod
    if dt        is not None: fname += "_"  + IcareHdfFileNameUTC(dt, utc=utc)
    if dt2       is not None: fname += "@"  + IcareHdfFileNameUTC(dt2)
    if duration  is not None: fname += "-"  + duration
    if spec      is not None: fname += "_"  + spec
    if isinstance(version, tuple): fname += "_V%d-%02d" % version[:2]
    elif version is not None     : fname += "_V" + version.replace('.', '-').strip('vV')
    if extension is not None: fname += "."  + extension
    return fname[1:]


def FromIcareHdfFileNameUTC(dtstr): return datetime.strptime(dtstr[:19], "%Y-%m-%dT%H-%M-%S")

def FromIcareHdfUTC(dtstr): return datetime.strptime(dtstr[:19], "%Y-%m-%dT%H:%M:%S")

def FromIcareHdfFullUTC(dtstr):
    try   : return datetime.strptime(dtstr, "%Y-%m-%dT%H:%M:%SZ%f")  # Format NetCDF SATMOS
    except:
        try   : return datetime.strptime(dtstr, "%Y-%m-%dT%H:%M:%S.%fZ") # Format ISO Icare
        except: return datetime.strptime(dtstr, "%Y-%m-%dT%H:%M:%S.%f") # Format ISO

def FromIcareHdfDate(dtstr): return datetime.strptime(dtstr, "%Y-%m-%d").date()


def IcarePath(path):
    abspth = os.path.abspath(path)
    idir, nam = os.path.split(abspth)
    ddir = os.path.split(idir)[1]
    try:
        dt = datetime.strptime(ddir, "%Y_%m_%d")
        return os.path.join(ddir,nam)
    except:
        return nam


def IcarePath2(path):
    idir, nam = os.path.split(os.path.abspath(path))
    ddir = os.path.split(idir)[1]
    try:
        dt = datetime.strptime(ddir, "%Y_%m_%d")
        return os.path.join(dt.strftime("%Y"), os.path.join(ddir,nam))
    except:
        return nam


def MakeNewDir(path):
    try:
        os.makedirs(path)
    except OSError as erreur:
        if erreur.errno != EEXIST:
            print("Error making %s dir." % path)
            return False
    return True


def MakeTmpDir(path, prefix):
    try:
        tmpdir = mkdtemp(prefix=prefix, dir=path)
        os.makedirs(tmpdir)
    except OSError as erreur:
        if erreur.errno != EEXIST:
            print("Error making temporary directory.")
            tmpdir = None
    return tmpdir


def JulToGregYearDay(year, month, day):
    return date(year, month, day).toordinal() + 1 - date(year, 1, 1).toordinal()

def GregToJulYearDay(year, jday):
    return date(year, 1, 1) + timedelta(days=jday-1)
