#!/usr/bin/env python3

# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import division
import os
from os import _exit, path

import gc

from numpy import *

try:
    from pyhdf.SD import *
except:
    print("ERROR loading pyhdf")
    pass

from datetime import *

import warnings
#warnings.filterwarnings("error", category=RuntimeWarning)

from ConfigFile import *

WithGeoSchedule = True
try:
    from GeoSchedule import *
except:
    WithGeoSchedule = False

C_RIDL =  6378.2064  ## IDL-like       earth radius (km)
#C_RPOL =  6356.7523  ## Polar          earth radius (km)  EUMETSAT
#C_REQU =  6378.1370  ## Equatorial     earth radius (km)  EUMETSAT
#C_HSAT = 42164.57    ## Geostationnary orbit radius (km)

C_RPOL =  6356.5838  ## Polar          earth radius (km)
C_REQU =  6378.1690  ## Equatorial     earth radius (km)
#C_HSAT = 42164.12    ## Geostationnary orbit radius (km)
C_HSAT = 42164       ## Geostationnary orbit radius (km)

C_MVD2 = C_HSAT*C_HSAT - C_REQU*C_REQU ## Squared maximum visible distance
C_MVD  = sqrt(C_MVD2)                  ## Maximum visible distance
C_RR   = C_RPOL/C_REQU                 ## Earth radii ratio
C_RR2  = C_RR*C_RR                     ## Squared earth radii ratio
C_1DRR2= 1./C_RR2                      ## C_RR2 inverse 

GEO_ECL_DT0 = datetime(1999,12,31) ## Base date for sun elements approximation

H_PLANCK = 6.6260688e-34   ## Planck constant   (J.s)
C_LIGHT  = 2.997924580e8   ## Light celerity    (m/s)
K_BOLTZ  = 1.3806503e-23   ## Boltzman constant (J/K)

PLANCK_C1 = 2.*H_PLANCK*C_LIGHT*C_LIGHT*1.e+11
PLANCK_C2 = 100.*H_PLANCK*C_LIGHT/K_BOLTZ


def GetConfig(cfg_in, path):
    cfg = Config(None, [""]).fget(cfg_in, path)
    return cfg


def SunApparentOrbitElements(utc):

##  Reference datetime for satellite longitude
    delta = utc - GEO_ECL_DT0
    d     = delta.days + delta.seconds/86400.

##  Sun's apparent Earth-orbit elements 
    peri_lon = (282.9404   + d*4.70935e-5  ) % 360  # Perihelion longitude
    ecc      =    0.016709 - d*1.151e-9             # Eccentricity
    mean_ano = (356.0470   + d*0.9856002585) % 360  # Mean anomaly
    ecl_incl = ( 23.4393   - d*3.563e-7    ) % 360  # Ecliptic inclination

    r_mean_ano = radians(mean_ano)
    r_ecl_incl = radians(ecl_incl)
    r_peri_lon = radians(peri_lon)

##  Eccentric anomaly approximation
    E = r_mean_ano + ecc*sin(r_mean_ano)*(1 + ecc*cos(r_mean_ano))

##  Sun coordinates in ecliptic plane (X axis towards perihelion)
    x, y = (cos(E) - ecc, sin(E)*sqrt(1 - ecc**2))

    return dict(Perihelion_Longitude=peri_lon, Eccentricity=ecc, Eccentric_Anomaly=E,
                Mean_Anomaly=mean_ano, Ecliptic_Inclination=ecl_incl, Ecliptic_Coordinates=(x,y))


def SunEarthDistance(utc):
    x, y =  SunApparentOrbitElements(utc)['Ecliptic_Coordinates']
    return double(sqrt(x**2 + y**2))


## This function computes dates et times of beginning and end of a geostationnary
## satellite eventual eclipse. The method is due to Paul Schlyter, Stockholm, Sweden
## and is descibed at http://www.stjarnhimlen.se/comp/tutorial.htm
##
def GeoEclipseDates(dt0, satlon):

##  Reference datetime for satellite longitude
    dt    = datetime.combine(dt0, time(0)) - timedelta(hours=satlon/15.)

##  Sun's apparent Earth-orbit elements 
    sun_elts = SunApparentOrbitElements(dt)

    r_peri_lon, r_ecl_incl = radians(sun_elts['Perihelion_Longitude']), radians(sun_elts['Ecliptic_Inclination'])
    x, y = sun_elts['Ecliptic_Coordinates']

##  Sun distance and apparent longitude
    r, sun_lon = (sqrt(x**2 + y**2), arctan2(y,x) + r_peri_lon)

##  Sun ecliptic coordinates
    xx, yy, zz = (r*cos(sun_lon), r*sin(sun_lon), 0.)

##  Sun equatorial coordinates
    xeq, yeq, zeq = (xx, yy*cos(r_ecl_incl) - zz*sin(r_ecl_incl), yy*sin(r_ecl_incl) - zz*cos(r_ecl_incl))

##  Sun declinaison
    Decl = arcsin(zeq/sqrt(xeq**2 + yeq**2 + zeq**2))
    
##  Eclipse duration
    try:
        ts = degrees(arccos(sqrt(1 - (C_REQU/C_HSAT)**2)/cos(Decl)))/7.5
        return (dt - timedelta(hours=ts/2.), dt + timedelta(hours=ts/2))
    except (ValueError, RuntimeWarning):
        # No eclipse...
        return None


## This function computes dates et times of beginning and end of a geostationnary
## satellite eventual eclipse. The method is due to Paul Schlyter, Stockholm, Sweden
## and is descibed at http://www.stjarnhimlen.se/comp/tutorial.htm
##
def GeoEclipseDatesOld(dt0, satlon):

##  Reference datetime for satellite longitude
    dt    = datetime.combine(dt0, time(0)) - timedelta(hours=satlon/15.)
    delta = dt - GEO_ECL_DT0
    d     = delta.days + delta.seconds/86400.

##  Sun's apparent Earth-orbit elements 
    peri_lon = (282.9404   + d*4.70935e-5  ) % 360  # Perihelion longitude
    ecc      =    0.016709 - d*1.151e-9             # Eccentricity
    mean_ano = (356.0470   + d*0.9856002585) % 360  # Mean anomaly
    ecl_incl = ( 23.4393   - d*3.563e-7    ) % 360  # Ecliptic inclination

    r_mean_ano = radians(mean_ano)
    r_ecl_incl = radians(ecl_incl)
    r_peri_lon = radians(peri_lon)
    
##  Eccentric anomaly approximation
    E = r_mean_ano + ecc*sin(r_mean_ano)*(1 + ecc*cos(r_mean_ano))

##  Sun coordinates in ecliptic plane (X axis towards perihelion)
    x, y = (cos(E) - ecc, sin(E)*sqrt(1 - ecc**2))

##  Sun distance and apparent longitude
    r, sun_lon = (sqrt(x**2 + y**2), arctan2(y,x) + r_peri_lon)

##  Sun ecliptic coordinates
    xx, yy, zz = (r*cos(sun_lon), r*sin(sun_lon), 0.)

##  Sun equatorial coordinates
    xeq, yeq, zeq = (xx, yy*cos(r_ecl_incl) - zz*sin(r_ecl_incl), yy*sin(r_ecl_incl) - zz*cos(r_ecl_incl))

##  Sun declinaison
    Decl = arcsin(zeq/sqrt(xeq**2 + yeq**2 + zeq**2))
    
##  Eclipse duration
    try:
        ts = degrees(arccos(sqrt(1 - (C_REQU/C_HSAT)**2)/cos(Decl)))/7.5
        return (dt - timedelta(hours=ts/2), dt + timedelta(hours=ts/2.))
    except ValueError :
        # No eclipse...
        return None


def LinCol2LatLon(lin, col, cfg, err=-999.):
    return LinCol2LatLon_(lin, col, cfg['LON'], cfg['CENTERED'], cfg['CFAC'], cfg['COFF'], \
                          cfg['LFAC'], cfg['LOFF'], cfg['SWEEP_XAXIS'], err=err)


def LinCol2LatLon_(lin, col, ssplon, centered, cfac, coff, lfac, loff, sweepx=0, err=-999.):
    x = radians(float64((array(col) + 0.5*centered - coff)*2**16/cfac))
    y = radians(float64((array(lin) + 0.5*centered - loff)*2**16/lfac))

    cosx, cosy, sinx, siny = (cos(x), cos(y), sin(x), sin(y))
    del x, y

    aux  = C_HSAT*cosx*cosy
    aux2 = cosy*cosy + C_1DRR2*siny*siny
    
    sd = sqrt(aux*aux - aux2*C_MVD2)
    sn = (aux - sd)/aux2
    del aux, aux2
    gc.collect()

    s0 = C_HSAT - sn*cosx*cosy
    if sweepx: sg, s1, s2, s3 = -1, C_HSAT - s0, -sn*sinx     , -sn*siny*cosx
    else     : sg, s1, s2, s3 =  1,          s0,  sn*sinx*cosy, -sn*siny
    del sn, cosx, cosy
    gc.collect()

#    sxy = sqrt(s1*s1 + s2*s2)
    sxy = sqrt(s0*s0 + s2*s2)

    lat = degrees(arctan(C_1DRR2*s3/sxy))
    del s3,sxy
#    lon = AbsoluteLongitude(degrees(arctan(s2/s1)), ssplon)
    lon = AbsoluteLongitude(degrees(arctan(sg*s2/s0)), ssplon)
    del s1,s2

#    if isinstance(lin, ndarray):
    if isinstance(sd, ndarray):
        bad = where(logical_not(isfinite(sd)))
        lat[bad] = err
        lon[bad] = err
    elif not isfinite(sd):
        lat, lon = ( err, err)
    return (lat, lon)


def LatLon2LinCol(lat, lon, cfg, err=-1):
    return LatLon2LinCol_(lat, lon, cfg['LON'], cfg['CENTERED'], cfg['CFAC'], \
                          cfg['COFF'], cfg['LFAC'], cfg['LOFF'], cfg['SWEEP_XAXIS'], err=err)

def LatLon2LinColBis(lat, lon, cfg, err=-1):
    return LatLon2LinColBis_(lat, lon, cfg['LON'], cfg['CENTERED'], cfg['CFAC'], \
                          cfg['COFF'], cfg['LFAC'], cfg['LOFF'], cfg['SWEEP_XAXIS'], err=err)

def LatLon2LinColBis_(lat, lon, ssplon, centered, cfac, coff, lfac, loff, sweepx=0, err=-1):
    warnings.filterwarnings("error", category=RuntimeWarning)

    C_2_16 = 2**(-16)
    
    rlat  = radians(array(lat))
#    rlon  = radians(array(lon))
#    rslon = radians(ssplon)
    rrlon  = radians(array(lon)) - radians(ssplon)
    
    c_lat = arctan(C_RR2*tan(rlat))
    cos_clat = cos(c_lat)
    sin_clat = sin(c_lat)
    del c_lat

    rl = C_RPOL/sqrt(1 - (1 - C_RR2)*cos_clat*cos_clat)
#    r1 = C_HSAT - rl*cos_clat*cos(rlon - rslon)
#    r2 = -rl*cos_clat*sin(rlon - rslon)
    r1 = C_HSAT - rl*cos_clat*cos(rrlon)
    r2 = -rl*cos_clat*sin(rrlon)
    r3 = rl*sin_clat
    del cos_clat, sin_clat

    rn = sqrt(r1*r1 + r2*r2 + r3*r3)
    ad2 = r1*r1 + r2*r2 + r3*r3/C_RR2
    bd = C_HSAT*r1;
    delta2 = bd*bd - ad2*C_MVD2;
    halfsom = bd*rn/ad2

    if not isinstance(lat, ndarray):
        if (delta2 < 0 or  rn > halfsom): return (-1, -1)
        if sweepx: x, y  = (degrees(arcsin(-r2/rn)), degrees(arctan(-r3/r1)))
        else     : x, y  = (degrees(arctan(-r2/r1)), degrees(arcsin(-r3/rn)))
    else:
        ok = where(logical_and(delta2 >= 0, rn <= halfsom))
        del delta2, halfsom
        if sweepx: x, y  = (degrees(arcsin(-r2[ok]/rn[ok])), degrees(arctan(-r3[ok]/r1[ok])))
        else     : x, y  = (degrees(arctan(-r2[ok]/r1[ok])), degrees(arcsin(-r3[ok]/rn[ok])))
        del r1, r2, r3, rn

    decal = 0
    if centered == 1: decal = -0.5
    xc, yl = (coff + x*C_2_16*cfac + decal, loff + y*C_2_16*lfac + decal)
    del x, y

    if isinstance(lat, ndarray):
        #lin, col = err*ones(rlat.shape, dtype=int), err*ones(rrlon.shape, dtype=int)
        lin, col = err*ones(rlat.shape, dtype=int16), err*ones(rrlon.shape, dtype=int16)
        lin[ok], col[ok] = asarray(rint(yl), dtype=int16), asarray(rint(xc), dtype=int16)
        return (lin, col)
    else: return (int(rint(yl)), int(rint(xc)))


def LatLon2LinCol_(lat, lon, ssplon, centered, cfac, coff, lfac, loff, sweepx=0, err=-1):
    C_2_16 = 2**(-16)
    
    rlat  = radians(array(lat))
    rlon  = radians(array(lon))
    rslon = radians(ssplon)
    
    c_lat = arctan(C_RR2*tan(rlat))
    cos_clat = cos(c_lat)
    sin_clat = sin(c_lat)
    rl = C_RPOL/sqrt(1 - (1 - C_RR2)*cos_clat*cos_clat)
    r1 = C_HSAT - rl*cos_clat*cos(rlon - rslon)
    r2 = -rl*cos_clat*sin(rlon - rslon)
    r3 = rl*sin_clat
    rn = sqrt(r1*r1 + r2*r2 + r3*r3)

    ad2 = r1*r1 + r2*r2 + r3*r3/C_RR2
    bd = C_HSAT*r1;
    delta2 = bd*bd - ad2*C_MVD2;
    halfsom = bd*rn/ad2

    lin, col = err*ones(rlat.shape, dtype=int), err*ones(rlon.shape, dtype=int)

    if not isinstance(lat, ndarray):
        if (delta2 < 0 or  rn > halfsom): return (-1, -1)
        if sweepx: x, y  = (degrees(arcsin(-r2/rn)), degrees(arctan(-r3/r1)))
        else     : x, y  = (degrees(arctan(-r2/r1)), degrees(arcsin(-r3/rn)))
    else:
        ok = where(logical_and(delta2 >= 0, rn <= halfsom))
        if sweepx: x, y  = (degrees(arcsin(-r2[ok]/rn[ok])), degrees(arctan(-r3[ok]/r1[ok])))
        else     : x, y  = (degrees(arctan(-r2[ok]/r1[ok])), degrees(arcsin(-r3[ok]/rn[ok])))

    decal = 0
    if centered == 1: decal = -0.5
    xc, yl = (coff + x*C_2_16*cfac + decal, loff + y*C_2_16*lfac + decal)

    if isinstance(lat, ndarray):
        lin[ok], col[ok] = asarray(rint(yl), dtype=int), asarray(rint(xc), dtype=int)
        return (lin, col)
    else: return (int(rint(yl)), int(rint(xc)))


def GeoSetLinColBounds(geo, zone):
    if   geo[zone] == -1 : return [1, 1, 0, 0]
    else: return geo[zone]


def EarthLatRadius(lat):
    return C_RPOL/ sqrt( 1 - cos(radians(float_(lat))**2 * (1 - (C_RPOL/C_REQU)**2)))


def GeodeticAngle(lon0, lat0, lon1, lat1):
    try:
        # Pour éviter les problèmes numériques éventuels
##         rlat0, rlon0, rlat1, rlon1 = (deg2rad(float_(lat0)), deg2rad(float_(lon0)),
##                                       deg2rad(float_(lat1)), deg2rad(float_(lon1)))
        rlat0, rlon0, rlat1, rlon1 = (radians(float_(lat0)), radians(float_(lon0)),
                                      radians(float_(lat1)), radians(float_(lon1)))
        cosc = minimum(1.0, maximum(-1.0, sin(rlat0)*sin(rlat1) + cos(rlat0)*cos(rlat1)*cos(rlon0 - rlon1)))
        return arccos(cosc)
    except:
        raise


def GeoDist(lon0, lat0, lon1, lat1):
    try: return GeodeticAngle(lon0, lat0, lon1, lat1)*(EarthLatRadius(lat0) + EarthLatRadius(lat1))/2.
    except: return -1.


# Pour avoir la même valeur que MAP_2POINTS d'IDL
def Map2Points(lon0, lat0, lon1, lat1):
    try: return GeodeticAngle(lon0, lat0, lon1, lat1)*C_RIDL
    except: return -1.


def GeoBt2Rad(T, a, b, v):
    return PLANCK_C1*v**3/(exp(PLANCK_C2*v/(b*T + a)) - 1.)


def GetGeostatConfig(geo_cfg, sched_sw=True, no_sat=False):
    geo_req_cfg = {
        'SAT'         : (str  , 's:', None),
        'LON'         : (float, 'l:', None),
        'SSP'         : (float, 'l:',-999.),
        'ALT'         : (float, 'a:', None),
        'COLS'        : (int  , 'w:', None),
        'LINES'       : (int  , 'h:', None),
        'CENTERED'    : (int  , 'g:', None),
        'SWEEP_XAXIS' : (int  , 'x:',    0),
        'CYCLE'       : (int  , 'f:', None),
        'SPIN_RATE'   : (int  , 'r:', None),
        'DIR'         : (str  , 'd:', None),
        'SCAN_DELTA_T': (float, 'D:', None),
        'SCAN_RATE'   : (str  , 'R:', None),
        'COFF'        : (float, 'c:', None),
        'LOFF'        : (float, 'l:', None),
        'CFAC'        : (float, 'C:', None),
        'LFAC'        : (float, 'L:', None),
        'GLOBE'       : (int  , 'G:',    1),
        'NORTH'       : (int  , 'N:',   -1),
        'SOUTH'       : (int  , 'S:',   -1),
        'NORTH_MIN'   : (int  ,   '',   -1),
        'SOUTH_MIN'   : (int  ,   '',   -1),
        'NORTH_MAX'   : (int  ,   '',   -1),
        'SOUTH_MAX'   : (int  ,   '',   -1),
        'SENSOR'      : (str  ,   '', None),
        'INSTR'       : (str  ,   '', None),
        'HRVIS'       : (str  ,   '',   ""),
        'VIS06'       : (str  ,   '',   ""),
        'VIS08'       : (str  ,   '',   ""),
        'IR016'       : (str  ,   '',   ""),
        'IR039'       : (str  ,   '',   ""),
        'WV062'       : (str  ,   '',   ""),
        'WV073'       : (str  ,   '',   ""),
        'IR087'       : (str  ,   '',   ""),
        'IR097'       : (str  ,   '',   ""),
        'IR108'       : (str  ,   '',   ""),
        'IR120'       : (str  ,   '',   ""),
        'IR134'       : (str  ,   '',   ""),
        'GEOCLD'      : (str  ,   '',   ""),
        'RES'         : (float,   '', None),
        'LRES'        : (float,   '',   0.),
        'LCENTERED'   : (int  ,   '',  -1 ),
        'HRES'        : (float,   '',   0.),
        'HCENTERED'   : (int  ,   '',  -1 ),
        'VRES'        : (float,   '',   0.),
        'VCENTERED'   : (int  ,   '',  -1 ),
        'LCOLS'       : (int  ,   '',   0 ),
        'LLINES'      : (int  ,   '',   0 ),
        'LCOFF'       : (float,   '',   0.),
        'LLOFF'       : (float,   '',   0.),
        'LCFAC'       : (float,   '',   0.),
        'LLFAC'       : (float,   '',   0.),
        'HCOLS'       : (int  ,   '',   0 ),
        'HLINES'      : (int  ,   '',   0 ),
        'HCOFF'       : (float,   '',   0.),
        'HLOFF'       : (float,   '',   0.),
        'HCFAC'       : (float,   '',   0.),
        'HLFAC'       : (float,   '',   0.),
        'VCOLS'       : (int  ,   '',   0 ),
        'VLINES'      : (int  ,   '',   0 ),
        'VCOFF'       : (float,   '',   0.),
        'VLOFF'       : (float,   '',   0.),
        'VCFAC'       : (float,   '',   0.),
        'VLFAC'       : (float,   '',   0.),
        }

    if no_sat: geo_req_cfg['SAT'] = geo_req_cfg['SAT'][:2] + ('',)
    cfg = Config(None, [""]).fget(geo_req_cfg, geo_cfg)
    
    if sched_sw and WithGeoSchedule: cfg['SCHED'] = GeoSchedule(cfg['SAT'])

    return cfg


def RelativeLongitude(lon, lon0): return (lon - lon0 + 180) % 360 - 180


def AbsoluteLongitude(lon, lon0): return RelativeLongitude(lon, -lon0)


def CheckInterval(interval, verb=True, sens=True):
    if   not isinstance(interval, list):
        if verb: print('This is not a sequence:' , interval)
    elif len(interval) != 2            :
        if verb: print('This is not an interval:', interval)
    elif sens and interval[0] > interval[1] :
        if verb: print('Invalid interval:'       , interval)
    else                               :
        return True
    return False


def CheckVal(value, vmin, vmax, verb=True):
    if value >= vmin and value <= vmax: return True
    if verb: print(value, 'is not in [', vmin, ',', vmax, ']')
    return False


def CheckValInRange(value, extr_vals, verb=True):
    if CheckInterval(extr_vals, verb): return CheckVal(value, extr_vals[0], extr_vals[1], verb)
    return False


def RangeIntersect(range1, range2, verb=True, check=False, closed=False):
    if CheckInterval(range1, verb) and CheckInterval(range2, verb):
        rngn, rngx = (max(range1[0], range2[0]), min(range1[1], range2[1]))
        if closed:flag = rngn <= rngx
        else: flag = rngn < rngx
        if check: return flag
        elif flag: return [rngn, rngx]
    if check: return False
    else: return None


def CheckLatitude (value, verb=True): return CheckVal(value,  -90,  90, verb)


def CheckLongitude(value, verb=True): return CheckVal(value, -180, 180, verb)


def CheckLatInterval(lats, verb=True):
    if CheckInterval(lats, verb):
        return (CheckLatitude(lats[0], verb) and CheckLatitude(lats[1], verb))
    return False


def CheckLonInterval(lons, verb=True, sens=False):
    if CheckInterval(lons, verb, sens):
        return (CheckLongitude(lons[0], verb) and CheckLongitude(lons[1], verb))
    return False


def CheckLatInLatRange(lat, extr_lats, verb=True):
    if CheckLatInterval(extr_lats, verb): return CheckValInRange(lat, extr_lats, False)
    if verb: print("Bad latitudes range")
    return False


def CheckLonInLonRange(lon, extr_lons, orig, verb=True):
    if CheckLonInterval(extr_lons, verb):
        rlon = RelativeLongitude(lon, orig)
        rlons = list(RelativeLongitude(array(extr_lons), orig))
        print(rlon, rlons)
        return CheckValInRange(rlon, rlons, False)
    return False


def LatRangeIntersect(lats1, lats2, closed=False, verb=True, check=False):
    if CheckLatInterval(lats1, verb) and CheckLatInterval(lats2, verb):
        return RangeIntersect(lats1, lats2, verb, check)
    if check: return False
    else: return None


def LonRangeIntersect(lons1, lons2, orig, closed=False, verb=True, check=False):
    print(lons1, lons2)
    if CheckLonInterval(lons1, verb) and CheckLonInterval(lons2, verb):
        rlons1, rlons2 = ( list(lons1), list(lons2) )
        while rlons1[0] > rlons1[1]: rlons1[0] -=360
        while rlons2[0] > rlons2[1]: rlons2[0] -=360
##         while rlons1[0] < 0 : rlons1 = map(lambda i:i+360,rlons1)
##         while rlons2[0] < 0 : rlons2 = map(lambda i:i+360,rlons2)
        while rlons1[0] < 0 or rlons2[0] < 0 :
            rlons1 = [i+360 for i in rlons1]
            rlons2 = [i+360 for i in rlons2]
        if rlons1[1] - rlons1[0] > 360: print("ERROR",lons1)
        if rlons2[1] - rlons2[0] > 360: print("ERROR",lons2)
        print(rlons1, rlons2)
        ret = RangeIntersect(rlons1, rlons2, verb, check)
        if check: return ret
        elif ret is not None:
            while ret[1] > 180 : ret = [i-360 for i in ret]
            while ret[0] < -180 : ret = [i+360 for i in ret]
##             return list(AbsoluteLongitude(array(ret), 0.))
            return ret
    if check: return False
    else: return None


def RangeRelativeLongitude(lon, lon0):
    rlon = RelativeLongitude(lon, lon0)
    if   rlon[0] ==  180: rlon[0] = -180
    elif rlon[1] == -180: rlon[1] =  180
    if rlon[0] > rlon[1]:
        if rlon[1] < 0:rlon[1] += 360
        else: rlon[0] -= 360
    if rlon[0] == rlon[1]:
        if rlon[1] < 0:rlon[1] += 360
        else: rlon[0] -= 360
    return rlon


def RangeAbsoluteLongitude(lon, lon0): return RangeRelativeLongitude(lon, -lon0)


def GeoViewAngle(lat, lon, lon0):
    R = (EarthLatRadius(lat) + EarthLatRadius(0))/2.
    beta = GeodeticAngle(0, 0, RelativeLongitude(lon, lon0), lat)
    return degrees(arctan2(R*sin(beta), C_HSAT - R*cos(beta)))


def GeoPixelArea(zen, res, bad=None):
    area = res*res/cos(radians(zen))
    if bad is not None : area[where(zen == bad)] = bad
    return area


## This function removes North and/or South Geostationnary slot file if Fulldisk file is available
## for the same slot to avoid repeated data.  This problem occurs with Goes-E for instance.
#def GeoListCorrection(geolist):
#    for fg in set(geolist):
#        if '_G' in path.split(fg)[1]:
#            for z in ('N', 'S'):
#                f = fg.replace('_G', '_' + z)
#                if f in geolist: geolist.remove(f)
#    return geolist


# This function removes North and/or South Geostationnary slot file if at least one Fulldisk file is
# present within the list to avoid repeated data.  This problem occurs with Goes-E for instance.
def GeoListGlobeOnly(geolist):
    if any(['_G' in x for x in geolist]):
        for f in set(geolist):
            for z in ('_N', '_S'):
                if z in path.split(f)[1]: geolist.remove(f)
    return geolist


# This function removes geostationnary slot files with unrequired zone if at least one file with the
# required zone is present within the list, to avoid repeated data. This problem occurs with Goes files.
def GeoListUnicZone(geolist, zone='G'):
    global z0, zones
    z0, zones = '_' + zone,  [ '_G', '_N', '_S' ]
    zones.remove(z0)
    if   any([z0       in x for x in geolist]): exclude = zones
    elif any([zones[0] in x for x in geolist]): exclude = zones[1]
    else                                           : exclude = []
    outlist = list(geolist)
    for z in exclude:
        for f in set(outlist):
            if z in path.split(f)[1]: outlist.remove(f)
    return outlist

