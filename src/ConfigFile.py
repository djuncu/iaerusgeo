#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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


from __future__ import print_function
#from builtins import str
#from builtins import map
#from builtins import range
#from builtins import object
import getopt
import sys
import os
from string import Template
#import parser

try               : from collections import OrderedDict
except ImportError: from ordereddict import OrderedDict

INFO__, WARN__, ERROR__ = (0,1,2)
DFT_VERBOSE__ = WARN__
BOOL_TRUE  = set([ 'true' , 'vrai', 'yes', 'oui' ])
BOOL_FALSE = set([ 'false', 'faux', 'no' , 'non' ])
BOOL_VALUES = BOOL_TRUE | BOOL_FALSE
#========================================================================

class ConfigError(Exception):
#----------------------------------------------------------------------
    def __init__(self, value): self.value = value
#----------------------------------------------------------------------
    def __str__(self): return repr(self.value)

#========================================================================

class Config(object):
#----------------------------------------------------------------------
    def __init__(self, verif=None, nones=None, verbose=DFT_VERBOSE__, comment='#',
                 include='@', sameline='\\', variable='$', listsep=',', getall=False):
        self.comment, self.include, self.sameline, self.variable, self.listsep = \
                      (comment, include, sameline, variable, listsep)
        self.vars = OrderedDict()
        self.verif, self.nones, self.verbose, self.getall = (verif, nones, verbose, getall)
        if verbose is None: self.verbose = DFT_VERBOSE__
#----------------------------------------------------------------------
    def nonevalues(self):
        if self.nones is None: return
        for v in self.nones:
            for k in list(self.cfg.keys()):
                if self.cfg[k] is v: self.cfg[k] = None
#----------------------------------------------------------------------
    def set(self, indict, path):
        dir, nam = os.path.split(path)

        try: f = open(path, "r")
        except IOError as e: self.err("Unable to open file %s" % path)

        d = OrderedDict()
        l, ll = ("", "")
        try:
            for l in f:
                l0 = l.strip()
                if len(l0) == 0 or l0[0] == self.comment: continue
                ll += l0.replace(self.sameline, ' ')
                if l0[-1] == self.sameline: continue
                if ll[0] == self.include:
                    self.set(indict, os.path.join(dir, self.evaluate(ll[1:], path)))
                else:
                    name, value = ll.split("=")
                    if name[0] == self.variable: self.vars[name[1:].strip()] = self.evaluate(value.strip(), path)
                    else: d[name.strip()] = value.strip()
                ll = ""
            f.close()
        except IOError as e: self.err("Unable to read file %s" % path)
        except Exception as e: self.err("Error in file %s: %s" % (path, l))

        self.update(indict, d, "file %s" % path)
#----------------------------------------------------------------------
    def evaluate(self, value, msg):
        try: return Template(value).substitute(self.vars).strip()
        except KeyError as e: self.err("%s: Undefined variable used in %s" % (value, msg))
#----------------------------------------------------------------------
    def typebool(self, value):
        if str(value).isdigit() : return int(value) != 0
        elif value.lower() in BOOL_VALUES : return value.lower() in BOOL_TRUE
        raise ValueError("Invalid boolean value : " + str(value))
#----------------------------------------------------------------------
    def typevar(self, value, dtype):
        if dtype != type(True): return dtype(value)
        return self.typebool(value)
#----------------------------------------------------------------------
    def addthiskey(self, k, value):
        val, values_str = None, value.split(self.listsep)
        try   : val = list(map(int, values_str))
        except:
            try   : val = list(map(float, values_str))
            except:
                try   : val = [self.typebool(x) for x in values_str]
                except: val = None
        if val is None   : self.cfg[k] = value
        elif len(val) > 1: self.cfg[k] = val
        else             : self.cfg[k] = val[0]
#----------------------------------------------------------------------
    def update(self, indict, dico, msg):
        for k in list(dico.keys()):
            if k not in indict:
                if self.getall:
                    # First evaluate variables
                    dicoval = self.evaluate(dico[k].strip(), msg)
                    self.addthiskey(k, dicoval)
                else: self.info("unused key %s in %s" % (k, msg))
        for k in list(indict.keys()):
            dtype, st, default = indict[k]
            if k in dico:
                # First evaluate variables
                dicoval = self.evaluate(dico[k].strip(), msg)
                # Then convert values with type
                try:
                    if len(dicoval) != 0 and dicoval[0] == '[':
                        # We have a list
                        if dicoval[-1] != ']': self.err("%s : Unbalanced delimitors in %s" % (k, msg))
                        elif len(dicoval) == 2: self.cfg[k] = []
                        else:
                            values = dicoval[1:-1].split(self.listsep)
                            if not isinstance(self.cfg[k], list):
                                self.cfg[k] = values
                            for i in range(len(values)):
                                val = values[i].strip()
                                if val == '': val = default
                                else: val = self.typevar(val, dtype)
                                if i < len(self.cfg[k]): self.cfg[k][i] = val
                                else: self.cfg[k].append(val)
                    else: self.cfg[k] = self.typevar(dicoval, dtype)
                except ValueError as e: self.err("%s : Unable to convert %s to %s" % (k, dico[k], dtype))
#----------------------------------------------------------------------
    def fget(self, indict, path):
        self.reinit(indict)
        self.set(indict, path)
        self.validate()
        return self.cfg
#----------------------------------------------------------------------
    def get(self, indict):
        self.reinit(indict)
        optstr = ""
        for k in list(indict.keys()): optstr += indict[k][1]
        optlist, args = getopt.getopt(sys.argv[1:], optstr)

        # Parse command files if some
        if len(args) != 0:
            for i in range(len(args)): self.set(indict, args[i])

        # Now parse command line
        dopt = OrderedDict()
        for opt in optlist:
            for k in list(indict.keys()):
                if len(indict[k][1]) == 0: continue
                if opt[0][1] == indict[k][1][0]:
                    if opt[1] == '' and indict[k][0] == 'int': dopt[k] = '1'
                    else: dopt[k] = opt[1]
                    break
        self.update(indict, dopt, "command line")

        # Dictionnary local validation
        self.validate()

        return self.cfg
#----------------------------------------------------------------------
    def reinit(self, indict):
        self.cfg = OrderedDict()
        self.errors, self.errors = (False, False)
        for k in list(indict.keys()): self.cfg[k] = indict[k][2]
#----------------------------------------------------------------------
    def validate(self):

        # Dictionnary local validation
        for k in list(self.cfg.keys()):
            if isinstance(self.cfg[k], list):
                for kk in self.cfg[k]:
                    if kk is None:
                        self.err("%s : at least one required value is not present in configuration" % k)
            elif self.cfg[k] is None:
                self.err("%s : required value not present in configuration" % k)

        self.nonevalues()

        # Dictionnary user validation if needed
        if self.verif is not None:
            ok, msg = self.verif(self.cfg)
            if not ok: self.err(msg)

        # Check if errors appeared
        if self.errors: raise ConfigError("Configuration Errors")

#----------------------------------------------------------------------
    def err(self, msg):
        self.errors = True
        if self.verbose <= ERROR__: self.information("Error", msg)
#----------------------------------------------------------------------
    def warn(self, msg):
        self.warns = True
        if self.verbose <= WARN__ : self.information("Warning", msg)
#----------------------------------------------------------------------
    def info(self, msg):
        if self.verbose <= INFO__ : self.information("Information", msg)
#----------------------------------------------------------------------
    def information(self, errtype, msg):
        if not isinstance(msg, list): msgs = [msg]
        else: msgs = msg
        for m in msgs:
            sys.stderr.write("***Config " + errtype + "*** " + m + "\n")
#----------------------------------------------------------------------
    def view(self, indict):
        tmpcfg = self.get(indict)
        for key, val in list(tmpcfg.items()): print(key, '=', val)
#----------------------------------------------------------------------
    def fview(self, indict, path):
        tmpcfg = self.fget(indict,path)
        for key, val in list(tmpcfg.items()): print(key, '=', val)
#----------------------------------------------------------------------
    def display(self, prefix):
        print("\n" + prefix + "Configuration:")
        keys = list(self.cfg.keys())
        keys.sort()
        for k in keys: print(prefix + "  " + k + " : ", self.cfg[k])
        print("")
#----------------------------------------------------------------------
def GetDefaultCfg(cfg, nones):
    d = OrderedDict()
    try:
        for k in list(cfg.keys()):
            t, o, v = cfg[k]
            if v in nones: v = None
            if v is None: d[k] = None
            else: d[k] = t(v)
        return d
    except:
        return None
#----------------------------------------------------------------------
