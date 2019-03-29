#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Definitions or Basic classes used in mol2genparams.py or Molecule type.

Author: Hahnbeom Park and Frank DiMaio 
'''
import sys
from optparse import OptionParser
from Types import *

class OptionClass:
    def __init__(self,argv):
        self.debug = False

        self.init_from_parser(argv)

        # basic options
        self.verbose = self.opt.debug
        self.resname_counter = 0

        # chi angle control
        self.opt.report_puckering_chi = False
        self.opt.report_Hapol_chi = False
        self.opt.report_nbonded_chi = False
        self.opt.report_amide_chi = False
        self.opt.report_ringring_chi = True

        # constants / unchanged
        self.opt.reassign_biaryl_aclass = False
        self.opt.report_as_atype = True
        self.opt.report_as_mmtype = False
        self.opt.define_aro_by_input_geometry = False #unused now
        self.opt.longest_puckering_ring = 6
        self.opt.ring_sampling_sp3_only = True
        
    def init_from_parser(self,argv):
        usagemsg = "USAGE: python mol2genparams.py [-s mol2file or -l mol2filelist] [options]\n"
        usagemsg += "       For more questions, please email hahnbeom@uw.edu or dimaio@uw.edu"
        parser = OptionParser(usage=usagemsg)
        parser.add_option("-s","--inputs",
                          default=[],
                          help="",
                          action="append"
                          )
        parser.add_option("-l",
                          default=False,
                          help="",
                          type="string",
                          )
        parser.add_option("--nm","--resname",
                          dest="resname",
                          default=None,
                          help="Residue name",
                          )
        parser.add_option("--auto_nm",
                          dest="auto_resprefix",
                          default=None,
                          help="Automatically rename resname starting with argument; default L[00-99]",
                          #type="string"
                          )
        parser.add_option("--am1bcc",
                          dest="do_am1bcc_calc",
                          default=False,
                          help="Calculate am1bcc charge",
                          action="store_true"
                          )
        parser.add_option("--prefix",
                          default=None,
                          help="Prefix of output names (prefix.params,prefix_0001.pdb), default as the prefix of input mol2 file",
                          )
        parser.add_option("--debug",
                          default=False,
                          help="Report verbose output for debugging",
                          action="store_true"
                          )
        parser.add_option("--no_output",
                          default=False,
                          help="Do not report params or pdb",
                          action="store_true"
                          )
        parser.add_option("--funcgrp",
                          dest="report_funcgrp",
                          default=False,
                          help="Report functional group assignment to stdout",
                          action="store_true"
                          )
        parser.add_option("--elec_cp_rep",
                          dest="write_elec_cp_rep",
                          default=False,
                          help="Report elec-countpair info to [prefix].elec_cp_ref",
                          action="store_true"
                          )
        parser.add_option("--elec_grpdef",
                          dest="write_elec_grpdef",
                          default=False,
                          help="Report elec-grp-definition info to [prefix].grpref",
                          action="store_true"
                          )
        parser.add_option("--puckering_chi",
                          default=False,
                          help="Define ring puckering torsions as rotatable CHI",
                          action="store_true"
                          )
        parser.add_option("--amide_chi",
                          default=False,
                          help="Define amide as rotatable CHI",
                          action="store_true"
                          )
        parser.add_option("--freeze_ringring",
                          default=False,
                          help="Define  as rotatable CHI",
                          action="store_true"
                          )

        if len(argv) < 2:
            parser.print_help()
            sys.exit()
            
        (self.opt, args) = parser.parse_args(args=argv[1:])
        if self.opt.l:
            self.opt.inputs = [l[:-1] for l in open(self.opt.l)]

        if self.opt.inputs == []:
            parser.print_help()
            sys.exit()
            
        if self.opt.resname != None and self.opt.auto_resprefix != None:
            sys.exit("--nm and --auto_nm cannot be used together!")

        if self.opt.resname != None:
            print ("Renaming ligand residue name as %s"%self.opt.resname)
        else:
            self.opt.resname = "LG1"
            
        if self.opt.auto_resprefix != None:
            self.opt.auto_resname = True
            if len(self.opt.auto_resprefix) > 1:
                print( "argument for --auto_nm should be 1-character! Using prefix 'L' instead (named as L00 ~ L99)" )
                self.opt.auto_resprefix = "L"
        else:
            self.opt.auto_resname = False
            self.opt.auto_resprefix = "L"

    def get_resname(self):
        if self.opt.auto_resname:
            if self.resname_counter > 99:
                print( "WARNING! Residue name will have > 3 characters which may cause issue with Rosetta. Try by splitting list to contain < 100 entries")
            resname = '%1s%02d'%(self.opt.auto_resprefix,
                                 self.resname_counter)
        else:
            resname = self.opt.resname
        return resname
        
class AtomClass:
    def __init__(self,name,atype,hyb,charge):
        self.hyb = hyb
        self.name = name
        self.atype = atype
        self.bonds = []
        self.connected_to_polar = False
        self.has_H = False
        self.is_H = (atype == 2)
        self.aclass = 0
        self.charge = charge
        self.icoord = []
        self.vrt_i = -1 #undefined
        self.root = -1 # undefined
        self.groot = -1 #undefined
        self.ring_index = -1 #undefined
    
    def add_bond(self,atm_connected,order):
        self.bonds.append((atm_connected,order))

    def report_self(self):
        return ' %2s %6s %3d'%(self.atype,self.name,self.hyb)

class BondClass:
    def __init__(self,atm1,atm2,order):
        self.atm1 = min(atm1,atm2)
        self.atm2 = max(atm1,atm2)
        self.order = order
        if self.order > 1:
            self.is_conjugated = True
        else:
            self.is_conjugated = False
        self.cut_bond = False

    def order_in_params(self):
        if self.order == 1 and self.is_conjugated: #use this until we get conjugated type working
            return 2
        else:
            return self.order

class RingClass:
    def __init__(self,atms,cut=None):
        #type:
        #1: puckering & ring-sampling
        #2: puckering & not ring-sampling
        #3: puckering & long (chi-sampling)
        #4: aromatic
        #5: sugar
        self.type = None 
        self.atms = atms
        self.natms = len(atms)
        self.cut_bond = cut

    def has(self,atms):
        for atm in atms:
            if atm not in self.atms:
                return False
        return True

class FunctionalGroupClass:
    def __init__(self,atms,grptype):
        self.atms = atms
        self.grptype = grptype

    def show(self,out):
        l = '%-15s '%self.grptype
        l += '.'.join([atm.name for atm in self.atms])
        l += ' '+'.'.join([ACLASS_ID[atm.aclass] for atm in self.atms])
        out.write(l+'\n')
        
