# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Jack Maguire

import pyrosetta
from pyrosetta.rosetta import core, numeric, protocols, basic
import re

def verbose():
    return False

class IsBad():
    def __init__(self):
        self.bad = False
        self.count = 0
    def register(self):
        self.bad = True
        self.count += 1

global_is_bad = IsBad()

class MyClass:
    def __init__(self):
        pass

class GreenListManager:
    def __init__(self):
        # Only add things here that are overloaded in a way that results in specification, not divergent behavior
        self.greenlist_regex = [
            "<class 'pyrosetta.rosetta.protocols.mean_field.jagged_array_*",
            "<class 'pyrosetta.rosetta.protocols.mean_field.*",
            "<class 'pyrosetta.rosetta.protocols.hbnet.HBNetScore'>$",
            "<class 'pyrosetta.rosetta.core.chemical.PoseResidueTypeSet'>$",
            "<class 'pyrosetta.rosetta.protocols.multistage_rosetta_scripts.cluster.metrics.*",
            "<class 'pyrosetta.rosetta.protocols.evaluation.SingleValuePoseEvaluator*",
            "<class 'pyrosetta.rosetta.core.pack.dunbrack.DunbrackRotamer_*",
        ]
    def bypass(self, name):
        for g in self.greenlist_regex:
            if re.match( g, name ) is not None:
                return True
        return False

greenlist = GreenListManager()

def extract_class_signature_from_overload_line( l ):
    #     1. (pose: pyrosetta.rosetta.core.pose.Pose, residue_positions: pyrosetta.rosetta.utility.vector1_bool) -> None
    l = l.split("->")[0]
    #    1. (pose: pyrosetta.rosetta.core.pose.Pose, residue_positions: pyrosetta.rosetta.utility.vector1_bool)
    s1 = l.find("(")+1
    s2 = l.rfind(")")
    if s1==s2: return []
    l = l[s1:s2]
    #pose: pyrosetta.rosetta.core.pose.Pose, residue_positions: pyrosetta.rosetta.utility.vector1_bool
    ls = [a.strip().rstrip() for a in l.split( "," )]
    #['pose: pyrosetta.rosetta.core.pose.Pose', 'residue_positions: pyrosetta.rosetta.utility.vector1_bool']
    ls = [ a[a.find(":")+1:].strip() for a in ls ]
    #['pyrosetta.rosetta.core.pose.Pose', 'pyrosetta.rosetta.utility.vector1_bool']
    ls = [ eval(a) for a in ls ]
    #[<class 'pyrosetta.rosetta.core.pose.Pose'>, <class 'pyrosetta.rosetta.utility.vector1_bool'>]
    # returns a list of classes
    return ls
    
def signatures_conflict( sig1, sig2, only_enforce_distinct_parents=True ):
    # if only_enforce_distinct_parents=True, we will not complain if foo(<class B>) shadows foo(<class B>), we will only complain if foo(<class B>) shadows foo(<class A>) where A is a parent of B
    if len(sig1) != len(sig2): return False
    
    if (sig1 == sig2) and only_enforce_distinct_parents: return False

    for i in range(len(sig1)):
        class1 = sig1[i]
        class2 = sig2[i]
        check12 = issubclass( class1, class2 ) and class1 is not bool
        check21 = issubclass( class2, class1 ) and class2 is not bool
        if not (check12 or check21):
            # if neither class clashes with the other one, these two signatures are safe
            return False

    return True
    

def test_function( F, custom_name = None ):

    if custom_name is None:
        custom_name = str(F)

    if greenlist.bypass( custom_name ):
        return
    
    failed = True
    try:
        F( MyClass() )
        failed = False
    except Exception as E:
        exception_str = str(E)

    assert failed, F
    assert "The following argument types are supported" in exception_str or "No module named 'numpy'" in exception_str, exception_str

    signatures = []
    for l in exception_str.split("\n"):
        if l[0:4] != "    " or l[5] != ".": continue
        #if "->" not in l: continue

        if "*" in l:
            if verbose(): print( "Skipping weird signature", l )
            # like     2. pyrosetta.rosetta.core.chemical.MutableChiRecord(atm_vec: pyrosetta.rosetta.utility.vector1_void_*)
            continue

        if "::" in l:
            if verbose(): print( "Skipping incompletely bound signature", l )
            continue

        if "Tuple[" in l or "Callable[" in l or "tuple[" in l:
            if verbose(): print( "Skipping complicated signature (for now?)", l )
            continue

        try:
            sig = extract_class_signature_from_overload_line( l )
        except Exception as e:
            if str(e).startswith("name '") and str(e).endswith("' is not defined"):
                if verbose(): print( "Skipping incompletely bound signature", l, ":", str(e) )
                continue
            print( "ERROR" )
            print( l )
            print( "-----" )
            print( e )
            exit( 0 )
        signatures.append( sig )

    for i in range(len(signatures)):
        for j in range(i+1,len(signatures)):
            if signatures_conflict( signatures[i], signatures[j] ):
                print( f"CONFLICT `{custom_name}` --- `{signatures[i]}` --- `{signatures[j]}`" )
                global_is_bad.register()

def test_class( C ):
    # 1. test __init__
    try:
        test_function( C )
    except AssertionError as E:
        if "No constructor defined" in str(E):
            if verbose(): print( "No constructor defined for", C )
            pass
        else:
            raise E

    # 2. test functions
    for cname in dir(C):
        if cname.startswith("__"): continue
        #if cname == "__init__": continue

        cattr = getattr(C,cname)
        type_str = str( type( cattr ) )
        if type_str == "<class 'instancemethod'>":
            try:
                test_function( cattr, custom_name=f"{C}::{cname}" )
            except AssertionError as e:
                if "takes no arguments" in str(e):
                    pass
                else:
                    raise e
        elif type_str == "<class 'module'>":
            print( cattr )
            exit( 0 ) #testing to see if this is possible
            test_module( cattr )
        else:
            if verbose(): print( "skipping type", type_str, "for", cname, "in", C )

def test_module( D ):
    for dname in dir(D):
        dattr = getattr(D,dname)
        type_str = str( type( dattr ) )
        if type_str == "<class 'pybind11_builtins.pybind11_type'>":
            if dname.endswith( "Creator" ):
                if verbose(): print( "skipping creator type", dname )
            else:
                test_class( dattr )
        elif type_str == "<class 'builtin_function_or_method'>":
            test_function( dattr )
        elif type_str == "<class 'module'>":
            test_module( dattr )
        else:
            if verbose(): print( "skipping type", type_str, "for", dname )

for module in [core, numeric, protocols, basic]:
    test_module( module )
if global_is_bad.bad:
    raise Exception(f"Found {global_is_bad.count} shadowing instances")
