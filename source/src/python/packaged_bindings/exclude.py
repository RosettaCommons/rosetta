# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   exclude.py
## @brief  functions to exclude varios elements from Python buidings
## @author Sergey Lyskov and William Sheffler

import os, os.path, re, time, commands, sys

_SconsFiles = []

def isFileInScons(fname):
    #print 'isFileInScons', fname, ' --> ', fname in _SconsFiles
    if not _SconsFiles:
        all_scons_files = [f for f in commands.getoutput('ls *.src.settings').split() if f not in ['apps.src.settings', 'devel.src.settings', 'pilot_apps.src.settings']]
        for scons_file in all_scons_files:
            f = file(scons_file).read();  exec(f)
            for k in sources:
                for f in sources[k]:
                    #all_sources.append( k + '/' + f + obj_suffix)
                    _SconsFiles.append( (k + '/' + f).replace('//', '/'))  # some people don't know the right syntax for scons...
        #print '_SconsFiles', _SconsFiles
    return fname in _SconsFiles


BannedFiles = [
    'utility/keys',
    'utility/options',
    'utility/options/keys',
    'utility/pointer',
    'utility/boinc',
    'utility/io',
    'utility/sql_database',
    'utility/factory',
    'utility/tools',
    'utility/signals',
    'utility/json_spirit',
    'utility/py/PyHelper.hh',
    'utility/exit.hh',
    'utility/vectorL.hh', # problem with friend 'swap' functions, definition need to me moved out of class

    'basic/options/keys',
    'numeric/xyzVector.hh', 'numeric/xyzMatrix.hh', 'numeric/xyzTransform.hh',

    'core/id/AtomID_Map.hh', # Stdlib incompatability on osx build, need to manually remap []/() operators to 'get' function in binding.

    'core/chemical/ElectronConfiguration.hh',
    'core/chemical/orbitals',

    'core/environment',

    'core/fragment/picking_old',
    'core/fragment/picking_old/concepts',


    'core/scoring/fiber_diffraction',
    'core/scoring/etable/BaseMembEtableEnergy.hh',  # abandoned?
    'core/scoring/etable/CoarseEtableEnergyCreator.hh',  # not in scons (.hh only)
    'core/scoring/memb_etable/BaseMembEtableEnergy.hh', # abandoned?
    'core/scoring/rna/RNA_FA_Stack.hh', # not in scons (.hh only)
    'core/scoring/methods/GaussianOverlapEnergyCreator.hh', # not in scons
    'protocols/forge/remodel/RemodelLoopMoverCreator.hh',
    #'protocols/toolbox/task_operations/RestrictToMoveMapChiOperationCreator.hh',

    #'core/pack/dunbrack/DunbrackRotamer.hh', # too many template args (more then 48)

    'protocols/abinitio/JumpingFoldConstraints.hh', # not in scons (.hh only)
    'protocols/boinc',
    'protocols/filters/RGFilterCreator.hh', # not in scons
    'protocols/flxbb/InterlockAromaCreator.hh', # not in scons
    'protocols/moves/ChangeFoldTreeMover.hh', # not in scons

    'protocols/jd2/MultiThreadingJob.hh',  # Linker errors, neeed somefunction to be explicitly defined

    #'protocols/jd2/JD2ResourceManagerInputterCreator.hh', # not in scons

    'protocols/canonical_sampling/mc_convergence_checks/MPIBPool_ConvergenceCheck.hh', # not in scons, MPI only
    'protocols/canonical_sampling/mc_convergence_checks/MPIHPool_ConvergenceCheck.hh', # not in scons, MPI only
    'protocols/canonical_sampling/mc_convergence_checks/MPIPool_ConvergenceCheck.hh', # not in scons, MPI only

    'protocols/moves/ReportToDB.hh',  # need sqlite3, not sure if we need it...

    'protocols/swa/rna/StepWiseRNA_Classes.hh', # not in scons (.hh only)

    'protocols/viewer',  # OpenGL
    'protocols/star',  # void* StarAbinitio_main(void*); etc

    #'protocols/wum/WorkUnitManager.hh', # strange linker errors, will deal with it later

    #'protocols/toolbox/task_operations/RestrictToInterfaceCreator.hh', # not in scons?
    #we need  -DBOOST_NO_INITIALIZER_LISTS or gccxml choke on protocols/genetic_algorithm/GeneticAlgorithm.hh

    'numeric/kdtree/WrappedType.hh', # Duplicate class name!

    'protocols/simple_moves/SymmetricFragmentMover.hh', # Somthing with linking/clang (gcc work fine)
    'protocols/nonlocal/DistributionSampler.hh',  # Seems to be abandoned

    # RELEASE exclude's. Do not remove them until got confirmation from the developers.
    'protocols/qsar', # requested by Sam DeLuca
    'protocols/qsar/scoring_grid', # requested by Sam DeLuca
    'protocols/noesy_assign', # requested by Oliver

    'protocols/nonlocal/Chunk.hh', # Problem with GCCXML on GCC 4.0
    'protocols/nonlocal',

    'protocols/medal/MedalMain.hh',  # void* (void*) declaration

    'protocols/wum2/EndPoint.hh',  # Function returning a function
]

def isBanned(fname):
    ''' Check if given path or file name is banned from PyRosetta
    '''
    if fname in BannedFiles:
        return True
    elif fname.endswith('.hh') and os.path.isfile(fname[:-3]+'.cc'):  # Aha! Cpp file present... let's check if it in scons...
        return not isFileInScons(fname[:-3])
    else:
        return False

def getIncludes(files):
    ''' return string with includes for finalize2
    '''
    files = files[:]

    output = ''

    output = set()
    for f in files[:]:
        cc_file = f.replace('.hh', '.cc')
        if os.path.isfile(cc_file):
            files.append(cc_file)

    for f in files:
        lines = commands.getoutput(  "cat %s | grep '^#include <core' | grep '\.hh\|\.hpp'" % f) + '\n'        \
                + commands.getoutput("cat %s | grep '^#include <numeric' | grep '\.hh\|\.hpp'" % f) + '\n'   \
                + commands.getoutput("cat %s | grep '^#include <utility' | grep '\.hh\|\.hpp'" % f) + '\n'   \
                + commands.getoutput("cat %s | grep '^#include <protocols' | grep '\.hh\|\.hpp'" % f) + '\n' \
                + commands.getoutput("cat %s | grep '^#include <basic' | grep '\.hh\|\.hpp'" % f) + '\n' \
                + commands.getoutput("cat %s | grep '^#include <ObjexxFCL' | grep '\.hh\|\.hpp'" % f) + '\n'
                #+ commands.getoutput("cat %s | grep '^#include <boost' | grep '\.hh\|\.hpp'" % f) + '\n'

        for i in lines.split('\n'):
            if i:
                #print 'i=', i
                iname = i[i.index('<')+1 : i.index('>')]
                if os.path.isfile( iname.replace('.fwd', '') ): i = i.replace('.fwd', '')
                output.add(i)

    output = '\n'.join(output)


    lines = list( set( output.split('\n') ) );  lines.sort()
    output = '\n'.join( lines ) + '\n'
    return output

def finalize2(fname, dest, path, module_name='_noname', includes=''):

    f = open(fname);  s = '#include <boost/python.hpp>\n' + includes + f.read();  f.close()

    f = open(fname,'w')
    f.write(s)
    f.close()
