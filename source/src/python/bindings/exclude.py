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


#if sys.platform != 'win32':
#    import pyplusplus, pygccxml
#    from pyplusplus.module_builder import call_policies


MAKE_SURE_IS_COPYABLE = [
  "core::conformation::Residue",
]

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
    # tmp debug only
    # 'protocols/features/AtomAtomPairFeatures.hh',
    # 'protocols/features/AtomInResidueAtomInResiduePairFeatures.hh',
    # 'protocols/features/HBondFeatures.hh',
    # 'protocols/features/PairFeatures.hh',
    # 'protocols/features/ProteinBackboneAtomAtomPairFeatures.hh',
    # 'protocols/features/ResidueBurialFeatures.hh',
    # 'protocols/features/ResidueScoresFeatures.hh',
    # 'protocols/features/RotamerFeaturesCreator.hh',
    # 'protocols/features/ScoreFunctionFeatures.hh',
    # 'protocols/features/ScoreTypeFeatures.hh',
    # 'protocols/features/ScoreTypeFeatures.hh',

    # 'protocols/features/RotamerFeatures.hh',  # array with negative size

    # Temporary, remove after transition to SP
    'utility/query', 'core/coarse',


    'utility/PyHelper.hh', 'utility/keys', 'utility/options', 'utility/options/keys',
    #'utility/pointer', 'utility/pointer/boost', 'utility/pointer/std',
    #'utility/pointer/std', 'utility/pointer/refcount',
    'basic/options/keys', 'utility/exit.hh',
    'numeric/xyzVector.hh',

    #'basic/Tracer.hh',
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

    'utility/vectorL.hh', # problem with friend 'swap' functions, definition need to me moved out of class
    'numeric/kdtree/WrappedType.hh', # Duplicate class name!

    'protocols/simple_moves/SymmetricFragmentMover.hh', # Somthing with linking/clang (gcc work fine)




    'utility/io', #/ozstream.hh',  # need bindings for std enum: std::_Ios_Openmode
    'utility/sql_database/DatabaseSessionManager.hh',  # SQLite have its own bindigs in Python

    'protocols/nonlocal/DistributionSampler.hh',  # Seems to be abandoned

    # RELEASE exclude's. Do not remove them until got confirmation from the developers.
    'protocols/qsar', # requested by Sam DeLuca
    'protocols/qsar/scoring_grid', # requested by Sam DeLuca
    'protocols/noesy_assign', # requested by Oliver

    'protocols/nonlocal/Chunk.hh', # Problem with GCCXML on GCC 4.0
    'protocols/nonlocal',

    'protocols/medal/MedalMain.hh',  # void* (void*) declaration

    'protocols/wum2/EndPoint.hh',  # Function returning a function

    #'protocols/rpc',  # Some problem with parser, but in general: do we need RPC calls in PyRosetta? (it easier to implement it in Python anyway)


#'numeric/xyzTransform.hh',

# 'protocols/sic_dock/RigidScore.hh',
# 'protocols/sic_dock/Rose.hh',
# 'protocols/sic_dock/SICFast.hh',
# 'protocols/sic_dock/designability_score.hh',
# 'protocols/sic_dock/loophash_util.hh',
# 'protocols/sic_dock/read_biounit.hh',
# 'protocols/sic_dock/types.hh',
# 'protocols/sic_dock/util.hh',
# 'protocols/sic_dock/xyzStripeHashPose.hh',
# 'protocols/sic_dock/xyzStripeHashPoseWithMeta.hh',

    #'protocols/frag_picker', # whole dir not in scons (ie been moved now), temporary

# Lion workarround
    # problem with gccxml using gcc 4.2.1 ???
    #'utility/signals',
    #'numeric/random',  # problem with interger been too large (gccxml)
'''
    'numeric/random/DistributionSampler.hh',  # problem with interger been too large (gccxml)  #include <boost/math/distributions.hpp>


    'core/io/silent/ProteinSilentStruct.tmpl.hh',

'core/io/silent/BinaryProteinSilentStruct.hh',
'core/io/silent/BinaryRNASilentStruct.hh',
'core/io/silent/ProteinSilentStruct.hh',
'''


]

def isBanned(fname):
    ''' Check if given path or file name is banned from PyRosetta
    '''
    if fname in BannedFiles: return True
    elif fname.endswith('.hh')  and  os.path.isfile(fname[:-3]+'.cc'):  # Aha! Cpp file present... let's check if it in scons...
        return not isFileInScons(fname[:-3])
    else: return False


def namespace(ns):
    ''' check if namespace 'ns' should be excluded from generating any bindings.
        This function also exclude all unrelated folders like .svn
    '''
    nslist = [
              'utility/excn',
              'numeric/deriv',
              'numeric/geometry',
              'numeric/internal',
              'numeric/interpolation',
              'numeric/interpolation/full',
              'numeric/interpolation/periodic_range/full',
              'numeric/interpolation/periodic_range/half',
              'numeric/interpolation/periodic_range/periodic_value/full',
              'numeric/interpolation/periodic_range/periodic_value/half',
              'numeric/SVD',
              'numeric/model_quality',
              'numeric/kdtree',
              'core/conformation/signals',
              'core/options',
              'core/options/keys',
              'core/fragment/io',
              'core/fragment/picking/concepts',
              'core/fragment/picking/vall/eval',
              'core/fragment/picking/vall/scores',
              #'core/fragment',
              'core/coarse', #??? new for 1 file scheme
              'core/scoring/NV',

              'core/scoring/disulfides',
              'core/scoring/carbon_hbonds',
              #'core/scoring/hbonds',
              'core/scoring/hbonds/hbtrie',
              'core/scoring/dunbrack',
              'core/scoring/elec',
              'core/scoring/etable',
              'core/scoring/etable/count_pair',
              'core/scoring/etable/etrie',
              'core/scoring/geometric_solvation',
              'core/scoring/packstat',
              'core/scoring/packing',
              'core/scoring/symE',
              #'core/scoring/constraints', # temp

              'core/io/raw_data',
              'basic/database',
              'core/io/sequence_comparation',
              'core/io/serialization',
              'core/io/silent',
              'core/util',
              'core/pack/interaction_graph',
              'core/pack/rotamer_set',
              #'core/pack/annealer',
              'core/pose/signals',
              'core/scoring/rna',

              #'protocols/abinitio',
              'protocols/boinc',
              'protocols/branch_angle',
              'protocols/checkpoint',
              'protocols/cluster',
              'protocols/comparative_modeling',
              'protocols/ddg',
              'protocols/enzdes',
              'protocols/evaluation',
              'protocols/filters',
              'protocols/flexpack',
              'protocols/flexpack/rotamer_set',
              'protocols/frags',
              'protocols/jd2',
              'protocols/jd2/archive',
              'protocols/geometry',
              'protocols/genetic_algorithm',
              'protocols/hotspot_hashing',
              'protocols/jobdist',
              'protocols/jumping',
              'protocols/looprelax',
              #'protocols/loops',
              'protocols/ligand_docking/ligand_options',
              'protocols/motifs',
              'protocols/moves/kinematic_closure',
              'protocols/multistate_design',
              'protocols/optimize_weights',
              'protocols/protein_interface_design',
              'protocols/rbsegment_Moves',
              'protocols/rna',
              'protocols/dna',
              'protocols/rna_denovo',
              'protocols/smanager',
              'protocols/toolbox',
              'protocols/toolbox/PoseMetricCalculators',
              'protocols/viewer',
              'protocols/topology_broker/weights',
              'protocols/mpi',
              'protocols/protein_interface_design/movers',

              'utility/boinc',
              'utility/io',
              'utility/query',
              'utility/factory',
              'utility/options',
              'utility/options/keys',
              'utility/keys',
              'utility/tools',
              'utility/signals',

    ]
    if ns.find('.svn') >= 0: return True
    return ns in nslist




exclude_header_list = []


def mb_exclude(path, mb, hfile):
    def E(*args):
        obj = mb
        fun = True
        try:
            for a in args:
                if fun:
                    obj = obj.__getattribute__(a)
                    fun = False
                else:
                    obj = obj(a)
                    fun = True
            if not fun:
                obj()

        except pygccxml.declarations.matcher.declaration_not_found_t: pass

    def E_(v):
        obj = mb
        try: eval(v)
        except pygccxml.declarations.matcher.declaration_not_found_t: pass


    if path == 'core/pose':
        try:
            pass
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        if hfile == 'core/pose/util.hh':
            mb.free_function("delete_comment").exclude()


    if path == 'core/pose/metrics':
        #mb.class_('PoseMetricCalculator').exclude()  # Pure virtual class
        pass


    if path == 'core/chemical':
        try: mb.class_("ResidueTypeSet").constructors().exclude() # default arg for ctor causes problems
        except pygccxml.declarations.matcher.declaration_not_found_t: pass


    if path == 'core/conformation':
        #t = mb.classes( lambda decl: decl.name.startswith( 'impl' ) )
        def fff(name):
            if name.find( 'std::map<core::id::StubID,core::kinematics::RT') >= 0: return True
            else: return False

        try:
            mb.class_("Conformation").mem_funs('insert_fragment').exclude()
            #mb.class_("Conformation").mem_funs('res_begin').exclude() -- apl note: removing this member function; so commenting out this line
            #mb.class_("Conformation").mem_funs('res_end').exclude() -- apl note: removing this member function; so commenting out this line
            #mb.class_("Conformation").mem_funs('get_stub_transform').exclude()

        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        if hfile=='core/conformation/Residue.hh':
            residue = mb.class_( 'Residue' )
            residue.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )

        if hfile=='core/conformation/Atom.hh':
            residue = mb.class_( 'Atom' )
            residue.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )



    if path == 'core/fragment':
        if hfile=='core/fragment/FragID_Iterator.hh':
            mb.class_("FragID_Iterator").exclude()

        if hfile=='core/fragment/util.hh':
            mb.free_function("flatten_list").exclude()  # 'undefined symbol' error on python import

        if hfile=='core/fragment/FragData.hh':
            def pr(x):
                print x
                print "  access_type:", x.access_type
                return True
            # FragData( Size nr_res ) <-- private
            #mb.class_("FragData").constructor(arg_types=['platform::Size']).exclude()
            for c in mb.class_("FragData").constructors(lambda x: pr(x)):
                if c.access_type == 'protected': c.exclude()
                #print c



        #E_(" mb.class_('SecstructSRFD').mem_funs('clone').exclude() ")

        #E('class_', "BBTorsionSRFD", 'mem_funs', "clone", 'exclude') # mb.class_("BBTorsionSRFD").mem_funs("clone").exclude()


    if path == 'core/graph':
        #mb.class_("PointGraphVertexData").var('NUM_EDGES_TO_RESERVE').exclude()  # static const
        E('class_', "PointGraphVertexData", 'var', 'NUM_EDGES_TO_RESERVE', 'exclude' )  # static const


    if path == 'core/kinematics':
        #mb.class_("Edge").var('PEPTIDE').exclude()
        E('class_', "Edge", 'var', 'PEPTIDE', 'exclude')
        #mb.class_("Edge").var('CHEMICAL').exclude()
        E('class_', "Edge", 'var', 'CHEMICAL', 'exclude')

        if hfile == 'core/kinematics/MoveMap.hh':
            #mb.class_('MoveMap').member_function('find').exclude()
            mb.class_('MoveMap').mem_funs('find').exclude()


        if hfile == 'core/kinematics/FoldTree.hh':
            cl = mb.class_( 'FoldTree' )
            cl.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )

        if hfile == 'core/kinematics/Jump.hh':
            #E('free_function', "distance", 'exclude')
            mb.free_functions("distance").exclude()
            #mb.class_('Jump').mem_funs('distance').exclude()


    if path == 'core/kinematics/tree':
        #mb.class_('BondedAtom').member_function('clone').exclude()
        E('class_', 'BondedAtom', 'member_function', 'clone', 'exclude')
        #mb.class_('JumpAtom').member_function('clone').exclude()
        E('class_', 'JumpAtom', 'member_function', 'clone', 'exclude')

        if hfile == 'core/kinematics/tree/Atom.hh':
            mb.free_functions("distance").exclude()
            mb.free_functions("distance_squared").exclude()


    if path == 'core/pack':
        if hfile == 'core/pack/pack_rotamers.hh':
            #E('free_function', "pack_rotamers_run", 'exclude')
            #fns = mb.free_function("pack_rotamers_run")
            #for f in fns: f.exclude()
            #mb.free_function("pack_rotamers_run", return_type="core::Real").exclude()
            mb.free_functions("pack_rotamers_run").exclude()
            #print 'QQQ', mb.free_functions("pack_rotamers_run")
            #print '______________________________'

            #f = mb.free_function('symmetric_pack_rotamers_run')
            #print f, dir(f)
            #print 'Retutn type:', f.return_type

            E('free_function', 'symmetric_pack_rotamers_run', 'exclude')


    if path == 'core/pack/task':
        try:
            #mb.class_("IGEdgeReweighter").exclude()
            mb.class_("IGEdgeReweightContainer").mem_funs("reweighters_begin").exclude()
            mb.class_("IGEdgeReweightContainer").mem_funs("reweighters_end").exclude()

        except pygccxml.declarations.matcher.declaration_not_found_t: pass


    if path == 'core/pack/task/operation':
        try:
            mb.class_("ResFilterFactory").mem_funs("newResFilter").exclude()

        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        E('class_', 'TaskOperationFactory', 'member_function', 'newTaskOperation', 'exclude')


    if path == 'core/import_pose':
        if hfile=='core/import_pose/import_pose.hh':
            #E('free_functions', 'pose_from_pose', 'exclude')
            pass



    if path == 'core/pose':
        pass
        #if hfile=='core/pose/Pose.hh':
        #    residue = mb.class_( 'Pose' )
        #    residue.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )


    if path == 'core/scoring':
        #mb.class_( "ResidueNeighborIterator" ).exclude()  # pure virtual class, can't be created
        E('class_', "ResidueNeighborIterator", 'exclude')
        #mb.class_( "ResidueNeighborConstIterator" ).exclude()  # pure virtual class, can't be created
        E('class_', "ResidueNeighborConstIterator", 'exclude')  # pure virtual class, can't be created

        try:
            mb.class_("ScoreFunction").mem_funs("long_range_energies_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("long_range_energies_end").exclude()
            mb.class_("ScoreFunction").mem_funs("cd_2b_intrares_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("cd_2b_intrares_end").exclude()
            mb.class_("ScoreFunction").mem_funs("ci_2b_intrares_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("ci_2b_intrares_end").exclude()
            mb.class_("ScoreFunction").mem_funs("all_energies_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("all_energies_end").exclude()

            mb.class_("ScoreFunction").mem_funs("ci_lr_2b_methods_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("ci_lr_2b_methods_end").exclude()
            mb.class_("ScoreFunction").mem_funs("cd_lr_2b_methods_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("cd_lr_2b_methods_end").exclude()

            mb.class_("ScoreFunction").mem_funs("ws_methods_begin").exclude()
            mb.class_("ScoreFunction").mem_funs("ws_methods_end").exclude()

        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try:
            mb.class_("EMapVector").mem_funs("begin").exclude() # these aren't real iterators... just Real*'s
            mb.class_("EMapVector").mem_funs("end").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try:
            mb.class_("TwoBodyEMapVector").mem_funs("begin").exclude() # these aren't real iterators... just Real*'s
            mb.class_("TwoBodyEMapVector").mem_funs("end").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        if hfile=='core/scoring/ScoreFunction.hh':
            scorefnx = mb.class_( 'ScoreFunction' )
            scorefnx.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )

        if hfile=='core/scoring/Energies.hh':
            Energies = mb.class_( 'Energies' )
            Energies.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )

        if hfile=='core/scoring/ScoringManager.hh':
            mb.class_("ScoringManager").mem_fun("get_DDPLookupTable").exclude()  # for some reason boost look up this function in core::conformation instead of core::scoring


    if path == 'core/scoring/constraints':
        E('class_', 'CircularHarmonicFunc', 'member_function', 'clone', 'exclude')  #mb.class_('CircularHarmonicFunc').member_function('clone').exclude()
        E('class_', 'CircularPowerFunc', 'member_function', 'clone', 'exclude')  #mb.class_('CircularPowerFunc').member_function('clone').exclude()
        E('class_', 'CharmmPeriodicFunc', 'member_function', 'clone', 'exclude')  #mb.class_('CharmmPeriodicFunc').member_function('clone').exclude()

        E('class_', 'GaussianFunc', 'member_function', 'clone', 'exclude') #mb.class_('GaussianFunc').member_function('clone').exclude()
        try: mb.class_('HarmonicFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('MixtureFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('ScalarWeightedFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_('BoundFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('ConstantFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('EtableFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('PeriodicBoundFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('PeriodicFunc').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_('MultiConstraint').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('Constraint').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass


    if path == 'core/scoring/electron_density':
        if hfile == 'core/scoring/electron_density/util.hh':
            #mb.free_function('alignVectorSets').exclude()
            pass


    if path == 'core/scoring/hbonds':
        if hfile == "__core/scoring/hbonds/HBondSet.hh":
            try:
                l = mb.decls("get_residue_residue_hbond_energy")
                for i in l: i.exclude()
            except pygccxml.declarations.matcher.declaration_not_found_t: pass

        if hfile == 'core/scoring/hbonds/HBondEnergy.hh':
            #core::scoring::hbonds::HBondEnergy::evaluate_rotamer_pair_energies
            E('class_', "HBondEnergy", 'member_function', "evaluate_rotamer_pair_energies", 'exclude')

        if hfile == 'core/scoring/hbonds/types.hh':
            mb.free_function('HBEval_lookup_initializer').exclude()
            mb.var('HBEval_lookup').exclude()




    if path == 'core/scoring/rna':
        E('class_', "RNA_TorsionEnergy", 'member_function', "indicate_required_context_graphs", 'exclude')  #mb.class_( "RNA_TorsionEnergy" ).member_function( "indicate_required_context_graphs" ).exclude()


    if path == 'core/scoring/methods':
        try: mb.class_('ChainbreakEnergy').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_('DistanceChainbreakEnergy').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_('LinearChainbreakEnergy').member_function('clone').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_('MMBondAngleEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('MMTorsionEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('OmegaTetherEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('P_AA_pp_Energy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('RMS_Energy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('DunbrackEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('RamachandranEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('ReferenceEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('WaterAdductHBondEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_('WaterAdductIntraEnergy').member_function('indicate_required_context_graphs').exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        #mb.class_('RG_Energy').member_function('clone').exclude()  # 'undefinded' error in Python, not in scons


    if path == 'core/scoring/trie':
        E('class_', 'TrieCollection', 'member_function', 'clone', 'exclude')  #mb.class_('TrieCollection').member_function('clone').exclude()


    if path == 'core/sequence':
        try: mb.class_("L1ScoringScheme").mem_funs("clone").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_("L1ScoringScheme").mem_funs("score").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_("MatrixScoringScheme").mem_funs("clone").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_("MatrixScoringScheme").mem_funs("score").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_("ProfSimScoringScheme").mem_funs("clone").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_("ProfSimScoringScheme").mem_funs("score").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_("SequenceProfile").mem_funs("clone").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass

        try: mb.class_("SimpleScoringScheme").mem_funs("clone").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass
        try: mb.class_("SimpleScoringScheme").mem_funs("score").exclude()
        except pygccxml.declarations.matcher.declaration_not_found_t: pass



    if path == 'protocols/abinitio':
        #mb.decl('::std::vector< core::pose::Pose >').exclude()
        #mb.decl('PoseList').exclude()
        if hfile == 'protocols/abinitio/ConstraintFragmentMover.hh':
            mb.class_('ConstraintFragmentMover').mem_funs("choose_fragment").exclude()

        if hfile == 'protocols/abinitio/KinematicAbinitio.hh':
            mb.class_('JumpingFoldConstraintsWrapper').mem_funs("register_options").exclude()

        if hfile == 'protocols/abinitio/KinematicControl.hh':
            mb.class_('KinematicControl').mem_funs("remove_chainbreak_variants").exclude()

        #if hfile == 'protocols/abinitio/ConstraintFragmentMover.hh':
        #    removeProtectedConstructos(mb, "ConstraintFragmentMover")


    if path == 'protocols/branch_angle':
        mb.free_function("bonded_neighbor_all_res").exclude()


    if path == 'protocols/docking':
        if hfile == 'protocols/docking/DockingProtocol.hh':
            mb.class_('DockingProtocol').mem_funs("show").exclude()
            #mb.class_('DockingProtocol').mem_funs("DockingProtocol").exclude()
            mb.class_('DockingProtocol').constructors().exclude()

            #c = mb.class_('DockingProtocol').constructors()
            #c[1].exclude()
            #c[0].exclude()
            #c[0].exclude()

        if hfile == 'protocols/docking/DockingLowRes.hh':
            mb.class_('DockingLowRes').mem_fun("show").exclude()


    if path == 'protocols/dna':
        mb.class_('PDBOutput').operator("()").exclude()


    if path == 'protocols/jobdist':
        if hfile == 'protocols/jobdist/standard_mains.hh':
            mb.free_function("main_plain_mover").exclude()
            mb.free_function("main_atomtree_diff_mover").exclude()
            mb.free_function("main_plain_pdb_mover").exclude()


    if path == 'protocols/ligand_docking':
        pass
        #mb.free_function("set_sphere").exclude()


    if path == 'protocols/moves':
        if hfile == 'protocols/moves/Mover.hh':
            #mb.class_('Mover').member_function('parse_tag').exclude()
            mb.class_('Mover').member_function('parse_my_tag').exclude()

        if hfile == 'protocols/simple_moves/PackRotamersMover.hh':
            mb.class_('PackRotamersMover').member_function('run').exclude()  # Template problem

        #if hfile == 'protocols/moves/WobbleMover.hh':
        #    mb.class_( "WobbleMover" ).exclude()  # Some constructors not defined

        if hfile == 'protocols/moves/SwitchResidueTypeSetMover.hh':
            c = mb.class_( 'SwitchResidueTypeSetMover' )
            c.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )

        if hfile == 'protocols/moves/ScoreMover.hh':
            #print dir(mb.class_("ScoreMover"))
            mb.class_("ScoreMover").member_function('register_options').exclude()  # static function


    #if path == 'protocols/rna':
    #    mb.class_('RNA_SecStructInfo').member_function('clone').exclude()


    #if path == 'protocols/loops':
    #    if hfile == 'protocols/loops/LoopClosure.hh':
    #        removeProtectedConstructos(mb, "LoopClosure")

            #for c in mb.class_("LoopClosure").constructors(lambda x: x.access_type == 'protected'):
            #    c.exclude()


def make_nonvirtual(c):
  # isabs = False
  # for b in c.bases:
  #   b = b.related_class
  #   if not b.is_abstract and not b.is_wrapper_needed(): continue
  #   print c
  #   print b
  #   assert namespace_of(b) is namespace_of(c)
  #   make_nonvirtual(b)
  for mf in c.calldefs(allow_empty=True):
    # if mf.virtuality
    mf.virtuality = pyplusplus.decl_wrappers.declarations.VIRTUALITY_TYPES.NOT_VIRTUAL
  c._redefined_funcs = []
  c.is_abstract = c.__WAS_ABSTRACT
  c.__NEEDS_NO_INIT = True
  c.__IS_HACKED_NONVIRTUAL = True
  if c.is_abstract:
    c.constructors(allow_empty=True,recursive=False).exclude()
  if c.is_wrapper_needed():
    print "NONVIRTUAL'ized still needs wrapper",c
    print c.is_wrapper_needed()
    assert False



def removeProtectedConstructos(mb):
    try:
        for c in mb.classes():
            for f in c.constructors(lambda x: x.access_type == 'protected', allow_empty=True):
                f.exclude()
    except RuntimeError: pass

    #for f in mb.class_(class_name).constructors(lambda x: x.access_type == 'protected'):
        #f.exclude()

def exclude_private_and_protected_members(mb, ns):
    try:
        for c in ns.classes(allow_empty=True,recursive=False):
            for f in c.protected_members: f.exclude()
            for f in c.private_members: f.exclude()
    except RuntimeError: pass


def add_print_operators(mb, ns):
    for c in ns.classes(allow_empty=True,recursive=False):
        for o in ns.operators(allow_empty=True,recursive=False):
            if o.name == 'operator<<'  and  len(o.argument_types) == 2 and \
                    str(o.argument_types[0]) == 'std::ostream &' and \
                    '::' + str(o.argument_types[1]) ==  c.decl_string + ' const &':
                print 'Adding print operator for:', c.decl_string
                c.add_registration_code( 'def( bp::self_ns::str( bp::self ) )' )

            #print 'str(o.argument_types[1])', str(o.argument_types[1])

            #print o, dir(o)
            #print 'name:', o.name
            #print 'argument_types', dir(o.argument_types[0]), ' ~~~ ', o.argument_types[0], str(o.argument_types[0]) == 'std::ostream &'
            #print 'decl_string', o.argument_types[0].decl_string(), ' partial_decl_string', o.argument_types[0].partial_decl_string()
            #print 'arguments', type(o.arguments[0]), ' ~~~ ', o.arguments[0]




def annotate(mb):
    try:
        for c in mb.classes():
            c.__WAS_ABSTRACT = c.is_abstract
    except RuntimeError: pass



def getNameSpace(path, mb):
    nspath = path.replace("/","::")
    nspath = nspath.replace("::::","::")
    #if nspath in CONFIG.NAMESPACE_EXCLUDE: return
    print "MAKING PACKAGE:",nspath
    if nspath: ns = namespace_from_full_path(mb, nspath)
    else: ns = mb.global_ns
    assert ns

    return ns


def exclude(path, mb, hfile):
    print "Manual excluding:", path
    #mb_exclude(path, mb)
    #return

    ns = getNameSpace(path, mb)
    #ns = '::'.join( path.split('/') ) + ' [namespace]'
    #print '~~~~~~~ Ns_old=%s Ns=%s Path=%s' % (ns_old, ns, path)

    set_default_call_policies(mb)

    mb_exclude(path, mb, hfile)
    #removeProtectedConstructos(mb)
    exclude_private_and_protected_members(mb, ns)
    add_print_operators(mb, ns)

    annotate(mb)
    for c in ns.classes(allow_empty=True,recursive=False):
        if c.is_wrapper_needed() or c.is_abstract:
            make_nonvirtual(c)


    #namespace_from_full_path(mb,"utility").classes(recursive=0,allow_empty=True),
    #namespace_from_full_path(mb,"utility").free_funs(recursive=0,allow_empty=True),


    exclude_impls_(mb,ns)
    setup_held_type(mb,ns)
    remove_tracer_default_args_and_vars(mb, ns)

    #make_sure_is_copyable(mb,ns)

    mb._module_builder_t__registrations_code_tail = []
    ic = make_implicitly_convertible(mb,ns)
    if ic:
        print '\nMaking Implicitly convertible:'
        for i in ic: print i
    mb.add_registration_code( "".join(ic) )



def exclude_impls_(mb,ns):
  for c in ns.classes(allow_empty=True,recursive=False):
    if c.decl_string.endswith("_"):
      print "EXCLUDE IMPL_:",c.decl_string
      c.exclude()


def setup_held_type(mb,ns):
  for c in ns.classes(allow_empty=1,recursive=False):
    if is_refcount(c):
        if c.decl_string == u'::utility::pointer::ReferenceCount': continue  # we cant hold our self...
        #print "c.decl_string:", c.decl_string
        # if c.held_type:
        #   if c.held_type.count("utility::pointer::owning_ptr"): continue
        if c.is_abstract:
            print 'Abstarct:', c.decl_string
            if hasattr(c,'__IS_HACKED_NONVIRTUAL'):
                c.held_type = "utility::pointer::owning_ptr< "+c.decl_string+" >"
        elif len(c.is_wrapper_needed()) > 0:
            #print "len(c.is_wrapper_needed()) > 0: c.decl_string:", c.decl_string
            #pass #
            #c.held_type = "utility::pointer::owning_ptr< "+c.wrapper_alias+" >"
            print "~~ SET HELD TYPE OP",c , ' exportable?:', c.exportable
            c.held_type = "utility::pointer::owning_ptr< "+c.decl_string+" >"

        else:
            print "SET HELD TYPE OP",c
            c.held_type = "utility::pointer::owning_ptr< "+c.decl_string+" >"


def make_implicitly_convertible(mb,ns):
  ic = []
  for c in ns.classes(allow_empty=True):
    if not c.held_type: continue
    if not isinc(c): continue
    print 'c.decl_string:', c.decl_string
    if c.held_type.count("utility::pointer::owning_ptr") and c.decl_string != u'::utility::pointer::ReferenceCount':
      tplt  = "bp::implicitly_convertible< utility::pointer::owning_ptr< %s >\n"
      tplt += "                          , utility::pointer::owning_ptr< %s > >();\n"
      ic.append(tplt % ( c.decl_string,c.decl_string+" const") )
      for b in get_class_bases(c):
          print 'c, b:', c, b
          if b is c: continue
          #if not isinc(b):
          #    print 'isinc...'
          #    continue
          if b.decl_string == u'::utility::pointer::ReferenceCount': continue
          ic.append(tplt % ( c.decl_string,b.decl_string) )
          #ic.append(tplt % ( c.decl_string,b.decl_string+" const") )
  return ic

def make_sure_is_copyable(mb,ns):
  for c in MAKE_SURE_IS_COPYABLE:
    c = class_from_full_path(mb,c)
    if c in ns.classes(allow_empty=True,recursive=True):
      c.noncopyable = False


def set_default_call_policies(mb):
  print 'set_default_call_policies'
  for f in mb.calldefs().declarations:
    f.create_with_signature = True
    # if f.call_policies is None:
    if not hasattr(f,'return_type') or not f.return_type: continue
    rt = f.return_type
    rtdecl = rt
    if hasattr(rtdecl,'base'):
      b = rt
      while hasattr(b,'base'): b = b.base
      if hasattr(b,'declaration'): rtdecl = b.declaration
    rtbase = rt.decl_string.strip().rstrip("&").rstrip("*").strip()
    if rtbase.endswith("const"): rtbase = rtbase[:-5].strip()
    # if f.name == "default_value":
    #   print "\nDEFAULT_VALUE",rtbase
    #   print rt.decl_string
    if rt.decl_string.endswith('&') and rtbase.split("::")[-1] in  ("Real","int","float","Length","long","unsigned","Size","unsigned int","long unsigned int","size_t","string","Energy","double","char","DistanceSquared","Distance","bool","TorsionType","PackerEnergy"):
      if rt.decl_string.count("const"):
        # if f.name == "default_value": print "const"
        f.call_policies = call_policies.return_value_policy( call_policies.copy_const_reference )
      else:
        # if f.name == "default_value": print "non-const"
        f.call_policies = call_policies.return_value_policy( call_policies.copy_non_const_reference )
    elif rt.decl_string.strip().endswith("&") or rt.decl_string.strip().endswith("*"):#if f.call_policies is None:
      if isinstance(rtdecl,pyplusplus.decl_wrappers.enumeration_wrapper.enumeration_t):
        if rt.decl_string.count('const'):
          f.call_policies = call_policies.return_value_policy( call_policies.copy_const_reference )
        else:
          f.call_policies = call_policies.return_value_policy( call_policies.copy_non_const_reference )
      else:
        # if f.name == "default_value": print "ref"
        f.call_policies = call_policies.return_value_policy( call_policies.reference_existing_object )


def exclude_private_and_protected_members(mb,ns):
  for c in ns.classes(allow_empty=True,recursive=False):
    for f in c.protected_members: f.exclude()
    for f in c.private_members: f.exclude()

def exclude_disallowed_virtuals(mb,ns):
  nowraper = map(lambda x: class_from_full_path(mb,x),CONFIG.DISALLOW_WRAPPERS)
  for c in ns.classes(allow_empty=True,recursive=False):
    if not c in nowraper: continue
    print "\nEXCLUDE ABSTRACT!!!!!",c,'\n'
    if c.is_abstract: c.exclude()
    # for mf in c.calldefs(allow_empty=True,recursive=False):
    #   if mf.virtuality == pyplusplus.decl_wrappers.declarations.VIRTUALITY_TYPES.VIRTUAL:
    #     # print "OVERRIDE VIRTUALITY:",mf
    #     mf.virtuality = pyplusplus.decl_wrappers.declarations.VIRTUALITY_TYPES.NOT_VIRTUAL
    # if c.is_wrapper_needed:
    #   c.is_wrapper_needed = []

def make_nonvirtual(c):
  # isabs = False
  # for b in c.bases:
  #   b = b.related_class
  #   if not b.is_abstract and not b.is_wrapper_needed(): continue
  #   print c
  #   print b
  #   assert namespace_of(b) is namespace_of(c)
  #   make_nonvirtual(b)
  for mf in c.calldefs(allow_empty=True):
    # if mf.virtuality
    mf.virtuality = pyplusplus.decl_wrappers.declarations.VIRTUALITY_TYPES.NOT_VIRTUAL
  c._redefined_funcs = []
  c.is_abstract = c.__WAS_ABSTRACT
  c.__NEEDS_NO_INIT = True
  c.__IS_HACKED_NONVIRTUAL = True
  if c.is_abstract:
    c.constructors(allow_empty=True,recursive=False).exclude()
  if c.is_wrapper_needed():
    print "NONVIRTUAL'ized still needs wrapper",c
    print c.is_wrapper_needed()
    #assert False
    # ^^^ nah, thats ok


def parse_full_name(full_name):
  pth = []#filter(len,full_class_name.split("::"))
  regex = re.compile("^([^<>, :]*)::")
  while regex.match(full_name):
    ns = regex.search(full_name).group(1)    # print ns
    if ns:
      pth.append( ns )
    full_name = regex.sub("",full_name)    # print full_name
  pth.append(full_name)
  pth = [x for x in pth if x]
  return pth

def class_from_full_path(mb,fullname):
  path = parse_full_name(fullname)
  nspath,cname = path[:-1],path[-1]
  ns = mb.global_ns
  for nsname in nspath:
    for tmp in ns.namespaces(nsname):
      if tmp.parent is ns:
        ns = tmp
  return ns.class_(cname)

def namespace_from_full_path(mb,fullname):
  nspath = parse_full_name(fullname)
  ns = mb.global_ns
  for nsname in nspath:
      #print 'namespace_from_full_path: ', nsname, ns
      ns = ns.namespace(nsname,recursive=False)
  return ns

def remove_tracer_default_args_and_vars(mb, ns):
  for cd in ns.calldefs(allow_empty=True):
    for arg in cd.arguments:
      if basetype(arg) == "basic::Tracer":
        if arg.default_value:
          print "REMOVING TRACER DEFAULT ARG",arg
          arg.default_value = None
  for fv in ns.vars(recursive=False,allow_empty=True):
    if basetype(fv) == "basic::Tracer":
      fv.exclude()
      print "REMOVING TRACER FREE VAR",fv,`fv`
      print "        ",fv.ignore

def basetype(t):
  """docstring for basetype"""
  # typestr = typestr.split("::")[-1].strip()
  # return re.sub("(?:const|[&]|<.*?>)","",typestr).strip()
  if hasattr(t,'type'): t = t.type
  while hasattr(t,"base"):
    t = t.base
  return str(t)


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
                try:
                    iname = i[i.index('<')+1 : i.index('>')]
                    if os.path.isfile( iname.replace('.fwd', '') ): i = i.replace('.fwd', '')
                    output.add(i)
                except ValueError: pass

    output = '\n'.join(output)


    lines = list( set( output.split('\n') ) );  lines.sort()
    output = '\n'.join( lines ) + '\n'
    return output


def finalize2(fname, dest, path, module_name='_noname', add_by_hand=False, includes=''):

    f = open(fname);  s = '#include <boost/python.hpp>\n' + includes + f.read();  f.close()

    namespace = os.path.basename(path)
    by_hand = dest + '/../src/' + path + '/_' + namespace + '__by_hand.cc'
    #print '!++ ', by_hand, add_by_hand, 'namespace=%s' % namespace, 'path=%s' % path
    if os.path.isfile(by_hand) and add_by_hand:
        print 'Adding by bindings writen by hand...'
        f = open(by_hand);
        s = f.read() + '\n\n' + s
        f.close()

        #s = re.sub('BOOST_PYTHON_MODULE\(.*\)\{',
        #           'BOOST_PYTHON_MODULE(_%(n)s){\n  wrap__%(n)s__by_hand();' % dict(n=namespace),
        #           s)

        if add_by_hand:
            s = re.sub('BOOST_PYTHON_MODULE\(.*\)\{',
                   'BOOST_PYTHON_MODULE(%(module_name)s){\n  wrap__%(n)s__by_hand();' % dict(n=namespace, module_name=module_name),
                   s)
        else:
            s = re.sub('BOOST_PYTHON_MODULE\(.*\)\{',
                   'BOOST_PYTHON_MODULE(%(module_name)s){' % dict(n=namespace, module_name=module_name),
                   s)


    f = open(fname,'w')
    f.write(s)
    f.close()




def finalize(fname, dest, path, mb, module_name='_noname', add_by_hand=False, files=[], add_includes=True):
    files = files[:]
    #print 'finalize... add_by_hand=%s' % add_by_hand, files
    #if len(files) != 1:  # Build mode 'All'
    if False:
        output = commands.getoutput("cd %s && cat *.cc *.hh | grep '#include <core'" % path)
        '''
        output += '\n' + commands.getoutput("cd %s && cat *.cc *.hh | grep '#include <protocols'" % path)
        output += '\n' + commands.getoutput("cd %s && cat *.cc *.hh | grep '#include <numeric'" % path)
        output += '\n' + commands.getoutput("cd %s && cat *.cc *.hh | grep '#include <utility'" % path)
        '''
        output += '\n' + commands.getoutput("cd %s && cat *.cc *.hh | grep '#include <protocols\|#include <numeric\|#include <utility' | grep '\.hh'" % path)

        # Some exceptional situation where we want to add additional includes
        if fname.find('AmbiguousMultiConstraint.cc') > 0: output += '\n#include <core/conformation/Conformation.hh>\n'
        if fname.find('RotamerTrieBase.cc') > 0:
            output += '\n#include <core/scoring/etable/EtableEnergy.hh>\n'
            output += '\n#include <core/scoring/etable/CoarseEtableEnergy.hh>\n'
            output += '\n#include <core/scoring/hbonds/HBondEnergy.hh>\n'
        if fname.find('TrieCountPairBase.cc') > 0:
            output += '\n#include <core/scoring/etable/EtableEnergy.hh>\n'
            output += '\n#include <core/scoring/etable/CoarseEtableEnergy.hh>\n'
            output += '\n#include <core/scoring/hbonds/HBondEnergy.hh>\n'
        if fname.find('EnergyMethod.cc') > 0:
            output += '\n#include <core/optimization/MinimizerMap.hh>\n'
            output += '\n#include <core/pack/task/PackerTask.hh>\n'
        if fname.find('RotamerTrieBase.cc') > 0:
            output += '\n#include <core/scoring/elec/FA_ElecEnergy.hh>\n'
        if fname.find('ScoreFunction.cc') > 0:
            output += '\n#include <core/optimization/MinimizerMap.hh>\n'
        #if fname.find('DNA_BaseEnergy.cc') > 0:
        #    output += '\n#include <core/pack/task/PackerTask.hh>\n'
        #if fname.find('GenBornEnergy.cc') > 0:
        #    output += '\n#include <core/pack/task/PackerTask.hh>\n'
        if fname.find('core/scoring/methods') > 0:
            output += '\n#include <core/pack/task/PackerTask.hh>\n'
            output += '\n#include <core/optimization/MinimizerMap.hh>\n'

        if fname.find('core/scoring/hbonds') > 0:
            output += '\n#include <core/pack/task/PackerTask.hh>\n'

        if fname.find('TrieCountPairBase.cc') > 0:
            output += '\n#include <core/scoring/elec/FA_ElecEnergy.hh>\n'

        if fname.find('FlexbbIGFactory.cc') > 0:
            output += '\n#include <core/pack/task/PackerTask.hh>\n'


        #if fname.find('AbrelaxApplication.hh') > 0:
        #    output += '\n#include <core/kinematics/MoveMap.hh>\n'

    output = ''
    if add_includes:

        if True:
            output = set()
            for f in files[:]:
                cc_file = f.replace('.hh', '.cc')
                if os.path.isfile(cc_file):
                    files.append(cc_file)

            for f in files:
                #print f
                lines = commands.getoutput(  "cat %s | grep '^#include <core' | grep '\.hh\|\.hpp'" % f) + '\n'        \
                        + commands.getoutput("cat %s | grep '^#include <numeric' | grep '\.hh\|\.hpp'" % f) + '\n'   \
                        + commands.getoutput("cat %s | grep '^#include <utility' | grep '\.hh\|\.hpp'" % f) + '\n'   \
                        + commands.getoutput("cat %s | grep '^#include <protocols' | grep '\.hh\|\.hpp'" % f) + '\n' \
                        + commands.getoutput("cat %s | grep '^#include <basic' | grep '\.hh\|\.hpp'" % f) + '\n' \
                        + commands.getoutput("cat %s | grep '^#include <ObjexxFCL' | grep '\.hh\|\.hpp'" % f) + '\n'
                        #+ commands.getoutput("cat %s | grep '^#include <boost' | grep '\.hh\|\.hpp'" % f) + '\n'

                #print f, lines

                for i in lines.split('\n'):
                    if i:
                        #print 'i=', i
                        iname = i[i.index('<')+1 : i.index('>')]
                        if os.path.isfile( iname.replace('.fwd', '') ): i = i.replace('.fwd', '')
                        output.add(i)

            output = '\n'.join(output)
            #print output
            #sys.exit(1)


            #if fname.find('ShortLoopClosure.cc') > 0:
            #    output += '\n#include <protocols/abinitio/FragmentMover.hh>\n'

            if fname.find('IterativeFullatom.cc') > 0:
                output += '\n#include <protocols/abinitio/PairingStatistics.hh>\n'

            if fname.find('FASelectSlidingWindowLoopClosure.cc') > 0:
                output += '\n#include <protocols/evaluation/ConstraintEvaluator.hh>\n'


        else:  # Build mode 'one file'
            print "finalize... Build mode 'one file'..."
            f = files[0].replace('.hh', '.cc')
            if os.path.isfile(f):
                output = commands.getoutput("cat %s | grep '#include <core'" % f)
                output += '\n' + commands.getoutput("cat %s | grep '#include <protocols'" % f)
                output += '\n' + commands.getoutput("cat %s | grep '#include <numeric'" % f)
            else: output = ''

        lines = list( set( output.split('\n') ) );  lines.sort()
        output = '\n'.join( lines ) + '\n'

        #output = re.sub('\.fwd\.', '.', output)

        # removing non existen header files
        #lines = output.split('\n')
        #while

    f = open(fname);  s = '#include <boost/python.hpp>\n' + output + f.read();  f.close()

    namespace = os.path.basename(path)
    by_hand = dest + '/../src/' + path + '/_' + namespace + '__by_hand.cc'
    #print '!++ ', by_hand, add_by_hand, 'namespace=%s' % namespace, 'path=%s' % path
    if os.path.isfile(by_hand) and add_by_hand:
        print 'Adding by bindings writen by hand...'
        f = open(by_hand);
        s = f.read() + '\n\n' + s
        f.close()

        #s = re.sub('BOOST_PYTHON_MODULE\(.*\)\{',
        #           'BOOST_PYTHON_MODULE(_%(n)s){\n  wrap__%(n)s__by_hand();' % dict(n=namespace),
        #           s)

        if add_by_hand:
            s = re.sub('BOOST_PYTHON_MODULE\(.*\)\{',
                   'BOOST_PYTHON_MODULE(%(module_name)s){\n  wrap__%(n)s__by_hand();' % dict(n=namespace, module_name=module_name),
                   s)
        else:
            s = re.sub('BOOST_PYTHON_MODULE\(.*\)\{',
                   'BOOST_PYTHON_MODULE(%(module_name)s){' % dict(n=namespace, module_name=module_name),
                   s)


    f = open(fname,'w')
    f.write(s)
    f.close()

    #print 'Creating __init__.py file...'
    #f = file( dest + '/' + path + '/__init__.py', 'a');
    #f.write('from %s import *\n' % os.path.basename(fname)[:-3]);
    #f.close()





def finalize_old(fname, path, mb):
    ns = getNameSpace(path, mb)

    print "FINALIZE"
    f = open(fname)
    s = f.read()
    f.close()
    f = open(fname,'w')
    s = re.sub("BOOST_PYTHON_MODULE\(.*\)","void wrap"+ns.decl_string.replace(":","_")+"()",s)
    # vec_unsigned    = "::utility::vector1<unsigned, std::allocator<unsigned> >"
    # vecvec_unsigned  = "::utility::vector1<utility::vector1<unsigned, std::allocator<unsigned> >, std::allocator<utility::vector1<unsigned, std::allocator<unsigned> > > >"
    # vec_size    = "::utility::vector1<Size>"
    # vecvec_size = "::utility::vector1<utility::vector1<Size> >"
    # if s.count(vec_unsigned) or s.count(vecvec_unsigned):
    #   print "REPLACING VEC<UNSIGNED> w/ VEC<SIZE_T>!!!!!!!!!!"
    #   s = s.replace(vec_unsigned,vec_size)
    #   s = s.replace(vecvec_unsigned,vecvec_size)
    # # n = fname.split("/")[-1][1:].replace(".cc","")
    # # s = s.replace('#include "all_headers.hh"','#include "include/all_forwards.hh"\n#include "include/%s.hh"'%n)
    incname = ns.decl_string.strip(":").replace("::","__")+".hh"
    # s = s.replace( '#include "include/%s/groups/all.hh"'%CONFIG.SVN_REV
    #              , '#include "include/%s/groups/'%CONFIG.SVN_REV+incname+'"')
    s = s.replace("protocols::moves::Real","core::Real")
    #
    pat1 = r"\(\(const std::allocator<([^<>]+)>\&\)\(\(const std::allocator<\1>\*\)\(\& std::allocator<\1>\(\)\)\)\)"
    pat2 = r" /* REMOVED ((const std::allocator<\1>&)((const std::allocator<\1>*)(& std::allocator<\1>())))*/ "
    s = re.sub(pat1,pat2,s)
    #
    f.write(s)
    f.close()
    #os.system('find %s -name \*~ | xargs rm'%ROOT_DIR)
    #os.system("find %s -name named_tuple.py | xargs rm"%ROOT_DIR)


def filter_code_creators(mb, path, wrappers):
  ns = getNameSpace(path, mb)

  print 'getting all creators...'
  creators = pyplusplus.code_creators.make_flatten(mb.code_creator.creators)
  print(len(creators))
  creators = filter( lambda x: hasattr(x,'declaration'), creators )
  print(len(creators))
  for cc in creators:
    nstmp = namespace_of(cc)
    if nstmp is not ns:
      print 'removing CC for',cc.declaration
      cc.parent.creators.remove(cc)
    if not wrappers:
      if str(type(cc)).count("_wrapper"):
        print cc
        print cc.declaration
        assert False



def isinc(x):
  return not x._ignore

def is_refcount(c):
  return "utility::pointer::ReferenceCount [class]" in [str(b) for b in get_class_bases(c)]

def postproc(code):
  code = re.sub(r'<unsigned([,>])',r'<std::size_t\1',code)
  return code

def ischildof(a,b):
  while a and a is not b: a = a.parent
  return a is b


def namespace_of(x):
  # t = time.time()
  if isinstance(x,pyplusplus.code_creators.code_creator_t): x = x.declaration
  # print "FIND NS:",time.time() - t
  while x:
    # t = time.time()
    if isinstance(x,pyplusplus.decl_wrappers.namespace_t): break
    # print 'isinstance',time.time()-t
    t = time.time()
    x = x.parent
    # print 'parent',time.time()-t
    # print "   ",x
  return x

def class_of(x):
  # t = time.time()
  if isinstance(x,pyplusplus.code_creators.code_creator_t): x = x.declaration
  # print "FIND NS:",time.time() - t
  while x:
    # t = time.time()
    if isinstance(x,pyplusplus.decl_wrappers.class_t): break
    # print 'isinstance',time.time()-t
    t = time.time()
    x = x.parent
    # print 'parent',time.time()-t
    # print "   ",x
  return x


def get_class_bases(argcls,baseslist=None):
  baseslist = []
  if hasattr(argcls,'bases'):
    for bi in argcls.bases:
      assert hasattr(bi,"related_class")
      baseslist.extend( get_class_bases(bi.related_class,baseslist) )
  if baseslist is None:
    baseslist = []
  baseslist.append(argcls)
  return baseslist
