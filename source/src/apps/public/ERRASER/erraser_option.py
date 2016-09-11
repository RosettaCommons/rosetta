# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

class erraser_option :

    def __init__( self ) :
        #General options
        self.input_pdb = ''
        self.map_file = ''
        self.out_pdb = ''
        self.out_prefix = ''
        self.map_reso = 2.5
        self.rosetta_folder = ''
        self.debug = False
        self.verbose = False
        self.kept_temp_folder = False
        self.use_existing_temp_folder = True
        self.new_torsional_potential = True
        self.corrected_geo = True
        self.rna_prot_erraser = False
        self.rosetta_folder = ""
        self.rosetta_bin = ""
        self.rosetta_database = ""
        self.log_out = ""
        self.log_err = ""
        self.nproc = 0
        self.multiproc_minimize = False
        self.guarantee_no_DNA = False

        #erraser options
        self.n_iterate = 1
        self.rebuild_rmsd = False
        self.rebuild_all = True
        self.fixed_res = []
        self.cutpoint = []
        self.extra_res = []
        self.fixed_res_rs = []
        self.cutpoint_rs = []
        self.rebuild_res_pdb = ""
        self.skip_single_res_minimize = False

        #Minimize
        self.vary_geometry = True
        self.skip_minimize = False
        self.attempt_pyrimidine_flip = False
        self.constrain_phosphate = True
        self.res_slice = []

        #SWA rebuilding
        self.native_screen_RMSD = 3.0
        self.use_native_edensity_cutoff = False
        self.native_edensity_cutoff = 0.8
        self.ideal_geometry = True
        self.allow_syn_pyrimidine = True
        self.include_native = False
        self.slice_nearby = True
        self.rebuild_res = 0
        self.scoring_file = ''
        self.cluster_RMSD = 0.1
        self.is_append = True
        self.native_screen_RMSD = 3.0
        self.num_pose_kept = 100
        self.num_pose_kept_cluster = 10
        self.rebuild_res_list = []
        self.search_syn_pyrimidine_only_when_native_syn =True
        self.constrain_chi = True
        
        #2015 fixes
        self.o2prime_legacy_mode = True
        self.use_2prime_OH_potential = False
        self.fcc2012_new_torsional_potential = True
        self.fcc2012_scoring_file = True
        
        self.enlarge_H_lj = False

    def read_cmdline_full( self, argv ) :
        #General options
        self.input_pdb = parse_options( argv, 'pdb', '' )
        self.map_file = parse_options( argv, 'map', '' )
        self.out_pdb = parse_options( argv, 'out_pdb', '' )
        self.out_prefix = parse_options( argv, 'out_prefix', '' )
        self.map_reso = parse_options( argv, 'map_reso', 2.5 )
        self.debug = parse_options( argv, "debug", "False" )
        self.verbose = parse_options( argv, "verbose", "False" )
        self.kept_temp_folder = parse_options( argv, "kept_temp_folder", "False" )
        self.use_existing_temp_folder = parse_options( argv, "use_existing_temp_folder", "True" )
        self.new_torsional_potential = parse_options( argv, "new_torsional_potential", "True" )
        self.corrected_geo = parse_options( argv, "corrected_geo", "True" )
        self.rna_prot_erraser = parse_options( argv, "rna_prot_erraser", "False" )
        self.rosetta_folder = parse_options( argv, 'rosetta_folder', '')
        self.rosetta_bin = parse_options( argv, 'rosetta_bin', '')
        self.rosetta_database = parse_options( argv, 'rosetta_database', '')
        self.nproc = parse_options( argv, 'nproc', 0 )
        self.multiproc_minimize = parse_options( argv, 'multiproc_minimize', "False" )
        self.guarantee_no_DNA = parse_options( argv, 'guarantee_no_DNA', "False" )


        #erraser options
        self.n_iterate = parse_options( argv, 'n_iterate', 1 )
        self.rebuild_rmsd = parse_options( argv, "rebuild_rmsd", "True" )
        self.rebuild_all = parse_options( argv, "rebuild_all", "False" )
        self.fixed_res = parse_option_chain_res_list ( argv, 'fixed_res' )
        self.cutpoint = parse_option_chain_res_list ( argv, 'cutpoint' )
        self.extra_res = parse_option_chain_res_list ( argv, 'rebuild_extra_res' )
        self.fixed_res_rs = parse_option_int_list ( argv, 'fixed_res_rs' )
        self.cutpoint_rs = parse_option_int_list ( argv, 'cutpoint_rs' )
        self.rebuild_res_pdb = parse_options( argv, "rebuild_res_pdb", '' )

        #Minimize
        self.vary_geometry = parse_options( argv, "vary_geometry", "True" )
        self.skip_minimize = parse_options( argv, "skip_minimize", "False" )
        self.attempt_pyrimidine_flip = parse_options( argv, "attempt_pyrimidine_flip", "False" )
        self.constrain_phosphate = parse_options( argv, "constrain_phosphate", "True" )
        self.res_slice = parse_option_int_list ( argv, 'res_slice' )

        #SWA rebuilding
        self.native_screen_RMSD = parse_options( argv, "native_screen_RMSD", 3.0 )
        self.native_edensity_cutoff = parse_options(argv, "native_edensity_cutoff", 0.8)
        self.use_native_edensity_cutoff = parse_options(argv, "use_native_edensity_cutoff", "False")
        self.ideal_geometry =  parse_options( argv, "ideal_geometry", "True" )
        self.include_native =  parse_options( argv, "include_native", "False" )
        self.allow_syn_pyrimidine =  parse_options( argv, "allow_syn_pyrimidine", "True" )
        self.search_syn_pyrimidine_only_when_native_syn = (
            parse_options( argv, "search_syn_pyrimidine_only_when_native_syn", "True" ) )
        self.slice_nearby =  parse_options( argv, "slice_nearby", "True" )
        self.rebuild_res = parse_options( argv, "rebuild_res", 0 )
        self.scoring_file = parse_options( argv, "scoring_file", "" )
        self.cluster_RMSD = parse_options( argv, "cluster_RMSD", 0.1 )
        self.is_append = parse_options( argv, "is_append", "True" )
        self.constrain_chi = parse_options( argv, "constrain_chi", "True" )
        self.native_screen_RMSD = parse_options( argv, "native_screen_RMSD", 3.0 )
        self.num_pose_kept =  parse_options( argv, "num_pose_kept", 100 )
        self.num_pose_kept_cluster =  parse_options( argv, "num_pose_kept_cluster", 10 )
        self.rebuild_res_list = parse_option_int_list ( argv, 'rebuild_res_list' )

        #2015 fixes
        self.o2prime_legacy_mode = parse_options( argv, "o2prime_legacy_mode", "True" )
        self.use_2prime_OH_potential = parse_options( argv, "use_2prime_OH_potential", "False" )
        self.fcc2012_new_torsional_potential = parse_options( argv, "fcc2012_new_torsional_potential", "True" )
        self.fcc2012_scoring_file = parse_options( argv, "fcc2012_scoring_file", "True" )
        
        self.enlarge_H_lj = parse_options( argv, "enlarge_H_lj", "True" )

        self.finalize()

    def read_cmdline_erraser( self, argv ) :
        #General options
        self.input_pdb = parse_options( argv, 'pdb', '' )
        self.map_file = parse_options( argv, 'map', '' )
        self.out_pdb = parse_options( argv, 'out_pdb', '' )
        self.map_reso = parse_options( argv, 'map_reso', 2.5 )
        self.rna_prot_erraser = parse_options( argv, "rna_prot_erraser", "False" )
        self.debug = parse_options( argv, "debug", "False" )
        self.rosetta_folder = parse_options( argv, 'rosetta_folder', '')
        self.rosetta_bin = parse_options( argv, 'rosetta_bin', '')
        self.rosetta_database = parse_options( argv, 'rosetta_database', '')
        self.nproc = parse_options( argv, 'nproc', 0 )
        self.multiproc_minimize = parse_options( argv, 'multiproc_minimize', "False" )
        self.guarantee_no_DNA = parse_options( argv, 'guarantee_no_DNA', "False" )

        #erraser_single_res options
        self.n_iterate = parse_options( argv, 'n_iterate', 1 )
        self.rebuild_rmsd = parse_options( argv, "rebuild_rmsd", "True" )
        self.rebuild_all = parse_options( argv, "rebuild_all", "False" )
        self.fixed_res = parse_option_chain_res_list ( argv, 'fixed_res' )
        self.extra_res = parse_option_chain_res_list ( argv, 'rebuild_extra_res' )

        #SWA rebuilding
        self.use_native_edensity_cutoff = parse_options(argv, "use_native_edensity_cutoff", "False")
        self.constrain_chi = parse_options( argv, "constrain_chi", "True" )
        self.native_screen_RMSD = parse_options( argv, "native_screen_RMSD", 3.0 )
        self.search_syn_pyrimidine_only_when_native_syn = (
                parse_options( argv, "search_syn_pyrimidine_only_when_native_syn", "True" ) )
        self.num_pose_kept =  parse_options( argv, "num_pose_kept", 100 )
        self.num_pose_kept_cluster =  parse_options( argv, "num_pose_kept_cluster", 10 )
        
        #2015 fixes
        self.o2prime_legacy_mode = parse_options( argv, "o2prime_legacy_mode", "True" )
        self.scoring_file = parse_options( argv, "scoring_file", "" )
        self.use_2prime_OH_potential = parse_options( argv, "use_2prime_OH_potential", "False" )
        self.fcc2012_new_torsional_potential = parse_options( argv, "fcc2012_new_torsional_potential", "True" )
        self.fcc2012_scoring_file = parse_options( argv, "fcc2012_scoring_file", "True" )

        self.enlarge_H_lj = parse_options(argv, "enlarge_H_lj", "True")

        self.finalize()

    def read_cmdline_erraser_single_res( self, argv ) :
        #General options
        self.input_pdb = parse_options( argv, 'pdb', '' )
        self.map_file = parse_options( argv, 'map', '' )
        self.out_prefix = parse_options( argv, 'out_prefix', '' )
        self.map_reso = parse_options( argv, 'map_reso', 2.5 )
        self.rna_prot_erraser = parse_options( argv, "rna_prot_erraser", "False" )
        self.debug = parse_options( argv, "debug", "False" )
        self.use_existing_temp_folder = parse_options( argv, "use_existing_temp_folder", "True" )
        self.rosetta_folder = parse_options( argv, 'rosetta_folder', '')
        self.rosetta_bin = parse_options( argv, 'rosetta_bin', '')
        self.rosetta_database = parse_options( argv, 'rosetta_database', '')
        self.nproc = parse_options( argv, 'nproc', 0 )
        self.multiproc_minimize = parse_options( argv, 'multiproc_minimize', "False" )

        #erraser options
        self.rebuild_res_pdb = parse_options( argv, "rebuild_res", '' )
        self.skip_single_res_minimize = parse_options( argv, "skip_single_res_minimize", "False" )

        #SWA rebuilding
        self.use_native_edensity_cutoff = parse_options(argv, "use_native_edensity_cutoff", "False")
        self.constrain_chi = parse_options( argv, "constrain_chi", "True" )
        self.native_screen_RMSD = parse_options( argv, "native_screen_RMSD", 3.0 )
        self.search_syn_pyrimidine_only_when_native_syn = (
                parse_options( argv, "search_syn_pyrimidine_only_when_native_syn", "True" ) )
        self.num_pose_kept =  parse_options( argv, "num_pose_kept", 100 )
        self.num_pose_kept_cluster =  parse_options( argv, "num_pose_kept_cluster", 10 )
        
        #2015 fixes
        self.o2prime_legacy_mode = parse_options( argv, "o2prime_legacy_mode", "True" )
        self.scoring_file = parse_options( argv, "scoring_file", "" )
        self.use_2prime_OH_potential = parse_options( argv, "use_2prime_OH_potential", "False" )
        self.fcc2012_new_torsional_potential = parse_options( argv, "fcc2012_new_torsional_potential", "True" )
        self.fcc2012_scoring_file = parse_options( argv, "fcc2012_scoring_file", "True" )

        self.enlarge_H_lj = parse_options(argv, "enlarge_H_lj", "True")

        self.finalize()

    def finalize( self ) :
        if self.input_pdb == '' :
            error_exit( 'input_pdb not specified for erraser_option!!!' )
        check_path_exist( self.input_pdb )
        self.input_pdb = abspath( self.input_pdb )

        if self.out_prefix == "" :
            self.out_prefix = basename( self.input_pdb ).replace('.pdb', '')

        if self.out_pdb == "" :
            self.out_pdb = basename( self.input_pdb ).replace('.pdb', '_erraser.pdb')
        self.out_pdb = abspath( self.out_pdb )
        if exists(self.out_pdb) :
            print "Output pdb file %s exists... Remove it..." % self.out_pdb
            remove(self.out_pdb)

        if self.map_file != "" :
            check_path_exist( self.map_file )
            self.map_file = abspath( self.map_file )

        if self.debug :
            self.verbose = True
            self.kept_temp_folder = True
            if self.num_pose_kept_cluster < 10 :
                self.num_pose_kept_cluster = 10

        if not self.use_native_edensity_cutoff :
            self.native_edensity_cutoff = -1

        if self.scoring_file == "" :
            if self.map_file == "" :
                self.scoring_file = "stepwise/rna/rna_loop_hires_04092010 "
            elif self.fcc2012_scoring_file is True:
                self.scoring_file = "stepwise/rna/rna_hires_elec_dens_FCC2012 "
            else :
                self.scoring_file = "stepwise/rna/rna_hires_elec_dens "

        if self.rosetta_bin != '' and self.rosetta_folder != '':
            error_exit("The options 'rosetta_folder' and 'rosetta_bin' are incompatible! Specify only one of them.")
        if self.rosetta_database != '' and self.rosetta_folder != '':
            error_exit("The options 'rosetta_folder' and 'rosetta_database' are incompatible! Specify only one of them.")
        if ( ( self.rosetta_database != '' and self.rosetta_bin == '' ) or
             ( self.rosetta_database == '' and self.rosetta_bin != '' ) ):
            error_exit("Both options 'rosetta_bin' and 'rosetta_database' need to be specified!")

        if self.rosetta_database != '' :
            self.rosetta_bin = self.rosetta_folder
            self.rosetta_database = self.rosetta_folder

        
        max_cpu = multiprocessing.cpu_count()
        if self.nproc > max_cpu:
            print "Number of processes requested exceeds CPU count (nproc=%d, max_cpu=%d) ... setting nproc to %d" % (self.nproc, max_cpu, max_cpu)
            self.nproc = max_cpu
