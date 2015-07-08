// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author ora


// Rosetta Headers
#include "options.h"
#include "after_opts.h"
#include "analyze_interface_ddg_ns.h"
#include "design.h"
#include "disulfides.h"
#include "dna.h"
#include "dna_ns.h"
#include "dna_motifs.h"
#include "dock_fab.h"
#include "dock_pivot_ns.h"
#include "docking_movement.h"
#include "docking_score.h"
#include "docking_ns.h"
#include "domins_score.h"
#include "favor_residue_ns.h"
#include "enzyme.h"
#include "enzyme_ns.h"
#include "files_paths.h"
#include "filters.h"
#include "fullatom.h"
#include "fullatom_setup.h"
#include "geometric_solvation_ns.h"
#include "hbonds.h" // jss set_fine_hb_categories()
#include "InteractionGraphBase.h"
#include "InteractionGraphSupport.h"
#include "ligand.h"
#include "loops.h"
#include "loops_ns.h"
#include "misc.h"
#include "MultiCoolAnnealer.h"
#include "namespace_assemble_options.h"
#include "namespace_options.h"
#include "param.h"
#include "param_pack.h"
#include "param_rotamer_trie.h"
#include "pdbstats_ns.h"
#include "pH_main.h"
#include "pKa_mode.h"
#include "pose.h"
#include "pose_ligand.h"
#include "read_aaproperties.h"
#include "relax.h"
#include "relax_structure.h"
#include "rotamer_trial_energies.h"
#include "rotamer_trie_calc_energies.h"
#include "runlevel.h"
#include "SimAnnealerBase.h"
#include "score.h"
#include "score_ns.h"
#include "symmetric_design.h"
#include "taboo_search.h"
#include "termini.h"
#include "version_rosetta.h"
#include "water.h"
#include "water_ns.h"
#include "weights_manager.h"
#include "vdw.h"
#include "add_pser.h"//KMa phospho_ser(SEP) 2006-01
#include "SurfaceMode.h"//KMa surface 2006-02

// ObjexxFCL Headers
#include <ObjexxFCL/DimensionExpressions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <algorithm>
#include <cstdlib>
#include <iostream>

//Utility Headers
#include <utility/basic_sys_util.hh>


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///car set command line and mode-specific options
///
/// @details
///
/// @param  mode - [in/out]? -
/// @param  number_of_output - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_rosetta_options(
	std::string & mode,
	int & number_of_output
)
{
	using namespace design_sym;
	using namespace param;

	if ( truefalseoption("version") ) write_version();
//car define mode
	mode = "abinitio";
	if ( truefalseoption("score") ) {
		mode = "score";
		get_score_options();
	} else if ( truefalseoption("abinitio") ) {
//js placed here to see if loops management can be ported to abinitio mode
		mode = "abinitio";
	} else if ( truefalseoption("refine") ) {
		mode = "refine";
		get_refine_options();
	} else if ( truefalseoption("assemble") ) {
		mode = "assemble";
		get_assemble_options();
	} else if ( truefalseoption("idealize") ) {
		mode = "idealize";
		get_idealize_options();
	} else if ( truefalseoption("relax") ) {
		mode = "relax";
		get_relax_options();
//dk abrelax_mode was used in boinc but is now deprecated
	} else if ( truefalseoption("abrelax") || truefalseoption("abrelax_mode") ) {
		mode = "abrelax";
		get_abrelax_options();
	} else if ( truefalseoption("design") ) {
		mode = "design";
		get_design_options();
	} else if ( truefalseoption("dock") ) {
		mode = "dock";
		get_dockmode_options();
	} else if ( truefalseoption("flexpepdock") ) {
		mode = "flexpepdock";
		get_flexpepdock_options();
	} else if (truefalseoption("pathways") ) {
		mode = "pathways";
		get_pathways_options();
	} else if ( truefalseoption("membrane") ) {
		mode = "membrane";
		get_membrane_options();
	} else if ( truefalseoption("loops") ) {
		mode = "loops";
		get_loopmode_options();
	} else if ( truefalseoption("pdbstats") ) {
		mode = "pdbstats";
		get_pdbstats_options();
	} else if ( truefalseoption("interface") ) {
		mode = "interface";
		get_interface_ddg_options();
	} else if ( truefalseoption("barcode_stats") ) {
		mode = "bc_stats";
		get_barcode_stats_options();
  } else if ( truefalseoption("featurize") ) {
		mode = "pose1";
		get_pose1_options();
		get_featurizer_options();
	} else if ( truefalseoption("pKa") ) {
		mode = "pKa";
		pKa_mode::get_commandline_options();
	} else if (truefalseoption("pose_looping")) {
			mode = "pose_looping";
			get_pose_looping_options();
	} else if (truefalseoption("domain_insertion")) {
			mode = "domain_insertion";
			get_domain_insertion_options();
	}	else if (truefalseoption("antibody_modeler")) { // aroop: Arvind's antibody modeling server
		mode = "antibody_modeler";
		get_antibody_modeler_options();
	} else if ( truefalseoption("pose1") ||
							truefalseoption("bk_min") ||
							truefalseoption("adna") ||
							truefalseoption("pdna") ||
							truefalseoption("prna") ||
							truefalseoption("jumping") ||
							truefalseoption("pose_abinitio") ||
							truefalseoption("pose_rhiju") ||
							truefalseoption("extract") ||
							truefalseoption("pose_sse") ||
							truefalseoption("extract_segment") ||
							truefalseoption("centroid_information") ||
							truefalseoption("pose_idealize") ||
							truefalseoption("close_chainbreaks") ||
							truefalseoption("pose_barcode_stats")  ||
							truefalseoption("evolution") ||
							truefalseoption("pose_stats") ||
							truefalseoption("fibril") ||
//							truefalseoption("jjh_loops") ||
							truefalseoption("jumpamer") ||
							truefalseoption("homolog_rescore") ||
							truefalseoption("homolog_score_list") ||
							truefalseoption("cluster_frags") ||
							truefalseoption("cluster_by_maxsub") ||
							truefalseoption("fragment_quality") ||
							truefalseoption("homolog_distances") ||
							truefalseoption("mj_min") ) {
		// calls pose routines after initialize_rosetta
		// no calls to initialize_query, initialize_start, or initialize_decoy
		// all options are false except for output_fa in order to trigger
		// a call to get_packing_options
		//
		// maybe will have "pose2/3/4..." modes that call
		// initialize_query/start/decoy?
		if ( truefalseoption("jump_relax") ) get_abrelax_options();

		mode = "pose1";
		get_pose1_options();
	}
	std::cout << "Rosetta mode: " << mode << std::endl;

	bool setup_fa; // does fullatom need to be setup?
	get_io_options( mode, setup_fa, number_of_output );

	initialize_score_reweights();

	// change max dis and size of etables if necessary
	if( setup_fa && truefalseoption("fa_max_dis") ) {
		using namespace pdbstatistics_pack;
		fa_max_dis = realafteroption("fa_max_dis");
		//std::cerr << "fullatom_setup.cc: CHANGING fa_max_dis to " << fa_max_dis << "!!!!!!!!!!!!!" << std::endl;
		fa_max_dis2 = fa_max_dis * fa_max_dis;
		safe_max_dis2 = fa_max_dis2 - epsilon ;
		etable_disbins = static_cast< int >( fa_max_dis2 * fa_bins_per_A2 ) + 1 ;
	}
	pdbstatistics_pack::MAX_ETABLE_DISBINS() = pdbstatistics_pack::etable_disbins;

	if( truefalseoption("hydrogen_interaction_cutoff") ){
		using namespace pdbstatistics_pack;
		hydrogen_interaction_cutoff = realafteroption("hydrogen_interaction_cutoff");
		hydrogen_interaction_cutoff = hydrogen_interaction_cutoff*hydrogen_interaction_cutoff;
		std::cerr << "setting hydrogen_interaction_cutoff to: "
							<< sqrt(hydrogen_interaction_cutoff) << std::endl;
	}

	if ( truefalseoption( "no_hb_env_dep" ) ){
		set_use_W_hb_env_dep_tk( false );
		set_use_W_hb_env_dep_lin( false );
	}

	//car other options that are required to control GLOBAL behavior
	disulfides::BOUNDARY::get_disulfide_options_and_init_module(); //chu initialize disulf flags
	bool const prna = truefalseoption("prna"); //pose_rna
	if( truefalseoption( "enable_dna" ) || prna) set_enable_dna( true );
	if( truefalseoption( "enable_rna" ) || prna) enable_rna( true );
	design::dna_interface = truefalseoption( "dna_interface" );
	//lin enable ligand as aa
	if ( truefalseoption("enable_ligand_aa") ) {
		set_enable_ligaa_flag( true );
		enzyme::read_lig2=truefalseoption( "lig2_file" );

	}
	//KMa phospho_ser(SEP) 2006-01
	if ( truefalseoption("phospho_ser") ) set_pser( true );
	if ( truefalseoption("surface") ) surface::mode.flag.enable(); //KMa surface
	if ( truefalseoption("loops") ) set_loop_flag(true);
	if ( truefalseoption("taboo") ) set_taboo_flag(true);
	if ( truefalseoption("vary_omega") ) score_set_vary_omega(true);

	//jss fine Hbond categories
	if ( truefalseoption("fine_hb_categories") ) 	set_fine_hb_categories(true);

	if ( truefalseoption("geometric_sol") ) {
		set_geometric_sol_flag(true);

		// jk Set the weight relative to the Hbond term
		float geo_weight;
		realafteroption("geometric_sol_weight", 0.4, geo_weight);
		set_geometric_sol_weight( geo_weight );

		// jk For now, turn off trie if geometric_sol is used
		param_rotamer_trie::use_rotamer_trie = false;
    rotamer_trial_energies::use_trials_trie_flag = false;

		if ( get_ligand_flag() ) {
			std::cout << "geometric solvation not yet supported with ligand mode" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}

	} else {
		set_geometric_sol_flag(false);
	}

	// Set MAX_AA
	int maxaa;
	if ( dna_variables::enable_dna ) {
		maxaa = 29;
	} else if ( add_pser() ) {
		maxaa = 21;//KMa phospho_ser(SEP) 2006-01
	} else {
		maxaa = 20;
	}

	//lin additioanl set for ligand as aa
	if( get_enable_ligaa_flag() ) {
		param_aa::ligand_aa_vector.push_back(++maxaa);
		param_aa::ligand_aa_vector.push_back(++maxaa);
		int lig_aa_num=2;
		if ( truefalseoption( "number_of_ligands" ) ) {
			int lig_aa_num1=0;
			intafteroption("number_of_ligands", 2, lig_aa_num1);
			if( lig_aa_num1 > 2 ){
				lig_aa_num=lig_aa_num1;
			}
			for( int lig_num = 3; lig_num <= lig_aa_num; lig_num++ ){
				param_aa::ligand_aa_vector.push_back(++maxaa);
			}
		}
		enable_ligaa_ns::MAX_LIGAA() =lig_aa_num;

//		param_aa::lig1 = ++maxaa;
//		param_aa::lig2 = ++maxaa;
		//std::cout<<"DEBUG:lig1 "<<param_aa::lig1<<" "<<param_aa::lig2<<" MAX_AA "<<maxaa<<std::endl;
	}

	MAX_AA() = maxaa;

	//KMa phospho_ser Set MAX_AA PLUS MODIFIED RESIDUES
	if ( add_pser() )	{
		MAX_AA_PLUS() = 21;
	} else {
		MAX_AA_PLUS() = 20;
	}

	if ( ( !files_paths::multi_chain ) && truefalseoption("multi_chain") )
		enable_multi_chain();

	if ( setup_fa ) get_packing_options();
	get_score_contact_options(); //mj
	setup_default_scorefxn(setup_fa);
	setup_runlevel();

	// jk If symmetric design, set num_clones, set MAX_CLONES
	intafteroption("sym_des_units",1,num_clones);
	MAX_CLONES() = num_clones;
	--num_clones;

	// Set MAX_AA_VARIANTS
	int max_aa_variants;
	if ( design::explicit_h2o ) {
		max_aa_variants = 120;
	} else if ( design::hydrate_dna ) {
		max_aa_variants = 100;
	} else if ( design::try_both_his_tautomers || get_pH_packing_flag() ) {
		max_aa_variants = 10; //rh Changed from 4
	} else {
		max_aa_variants = 1;
	}
	int more_var = 1;
	if ( termini_ns::use_N_terminus ) ++more_var;
	if ( termini_ns::use_C_terminus ) ++more_var;
	MAX_AA_VARIANTS() = max_aa_variants * more_var;

	// cant do this earlier b/c of some code in packer_weights setup
	// triggered by get_ligand_flag()
	if ( truefalseoption("mj_min") && !get_ligand_flag() ) {
		set_ligand_flag( true );
		std::cout << "temporarily setting ligand flag from FALSE to TRUE" << std::endl;
	}

	// Set MAX_ATOM - bkidd increased ligand number for GTP/GDP ligands
	int max_atom;
	if ( get_enable_ligaa_flag() || get_ligand_flag() ) {
		max_atom = 200;
	} else if ( dna_variables::enable_dna ) {
		max_atom = 34;
		if( design::hydrate_dna ) {
			max_atom = 36;
		}
	} else if ( design::explicit_h2o ) {
		max_atom = 29;
	} else {
		max_atom = 24;
	}
	if ( termini_ns::use_N_terminus || termini_ns::use_C_terminus )
		max_atom += 2;
	MAX_ATOM() = max_atom;


	// Set MAX_ATOMTYPES  --  changed max_atomtypes to account for 5 new sugar atoms, bkidd 12/06
	if ( get_ligand_flag() || get_enable_ligaa_flag() ) {
		MAX_ATOMTYPES() = 71;
	} else if ( design::explicit_h2o || water::read_hetero_h2o ) {
		MAX_ATOMTYPES() = 26;
	} else if ( dna_variables::enable_dna ) {
		MAX_ATOMTYPES() = 71;
	} else if ( add_pser() ) {
		MAX_ATOMTYPES() = 71;//KMa phospho_ser proteins
	} else {
		MAX_ATOMTYPES() = 25;
	}

	// Set MAX_REALTYPES  --  changed max_realtypes to account for 5 new sugar atoms, bkidd 12/06
	if ( get_ligand_flag() || get_enable_ligaa_flag() ) {
		MAX_REALTYPES = 54;
	} else if ( add_pser() ) {
		MAX_REALTYPES = 55;//KMa phospho_ser proteins
	} else {
		MAX_REALTYPES = 26;
	}

	// Set HETERO_ATOM_MAX
//	if ( ( get_ligand_flag() ) || ( water::read_hetero_h2o ) ) {
		HETERO_ATOM_MAX() = 5000;
//	} else {
//		HETERO_ATOM_MAX() = 0;
//	}

	// Set MAX_NEIGH
	MAX_NEIGH() = files_paths::max_frags;

	// Activate (allocate) trim_hijack arrays
	if ( loops_ns::loop_bool::trim ) loops_ns::trim_hijack::active = 1;

	// modify_pot stuff (sheffler)
	// modify_pot stuff (sheffler)
	realafteroption("mod_hhrep_height",   fullatom_setup_ns::mod_hhrep_height  ,
									                      fullatom_setup_ns::mod_hhrep_height   );
	realafteroption("mod_hhrep_width",    fullatom_setup_ns::mod_hhrep_width   ,
									                      fullatom_setup_ns::mod_hhrep_width    );
	realafteroption("mod_hhrep_center",   fullatom_setup_ns::mod_hhrep_center  ,
 									                      fullatom_setup_ns::mod_hhrep_center   );
	realafteroption("mod_hhrep_exponent", fullatom_setup_ns::mod_hhrep_exponent,
									                      fullatom_setup_ns::mod_hhrep_exponent );
	realafteroption("smooth_etable_ljweight", fullatom_setup_ns::smooth_etable_ljweight,
									                          fullatom_setup_ns::smooth_etable_ljweight );
	realafteroption("smooth_etable_solvweight", fullatom_setup_ns::smooth_etable_solvweight,
									                            fullatom_setup_ns::smooth_etable_solvweight  );


}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param  mode - [in/out]? -
/// @param[out]   setup_fa - out - will fullatom functions be used?
/// @param  number_of_output - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///   The only options that belong in this file are those that affect
///   global i/o.  If you have mode-specific i/o DO NOT SET IT UP HERE!!
///
///   These are options that define the behavior of initialize, input_pdb &
///   output_decoy only!!
///
///   IF ITS NOT DECLARED IN FILES_PATHS.H, DO NOT SET ITS VALUE HERE!!
///
/// @references
///
/// @author car 10/14/2003
///
/////////////////////////////////////////////////////////////////////////////////
void
get_io_options(
	std::string const & mode,
	bool & setup_fa,
	int & number_of_output
)
{
	using namespace files_paths;
	using namespace basic::options;
	using namespace design;

//--------------------------
//car define query: code,protein_name,protein_chain
//--------------------------
	int i = arg_count - 1;
	if ( i >= 3 ) {
		code = sized( arg_vector[ 1 ], 2 ); // Enforce length = 2
		protein_name = sized( arg_vector[ 2 ], NAME_LENGTH ); // Enforce length = NAME_LENGTH
		protein_chain = arg_vector[ 3 ][ 0 ];
	} else {
		code = "--";
		protein_name = "----";
		protein_chain = '-';
	}

//car allow explicit command line specification (overwrite 3 args if flags used)
	if ( truefalseoption("chain") ) {
		stringafteroption( "chain", '-', protein_chain );
		if ( protein_chain == '-' ) {
			std::cout << "must specify chain id after -chain" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
	if ( truefalseoption("protein") ) {
		stringafteroption( "protein", "----", protein_name );
		if ( protein_name[0] == '-' ) {
			std::cout << "must specify 4 letter pdb code after -protein" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}
	if ( truefalseoption("series") ) {
		stringafteroption( "series", "--", code );
		if ( code[0] == '-' ) {
			std::cout << "must specify 2 letter series code after -series" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}

//car ensure that code, protein_name and chain are known if needed
	if ( code[0] == '-' || protein_name[0] == '-' || protein_chain == '-' ) {
		if ( require_3args ) {
			std::cout << "This mode requires that the series_code " <<
			 "protein_name and chain_id be defined (no hyphens allowed)" << std::endl;
			std::cout << std::endl;
			std::cout << "EITHER:" << std::endl;
			std::cout << "The first three arguments must be:" << std::endl;
			std::cout << "series_code protein_name chain_id" << std::endl;
			std::cout << "OR:" << std::endl;
			std::cout << "use command line options:" << std::endl;
			std::cout << "-series <series_code> -protein <protein_name> " <<
			 "-chain <chain_id>" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		query_defined = false;
		code = "--"; // set these back to undefined values
		protein_name = "----";
		protein_chain = '-';
	} else {
		std::cout << "series_code = " << code << " :: protein_name is " <<
		 protein_name << ":: chain_id is " << protein_chain << '.' << std::endl;
		query_defined = true;
	}
//--------------------------
//  command line options
//--------------------------
//car options which change the default mode-specific behavior if specified
//car (ie rosetta options can be changed if these flags appear on
//car the command line, but otherwise are determined by the default
//car or mode-specific setting)
//car these option flags have previously set values;  those values should
//car NOT be altered UNLESS a flag appears on the command line;  ie use
//car syntax like: if ( command line flag detected ) set option value


	intafteroption("nstruct",default_nstruct,number_of_output);

	if ( use_fasta ) use_fasta = !truefalseoption("use_pdbseq");

	if ( truefalseoption("read_all_chains") ) read_all_chains = true;

	if ( truefalseoption("preserve_header") ) preserve_header = true;

	if ( truefalseoption("use_pdb_numbering") ) use_pdb_numbering = true;

	if ( truefalseoption("flip_symmetric_sidechains") ) {
		flip_symmetric_sidechains = true;
	}

	if ( !input_fa ) input_fa = truefalseoption("fa_input");
	if ( input_fa ) {
		output_fa = true;
		if ( !repack_input ) repack_input = truefalseoption("repack");
		if ( truefalseoption("no_optH") ) no_optH = true;
	}
	if ( !output_fa ) output_fa = truefalseoption("fa_output");

	if ( truefalseoption("overwrite") ) overwrite_pdbs = true;

	if ( !disable_filters ) disable_filters = truefalseoption("no_filters");

  if ( !output_pdb_gz ) output_pdb_gz = truefalseoption("output_pdb_gz");

  if ( !output_silent_gz ) output_silent_gz = truefalseoption("output_silent_gz");

  if ( !output_scorefile_gz ) output_scorefile_gz = truefalseoption("output_scorefile_gz");

	if ( !pose_silent_out ) pose_silent_out = truefalseoption("pose_silent_out");
  if ( !sasapack_bvalues ) sasapack_bvalues = truefalseoption("sasapack_bvalues");

	if ( truefalseoption("cenlist_values")) output_cen_vals = true;


	if ( truefalseoption( "termini" ) ){
		termini_ns::use_N_terminus = true;
		termini_ns::use_C_terminus = true;
	} else {
		termini_ns::use_N_terminus = truefalseoption( "Nterminus" );
		termini_ns::use_C_terminus = truefalseoption( "Cterminus" );
	}

  //apl trie vs trie algorithm default for design mode as of 6/30/2005
  //apl and relax mode as of 8/19/2005
  bool use_trie_tf = truefalseoption("use_trie");
  bool no_trie_tf  = truefalseoption("no_trie");

	//psh insert conformer switches in between the trie statements because
	//psh conformer has to turn trie off
	use_conformer = truefalseoption("use_conformer");
	if (use_conformer) {
		stringafteroption("use_conformer", "bbdepConfLib1.txt", confLibChoice);
		no_trie_tf = true; //apl unnecessary; no intrinsic trie/conformer conflict
	}

	use_bbind_conformer = truefalseoption("use_bbind_conformer");
	if (use_bbind_conformer){
		use_conformer = true;
		// default to the smallest library 1.8
		stringafteroption("use_bbind_conformer", "bbindConfLib1.txt", confLibChoice);
		no_trie_tf = true; //apl unnecessary; no intrinsic trie/conformer conflict
	}


  if ( ( ( use_trie_tf && ! no_trie_tf ) ||
			   ( ! no_trie_tf && mode == "design" ) ||
			   ( ! no_trie_tf && mode == "relax") ||
			   ( ! no_trie_tf && mode == "abrelax") ||
				 ( ! no_trie_tf && mode == "domain_insertion") ) &&
			 ! get_geometric_sol_flag() )
	{
		param_rotamer_trie::use_rotamer_trie = true;
	}

  //apl trie vs background algorithm default for relax mode as of 8/19/2005
  bool trials_trie_tf = truefalseoption("trials_trie");
  bool no_trials_trie_tf = truefalseoption("no_trials_trie");
  if ( ( (trials_trie_tf && !no_trials_trie_tf) ||
			   (mode == "relax" && ! no_trials_trie_tf) ||
			   (mode == "abrelax" && ! no_trials_trie_tf) ) &&
			 ! get_geometric_sol_flag() )
	{
    rotamer_trial_energies::use_trials_trie_flag = true;
  }

	if ( truefalseoption("pose_relax")) {
		if (!no_trie_tf) param_rotamer_trie::use_rotamer_trie = true;
		if (!no_trials_trie_tf) rotamer_trial_energies::use_trials_trie_flag = true;
	}

	//apl Interaction graph memory usage flag
if (truefalseoption("output_interaction_graph_memory_usage")) {
		pack::output_interaction_graph_memory_usage = true;
	}

  //apl read from interaction graph file
  bool read_from_ig_tf = truefalseoption("read_interaction_graph");
  bool write_to_ig_tf = truefalseoption("write_interaction_graph");
  if (read_from_ig_tf && ! write_to_ig_tf)
  {
    pack::read_interaction_graph_file = true;
		pack::write_interaction_graph_file = false;
  }
  //else if (read_from_ig_tf && write_to_ig_tf )
  //{
  //pack::read_interaction_graph_file = true;
  //pack::write_interaction_graph_file = true;
  //}
  else if ( write_to_ig_tf )
  {
    pack::read_interaction_graph_file = false;
    pack::write_interaction_graph_file = true;
  }
  else
  {
    pack::read_interaction_graph_file = false;
    pack::write_interaction_graph_file = false;
  }
  stringafteroption( "ig_file", "", pack::ig_file_name );

	if ( truefalseoption("packer_precompute_only"))
	{
		//apl 10/24/2006, by default use Lazy IG whenever faster
		//apl turn off LazyIG completely using packer_precompute_only flag
		pack::dynamically_decide_between_lazy_and_precomputed = false;
	}
	else if ( truefalseoption("tight_memory_restrictions") )
	{
		pack::tight_memory_restrictions = true;
		intafteroption( "MB_limit_for_rpes", 256, pack::megs_for_rpes_limit );
	}
	else if ( truefalseoption("lazy_ig"))
	{
		pack::lazy_ig = true;
		pack::dynamically_decide_between_lazy_and_precomputed = false;
	}
	else if ( truefalseoption("linmem_ig"))
	{
		pack::linmem_ig = true;
		pack::dynamically_decide_between_lazy_and_precomputed = false;
	}
	else if ( truefalseoption("minimalist_ig"))
	{
		pack::minimalist_ig = true;
		pack::dynamically_decide_between_lazy_and_precomputed = false;
	}

  //apl use sasa pack
  if ( truefalseoption("use_sasa_pack_score") ){
		pack::use_sasa_pack_score = true;
	}
	if (truefalseoption("output_dot_kinemage"))
	{
		pack::output_dot_kinemage = true;
	}

	if (truefalseoption( "pack_low_temp_annealing"))
	{
		pack::annealing_starts_at_low_temperature = true;
	}

	if ( truefalseoption( "multi_cool_annealer") )
	{
		pack::use_multi_cool_annealer = true;
		int command_line_top_to_keep;
		optional_positive_intafteroption("multi_cool_annealer", 10, command_line_top_to_keep );
		pack::MultiCoolAnnealer::set_top_to_keep( command_line_top_to_keep );
	}

	if ( truefalseoption( "debug_annealer_design" ))
	{
		design::debug_annealer_design = true;
	}

  if ( truefalseoption("no_his_his_pairE") ){
		design::no_his_his_pairE = true;
	}

//------------------------------------------------------------------------------
//car options set solely by command line
//car ie, values determined only by presence or absence of command line flag
//car these option flags do NOT have previously set values;
//car Their value MUST BE SET HERE!!!
//car ie use syntax like: option_value = command line flag exist/not-exists
//car                 or  set_option_value(command line flag exists/not-exists)
//------------------------------------------------------------------------------

	silent_input = truefalseoption("silent_input");
	if ( silent_input ) {
		by_index = truefalseoption("by_index");
		new_silent_reader = truefalseoption("new_reader") ||
			truefalseoption("fa_input");
	}
  skip_scorefile_check = truefalseoption("skip_scorefile_check");

	use_timer = truefalseoption("timer");
	count_attempts = truefalseoption("count_attempts");
	use_status = truefalseoption("status");
	use_decoy_status = truefalseoption("decoy_status");
	make_ise_movie = truefalseoption("ise_movie");
	if (!output_all) output_all = truefalseoption("output_all");
	output_chi_silent = truefalseoption("output_chi_silent");
	if (!accept_all) accept_all = truefalseoption("accept_all");

	skip_missing = truefalseoption("skip_missing_residues");
	if ( skip_missing ) allow_missing = true;

	if ( input_fa ) {
		active_rotamer_options.use_input_sc = truefalseoption("use_input_sc");
		active_rotamer_options.use_input_cb = truefalseoption("use_input_cb");
		include_inputchi = active_rotamer_options.use_input_sc;
	}

	// jk read packer weights from a file
	stringafteroption( "weightfile", "none", weightfile::weight_fname );
	if ( weightfile::weight_fname != "none" ) {
		weightfile::use_weightfile = true;
	}

	stringafteroption( "cst", "cst", cst_ext );
	stringafteroption( "dpl", "dpl", dpl_ext );
	stringafteroption( "resfile", "none", resfile ); // which resides to repack and redesign
	stringafteroption( "equiv_resfile", "none", equiv_resfile); // which residues to make equivalent
	                                                            // when packing multiple structures

//lin
	auto_resfile = truefalseoption("auto_resfile");
	if ( auto_resfile ) resfile = "auto";

	chain_last_char = truefalseoption("chain_inc");
	check_options(mode,number_of_output);

	scorefile_output_full_filename = truefalseoption("full_filename"); // no cutting

	//pb: allow re-mapping of the start sequence
	map_start_sequence = truefalseoption("map_sequence");
	if ( map_start_sequence ) {
		if ( !query_defined ) {
			std::cout << "STOP: need query defined for map_sequence option.";
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		require_frags = true;
	}

	intafteroption( "max_frags", 200, max_frags );

//car set return values
	setup_fa = output_fa;

	//rhiju
	output_centroids = truefalseoption( "output_centroids" );

//dek for BOINC to allow addition of a prefix to input file names
// that start with protein_name by default
  stringafteroption( "protein_name_prefix", "", protein_name_prefix );
// and for fragment files
  stringafteroption( "frags_name_prefix", "", frags_name_prefix );

	//rhiju A new mode that allows ab initio to be carried out
	//rhiju with a homolog sequence, and relax to be carried out with
	//rhiju a query sequence. There's a mapping in between, hence the
	//rhiju option is called "homolog_to_query_mapping".
	if (truefalseoption("protein_name_prefix_homolog")) {
		homolog_to_query_mapping = true;
		adjust_options_protein_prefix_homolog();
	}

	if ( truefalseoption("use_homolog_env_stats") ) use_homolog_env_stats = true;
	if ( truefalseoption("use_homolog_pair_stats") ) use_homolog_pair_stats = true;
	if ( truefalseoption("use_homolog_cendist_stats") ) use_homolog_cendist_stats = true;

	if ( truefalseoption("use_homolog_env_stats_in_farlx") ) use_homolog_env_stats_in_farlx = true;
	if ( truefalseoption("use_homolog_pair_stats_in_farlx") ) use_homolog_pair_stats_in_farlx = true;
	if ( truefalseoption("use_homolog_cendist_stats_in_farlx") ) use_homolog_cendist_stats_in_farlx = true;

	//rj option to turn on printing hbond info in the pdb stats when writing out pdb
	//rj files. currently, hbond info is only printed when nucleic acids are present.
	//rj making this a mode-independent flag since multiple modes can write out
	//rj pdb files (I think).
	output_hbond_info = truefalseoption("output_hbond_info");

	//For docking_pose_symmetry stuff.
	docking::monomer_input =  truefalseoption("monomer_input");

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_refine_options()
{
	using namespace files_paths;

	require_start = true;
	disable_filters = true;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_score_options()
{
	using namespace files_paths;

	require_3args    = false;
	require_start    = true;
	use_fasta        = false;
	require_frags    = false;
	disable_filters  = true;
	output_coord     = false;
	refold_input     = false;
	idealized_structure = false;
	allow_missing = true;
	default_nstruct = 0;

	if ( truefalseoption("dock") ) {
		get_docking_options();
	} else if ( truefalseoption("int_score") ) {
		get_interface_options();
	}

	if ( truefalseoption("refold") || truefalseoption("use_native_centroid") ) {
		refold_input     = true;
		idealized_structure = true;
		allow_missing = false;
	}

	if ( truefalseoption("silent_input") ) {
		idealized_structure = true;
	}

	if ( truefalseoption("domain_insertion") ) {
		get_domain_insertion_options();
	}

	if ( truefalseoption("antibody_modeler") ) {
		get_antibody_modeler_options();
	}

	if ( truefalseoption("apply_filters") ) {
		disable_filters = false;
		output_all = true;
		accept_all = true;
	}
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
void
get_barcode_stats_options()
{
	using namespace files_paths;

	require_3args    = true;
	require_start    = false;
	use_fasta        = true;
	require_frags    = true;
	disable_filters  = true;
	output_coord     = false;
	refold_input     = false;
	idealized_structure = true;
	allow_missing = false;
	default_nstruct = 1;
}

void
get_featurizer_options()
{
	using namespace files_paths;

	require_3args    = true;
	require_start    = false;
	use_fasta        = false;
	require_frags    = false;
	disable_filters  = true;
	output_coord     = false;
	refold_input     = false;
	idealized_structure = true;
	allow_missing = false;
	default_nstruct = 1;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_idealize_options()
{
	using namespace files_paths;

	require_3args = false;
	require_frags = false;
	require_start = true;
	disable_filters = true;
	use_fasta = false;
	default_nstruct = 1;
	use_constraints = false;
	use_scorefile = false;
	refold_input = false;
	use_pdb_numbering = true;
	overwrite_pdbs = true;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///az private relax_options have been gathered in a separate header file
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_relax_options()
{
	using namespace files_paths;
	using namespace relax_options;

	default_nstruct = 1;
	require_start = true;
	disable_filters = true;
	set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score6","score12");

	sim_aneal = truefalseoption("sim_aneal");
	cenrlx = truefalseoption("cenrlx");
	farlx = truefalseoption("farlx");
	record_rms_before_relax = truefalseoption("record_irms_before_relax");
	force_expand = truefalseoption("force_expand");
	looprlx = truefalseoption("looprlx");
	rb_rlx = truefalseoption("rb_relax");
	//psh temporary
	//looprlx_nonideal_silent_output = truefalseoption("looprlx_nonideal_silent_output");
	use_pose_relax = truefalseoption("pose_relax");
	if ( !cenrlx && !farlx && !looprlx ) {
		cenrlx = true;
		farlx = true;
	}
	if ( looprlx ){
		idealized_structure = false;
		allow_missing = true;
		randomize_missing = true;
	}
	if ( farlx ) {
		score_set_vary_omega( true );
		minimize = truefalseoption("minimize");
		if ( minimize ) {
			sc_only = truefalseoption("sc_only");
			bb_only = truefalseoption("bb_only");
			if ( bb_only && sc_only ) {
				std::cout << "WARNING: bb_only and sc_only are both TRUE" << std::endl;
				std::cout << "Both backbone and sidechain will be minimized!!!" << std::endl;
				sc_only = false; // "
				bb_only = false;
			}
			if ( !cenrlx && !looprlx ) {
				use_fasta = false;
				require_frags = false;
				require_3args = false;
			}
		}
	}
	if ( sc_only && !cenrlx ) idealized_structure = false;
	if( farlx || truefalseoption("loop_farlx") || truefalseoption("loop_farlx_only") ) output_fa = true;
	//ja don't need fragments or secondary structure pred. for relax
	if( truefalseoption("loop_farlx_only") ) {
		require_frags = false;
		use_fasta = false;
	}
	if ( use_pose_relax ) {
		if ( !truefalseoption("pose_relax_fragment_moves") && !looprlx && !cenrlx )	require_frags = false;
		vary_sidechain_bond_angles = truefalseoption("vary_bonds");
	}
	if ( truefalseoption("use_input_bond") ) {
		use_fasta = false;
	}

	if ( get_skip_fragment_moves() ){
		require_frags = false;
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief abrelax is a mode where ab initio is run first, then relax.
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_abrelax_options()
{
	using namespace files_paths;
	using namespace relax_options;

	//default_nstruct = 1;
	//require_start = true;
	disable_filters = false; // allows filters to be applied after abinitio
	disable_output_decoy_filters = true; // always output a decoy
	//set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score4","score12");

	sim_aneal = truefalseoption("sim_aneal");
	cenrlx = truefalseoption("cenrlx");
	farlx = true;
	record_rms_before_relax = truefalseoption("record_irms_before_relax");
  force_expand = truefalseoption("force_expand");
	if ( farlx ) {
		score_set_vary_omega(true);
		minimize = truefalseoption("minimize");
		if ( minimize ) {
			sc_only = truefalseoption("sc_only");
			bb_only = truefalseoption("bb_only");
			if ( bb_only && sc_only ) {
				std::cout << "WARNING: bb_only and sc_only are both TRUE" << std::endl;
				std::cout << "Both backbone and sidechain will be minimized!!!" << std::endl;
				sc_only = false; // "
				bb_only = false;
			}
// 			if ( !cenrlx ) {
// 				use_fasta = false;
// 				require_frags = false;
// 			}
		}
	}
	//if ( sc_only && !cenrlx ) idealized_structure = false;
	idealized_structure = true;
	output_fa = farlx;

}

///////////////////////////////////////////////////////////////////////////////
void
get_pose1_options()
{
	using namespace files_paths;
	using namespace relax_options;

	require_3args = false;
	require_frags = false;
	require_start = false;
	use_fasta = false;
	use_constraints = false;
	input_fa = false;
	output_fa = true; // trigger fa initialization?
	idealized_structure = false;
	allow_missing = false;
	skip_missing = false;
	set_pose_flag(true);

	// mode-specific stuff
	if( truefalseoption("cst_mode") ) {
		input_fa = true;
		require_start = true;
	  use_scorefile = false;
  	output_coord = false;
		//multi_chain = true;
	}
	if ( truefalseoption("mj_min") ) {
		dna_variables::enable_dna = true; // make space for ligand aa
	}
}

///////////////////////////////////////////////////////////////////////////////
///
/// @brief
//  options for pose_looping mode
///
/// @details
//  Things to note:
//  disable filters is true
//  pose flag is not set!!
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author Monica Berrondo
///
/////////////////////////////////////////////////////////////////////////////////
void
get_domain_insertion_options()
{
	using namespace files_paths;
	using namespace relax_options;

	require_3args = false;
	require_frags = true;
	require_start = false;
	use_fasta = false;
	use_constraints = false;
	input_fa = false;
	output_fa = false;
	idealized_structure = false;
	allow_missing = false;
	skip_missing = false;
	disable_filters = true;
	domain_insertion = true;
//	param_rotamer_trie::use_rotamer_trie = true;
	if (truefalseoption("create_fasta")) create_domins_fasta = true;
	if (truefalseoption("score")) {
		set_pose_flag(true);
		select_default_scorefxns("score8di", "score12di");
	}
}

///////////////////////////////////////////////////////////////////////////////
///
/// @brief
//  options for antibody modelling mode
///
/// @details
//  Things to note:
//  disable filters is true
//  pose flag is not set!!
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author Aroop Sircar
///
/////////////////////////////////////////////////////////////////////////////////
void
get_antibody_modeler_options()
{
	using namespace files_paths;
	using namespace relax_options;
	using namespace docking;

	idealized_structure = false;
	require_frags = false;
	require_start = true;
	use_fasta = false;
	default_nstruct = 1;
	use_constraints = false;
	disable_all_filters();
	disable_filters = true;
	input_fa = true;
	no_optH = true; // prevent Hydrogen optimization in input structure
	output_fa = true;
	use_pdb_numbering = true;
	read_all_chains = true;
	refold_input = false;
	repack_input = false;
	multi_chain = true;
	select_default_scorefxns("score4L","score12");
	antibody_modeler = true;
	return;
}

///////////////////////////////////////////////////////////////////////////////
///
/// @brief
//  options for pose_looping mode
///
/// @details
//  Things to note:
//  disable filters is true
//  pose flag is not set!!
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author Monica Berrondo
///
/////////////////////////////////////////////////////////////////////////////////
void
get_pose_looping_options()
{
	using namespace files_paths;
	using namespace relax_options;

	require_3args = false;
	require_frags = true;
	require_start = false;
	use_fasta = false;
	use_constraints = false;
	input_fa = false;
	output_fa = false;
	idealized_structure = false;
	allow_missing = false;
	skip_missing = false;
	disable_filters = true;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_loopmode_options()
{
	using namespace files_paths;
	using namespace loops_ns;

	set_loop_flag(true);

	require_start = true;
	allow_missing = true;
	use_filter(knot_type) = true;
	use_filter(co_type) = false;
	use_filter(rg_type) = false;
	use_filter(sheet_type) = false;
	refold_input = true;
	idealized_structure = false;

	grow = truefalseoption("grow");
	permute = truefalseoption("permute");
	screen = truefalseoption("screen");
	sequential = truefalseoption("sequential");
	loop_bool::trim = truefalseoption("trim"); // Disambiguate from ObjexxFCL function trim
	fa_refine = truefalseoption("fa_refine");
	idl_breaks = truefalseoption("idl_breaks");
	fold_with_dunbrack = truefalseoption("fold_with_dunbrack");

	fold = truefalseoption("fold");
	if ( !grow && !permute && !screen && !fa_refine && !idl_breaks && !fold_with_dunbrack)
		fold = true;

	fast = truefalseoption("fast");
	dunbrack_close = truefalseoption("dunbrack_close");
	set_fa_lj_rep_slope("lowres");

	require_frags = false;
	use_fasta = false;
	if ( screen ) {
		default_nstruct = 1;
		output_coord = false;
		use_scorefile = false;
		use_fasta = true;
		select_default_scorefxns("score0L","score12");
	}
	if ( grow ) {
		default_nstruct = 1;
		select_default_scorefxns("score3L","score12");
		output_coord = true;
		use_fasta = true;
		vdw_max_on = truefalseoption("vdw_max");
		if (vdw_max_on) {
			realafteroption("vdw_max",999., vdw_max);
		}
		rg_max_on = truefalseoption("rg_max");
		if (rg_max_on) {
			realafteroption("rg_max",999., rg_max);
		}
		wiggle_jxn = truefalseoption("wiggle_jxn");
		use_filter(knot_type) = false;
	}
	if ( permute ) {
		if ( fast ) {
			default_nstruct = 1;
			select_default_scorefxns("score3L","score12");
		} else {
#ifdef BOINC
			default_nstruct = 10;
#else
			default_nstruct = 1000;
#endif
			select_default_scorefxns("score4L","score12");
		}
		output_coord = true;
		use_fasta = true;
	}
	if ( idl_breaks ) {
		output_coord = true;
		select_default_scorefxns("score4L","score12");
		default_nstruct = 1;
	}

	if ( fold || fold_with_dunbrack ) {
#ifdef BOINC
		default_nstruct = 10;
#else
		default_nstruct = 1000;
#endif
		output_coord = true;
		if ( fold_with_dunbrack && truefalseoption("refine_only") ) {
			require_frags = false;
			use_fasta = false;
		} else {
			require_frags = true;
			use_fasta = true;
		};

		select_default_scorefxns("score4L","score12");
	}

	if ( fa_refine ) {
#ifdef BOINC
		default_nstruct = 10;
#else
		default_nstruct = 1000;
#endif
		output_fa = true;
		select_default_scorefxns("score4L","score12");
		set_fa_lj_rep_slope("highres");
	}

	if ( default_nstruct == 1 ) overwrite_pdbs = true;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief setup options for building loops for general use
///     (not necessarily loop mode)
///
/// @details
///
/// @global_read flags from command line
///
/// @global_write loops.h flags
///
/// @remarks written for use in dock mode
///     copied get_loopmode_options and removed mode-specific flags
///
/// @references
///
/// @author Jeff Gray 3/11/4
///
/////////////////////////////////////////////////////////////////////////////////
void
get_loop_options()
{
	using namespace files_paths;
	using namespace loops_ns;

	fast = truefalseoption("fast");

	require_frags = true;
	use_fasta = true;
//	select_default_scorefxns("score4L","score10L");

	use_filter(knot_type) = true;
	use_filter(co_type) = false;
	use_filter(rg_type) = false;
	use_filter(sheet_type) = false;
	disable_filters = false;
	idealized_structure = false;
	refold_input = false;
	set_pose_loop_flag( truefalseoption("pose") );
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_assemble_options()
{
	using namespace assemble_options;
	using namespace files_paths;

	if ( truefalseoption( "extend" )  || truefalseoption( "perturbe_starting" ) || truefalseoption("thread")) {
		require_start = true;
		disable_filters = true;
		default_nstruct = 1;

		input_fa = truefalseoption( "fa_input" );
		if ( input_fa ) {
			output_fa = true;
		}

//az no frag required => don't generate errors if
//az they are not present
		require_frags = false;
		use_scorefile = false;
		use_fasta = false;
		require_3args = false;
		read_all_chains = true;
		idealized_structure = false;

//az "true" assemble domain ...
	} else {

		require_start = true;
		disable_filters = true;
#ifdef BOINC
		default_nstruct = 10;
#else
		default_nstruct = 100;
#endif

//		select_default_scorefxns("score4","score11"); // why  ?
		if ( truefalseoption( "fa_input" ) ) {
			input_fa = true;
			output_fa = true;
		} else {
			input_fa = false;
			output_fa = false;
		}

		require_frags = true;
		use_scorefile = true;
		use_fasta = true;
		require_3args = true;
		read_all_chains = true;
		idealized_structure = false;


	}
	allow_missing = true;
	assemble_flag = true;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_design_options()
{
	using namespace design;
	using namespace files_paths;
	using namespace docking;
	using namespace favor_residue;
	using namespace analyze_interface_ddg_ns;
	using namespace design_loops_flags;

	disable_filters = true;
	require_start = true;
	refold_input = false;
	use_fasta = false;
	default_nstruct = 1;
	input_fa = true;
	output_fa = true;
	repack_input = false;
	select_default_scorefxns("score4","score12");

//ds   use pdb numbering
	use_pdb_numbering = true;

//bk   only pack rotamers, no backbone movement
	onlypack = false;
	if ( truefalseoption("onlypack") ) onlypack = true;

//bk   only design sequence, no backbone movement
	fixbb = false;
	if ( truefalseoption("fixbb") ) fixbb = true;

//xa   use parallel packing during design
	pack_in_parallel = false;
	if ( truefalseoption("pack_in_parallel") ) pack_in_parallel = true;
	if( pack_in_parallel ){
		fixbb = true;
		enable_multi_chain();
	}

// murphp use directed design
	if( truefalseoption("directed_design") ) {
		use_directed_design = true;
		fixbb = true;
		stringafteroption("directed_design","none",directed_design_filename);
	}

	// jk apply an offset to Eref of charged residues, with sign dependent on chain
	// jk (ie. to get electrostatic complementarity for interface design)
	design_interface_electrostatic_complementarity = false;
	realafteroption("interface_electrostatic_complementarity_shift", 0.,
									interface_electrostatic_complementarity_shift);
	if ( std::abs(interface_electrostatic_complementarity_shift) > 0.0001 ) {
		design_interface_electrostatic_complementarity = true;
	}

	///////////
//oe  reverse complement peptides during design
	reverse_comp = false;
	if ( truefalseoption("reverse_comp") ) reverse_comp = true;
	if( reverse_comp ){
		fixbb = true;
		enable_multi_chain();
	}

//xa   integer multiple of five annealing loops in outer (mutating)
//xa   loops of ParallelSimAnneal. Default of 3 results in 15 loops
	intafteroption("conv_limit_mod", 3, conv_limit_mod);

//xa   if parallel packing is being used for negative design,
//xa   the preferred state is specified by naming one chain in the structure/complex

	stringafteroption("chain_in_preferred_state", "none", chain_in_preferred_state);

	stringafteroption("chain_in_ref_state", "none", chain_in_ref_state);

//jjh option to select closest rotamer to native, skipping energies and mc
	rotamerize = false;
	if ( truefalseoption("rotamerize") ) {
		rotamerize = true;
		///jjh Rotamerize implies that you are limited to the native sequence
		fixbb = false;
		onlypack = true;
	}

	if ( truefalseoption("point_mutation") ) point_mutation = true;

//bk use monte carlo minimization as the last step during rotamer optimization
	mcmin_trials = false;
	if ( truefalseoption("mcmin_trials") ) mcmin_trials = true;

//jk don't allow Cys or Gly unless they're present in the starting struct
  if ( truefalseoption("no_new_CG") ) no_new_CG = true;

//jjh option to use multistate design techniques
	multistate = false;
	if ( truefalseoption("multistate") ) multistate = true;

//jjh option to use dna-protein structural motifs
	dna_motifs = false;
	if ( truefalseoption("dna_motifs") ) dna_motifs = true;

	float fdummy;
	realafteroption("motifs_close_enough", 1.2, fdummy);
	DNA_Motifs::Motifs_close_enough = fdummy;

	realafteroption("motifs_far_dist", 5.0, fdummy);
	DNA_Motifs::Motifs_far_dist = fdummy;

	realafteroption("Motifs_bump_value", 100.0, fdummy);
	DNA_Motifs::Motifs_bump_value = fdummy;

	realafteroption("motifs_match_weight", 1.0, fdummy);
	DNA_Motifs::motifs_match_weight = fdummy;

	int idummy;
	intafteroption("motifs_rots_required", 2, idummy);
	DNA_Motifs::Motifs_num_rots_required = idummy;

	design_loops = false;
	if ( truefalseoption("design_loops") ) {
		design_loops = true;
		select_default_scorefxns("score4","score12");
		output_coord = true;
		use_scorefile = true;
		//design_loops_dock, turn on multiple chain option
		stringafteroption( "design_loops", "hold", design_loops_option );
		if ( design_loops_option == "dock"){
			enable_multi_chain();
			design::design_loops_flags::design_loops_dock = true;
		}else{
			design::design_loops_flags::design_loops_hold = true;
		}
	}

	design_min_inter = false;
	if ( truefalseoption("design_min_inter") ) {
		design_min_inter = true;
		enable_multi_chain();
		select_default_scorefxns("score4","score12");
		output_coord = true;
		use_scorefile = true;
	}

	dock_des_min_inter = false;
	if ( truefalseoption("dock_des_min_inter") )
	{
		dock_des_min_inter = true;
		enable_multi_chain();
		select_default_scorefxns("score4","score12");
		output_coord = false;
		use_scorefile = true;
		if (truefalseoption("timelimit")) timelimit = intafteroption("timelimit");
		n_partner = 1;
	}

	flex_peptide = false;
	if ( truefalseoption("flex_peptide") )
	{
		flex_peptide = true;
		enable_multi_chain();
		select_default_scorefxns("score4","score12");
//		output_coord = true;
		use_scorefile = true;
//    docking_score_norepack = true; // structures are expected to have been minimized beforehand
	}


//bk   designs based on a mutation list
	mut_list = false;
	if ( truefalseoption("mut_list") ) mut_list = true;
	if ( mut_list ) {
		read_all_chains = true;
	}

//lin   cluster designs
	cluster_design = false;
	if ( truefalseoption("cluster_design") ) {
		cluster_design = true;
		use_cluster_size = truefalseoption("cluster_size");
		if ( use_cluster_size ) {
			intafteroption("cluster_size",1,cluster_size);
			use_real_cluster_only=truefalseoption("use_spatial_cluster");
			if ( use_real_cluster_only ) {
				real3afteroption("use_spatial_cluster",5.0,cluster_dis_cut1,8.0,
												 cluster_dis_cut2,12.0,cluster_dis_cut3);
			}
			design_neighbors=truefalseoption("design_neighbors");
			if ( design_neighbors ) {
				real2afteroption("design_neighbors",5.0,design_dis_cut1,8.0,
												 design_dis_cut2);
			}
			check_point=truefalseoption("check_point");
			if ( check_point ) {
				intafteroption("check_point",1,last_cpts_size);
			}
			if ( last_cpts_size < 1 || last_cpts_size > cluster_size ) {
				std::cout << "ERROR: check_point size should be in [1,cluster_size) " << std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
		}

		//the starting number of partners  3: single protein, 1: protein complex
		if ( truefalseoption("protein_complex") ) {
			n_partner = 1;
		} else {
			n_partner = 3;
		}

		repack_neighbors=truefalseoption("repack_neighbors");
		if ( repack_neighbors ) {
			real2afteroption("repack_neighbors",5.0,repack_dis_cut1,8.0,repack_dis_cut2);
		}
		minimize_structures=truefalseoption("minimize_structures");
		minimize_wt_struct=truefalseoption("minimize_wt_struct_only");
		if ( minimize_wt_struct || minimize_structures ) {
			minimize_bb=truefalseoption("min_bb");
			if ( minimize_bb && !truefalseoption("use_input_bond") ) {
				std::cout<<"ERROR: min_bb T, use_input_bond F, try use use_input_bond flag when use min_bb"<<std::endl;
				utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
			}
		}
		debug_output=truefalseoption("debug_output");
		intafteroption("n_mutatpos",1,n_mutatpos);
		if ( n_mutatpos < 1 ) {
			std::cout<<"WARNING: can not set n_mutatpos<1 "<<std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		fix_mut=truefalseoption("fix_mut");
		if ( fix_mut ) stringafteroption("mutin","mut",mutenergy_in);
		output_structures=truefalseoption("output_structures");
		if ( output_structures ) n_mutatpos = 1;
		stringafteroption("mutout","mut",mutenergy_out);
// 		linmin_trials=truefalseoption("linmin_trials");
// 		realafteroption("high_nb",100.0,high_nb);
// 		realafteroption("low_nb",0.0,low_nb);
		if ( minimize_structures && mcmin_trials ) {
			std::cout << "ERROR: can not use both flags linmin_trials and mcmin_trials " << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
	}

//jjh optimize rotamers prior to repack
	active_rotamer_options.rot_opt = truefalseoption("rot_opt");

//jjh Call Chu's rotamer_trials after design
	design_trials = false;
	if ( truefalseoption("design_trials") ) design_trials = true;

//bk   if you want the backbone to move
	mvbb = false;
	if ( truefalseoption("mvbb") ) mvbb = true;

//jk   design and dock (dock_mcm only)
	desock = false;
	if ( truefalseoption("desock") ) {
		desock = true;
		enable_multi_chain();
		select_default_scorefxns("score4","score12");
		output_coord = true;
	}

//bk   don't vary the sequence
	fixseq = false;
	if ( truefalseoption("fixseq") )  fixseq = true;

//glb flexible tail
	tail = false;
	if ( truefalseoption("tail") ) tail = true;
	if ( truefalseoption("tail_fix_helix") ) tail_fix_helix = true;
	if ( tail || tail_fix_helix ) {
		enable_multi_chain();
	}

//bk   redesign an interface (or just one side of an interface)
	design_inter = false;
	repack_inter = false;
	fix_target_seq = 0;
	if ( truefalseoption("repack_inter") ) {
		repack_inter = true;
		design_inter = true;
	} else {
		fix_target_seq = intafteroption( "fix_target_seq", 0 );
		if ( fix_target_seq > 0 ) {
			design_inter = true;
		}
	}
	if ( ! design_inter ) {
		design_inter = truefalseoption("design_inter");
	}

	if ( design_inter ) {
		enable_multi_chain();
	}

//bk  used to alter binding specificity, Tanja's second site suppressor stategy
	alter_spec = false;
	if ( truefalseoption("alter_spec") ) {
		alter_spec = true;
		enable_multi_chain();
	}

//xh used to mutate clusters in a protein
	mutate_cluster = false;
	if (truefalseoption("mutate_cluster") ){
		mutate_cluster = true;
	}

//xh used to output structures
	if ( truefalseoption("output_pdb") ) output_pdb = true;

//bk   design using small cluster - good for large rotamer libraries
	if ( truefalseoption("design_in_pieces") ) design_in_pieces = true;

	if ( truefalseoption( "chain_limit" ) ) {
		chain_limit = true;
		charafteroption( "chain_limit", ' ', chain_choose );
	}

	/// DNA interface design options
	if ( truefalseoption( "dna_interface" ) ) get_dna_interface_options();
/// End DNA interface options
	//////////////////////////////////////////////////////////////////////////////

//bk   used for scoring a structure
	natrot = false;
	optE_inter = false;
	if ( truefalseoption("natrot") ) {
		natrot = true;
		fixbb = true;
	}

	if ( onlypack || fixbb || desock || design_inter || alter_spec ||
			 cluster_design || mutate_cluster || mut_list || tail || tail_fix_helix ||
			 design_loops || design_min_inter || dock_des_min_inter || point_mutation|| flex_peptide
	) {
		idealized_structure = false;
		repack_input = false;
	}
	if ( truefalseoption("optE") || truefalseoption("optE_ctsa") ||
	 truefalseoption("optE_inter") ) {
		optE = true; // output E for opt Weights
		use_scorefile = false;
		output_coord = false;
		repack_input = false;
		if ( truefalseoption("optE_rep_cut") ) {
			optE_rep_cut = true;
			realafteroption( "optE_rep_cut", 500., optE_rep_cut_val );
		}
	}
	if ( truefalseoption("optE_ctsa") ) {
		optE_ctsa = true; // add extra custom junk w/o messing up brian's code
	}
	if ( truefalseoption("optE_jk") ) {
		optE_jk = true; // add extra custom junk w/o messing up brian's code
	}
	if ( truefalseoption("optE_inter") ) {
		optE_inter = true; // to use optE with interface residues.
		enable_multi_chain();
	}

	if ( fixbb || onlypack || optE || desock  || design_inter || alter_spec ||
			 mut_list || cluster_design || mutate_cluster || design_min_inter ||
			 dock_des_min_inter || design_loops || point_mutation || flex_peptide
	) {
		require_frags = false;
		require_3args = false;
	} else {
		require_frags = true;
		require_3args = true;
	}

//jjh  a file used to define competing states
	stringafteroption( "ms_comp_list", "default", ms_comp_list );

//jjh  filename for multistate output of best fitness sequence on target struct
	stringafteroption( "ms_out_name", "best_fit.pdb", ms_out_name );

//bk   a file used to calibrate weights
	stringafteroption( "Eout", "default", Eout );

//bk   a file used to compare sequence to the starting sequence
	stringafteroption( "sqc", "none", sqc );

//bk   beginning of output file name, used when fixbb = true
	stringafteroption( "pdbout", "des", pdbout );

//bk   if doing only design, specifies number of design runs
	intafteroption("ndruns",1,ndruns);
	if ( !( fixbb || onlypack || design_inter ) ) ndruns = 1;

//glb  options for tail
	intafteroption("begin",0,tail_begin);
	intafteroption("end",0,tail_end);

//ds   a file which can be used to calculate ddG binding using analyze_interface
	stringafteroption( "alter_spec_mutlist", "mutlist", alter_spec_mutfile );

//ds   a file which shows sc-sc and sc-bb energies for interface point mutations
	stringafteroption( "point_mut_energies", "none", point_mut_energies );

//ds   a file which shows "recovery" energies for interface with point mutation
	stringafteroption( "mut_mut_energies", "none", mut_mut_energies );

//ds   a file which contains a list of mutations
	stringafteroption( "design_mutlist", "none", design_mutlist );

	safety = false;
//ds   enable safety checks for design_mutlist
	if ( truefalseoption("safety") ) safety = true;

//lin  small perturbation for partner 2?
	design_dock_pert = truefalseoption("design_dock");
	if ( design_dock_pert || dock_des_min_inter ) {
			std::cout << "***apply docking perturb decoy rountin*****" << std::endl;
			real3afteroption("dock_pert",1.0,normal_perturbation,1.0,
	    parallel_perturbation,3.0,rotational_perturbation);
			enable_multi_chain();
	}

//lin give native residues a bonus so that less residues are mutated
	favor_native_residue = truefalseoption("favor_native_residue");
	if (favor_native_residue) {
		realafteroption("favor_native_residue",-1.5,native_bonus);
		std::cout << "favor_native_residue set to true and native_bonus set to " << native_bonus
		 << " as requested by the -favor_native_residue flag " << std::endl;
	}

	//mj give bonus to selected amino acids
	if ( truefalseoption("favor_polar"   ) ) { realafteroption( "favor_polar" ,    0.0, favor_polar );    favor_property_residue = true; }
	if ( truefalseoption("favor_nonpolar") ) { realafteroption( "favor_nonpolar" , 0.0, favor_nonpolar ); favor_property_residue = true; }
	if ( truefalseoption("favor_charged" ) ) { realafteroption( "favor_charged" ,  0.0, favor_charged );  favor_property_residue = true; }
	if ( truefalseoption("favor_aromatic") ) { realafteroption( "favor_aromatic" , 0.0, favor_aromatic ); favor_property_residue = true; }

	//mj in ligand small perturbation design mode allow to output only decoys that pass a given energy filter
	if ( design_dock_pert && truefalseoption( "ligand")
	 && truefalseoption( "ligand_filter_score" )
	 && truefalseoption( "ligand_filter_count" ) ) {
		ligand_filter = true;
		realafteroption( "ligand_filter_score" , 0.0, ligand_filter_score );
		intafteroption(  "ligand_filter_count" ,   0, ligand_filter_count );
	}

//bk   if designing with a fixed backbone don't use call to output in
//bk   main_rosetta, use an output function built into the packer.
	if ( fixbb || onlypack || desock || design_inter || alter_spec || mut_list || cluster_design || mutate_cluster ) {
		output_coord = false;
		use_scorefile = false;
	}

	//ja output rotamers and energies for a particular design residue to PDB file
	// for visual/energy analysis
	if ( truefalseoption( "dump_rotamers_pdb" ) ) {
		active_rotamer_options.dump_rotamers_pdb_flag = true;
		intafteroption( "dump_rotamers_pdb", 0, active_rotamer_options.dump_rotamers_pdb_res );
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @author
/// ashworth
///
////////////////////////////////////////////////////////////////////////////////
void
get_dna_interface_options()
{
	using namespace dna_variables;

	std::cout << "\nDNA-interface design mode:\n";
	enable_dna = true;

	dna_verbose = truefalseoption( "dna_verbose" );
	design_by_base = truefalseoption( "design_by_base" );

	std::cout << "\n";

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief set options required whenever -int_score is specified, regardless of
/// mode
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_interface_options()
{
// these all true in get_interface_ddg_options  -should call this function
// from get_interface_ddg_options
// (except use_pdb_numbering -- and this maybe should be true??)

	using namespace files_paths;
	using namespace water;

	enable_multi_chain();

	Wint_score_only    = true;
	Wint_repack_only   = true;
}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief set options required for interface_ddg mode
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_interface_ddg_options()
{
//tk   from design_options
//tk   needed as interface mode uses pack_rotamers
//car files_paths.h options

	using namespace analyze_interface_ddg_ns;
	using namespace design;
	using namespace files_paths;
	using namespace water;

	disable_filters       = true;
	require_start         = true;
	refold_input          = false;
	use_fasta             = false;
	default_nstruct       = 1;
	input_fa              = true;
	output_fa             = true;
	repack_input          = false;
	idealized_structure   = false;
	require_frags         = false;
	require_3args         = false;
	output_coord          = false;
	use_scorefile         = false;
	alter_spec_format     = false;
	chain_energies        = false;
	repack_neighbors      = false;
	affin_incr            = false;

//ds   use pdb numbering
	use_pdb_numbering = true;

	select_default_scorefxns("score4","score12");

//car design.h settings
	ndruns                = 1;

//lin  set packer weights to int weights and use env dep HB weights
//lin  default weight in ddg mode is pack weights for repacking
//lin  and int weights for scoring
	set_use_W_int_repack(false);
	set_use_W_int_score(true);

//tk   needed as interface mode has multiple chains
//tk   and uses docking convention:
//tk   the two protein partner are separated by "TER"
	enable_multi_chain();

//tk   ***********************************************
//tk   command line options specific to interface mode
//tk   ***********************************************
//tk   input:
//tk   ------
//tk   calculate binding energy only
	ddg_bind_only = false;
	if ( truefalseoption("ddg_bind_only") ) ddg_bind_only = true;
//tk   otherwise you need a
//tk   file which specifies mutations
	stringafteroption( "mutlist", "none", mutlist );
//tk   output:
//tk   -------
//tk   file name for containing binding energy output
	stringafteroption( "intout", "INT.OUT", intout );
//tk   run:
//tk   ----
//tk   energy function  during pack and scoring
//tk   default: pack_db.cc weights for packing
//tk   int weights for scoring (set in reset_weights_to_int_weights)
	Wpack_only = false;
	Wint_only  = false;
//tk   packing weights for packing and scoring:
	if ( truefalseoption("Wpack_only") ) {
		Wpack_only = true;
		Wint_repack_only = false;
		Wint_score_only = false;
	}
//tk   int weights for packing and scoring:
	if ( truefalseoption("Wint_only") ) {
		Wint_only = true;
		Wint_repack_only = true;
		Wint_score_only = true;
	}

//bk scan all point mutants at the interface
	pmut_scan = truefalseoption("pmut_scan");

	//ds this flag adds a filter to pmut scan so that the only mutations
	//ds allowed increase buried hydrophobic surface area across the
	//ds interface.  The number of neighbors for a wild type amino acid is
	//ds also output in the INTOUT energy file.
	if ( truefalseoption("affin_incr") ) {
		affin_incr = true;
		pmut_scan = true;
	}

//lin  calculate delta delta binding energy for alanine scan
  ala_scan = truefalseoption("ala_scan_only");
  output_mabo_input = false;
	if ( ala_scan ) output_mabo_input = truefalseoption("output_mabo_input");

//lin  the starting number of partners  3: single protein, 1: protein complex
//lin  Currenlty the monomeric ddg is using the same weight as interface ddg
//lin  the weight for monomeric ddg may need to be reparametrized and
//lin  the weight for dunbrack term and the reference energies may need recomputing
  if ( truefalseoption("monomeric_protein") ) {
    n_partner = 3;
    single_protein = true;
   } else {
     n_partner = 1;
     single_protein = false;
   }

//tk   enable safety checks for binding energy calculations
	if ( truefalseoption("safe") ) safety_check = true;
//lin  output the structure or not
	if ( truefalseoption("output_structure") ) output_structure = true;

//ds   allow analyze_interface to repack neighbors to each mutation
	if ( truefalseoption("repack_neighbors") ) repack_neighbors = true;

//bk  relax structures before and after mutation with gradient based minimization
//bk  the default is to vary all torsion angles at the interface as well
//bk  rigid body optimization
	if ( truefalseoption("min_interface") ) {
		min_interface = true;
		// constraints: default is to use them
		use_cst = !truefalseoption( "no_cst" );
		cstW = realafteroption("cstW",2.0);
		// what can move
		bb = truefalseoption("bb");     // backbone torsion
	  int_bb = truefalseoption("int_bb");  // interface backbone torsion
		chi = truefalseoption("chi");  // chi angles
		int_chi = truefalseoption("int_chi"); // interface chi angles
		rb = truefalseoption("rb"); // rigid body motion
		// default if no move flags were specified
		if (!bb && !int_bb && !chi && !int_chi && !rb ) {
			int_bb = true;
			int_chi = true;
			rb = true;
		}
	}

//bk  when repacking or minimizing, also relax the unbound state
	if ( truefalseoption("relax_unbound") ) relax_unbound = true;

//ja extra conditions moved to the appropriate place in analyze_interface_ddg.cc
//if ( truefalseoption("relax_unbound") &&
//		 (repack_neighbors || min_interface) ) relax_unbound = true;

//ds   if using analyze_interface_ddg with a mutlist file generated from
//ds   altered_specificity (design mode) then use a special format designed
//ds   to accommodate the large number of mutated structures.
	if ( truefalseoption("alter_spec_format") ) alter_spec_format = true;

//ds   include energies for each chain in binding energies table.
	if ( truefalseoption("chain_energies") ) chain_energies = true;


}


////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_pdbstats_options()
{

	using namespace files_paths;
	using namespace pdbstats;

	std::string whatkind;

	disable_filters = true;
	require_start = true;
	refold_input = false;
	idealized_structure = false;
	use_fasta = false;
	default_nstruct = 1;
	input_fa = true;
	output_coord = false;
	output_fa = false;
	repack_input = false;
	use_scorefile = false;
	allow_missing = true;
	require_frags = false;
	require_3args = false;

	select_default_scorefxns( "score4", "score12" );

//bk   different options for compiling statistics
	stringafteroption( "pdbstats", "none", whatkind );
	if ( whatkind == "none" ) {
		std::cout << "You must specify the desired kind of statistics" << std::endl;
		std::cout << "format: -pdbstats [charid]" << std::endl;
		std::cout << "the current options are:" << std::endl;
		std::cout << "distances" << std::endl;
		std::cout << "hbonds" << std::endl;
		std::cout << "avg_energies" << std::endl;
		std::cout << "phipsichi" << std::endl;
		std::cout << "dun_rot" << std::endl;
		std::cout << "rep_atm_pair" << std::endl;
		std::cout << "planes" << std::endl;
		std::cout << "cation" << std::endl;
		std::cout << "carbonyl" << std::endl;
		std::cout << "ramachandran" << std::endl;
		std::cout << "solvation" << std::endl;
		std::cout << "born_placeholder" << std::endl;
		std::cout << "pair_energies" << std::endl;
		std::cout << "env_stats" << std::endl;
    std::cout << "arom_hbonds" << std::endl; // SJF
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	avg_dis = false;
	hbond_stats = false;
	avg_E = false;
	phipsichi = false;
	dun_rot = false;
	pdbstat_dotcutoff = false;
	rep_atm_pair = false;
	phe_planes = false;
	cation_pi = false;
	carbonyl_o = false;
	rama_stats = false;
	solv_stats = false;
	born_placeholder = false;
	sasa_stats = false;
	pair_energies = false;
	env_stats = false;
	if ( whatkind == "distances" ) {
		avg_dis = true;
	} else if ( whatkind == "hbonds" ) {
		hbond_stats = true;
	} else if ( whatkind == "avg_energies" ) {
		avg_E = true;
	} else if ( whatkind == "phipsichi" ) {
		phipsichi = true;
	} else if ( whatkind == "dun_rot" ) {
		dun_rot = true;
	} else if ( whatkind == "dot_cutoff" ) {
		pdbstat_dotcutoff = true;
	} else if ( whatkind == "rep_atm_pair" ) {
		rep_atm_pair = true;
	} else if ( whatkind == "planes" ) {
		phe_planes = true;
	} else if ( whatkind == "cation" ) {
		cation_pi = true;
	} else if ( whatkind == "carbonyl" ) {
		carbonyl_o = true;
	} else if ( whatkind == "ramachandran" ) {
		rama_stats = true;
	} else if ( whatkind == "solvation" ) {
		solv_stats = true;
	} else if ( whatkind == "born_placeholder" ) {
		born_placeholder = true;
	} else if ( whatkind == "sasa_stats" ) {
		sasa_stats = true;
	} else if ( whatkind == "pair_energies" ) {
		pair_energies = true;
	} else if ( whatkind == "env_stats" ){
		env_stats = true;
		input_fa = false;
  } else if ( whatkind == "arom_hbonds" ) {
    arom_hbonds = true;
	} else {
		std::cout << "your pdbstats option was not recognized" << std::endl;
		std::cout << "the current options are:" << std::endl;
		std::cout << "distances" << std::endl;
		std::cout << "hbonds" << std::endl;
		std::cout << "avg_energies" << std::endl;
		std::cout << "phipsichi" << std::endl;
		std::cout << "dun_rot" << std::endl;
		std::cout << "sasa_stats" << std::endl;
		std::cout << "pair_energies" << std::endl;
		std::cout << "env_stats" << std::endl;
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///jk options related to docking, independent of mode
///jk BEWARE: changes scorefxn, atomvdw, etc.
///jk to setup minimal amount for a multi_chain system (gaining
///jk access to most docking subroutines), use enable_multi_chain
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_docking_options()
{

	using namespace design;
	using namespace dock_fab;
	using namespace docking;
	using namespace files_paths;

	docking_flag = true;

//     packing/fullatom
	active_rotamer_options.extrachi_cutoff = 0;
	set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score4d","score10d");

	enable_multi_chain();
	// jk return flags below to docking defaults
	// jk (changed in enable_multi_chain)
	docking_fullatom_flag = false;
	norepack1 = false;
	norepack2 = false;
	docking_score_norepack = false;
	dock_loop_ensemble_flag = false; // aroop: fragment insert in dock
	dock_loop_ensemble_loops_flag = false; // aroop: loop building in dock

//     details
	fab1 = truefalseoption("fab1"); // is either partner an antibody?
	fab2 = truefalseoption("fab2");

	norepack1 = truefalseoption("norepack1"); // should anything not be packed?
	norepack2 = truefalseoption("norepack2");
	if ( norepack1 ) std::cout << "Not repacking partner 1" << std::endl;
	if ( norepack2 ) std::cout << "Not repacking partner 2" << std::endl;

	//ia pose symmetry options
	docking_pose_symmetry = truefalseoption("docking_pose_symmetry");
	if ( docking_pose_symmetry ) {
		docking_pose_symm_full = truefalseoption("docking_pose_symm_full");
		intafteroption("pose_symm_n_monomers",2,pose_symm_n_monomers);
		stringafteroption("symm_type", "cn", symm_type);
		std::cout << " Pose symmetry mode with " << pose_symm_n_monomers << " monomers"
			<< std::endl;
		docking_pose_symm_subsystem = truefalseoption("docking_pose_symm_subsystem");
		docking_pose_symm_loops =  truefalseoption("docking_pose_symm_loops");
		docking_pose_symm_looprlx = truefalseoption("docking_pose_symm_looprlx");
	}

//	std::cout << "fun: n_monomer " << SS( n_monomers ) << std::endl;
//chu   if use unbound_backbone as starting structure
	unbound_start = truefalseoption("unbound_start");
//chu   if include native rotamers in unbound structures when packing
	unboundrot = truefalseoption("unboundrot");
	dock_rtmin = truefalseoption("dock_rtmin");

	// aroop docking with backbone flexibility
	dock_loop_ensemble_flag = truefalseoption("dock_flex");
	dock_loop_ensemble_loops_flag = truefalseoption("dock_loops");

//     scoring options/protocols
	docking_jumpout_mode         = truefalseoption("jumpout");
	docking_fullatom_flag        = truefalseoption("dockFA");
	docking_mcm_flag             = truefalseoption("dock_mcm");
	dock_mcm_loop_min_flag       = truefalseoption("dock_mcm_loopmin");
	docking_minimize_flag        = truefalseoption("dock_min") ||
	 docking_mcm_flag;
	docking_fullatom_search_flag = truefalseoption("FAsearch");
	prepack_full                 = truefalseoption("prepack_full");
	prepack_rtmin                = truefalseoption("prepack_rtmin");
	prepack_mode                 = prepack_full || prepack_rtmin;
	if ( truefalseoption("dockFA") ) input_fa = true;
	docking_score_norepack       = truefalseoption("dock_score_norepack");
	set_docking_interface_pack(!prepack_mode);
	docking_local_refine         = truefalseoption("docking_local_refine");
	docking_loop_trials          = truefalseoption("loop_trials");

	// pose docking
	if ( truefalseoption("pose") ) {
		pose_docking_flag = true;
		score_only = truefalseoption("score_only");
		minimize_backbone= truefalseoption("minimize");
		if (minimize_backbone) flexbb_docking_flag = true;
	}

	// loop modelling if necessary
	if ( truefalseoption("loop") ) {
		get_loop_options();
		if ( truefalseoption("minimize_loop")) {
			require_frags = false;
			use_fasta = false;
			minimize_backbone = true;
		}
		set_pose_loop_flag( pose_docking_flag );
		flexbb_docking_flag = true;
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_dockmode_options()
{
	using namespace dock_pivot;
	using namespace docking;
	using namespace files_paths;

	float scorefilter_val;

	docking_flag          = true;
	docking_fullatom_flag = false;

	idealized_structure = false;
	require_frags    = false;
	require_start    = true;
	use_fasta        = false;
	default_nstruct  = 1;
	use_constraints  = false;
	input_fa         = false;
	repack_input     = false; // input is repacked in special docking sub

	get_docking_options();

	//car setup default docking filtering-- only docking filters
	disable_all_filters();
	disable_filters = false;
	use_filter(dock_type) = true;

	docking_silent_input = truefalseoption("docking_silent_input");

	if ( (prepack_mode && !pose_docking_flag)
		|| (score_only && !docking_silent_input) )
		output_coord = false;

	// start conditions (docking flags)
	docking_randomize1     = truefalseoption("randomize1");
	docking_randomize2     = truefalseoption("randomize2");

	docking_axis_spin_flag = truefalseoption("spin") || docking_randomize1 ||
		docking_randomize2;

	docking_small_perturbation = truefalseoption("dock_pert") ||
		(!docking_randomize1 && !docking_randomize2);

	if ( docking_small_perturbation ) real3afteroption("dock_pert",3.0,
		normal_perturbation,8.0,parallel_perturbation,8.0,rotational_perturbation);

	// post-filter
	realafteroption("scorefilter",9999.,scorefilter_val);
	set_score_filter(scorefilter_val);
	realafteroption("I_sc_filter",5.0,scorefilter_val);
	set_docking_interf_energy_filter( scorefilter_val );
	realafteroption("chbrk_filter",1.0,scorefilter_val);
	set_chainbreak_score_filter( scorefilter_val );
	realafteroption("score_delta0",500.0, delta_before_mcm);
	realafteroption("score_delta1", 10.0, delta_after_one_min);
	realafteroption("score_delta5",  5.0, delta_after_five_mcm);
	smart_scorefilter_flag = truefalseoption("smart_scorefilter");
	realafteroption("smart_scorefilter",0.1,scorefilter_fraction);

	// fullatom modes
	if ( prepack_mode ) { // prepacking protocol to make start
		docking_fullatom_flag = true;
		input_fa              = true;
		output_fa             = true;
		disable_filters       = true;
		set_docking_interface_pack(false); // repack everything
	} else if ( docking_fullatom_flag || docking_fullatom_search_flag ||
		docking_minimize_flag || docking_mcm_flag) { // fullatom docking
		docking_fullatom_flag = true;
		input_fa              = true;
		output_fa             = true;
		disable_filters       = score_only;
		set_docking_interface_pack(true); // repack only the interface
	}

// pivot protocol for special flexible backbone case
	using_pivot_residues   = truefalseoption("pivot_residues");
	pivot_residue_filename = "PIVOT_FILE";
	if ( using_pivot_residues ) {
		stringafteroption( "pivot_residues", pivot_residue_filename,
			pivot_residue_filename );
		flexbb_docking_flag = true;
	}

}


///////////////////////////////////////////////////////////////////////////////
///
/// @brief options for docking of a flexible peptide onto a globular protein
///        Note: Perhaps a lot of these flags should be obsolete, examine this later
///        Note2: This is the place to add future flags to FlexPepDock
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author barak
///
/////////////////////////////////////////////////////////////////////////////////
void
get_flexpepdock_options(){
	// TODO: rethink about the need for all of these flags
	//using namespace dock_pivot;
	using namespace docking;
	using namespace files_paths;
	using namespace design;
	using namespace dock_fab;

	docking_flag          = true;
	docking_fullatom_flag = false;

	idealized_structure = false;
	require_frags    = false;
	require_start    = true;
	use_fasta        = false;
	default_nstruct  = 1;
	use_constraints  = false;
	input_fa         = true; // default of whether to read coordinates of side-chain atoms
	repack_input     = false; // input is repacked in special docking sub

	// TODO: these are options borrowed from get_docking_options():
	// ============================================================


//     packing/fullatom
	active_rotamer_options.extrachi_cutoff = 0; // only positions with [GE (?equal?) cutoff] neighbours, get extra rotamers
	set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score4d","score10d");

	enable_multi_chain();
	// jk return flags below to docking defaults
	// jk (changed in enable_multi_chain)
	norepack1 = false;
	norepack2 = false;
	docking_score_norepack = false;
	dock_loop_ensemble_flag = false; // aroop: fragment insert in dock
	dock_loop_ensemble_loops_flag = false; // aroop: loop building in dock

//     details
	fab1 = truefalseoption("fab1"); // is either partner an antibody?
	fab2 = truefalseoption("fab2");

	norepack1 = truefalseoption("norepack1"); // should anything not be packed?
	norepack2 = truefalseoption("norepack2");
	if ( norepack1 ) std::cout << "Not repacking partner 1" << std::endl;
	if ( norepack2 ) std::cout << "Not repacking partner 2" << std::endl;

//ora symmetry_mode: symmetry defined by rotation around z axis and translation
//ora along x-axis by size "init_docking_T_size"
//	docking_symmetry       = truefalseoption("symmetry");
//	docking_symmetry_save  = docking_symmetry;
	// allows to switch symmetry temporarily off
	//	dock_sym_multimer_start = truefalseoption("multimer");
	// indicates multimer has been read in and is to
	// be used to derive symmetry frame


//	std::cout << "fun: n_monomer " << SS( n_monomers ) << std::endl;
//chu   if use unbound_backbone as starting structure
	unbound_start = truefalseoption("unbound_start");
//chu   if include native rotamers in unbound structures when packing
	unboundrot = truefalseoption("unboundrot");
	dock_rtmin = truefalseoption("dock_rtmin");

	// aroop docking with backbone flexibility
	dock_loop_ensemble_flag = truefalseoption("dock_flex");
	dock_loop_ensemble_loops_flag = truefalseoption("dock_loops");

//     scoring options/protocols
	docking_jumpout_mode         = truefalseoption("jumpout");
	docking_fullatom_flag        = truefalseoption("flexpepdock_FA");
	docking_mcm_flag             = truefalseoption("dock_mcm");
	dock_mcm_loop_min_flag       = truefalseoption("dock_mcm_loopmin");
	docking_minimize_flag        = truefalseoption("dock_min") ||
	 docking_mcm_flag;
	docking_fullatom_search_flag = truefalseoption("FAsearch");
	prepack_full                 = truefalseoption("prepack_full");
	prepack_rtmin                = truefalseoption("prepack_rtmin");
	prepack_mode                 = prepack_full || prepack_rtmin;
	if ( truefalseoption("dockFA") ) input_fa = true;
	docking_score_norepack       = truefalseoption("dock_score_norepack");
	set_docking_interface_pack(!prepack_mode);
	docking_local_refine         = truefalseoption("docking_local_refine");
	docking_loop_trials          = truefalseoption("loop_trials");

	// pose docking
	if ( truefalseoption("pose") ) {
		pose_docking_flag = true;
		score_only = truefalseoption("score_only");
		minimize_backbone= truefalseoption("minimize");
		if (minimize_backbone) flexbb_docking_flag = true;
	}

	// loop modelling if necessary
	if ( truefalseoption("loop") ) {
		get_loop_options();
		set_pose_loop_flag( pose_docking_flag );
		flexbb_docking_flag = true;
	}
	// =======================================================
	// end of options from get_docking_options()

	//car setup default docking filtering-- only docking filters
	disable_all_filters();
	disable_filters = false;
	use_filter(dock_type) = true;

//	something that Chu wrote, see later if we need this
//	docking_silent_input = truefalseoption("docking_silent_input");

	if ( (prepack_mode && !pose_docking_flag)
		|| (score_only && !docking_silent_input) )
		output_coord = false;

	// start conditions (docking flags)
//	docking_randomize1     = truefalseoption("randomize1");
//	docking_randomize2     = truefalseoption("randomize2");

//	docking_axis_spin_flag = truefalseoption("spin") || docking_randomize1 ||
//		docking_randomize2;

//	docking_small_perturbation = truefalseoption("dock_pert") ||
//		(!docking_randomize1 && !docking_randomize2);

//	if ( docking_small_perturbation ) real3afteroption("dock_pert",3.0,
//		normal_perturbation,8.0,parallel_perturbation,8.0,rotational_perturbation);

	// post-filter
//	realafteroption("scorefilter",9999.,scorefilter_val);
//	set_score_filter(scorefilter_val);
//	realafteroption("I_sc_filter",5.0,scorefilter_val);
//	set_docking_interf_energy_filter( scorefilter_val );
//	realafteroption("chbrk_filter",1.0,scorefilter_val);
//	set_chainbreak_score_filter( scorefilter_val );
//	realafteroption("score_delta0",500.0, delta_before_mcm);
//	realafteroption("score_delta1", 10.0, delta_after_one_min);
//	realafteroption("score_delta5",  5.0, delta_after_five_mcm);
//	smart_scorefilter_flag = truefalseoption("smart_scorefilter");
//	realafteroption("smart_scorefilter",0.1,scorefilter_fraction);

	// fullatom modes
	if ( prepack_mode ) { // prepacking protocol to make start
		docking_fullatom_flag = true;
		input_fa              = true;
		output_fa             = true;
		disable_filters       = true;
		set_docking_interface_pack(false); // repack everything
	} else if ( docking_fullatom_flag || docking_fullatom_search_flag ||
		docking_minimize_flag || docking_mcm_flag) { // fullatom docking
		docking_fullatom_flag = true;
		input_fa              = true;
		output_fa             = true;
		disable_filters       = score_only;
		set_docking_interface_pack(true); // repack only the interface
	}

	// we want to work with pose mode only
	set_pose_flag(true);


//     packing/fullatom
	active_rotamer_options.extrachi_cutoff = 0;
	set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score4d","score10d");

	enable_multi_chain();
	// jk return flags below to docking defaults
	// jk (changed in enable_multi_chain)
	docking_fullatom_flag = false;
	norepack1 = false;
	norepack2 = false;
	docking_score_norepack = false;
	dock_loop_ensemble_flag = false; // aroop: fragment insert in dock
	dock_loop_ensemble_loops_flag = false; // aroop: loop building in dock

//     details
	fab1 = truefalseoption("fab1"); // is either partner an antibody?
	fab2 = truefalseoption("fab2");

	norepack1 = truefalseoption("norepack1"); // should anything not be packed?
	norepack2 = truefalseoption("norepack2");
	if ( norepack1 ) std::cout << "Not repacking partner 1" << std::endl;
	if ( norepack2 ) std::cout << "Not repacking partner 2" << std::endl;

	//ia pose symmetry options
	docking_pose_symmetry = truefalseoption("docking_pose_symmetry");
	if ( docking_pose_symmetry ) {
		docking_pose_symm_full = truefalseoption("docking_pose_symm_full");
		intafteroption("pose_symm_n_monomers",2,pose_symm_n_monomers);
		stringafteroption("symm_type", "cn", symm_type);
		std::cout << " Pose symmetry mode with " << pose_symm_n_monomers << " monomers"
			<< std::endl;
		docking_pose_symm_subsystem = truefalseoption("docking_pose_symm_subsystem");
		docking_pose_symm_loops =  truefalseoption("docking_pose_symm_loops");
		docking_pose_symm_looprlx = truefalseoption("docking_pose_symm_looprlx");
	}

//	std::cout << "fun: n_monomer " << SS( n_monomers ) << std::endl;
//chu   if use unbound_backbone as starting structure
	unbound_start = truefalseoption("unbound_start");
//chu   if include native rotamers in unbound structures when packing
	unboundrot = truefalseoption("unboundrot");
	dock_rtmin = truefalseoption("dock_rtmin");

	// aroop docking with backbone flexibility
	dock_loop_ensemble_flag = truefalseoption("dock_flex");
	dock_loop_ensemble_loops_flag = truefalseoption("dock_loops");

//     scoring options/protocols
	docking_jumpout_mode         = truefalseoption("jumpout");
	docking_fullatom_flag        = truefalseoption("flexpepdock_FA");
	docking_mcm_flag             = truefalseoption("dock_mcm");
	dock_mcm_loop_min_flag       = truefalseoption("dock_mcm_loopmin");
	docking_minimize_flag        = truefalseoption("dock_min") ||
	 docking_mcm_flag;
	docking_fullatom_search_flag = truefalseoption("FAsearch");
	prepack_full                 = truefalseoption("prepack_full");
	prepack_rtmin                = truefalseoption("prepack_rtmin");
	prepack_mode                 = prepack_full || prepack_rtmin;
	if ( truefalseoption("flexpepdock_FA") ) input_fa = true;
	docking_score_norepack       = truefalseoption("dock_score_norepack");
	set_docking_interface_pack(!prepack_mode);
	docking_local_refine         = truefalseoption("docking_local_refine");
	docking_loop_trials          = truefalseoption("loop_trials");

	// pose docking
	if ( truefalseoption("pose") ) {
		pose_docking_flag = true;
		score_only = truefalseoption("score_only");
		minimize_backbone= truefalseoption("minimize");
		if (minimize_backbone) flexbb_docking_flag = true;
	}

	// loop modelling if necessary
	if ( truefalseoption("loop") ) {
		get_loop_options();
		if ( truefalseoption("minimize_loop")) {
			require_frags = false;
			use_fasta = false;
			minimize_backbone = true;
		}
		set_pose_loop_flag( pose_docking_flag );
		flexbb_docking_flag = true;
	}


}


///////////////////////////////////////////////////////////////////////////////
///
/// @brief options for conformation pathway predictions protocol
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author Barak Raveh & Angela Enosh
///
/// @created 29/11/2007
///
/////////////////////////////////////////////////////////////////////////////////
void
get_pathways_options(){
	// TODO: rethink about the need for all of these flags
	//using namespace dock_pivot;
	using namespace docking;
	using namespace files_paths;
	using namespace design;
	using namespace dock_fab;

	docking_flag          = true;
	docking_fullatom_flag = false; // default

	idealized_structure = false;
	require_frags    = false;
	require_start    = true;
	use_fasta        = false;
	default_nstruct  = 1;
	use_constraints  = false;
	input_fa         = true; // default of whether to read coordinates of side-chain atoms
	repack_input     = false; // input is repacked in special docking sub

	// TODO: these are options borrowed from get_docking_options():
	// ============================================================


//     packing/fullatom
	active_rotamer_options.extrachi_cutoff = 0; // only positions with [GE (?equal?) cutoff] neighbours, get extra rotamers
	set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score4d","score10d");

	enable_multi_chain();
	// jk return flags below to docking defaults
	// jk (changed in enable_multi_chain)
	norepack1 = false;
	norepack2 = false;
	docking_score_norepack = false;
	dock_loop_ensemble_flag = false; // aroop: fragment insert in dock
	dock_loop_ensemble_loops_flag = false; // aroop: loop building in dock

//     details
	fab1 = truefalseoption("fab1"); // is either partner an antibody?
	fab2 = truefalseoption("fab2");

	norepack1 = truefalseoption("norepack1"); // should anything not be packed?
	norepack2 = truefalseoption("norepack2");
	if ( norepack1 ) std::cout << "Not repacking partner 1" << std::endl;
	if ( norepack2 ) std::cout << "Not repacking partner 2" << std::endl;

//ora symmetry_mode: symmetry defined by rotation around z axis and translation
//ora along x-axis by size "init_docking_T_size"
//	docking_symmetry       = truefalseoption("symmetry");
//	docking_symmetry_save  = docking_symmetry;
	// allows to switch symmetry temporarily off
	//	dock_sym_multimer_start = truefalseoption("multimer");
	// indicates multimer has been read in and is to
	// be used to derive symmetry frame


//	std::cout << "fun: n_monomer " << SS( n_monomers ) << std::endl;
//chu   if use unbound_backbone as starting structure
	unbound_start = truefalseoption("unbound_start");
//chu   if include native rotamers in unbound structures when packing
	unboundrot = truefalseoption("unboundrot");
	dock_rtmin = truefalseoption("dock_rtmin");

	// aroop docking with backbone flexibility
	dock_loop_ensemble_flag = truefalseoption("dock_flex");
	dock_loop_ensemble_loops_flag = truefalseoption("dock_loops");

//     scoring options/protocols
	docking_jumpout_mode         = truefalseoption("jumpout");
	docking_fullatom_flag        = truefalseoption("pathways_FA");
	docking_mcm_flag             = truefalseoption("dock_mcm");
	dock_mcm_loop_min_flag       = truefalseoption("dock_mcm_loopmin");
	docking_minimize_flag        = truefalseoption("dock_min") ||
	 docking_mcm_flag;
	docking_fullatom_search_flag = truefalseoption("FAsearch");
	prepack_full                 = truefalseoption("prepack_full");
	prepack_rtmin                = truefalseoption("prepack_rtmin");
	prepack_mode                 = prepack_full || prepack_rtmin;
	if ( truefalseoption("pathways_FA") ) input_fa = true;
	docking_score_norepack       = truefalseoption("dock_score_norepack");
	set_docking_interface_pack(!prepack_mode);
	docking_local_refine         = truefalseoption("docking_local_refine");
	docking_loop_trials          = truefalseoption("loop_trials");

	// pose docking
	if ( truefalseoption("pose") ) {
		pose_docking_flag = true;
		score_only = truefalseoption("score_only");
		minimize_backbone= truefalseoption("minimize");
		if (minimize_backbone) flexbb_docking_flag = true;
	}

	// loop modelling if necessary
	if ( truefalseoption("loop") ) {
		get_loop_options();
		set_pose_loop_flag( pose_docking_flag );
		flexbb_docking_flag = true;
	}
	// =======================================================
	// end of options from get_docking_options()

	//car setup default docking filtering-- only docking filters
	disable_all_filters();
	disable_filters = false;
	use_filter(dock_type) = true;

//	something that Chu wrote, see later if we need this
//	docking_silent_input = truefalseoption("docking_silent_input");

	if ( (prepack_mode && !pose_docking_flag)
		|| (score_only && !docking_silent_input) )
		output_coord = false;

	// start conditions (docking flags)
//	docking_randomize1     = truefalseoption("randomize1");
//	docking_randomize2     = truefalseoption("randomize2");

//	docking_axis_spin_flag = truefalseoption("spin") || docking_randomize1 ||
//		docking_randomize2;

//	docking_small_perturbation = truefalseoption("dock_pert") ||
//		(!docking_randomize1 && !docking_randomize2);

//	if ( docking_small_perturbation ) real3afteroption("dock_pert",3.0,
//		normal_perturbation,8.0,parallel_perturbation,8.0,rotational_perturbation);

	// post-filter
//	realafteroption("scorefilter",9999.,scorefilter_val);
//	set_score_filter(scorefilter_val);
//	realafteroption("I_sc_filter",5.0,scorefilter_val);
//	set_docking_interf_energy_filter( scorefilter_val );
//	realafteroption("chbrk_filter",1.0,scorefilter_val);
//	set_chainbreak_score_filter( scorefilter_val );
//	realafteroption("score_delta0",500.0, delta_before_mcm);
//	realafteroption("score_delta1", 10.0, delta_after_one_min);
//	realafteroption("score_delta5",  5.0, delta_after_five_mcm);
//	smart_scorefilter_flag = truefalseoption("smart_scorefilter");
//	realafteroption("smart_scorefilter",0.1,scorefilter_fraction);

	// fullatom modes
	if ( prepack_mode ) { // prepacking protocol to make start
		docking_fullatom_flag = true;
		input_fa              = true;
		output_fa             = true;
		disable_filters       = true;
		set_docking_interface_pack(false); // repack everything
	} else if ( docking_fullatom_flag || docking_fullatom_search_flag ||
		docking_minimize_flag || docking_mcm_flag) { // fullatom docking
		docking_fullatom_flag = true;
		input_fa              = true;
		output_fa             = true;
		disable_filters       = score_only;
		set_docking_interface_pack(true); // repack only the interface
	}

	// we want to work with pose mode only
	set_pose_flag(true);


//     packing/fullatom
	active_rotamer_options.extrachi_cutoff = 0;
	set_default_atomvdw( "hybrid" );

	select_default_scorefxns("score4d","score10d");

	enable_multi_chain();
	// jk return flags below to docking defaults
	// jk (changed in enable_multi_chain)
	docking_fullatom_flag = false;
	norepack1 = false;
	norepack2 = false;
	docking_score_norepack = false;
	dock_loop_ensemble_flag = false; // aroop: fragment insert in dock
	dock_loop_ensemble_loops_flag = false; // aroop: loop building in dock

//     details
	fab1 = truefalseoption("fab1"); // is either partner an antibody?
	fab2 = truefalseoption("fab2");

	norepack1 = truefalseoption("norepack1"); // should anything not be packed?
	norepack2 = truefalseoption("norepack2");
	if ( norepack1 ) std::cout << "Not repacking partner 1" << std::endl;
	if ( norepack2 ) std::cout << "Not repacking partner 2" << std::endl;

	//ia pose symmetry options
	docking_pose_symmetry = truefalseoption("docking_pose_symmetry");
	if ( docking_pose_symmetry ) {
		docking_pose_symm_full = truefalseoption("docking_pose_symm_full");
		intafteroption("pose_symm_n_monomers",2,pose_symm_n_monomers);
		stringafteroption("symm_type", "cn", symm_type);
		std::cout << " Pose symmetry mode with " << pose_symm_n_monomers << " monomers"
			<< std::endl;
		docking_pose_symm_subsystem = truefalseoption("docking_pose_symm_subsystem");
		docking_pose_symm_loops =  truefalseoption("docking_pose_symm_loops");
		docking_pose_symm_looprlx = truefalseoption("docking_pose_symm_looprlx");
	}

//	std::cout << "fun: n_monomer " << SS( n_monomers ) << std::endl;
//chu   if use unbound_backbone as starting structure
	unbound_start = truefalseoption("unbound_start");
//chu   if include native rotamers in unbound structures when packing
	unboundrot = truefalseoption("unboundrot");
	dock_rtmin = truefalseoption("dock_rtmin");

	// aroop docking with backbone flexibility
	dock_loop_ensemble_flag = truefalseoption("dock_flex");
	dock_loop_ensemble_loops_flag = truefalseoption("dock_loops");

//     scoring options/protocols
	docking_jumpout_mode         = truefalseoption("jumpout");
	docking_fullatom_flag        = truefalseoption("pathways_FA");
	docking_mcm_flag             = truefalseoption("dock_mcm");
	dock_mcm_loop_min_flag       = truefalseoption("dock_mcm_loopmin");
	docking_minimize_flag        = truefalseoption("dock_min") ||
	 docking_mcm_flag;
	docking_fullatom_search_flag = truefalseoption("FAsearch");
	prepack_full                 = truefalseoption("prepack_full");
	prepack_rtmin                = truefalseoption("prepack_rtmin");
	prepack_mode                 = prepack_full || prepack_rtmin;
	if ( truefalseoption("pathways_FA") ) input_fa = true;
	docking_score_norepack       = truefalseoption("dock_score_norepack");
	set_docking_interface_pack(!prepack_mode);
	docking_local_refine         = truefalseoption("docking_local_refine");
	docking_loop_trials          = truefalseoption("loop_trials");

	// pose docking
	if ( truefalseoption("pose") ) {
		pose_docking_flag = true;
		score_only = truefalseoption("score_only");
		minimize_backbone= truefalseoption("minimize");
		if (minimize_backbone) flexbb_docking_flag = true;
	}

	// loop modelling if necessary
	if ( truefalseoption("loop") ) {
		get_loop_options();
		if ( truefalseoption("minimize_loop")) {
			require_frags = false;
			use_fasta = false;
			minimize_backbone = true;
		}
		set_pose_loop_flag( pose_docking_flag );
		flexbb_docking_flag = true;
	}

	require_frags = false;


}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
get_membrane_options()
{
	using namespace files_paths;

	require_start = false;
	input_fa  = true;
	output_fa = true;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///car higher run levels give increasing levels of detailed output
///car current keywords:
///car verbose > gush > yap > chat > inform > standard > quiet > silent
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
setup_runlevel()
{
	using namespace files_paths;
	using namespace runlevel_ns;

	intafteroption("run_level",standard,runlevel);

	if ( truefalseoption("silent") || truefalseoption("output_silent_gz")
		|| pose_silent_out ) {
		runlevel = silent;
	} else if ( truefalseoption("verbose") ) {
		runlevel = verbose;
	} else if ( truefalseoption("gush") ) {
		runlevel = gush;
	} else if ( truefalseoption("yap") ) {
		runlevel = yap;
	} else if ( truefalseoption("chat") ) {
		runlevel = chat;
	} else if ( truefalseoption("inform") ) {
		runlevel = inform;
	} else if ( truefalseoption("quiet") ) {
		runlevel = quiet;
	}

	std::cout << "run level:  " << runlevel << std::endl;
	if ( runlevel == silent ) {
		use_status = false;
		use_decoy_status = false;
		use_scorefile = false;
	}

	benchmark = truefalseoption("benchmark");
	debug = truefalseoption("debug");

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief select scorefxns to be used as the defaults
///
/// @details
///
/// @param[in] fxn   in -
/// @param[in] fa_fxn   in -
///
/// @global_read
///
/// @global_write
///
/// @remarks an invalide selection sets the scorefxn to undefined
///
/// @references
///
/// @author car 10/22/03
///
/////////////////////////////////////////////////////////////////////////////////
void
select_default_scorefxns(
	std::string const & fxn,
	std::string const & fa_fxn
)
{
	using namespace scorefxns;

	default_scorefxn = score_undef;
	default_fascorefxn = score_undef;

//car allowed centroid values
	if ( fxn == "score0" ) default_scorefxn = score_0;
	if ( fxn == "score1" ) default_scorefxn = score_1;
	if ( fxn == "score2" ) default_scorefxn = score_2;
	if ( fxn == "score3" ) default_scorefxn = score_3;
	if ( fxn == "score4" ) default_scorefxn = score_4;
	if ( fxn == "score5" ) default_scorefxn = score_5;
	if ( fxn == "score6" ) default_scorefxn = score_6;
	if ( fxn == "score7" ) default_scorefxn = score_7;
	if ( fxn == "score8" ) default_scorefxn = score_8;

	if ( fxn == "score8di" ) default_scorefxn = score_8di;

	if ( fxn == "score4Lc" ) default_scorefxn = score_4Lc;

	if ( fxn == "score4d" ) default_scorefxn = score_4d;

	if ( fxn == "score3L" ) default_scorefxn = score_3L;
	if ( fxn == "score4L" ) default_scorefxn = score_4L;

//car allowed fullatom values
	if ( fa_fxn == "score9" ) default_fascorefxn = score_9;
	if ( fa_fxn == "score10" ) default_fascorefxn = score_10;
	if ( fa_fxn == "score11" ) default_fascorefxn = score_11;
	if ( fa_fxn == "score12" ) default_fascorefxn = score_12;

	if ( fa_fxn == "score10d" ) default_fascorefxn = score_10d;

	if ( fa_fxn == "score9L" ) default_fascorefxn = score_9L;
	if ( fa_fxn == "score10L" ) default_fascorefxn = score_10L;
	if ( fa_fxn == "score11L" ) default_fascorefxn = score_11L;
	if ( fa_fxn == "score12L" ) default_fascorefxn = score_12L;
	if ( fa_fxn == "score13" ) default_fascorefxn = score_13;
	if ( fa_fxn == "score14" ) default_fascorefxn = score_14;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief  initializes default_scorefxn using values either set on
///        command line or by method in options.cc
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
setup_default_scorefxn( bool & output_fa )
{
	using namespace scorefxns;


//car set default scorefxns
	if ( truefalseoption("scorefxn") ) {
		//car local
		int value;
		intafteroption("scorefxn",score_undef,value);

//car set the default_scorefxn
		if ( output_fa ) {
			default_fascorefxn = value;
		} else {
			default_scorefxn = value;
		}

//car check to see if fullatom output/fullatom function agree
		score_set_default_function(output_fa);

		if ( output_fa ) {
			if ( default_fascorefxn == score_undef ) {
				std::cout << A( 22, " WARNING!  scorefxn: " ) <<
				 I( 4, value ) <<
				 A( 60, " is not an allowed scorefxn for fullatom output" ) <<
				 std::endl;
				std::cout << "Using default fullatom scorefxn instead" << std::endl;
			} else if ( !fa_scorefxn ) {
				std::cout << A( 22, " WARNING!  scorefxn: " ) <<
				 I( 4, value ) <<
				 A( 60, "is a centroid scorefxn and fullatom output requested" ) <<
				 std::endl;
			}
		} else if ( default_scorefxn == score_undef ) {
			std::cout << A( 22, " WARNING!  scorefxn: " ) <<
			 I( 4, value ) <<
			 A( 60, " is not an allowed scorefxn for centroid output" ) <<
			 std::endl;
			std::cout << "Turning on fullatom output" << std::endl;
			output_fa = true;
			default_fascorefxn = value;
		}
	}

	if ( default_scorefxn == score_undef ) default_scorefxn = score_4;
	if ( default_fascorefxn == score_undef ) default_fascorefxn = score_12;

	std::cout << "default centroid scorefxn: " << SS( default_scorefxn ) << std::endl;
	std::cout << "default fullatom scorefxn: " << SS( default_fascorefxn ) << std::endl;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///car check for and fix  inconsistent option settings
///
/// @details
///
/// @param  mode - [in/out]? -
/// @param  number_of_output - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
check_options(
	std::string const & mode,
	int & number_of_output
)
{

	using namespace files_paths;

	if ( !query_defined ) {
		if ( !require_start && !pose_flag() ) {
			std::cout << "ERROR:: starting structures must be given if " << std::endl;
			std::cout << "query is not defined by series_code, " << std::endl;
			std::cout << "protein_name,and chain_id" << std::endl;
			std::cout << "stopping..." << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		if ( require_frags ) {
			std::cout << "WARNING::" << std::endl;
			std::cout << "fragments cannot be used: query not defined " << std::endl;
			std::cout << "by series_code,protein_name,and chain_id" << std::endl;
			std::cout << "THIS MAY CAUSE SERIOUS PROBLEMS!!" << std::endl;
			require_frags = false;
		}
		if ( use_constraints ) {
			std::cout << "WARNING::" << std::endl;
			std::cout << "constraints will not be used: query not defined" << std::endl;
			std::cout << "by series_code,protein_name,and chain_id" << std::endl;
			use_constraints = false;
		}
		if ( use_fasta ) {
			std::cout << "WARNING::" << std::endl;
			std::cout << "fasta/dat cannot be used: query not defined " << std::endl;
			std::cout << "by series_code,protein_name,and chain_id" << std::endl;
			std::cout << "sequence will be taken from starting pdb" << std::endl;
			use_fasta = false;
		}
	}

	if ( silent_input ) {
		if ( !idealized_structure ) {
			std::cout << "ERROR::  silent_input T, idealized_structure F" << std::endl;
			std::cout << "         silent input not compatible with" <<
			 "non-ideal structures" << std::endl;
			utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		}
		if ( !refold_input ) {
			std::cout << "WARNING:: silent_input T,  refold_input F" << std::endl;
			std::cout << " silent input structures must be refolded" << std::endl;
			std::cout << " setting refold_input T" << std::endl;
			refold_input = true;
		}
	}

	if ( refold_input && input_fa ) {
		std::cout << "WARNING:: refold_input T, input_fa T" << std::endl;
		std::cout << " refolding of fullatom input not suppported" << std::endl;
		std::cout << " setting refold_input F" << std::endl;
		refold_input = false;
	}

	if ( refold_input && allow_missing ) {
		std::cout << "WARNING::refold_input T, allow_missing T" << std::endl;
		std::cout << "         setting refold_input F" << std::endl;
		refold_input = false;
	}
	// chu DO NOT CHANGE THE ORDER OF THE FOLLOWING THREE CHECKS
	// -accept_all: even if the decoy does not pass the filter, we accept it,
	//              increase decoy counter and output it. So imply filters is on
	//              (no_filters is false ) and output_all is true.
	//
	// -no_filters: do not use filters and all decoys will be accepted and output
  //
	// -output_all: for decoys which fail to pass filters, output them but do not
	//              increase the decoy counter.
	// override order: accept_all > no_filters > output_all.
	if ( disable_filters && accept_all ) {
		std::cout << "WARNING:: accept_all T, no_filters T" << std::endl;
		std::cout << "if no filter is used, all structures will be accepted. "
							<< "so you probably want to use filters !" << std::endl;
		std::cout << "setting no_filter F" << std::endl;
		disable_filters = false;
	}

	if ( accept_all && ! output_all ) {
		std::cout << "WARNING:: accept_all T, output_all F" << std::endl;
		std::cout << "if all structures are accepted, they should be all output!"
							<< std::endl;
		std::cout << "setting output_all T" << std::endl;
		output_all = true;
	}

	if ( disable_filters && output_all ) {
		std::cout << "WARNING:: disable_filters T, output_all T" << std::endl;
		std::cout << "output_all option is not applicable because filters " <<
		 "are not in use" << std::endl;
		output_all = false;
	}
	// chu END

	if ( !output_coord && !use_scorefile && !disable_filters ) {
		std::cout << "WARNING::" << std::endl;
		std::cout << "no coordinate or scorefile output requested" << std::endl;
		std::cout << "filters are irrelevant and are being disabled" << std::endl;
	}

	if ( output_all && ! output_coord ) {
		std::cout << "WARNING::" << std::endl;
		std::cout << "output_all option selected" << std::endl;
		std::cout << "turning on coordinate output" << std::endl;
		output_coord = true;
	}

	if ( number_of_output > 0 && !output_coord ) {
		std::cout << "WARNING:: nstruct > 0 but output_coord F" << std::endl;
		if ( mode == "score" ) {
			output_coord = true;
			std::cout << "turning on coordinate output" << std::endl;
		} else {
			std::cout << "no pdb files will be output" << std::endl;
		}
	}

	if ( allow_missing && idealized_structure ) {
		std::cout << "WARNING:: allow_missing T, idealized_structure T" << std::endl;
		std::cout << "missing backbone atoms not allowed in idealized structures" << std::endl;
		std::cout << "setting allow_missing F" << std::endl;
		allow_missing = false;
	}

	if ( !require_start ) {
		if ( input_fa ) {
			std::cout << "WARNING:: input_fa T, require_start F" << std::endl;
			std::cout << "fullatom input not allowed when no starting structures are being input " << std::endl;
			std::cout << "setting input_fa F" << std::endl;
			input_fa = false;
		}
		if ( repack_input ) {
			std::cout << "WARNING:: repack_input T, require_start F" << std::endl;
			std::cout << "fullatom input not allowed when no starting structures are being input " << std::endl;
			std::cout << "setting repack_input F" << std::endl;
			repack_input = false;
		}
		if ( include_inputchi ) {
			std::cout << "WARNING:: include_inputchi T, require_start F" << std::endl;
			std::cout << "fullatom input not allowed when no starting structures are being input " << std::endl;
			std::cout << "setting include_inputchi F" << std::endl;
			include_inputchi = false;
		}
	}

	if ( !input_fa ) {
		if ( repack_input) {
			std::cout << "WARNING:: repack_input T, input_fa F" << std::endl;
			std::cout << "repacking only allowed for fullatom input" << std::endl;
			std::cout << "setting repack_input F" << std::endl;
			repack_input = false;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// jk setup to use (fullatom) docking flags in modes other than docking
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
enable_multi_chain()
{

	using namespace docking;
	using namespace files_paths;

	// jk This should be the only way to turn on multi_chain
	// jk If we've already been through here, don't do it again
	if ( multi_chain ) return;

	use_pdb_numbering = true;
	read_all_chains = true;
	refold_input = false;
	n_monomers = 2;
	norepack1 = true;
	norepack2 = true;
	docking_score_norepack = true;

	docking_fake_native = truefalseoption("fake_native"); // is native real?

//bs   option to output stats on interfacial hbonds
	docking_hb_stats = truefalseoption("docking_hb_stats");

//jk   option to suppress writing interface centroid in output PDBs
	docking_output_position_hetatm = truefalseoption("output_position_hetatm");

	docking_fullatom_flag = true;

	multi_chain = true;
	return;

}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///
/// @brief
///  rosetta will read in fragments, sequence, and native from
///   files with prefix specified after protein_name_prefix_homolog and frags_name_prefix_homolog.
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
/// rhiju
///////////////////////////////////////////////////////////////////////////////
void adjust_options_protein_prefix_homolog(){
	using namespace files_paths;
	using namespace basic::options;

	stringafteroption( "protein_name_prefix_homolog","homolog_",protein_name_prefix);
	stringafteroption( "frags_name_prefix_homolog", "homolog_", frags_name_prefix );
	map_start_sequence = true;
}

/////////////////////////////////////////////////////////////////////////////
///
/// @brief
///  rosetta will read in fragments, sequence, and native from
///   files with prefix specified after protein_name_prefix_query and frags_name_prefix_query.
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
/// rhiju
///////////////////////////////////////////////////////////////////////////////
void adjust_options_protein_prefix_query(){
	using namespace files_paths;
	using namespace basic::options;

	stringafteroption( "protein_name_prefix_query","query_",protein_name_prefix);
	stringafteroption( "frags_name_prefix_query", "query_", frags_name_prefix );
}

void initialize_score_reweights(){
	using namespace scorefxns::score_reweights;

	static bool init = {false};

	if (!init){
		realafteroption("vdw_reweight", 1.0, vdw_reweight ) ;
		realafteroption("env_reweight", 1.0, env_reweight ) ;
		realafteroption("pair_reweight", 1.0, pair_reweight ) ;
		realafteroption("cb_reweight", 1.0, cb_reweight ) ;
		realafteroption("sheet_reweight", 1.0, sheet_reweight ) ;
		realafteroption("ss_reweight", 1.0, ss_reweight ) ;
		realafteroption("hs_reweight", 1.0, hs_reweight ) ;
		realafteroption("rsigma_reweight", 1.0, rsigma_reweight ) ;
		realafteroption("rg_reweight", 1.0, rg_reweight ) ;
		realafteroption("pc_reweight", 1.0, pc_reweight ) ;
		realafteroption("fa_atr_reweight", 1.0, fa_atr_reweight ) ;
		realafteroption("fa_rep_reweight", 1.0, fa_rep_reweight ) ;
		realafteroption("fa_dun_reweight ", 1.0, fa_dun_reweight  ) ;
		realafteroption("fa_pair_reweight", 1.0, fa_pair_reweight ) ;
		realafteroption("fa_plane_reweight", 1.0, fa_plane_reweight ) ;
		realafteroption("fa_solv_reweight", 1.0, fa_solv_reweight ) ;
		realafteroption("fa_ref_reweight ", 1.0, fa_ref_reweight  ) ;
		realafteroption("fa_pH_reweight", 1.0, fa_pH_reweight ) ;
		realafteroption("fa_h2o_reweight", 1.0, fa_h2o_reweight ) ;
		realafteroption("fa_prob1b_reweight", 1.0, fa_prob1b_reweight ) ;
		realafteroption("fa_gb_elec_reweigh", 1.0, fa_gb_elec_reweight ) ;
		realafteroption("hb_srbb_reweight", 1.0, hb_srbb_reweight ) ;
		realafteroption("hb_lrbb_reweight", 1.0, hb_lrbb_reweight ) ;
		realafteroption("hb_sc_reweight", 1.0, hb_sc_reweight ) ;
		realafteroption("chainbreak_reweight", 1.0, chainbreak_reweight ) ;

		realafteroption("electron_density_reweight", 1.0, electron_density_reweight ) ;
		realafteroption("dummy_model_reweight", 1.0, dummy_model_reweight ) ;
		realafteroption("saxs_model_reweight", 1.0, saxs_model_reweight ) ;
		param_pack::pack_wts.set_Wplane_total(realafteroption( "Wplane_total", 0.0));

		realafteroption("barcode_reweight", 1.0, barcode_reweight ) ;
		realafteroption("barcode_energy_reweight", 1.0, barcode_energy_reweight ) ;

		init = true;
	}
	return;
}
