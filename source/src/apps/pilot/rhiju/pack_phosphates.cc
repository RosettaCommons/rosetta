// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file swa_monte_carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/id/NamedAtomID.hh>
#include <core/init/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.hh>
#include <core/pack/rotamer_trials.hh>
#include <protocols/farna/util.hh>
#include <protocols/farna/RNA_Minimizer.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarlo.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/rna/util.hh>
#include <protocols/viewer/viewers.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <numeric/xyz.functions.hh>

//////////////////////////////////////////////////////////
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static thread_local basic::Tracer TR( "apps.pilot.rhiju.pack_phosphates" );

OPT_KEY( Boolean, icoor_test )
OPT_KEY( Boolean, pack_all_phosphates )
OPT_KEY( IntegerVector, five_prime_phosphate_pack_res )
OPT_KEY( IntegerVector, three_prime_phosphate_pack_res )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fix_alpha_for_pack_phosphate( pose::Pose & pose, Size const res ){

	using namespace numeric;
	using namespace core::id;

	Real const current_alpha = dihedral(pose.xyz( NamedAtomID( "XO3'", res ) ),
																			pose.xyz( NamedAtomID( "XP  ", res ) ),
																			pose.xyz( NamedAtomID( "XO5'", res ) ),
																			pose.xyz( NamedAtomID( " C5'", res ) ) );

	Real const current_OP2_torsion = dihedral(pose.xyz( NamedAtomID( "XOP2", res ) ),
																						pose.xyz( NamedAtomID( "XP  ", res ) ),
																						pose.xyz( NamedAtomID( "XO5'", res ) ),
																						pose.xyz( NamedAtomID( " C5'", res ) ) );

	Real const desired_OP2_torsion = dihedral(pose.xyz( NamedAtomID( " OP2", res ) ),
																						pose.xyz( NamedAtomID( " P  ", res ) ),
																						pose.xyz( NamedAtomID( " O5'", res ) ),
																						pose.xyz( NamedAtomID( " C5'", res ) ) );

	Real const new_alpha = current_alpha + (desired_OP2_torsion - current_OP2_torsion);
	Size const chi_number_pseudoalpha = pose.residue_type( res ).RNA_type().chi_number_pseudoalpha();
	pose.set_chi( chi_number_pseudoalpha /*pseudo-alpha for packable phosphate*/, res,
								principal_angle_degrees( new_alpha ) );

	Real const new_alpha_check = dihedral(pose.xyz( NamedAtomID( "XO3'", res ) ),
																				pose.xyz( NamedAtomID( "XP  ", res ) ),
																				pose.xyz( NamedAtomID( "XO5'", res ) ),
																				pose.xyz( NamedAtomID( " C5'", res ) ) );

	TR << current_alpha << " " << current_OP2_torsion << " " << desired_OP2_torsion << " " << new_alpha << " " << new_alpha_check <<  std::endl;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
pack_phosphates()
{
  using namespace core::pose;
  using namespace core::scoring;
  using namespace core::chemical;
  using namespace core::chemical::rna;
  using namespace core::pose::full_model_info;
  using namespace protocols::stepwise;
  using namespace protocols::stepwise::monte_carlo;
  using namespace protocols::stepwise::monte_carlo::rna;
  using namespace utility::file;

	// Following could be generalized to fa_standard, after recent unification, but
	// probably should wait for on-the-fly residue type generation.
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	PoseOP native_pose;
	if ( option[ in::file::native ].user() ) native_pose = get_pdb_and_cleanup( option[ in::file::native ](), rsd_set );

	// Following could go to a FullModelSetup class.
	// read starting pose(s) from disk
	utility::vector1< std::string > const & input_files = option[ in::file::s ]();
	utility::vector1< pose::PoseOP > input_poses;
	if ( input_files.size() == 0 ) input_poses.push_back( new Pose ); // just a blank pose for now.
	for ( Size n = 1; n <= input_files.size(); n++ ) 	input_poses.push_back( get_pdb_and_cleanup( input_files[ n ], rsd_set ) );
	if ( option[ full_model::other_poses ].user() ) get_other_poses( input_poses, option[ full_model::other_poses ](), rsd_set );
	fill_full_model_info_from_command_line( input_poses ); 	//FullModelInfo (minimal object needed for add/delete)


	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	else  scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy.wts" );

	// actual pose to be sampled...
	pose::Pose & pose = *input_poses[ 1 ];
	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 500, 500 );


	utility::vector1< Size > const & cutpoint_open_in_full_model = const_full_model_info( pose ).cutpoint_open_in_full_model();
	utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
	TR << "CUTPOINT_OPEN_IN_FULL_MODEL " << cutpoint_open_in_full_model << std::endl;
	TR << "RES_LIST " << res_list << std::endl;
	TR << pose.annotated_sequence() << std::endl;
	fix_up_residue_type_variants( pose );
	TR << std::endl;

	scorefxn->show( pose );

	bool pack_all_ = option[ pack_all_phosphates ]();
	utility::vector1< Size > five_prime_phosphate_pack_res_list, three_prime_phosphate_pack_res_list;

	if ( pack_all_ ){
		Size const nres_full_model = const_full_model_info( pose ).size();
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			Size const seqpos_in_full_model = res_list[ n ];
			if ( ( n == 1 || pose.fold_tree().is_cutpoint( n - 1 ) ) &&
					 seqpos_in_full_model > 1 &&
					 !cutpoint_open_in_full_model.has_value( seqpos_in_full_model - 1 ) ){
				five_prime_phosphate_pack_res_list.push_back( res_list[ n ] );
			}
			if ( ( n == pose.total_residue() || pose.fold_tree().is_cutpoint( n ) ) &&
						 seqpos_in_full_model < nres_full_model &&
						 !cutpoint_open_in_full_model.has_value( seqpos_in_full_model ) ){
						 three_prime_phosphate_pack_res_list.push_back( res_list[ n ] );
			}
		}
	} else {
		five_prime_phosphate_pack_res_list  = option[ five_prime_phosphate_pack_res ]();
		three_prime_phosphate_pack_res_list = option[ three_prime_phosphate_pack_res ]();
	}

	Pose pose_copy = pose;

	TR << "CONTROL " << std::endl;
	protocols::farna::print_hbonds( pose_copy );

	TR << "FIVE_PRIME_PHOSPHATE_PACK_RES_LIST: " << five_prime_phosphate_pack_res_list << std::endl;
	TR << "THREE_PRIME_PHOSPHATE_PACK_RES_LIST: " << three_prime_phosphate_pack_res_list << std::endl;

	for ( Size n = 1; n <= five_prime_phosphate_pack_res_list.size(); n++ ){
		Size const sample_res = five_prime_phosphate_pack_res_list[ n ];
		Size const working_sample_res = get_res_list_from_full_model_info( pose ).index( sample_res );

		remove_variant_type_from_pose_residue( pose_copy,  "VIRTUAL_PHOSPHATE", working_sample_res );

		add_variant_type_to_pose_residue( pose, "FIVE_PRIME_PACKABLE_PHOSPHATE", working_sample_res);

		RNA_ResidueType const & rna_type = pose.residue_type( working_sample_res ).RNA_type();
		TR << rna_type.chi_number_pseudogamma() << " " << rna_type.chi_number_pseudobeta() << " " << rna_type.chi_number_pseudoalpha() << std::endl;
		pose.set_chi( rna_type.chi_number_pseudogamma(), working_sample_res, pose.residue(working_sample_res).mainchain_torsion(3) );
		pose.set_chi( rna_type.chi_number_pseudobeta() , working_sample_res, pose.residue(working_sample_res).mainchain_torsion(2) );
		pose.set_chi( rna_type.chi_number_pseudoalpha(), working_sample_res, pose.residue(working_sample_res).mainchain_torsion(1) );

		pose.dump_pdb( "BEFORE_ALPHA_FIX.pdb" );

		fix_alpha_for_pack_phosphate( pose, working_sample_res );

		if ( native_pose ){
			TR << "Copying in conformation from native" << std::endl;
			TR << "alpha?" << native_pose->residue(sample_res).mainchain_torsion(1) << std::endl;
			pose.set_chi( 4 + 1, working_sample_res, native_pose->residue(sample_res).mainchain_torsion(3) );
			pose.set_chi( 4 + 2, working_sample_res, native_pose->residue(sample_res).mainchain_torsion(2) );
			pose.set_chi( 4 + 3, working_sample_res, native_pose->residue(sample_res).mainchain_torsion(1) );
		}
	}

	for ( Size n = 1; n <= three_prime_phosphate_pack_res_list.size(); n++ ){
		Size const sample_res = three_prime_phosphate_pack_res_list[ n ];
		Size const working_sample_res = get_res_list_from_full_model_info( pose ).index( sample_res );
		add_variant_type_to_pose_residue( pose, "THREE_PRIME_PACKABLE_PHOSPHATE", working_sample_res);
		RNA_ResidueType const & rna_type =  pose.residue_type( working_sample_res ).RNA_type();

		if ( native_pose ){
			TR << "Copying in conformation from native" << std::endl;
			pose.set_chi( rna_type.chi_number_pseudoepsilon(), working_sample_res, native_pose->residue(sample_res).mainchain_torsion(5) );
			pose.set_chi( rna_type.chi_number_pseudozeta(), working_sample_res, native_pose->residue(sample_res).mainchain_torsion(6) );
		}
	}

	pose.dump_pdb( "START_POSE.pdb" );
	pose_copy.dump_pdb( "CONTROL_POSE.pdb" );

	TR << "WITHOUT PACK PHOS (CONTROL) " << std::endl;
	scorefxn->show( pose_copy );

	TR << "WITH PACK PHOS " << std::endl;
	scorefxn->show( pose );
	protocols::farna::print_hbonds( pose );

	if ( false ){
		// do minimizing
		protocols::farna::RNA_Minimizer rna_minimizer;
		rna_minimizer.deriv_check( option[ OptionKeys::rna::deriv_check ]() );
		rna_minimizer.use_coordinate_constraints( !option[ OptionKeys::rna::skip_coord_constraints]() );
		rna_minimizer.skip_o2prime_trials( option[ OptionKeys::rna::skip_o2prime_trials] );
		rna_minimizer.vary_bond_geometry( option[ OptionKeys::rna::vary_geometry ] );
		rna_minimizer.apply( pose );
	}

	pack::task::PackerTaskOP pack_task = pack::task::TaskFactory::create_packer_task( pose );
	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		pack_task->nonconst_residue_task( i ).and_extrachi_cutoff( 0 );
		//		pack_task->nonconst_residue_task( i ).or_ex4( true ); //extra O2prime modeler
		//		pack_task->nonconst_residue_task( i ).or_include_current( true );
		if ( pose.residue_type( i ).has_variant_type( "FIVE_PRIME_PACKABLE_PHOSPHATE" ) ) {
			pack_task->nonconst_residue_task( i ).nonconst_rna_task().set_sample_five_prime_phosphate( true );
		}
		if ( pose.residue_type( i ).has_variant_type( "THREE_PRIME_PACKABLE_PHOSPHATE" ) ) {
			TR << "GOING TO PACK 3' PHOSPHATE AT: " << i << std::endl;
			pack_task->nonconst_residue_task( i ).nonconst_rna_task().set_sample_three_prime_phosphate( true );
		}
		pack_task->nonconst_residue_task( i ).nonconst_rna_task().set_allow_phosphate_virtualization( true );
	}

	pack::rotamer_trials( pose,	*scorefxn,	pack_task );

	TR << "AFTER PACK " << std::endl;
	scorefxn->show( pose );
	protocols::farna::print_hbonds( pose );

	pose.dump_pdb( "PACK_POSE.pdb" );

}

void
get_icoor(){
	using namespace core::pose;
  using namespace core::scoring;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace protocols::farna;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna_phenix" );
	Pose pose;
	core::pose::make_pose_from_sequence( pose, "aaa", *rsd_set );

	// for 3' packable phosphate variant
	print_internal_coords( pose );

	// for 5' packable phosphate variant
	FoldTree f = pose.fold_tree();
	f.reorder(3);
	pose.fold_tree(f);
	print_internal_coords( pose );
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	if ( option[ icoor_test ]() ){
		get_icoor();
	} else {
		pack_phosphates();
	}

	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
  exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence> -s <start pdb> -input_res <input pdb1> [ -native <native pdb file> ] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		utility::vector1< core::Size > blank_size_vector;
		NEW_OPT( icoor_test, "Silly code snippet to get ICOOR for packable phosphate variants", false );
		NEW_OPT( pack_all_phosphates, "Try packing phosphates at all residues where they could go", false );
		NEW_OPT( five_prime_phosphate_pack_res, "Which residues (in full model numbering) to add 5' packable phosphates to", blank_size_vector );
		NEW_OPT( three_prime_phosphate_pack_res, "Which residues (in full model numbering) to add 3' packable phosphates to", blank_size_vector );

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( in::file::native );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( score::weights );
		option.add_relevant( basic::options::OptionKeys::stepwise::monte_carlo::cycles );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::verbose_scores );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::skip_deletions );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::add_delete_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::minimize_single_res_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::switch_focus_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::just_min_after_mutation_frequency );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_internal_hinge_moves );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_internal_local_moves );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_skip_bulge );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::temperature );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::extra_min_res );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::allow_variable_bond_geometry );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::constraint_x0 );
		option.add_relevant( OptionKeys::stepwise::monte_carlo::constraint_tol );
		option.add_relevant( OptionKeys::stepwise::rna::num_random_samples );
		option.add_relevant( OptionKeys::stepwise::rna::erraser );
		option.add_relevant( OptionKeys::stepwise::rna::sample_res );
		option.add_relevant( OptionKeys::full_model::rna::force_syn_chi_res_list );
		option.add_relevant( OptionKeys::stepwise::rna::virtual_sugar_keep_base_fixed );
		option.add_relevant( OptionKeys::stepwise::rna::force_centroid_interaction );
		option.add_relevant( basic::options::OptionKeys::stepwise::rna::bulge_res );
		option.add_relevant( basic::options::OptionKeys::full_model::rna::terminal_res );
		option.add_relevant( OptionKeys::rna::corrected_geo );

		core::init::init(argc, argv);
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


