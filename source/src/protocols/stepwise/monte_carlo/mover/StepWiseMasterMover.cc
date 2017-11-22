// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/ResampleMover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/mover/VaryLoopLengthMover.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/modeler/precomputed/PrecomputedLibraryMover.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/magnesium/MgMonteCarlo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>

#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <fstream>

static basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.StepWiseMasterMover" );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// StepWiseMasterMover for stepwise monte carlo:
//
//  1. makes decisions on add/move/resample, based on random number generator
//  2. will handle proposal ratios for detailed balance
//  3. can make specific moves, e.g. user_inputted or preminimizes.
//  4. has a 'build_full_model' mode used to put in 'garbage' residues to create full pose.
//
// Factored out of StepWiseMonteCarlo.
//
//   -- rhiju, 2015
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace core::pose;
using namespace protocols::stepwise::monte_carlo::rna;
using namespace protocols::stepwise::modeler;
using namespace core::pose::full_model_info;
using namespace protocols::magnesium;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

//Constructor
StepWiseMasterMover::StepWiseMasterMover( core::scoring::ScoreFunctionCOP scorefxn,
	protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options ):
	minimize_single_res_( false ), // changes during run
	success_( false ),
	proposal_density_ratio_( 1.0 ),
	num_tested_moves_( 0 )
{
	initialize( scorefxn, options );
}

//Destructor
StepWiseMasterMover::~StepWiseMasterMover()
{}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::apply( pose::Pose & pose ) {

	move_type_string_ = "";
	success_ = false;
	proposal_density_ratio_ = 1;

	if ( apply_legacy( pose ) ) return;

	success_ = stepwise_move_selector_->figure_out_all_possible_moves( pose );
	if ( !success_ ) return;

	StepWiseMove const stepwise_move = ( options_->enumerate() ) ? stepwise_move_selector_->swa_moves()[ 1 ] :
		stepwise_move_selector_->select_random_move( pose ) ;

	apply( pose, stepwise_move, false /* figure_out_all_possible_moves */ );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::apply( pose::Pose & pose,
	StepWiseMove const & stepwise_move,
	bool const figure_out_all_possible_moves /* = true */ ) {

	if ( figure_out_all_possible_moves ) stepwise_move_selector_->figure_out_all_possible_moves( pose );

	Real const forward_proposal_probability = stepwise_move_selector_->proposal_probability( stepwise_move );
	move_type_string_ = get_move_type_string( stepwise_move );
	pose::Pose pose_original = pose; // used for figuring out how to go backwards...

	if ( stepwise_move.attachments().size() > 0 ) {
		switch_focus_to_other_pose( pose, const_full_model_info( pose ).get_idx_for_other_pose_with_residue( stepwise_move.attachments()[1].attached_res() ),
			scorefxn_ /* allows score-based check on other poses */  );
	}

	do_the_move( stepwise_move, pose );

	// check reverse -- to get detailed balance factor:
	StepWiseMove const reverse_move = stepwise_move_selector_->reverse_move( stepwise_move, pose_original, pose );

	// silly hack -- reverse moves are not defined yet for ADD_LOOP_RES and DELETE_LOOP_RES
	if ( reverse_move != StepWiseMove() || ( stepwise_move.move_type() != ADD_LOOP_RES && stepwise_move.move_type() != DELETE_LOOP_RES ) ) {
		StepWiseMoveSelectorOP reverse_stepwise_move_selector = stepwise_move_selector_->clone();
		reverse_stepwise_move_selector->figure_out_all_possible_moves( pose );
		TR << TR.Magenta << "Reverse of " << stepwise_move << " is " << reverse_move << TR.Reset << std::endl;
		// reverse_stepwise_move_selector->output_moves();
		Real const reverse_proposal_probability = reverse_stepwise_move_selector->proposal_probability( reverse_move );
		proposal_density_ratio_ = ( reverse_proposal_probability / forward_proposal_probability );
	}

	// bias the acceptance of add moves, default factor is 1.0
	if ( stepwise_move.move_type() == ADD ) proposal_density_ratio_ *= options_->add_proposal_density_factor();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::do_the_move( StepWiseMove const & move, core::pose::Pose & pose ) {
	MoveType const & move_type = move.move_type();
	if ( move_type == ADD || move_type == DELETE || move_type == FROM_SCRATCH || move_type == ADD_SUBMOTIF ) {
		add_or_delete_mover_->apply( pose, move );
	} else if ( move_type == ADD_LOOP_RES || move_type == DELETE_LOOP_RES ) {
		vary_loop_length_mover_->apply( pose, move );
	} else {
		runtime_assert( move_type == RESAMPLE ||
			move_type == RESAMPLE_INTERNAL_LOCAL );
		resample_mover_->apply( pose, move, move_type_string_ );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseMasterMover::test_all_moves( pose::Pose & pose ) {
	// Pose working_pose = pose; // this is to prevent graphics threads from having a cow. [really?]
	num_tested_moves_ = 0;
	test_all_moves_recursively( pose );
	return true;
}

void
StepWiseMasterMover::test_all_moves_recursively( pose::Pose & pose ) {
	pose::Pose pose_save = pose;
	runtime_assert( options_->num_random_samples() == 0 ); // silly hack to prevent computation beyond rerooting...
	add_mover_->set_start_added_residue_in_aform( true  );

	success_ = stepwise_move_selector_->figure_out_all_possible_moves( pose, true /*verbose*/ );
	utility::vector1< StepWiseMove > const stepwise_moves = stepwise_move_selector_->swa_moves();
	for ( Size n = 1; n <= stepwise_moves.size(); n++ ) {
		StepWiseMove const & stepwise_move = stepwise_moves[ n ];
		pose::Pose pose_test = pose_save;
		TR << TR.Blue << "move " << n << " out of " << stepwise_moves.size() << " for " << get_all_res_list( pose_test ) << ". Total moves so far: " << ++num_tested_moves_ << TR.Reset << std::endl;
		apply( pose_test, stepwise_move, true );
		if ( stepwise_move.move_type() == ADD ) {
			pose = pose_test; // show in graphics.
			test_all_moves_recursively( pose ); // recurse through all possibilities
		}
	}

	pose = pose_save;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Original move selection had a cascade of silly random number choices.
// Some ad hoc choices when no moves were found.
// Did not allow for easy computation of proposal ratios.
//
//  Deprecate after May 2015 if not in use! -- rhiju, jan. 2015.
//
bool
StepWiseMasterMover::apply_legacy( pose::Pose & pose ) {

	if ( options_->new_move_selector() ) return false;

	if ( numeric::random::rg().uniform() < options_->switch_focus_frequency() ) switch_focus_among_poses_randomly( pose );

	set_minimize_single_res( numeric::random::rg().uniform() <= options_->minimize_single_res_frequency() );
	if ( numeric::random::rg().uniform() < options_->add_delete_frequency() ) {
		success_ = add_or_delete_mover_->apply( pose, move_type_string_ );
	} else {
		success_ = resample_mover_->apply( pose, move_type_string_ );
	}

	// bias the acceptance of add moves, default factor is 1.0
	if ( move_type_string_ == "add" ) proposal_density_ratio_ *= options_->add_proposal_density_factor();

	// this is probably deprecated now -- we may try to optimize minimization in the future,
	// but that will almost certainly require minimizing more than single residues.
	if ( minimize_single_res_ ) move_type_string_ += "-minsngl";

	return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::initialize( core::scoring::ScoreFunctionCOP scorefxn,
	protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options ){
	scorefxn_ = scorefxn;
	options_ = options;
	initialize();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::initialize(){

	using namespace protocols::stepwise::modeler::rna;
	runtime_assert( scorefxn_ != 0 );
	runtime_assert( options_  != 0 );

	stepwise_modeler_ = setup_unified_stepwise_modeler( options_, scorefxn_ );
	stepwise_modeler_->set_native_pose( get_native_pose() );

	// maybe AddMover could just hold a copy of ResampleMover...
	//  [ will revisit this when testing detailed balance via
	//    a 'sample' move that can choose whether to add or resample inside
	//    the Modeler. ]
	add_mover_ = AddMoverOP( new AddMover( scorefxn_ ) );
	add_mover_->set_native_pose( get_native_pose() );
	add_mover_->set_start_added_residue_in_aform( false );
	add_mover_->set_presample_added_residue(  true );
	add_mover_->set_presample_by_swa(  true );
	add_mover_->set_stepwise_modeler( stepwise_modeler_->clone_modeler() );
	add_mover_->set_designing_with_noncanonicals( options_->designing_with_noncanonicals() );
	add_mover_->set_sample_pH( options_->sample_pH() );

	delete_mover_ = DeleteMoverOP( new DeleteMover );
	delete_mover_->set_native_pose( get_native_pose() );
	delete_mover_->set_stepwise_modeler( stepwise_modeler_->clone_modeler() );
	delete_mover_->set_options( options_ );

	from_scratch_mover_ = FromScratchMoverOP( new FromScratchMover );
	from_scratch_mover_->set_native_pose( get_native_pose() );
	from_scratch_mover_->set_stepwise_modeler( stepwise_modeler_->clone_modeler() );

	add_or_delete_mover_ = AddOrDeleteMoverOP( new AddOrDeleteMover( add_mover_, delete_mover_, from_scratch_mover_ ) );
	add_or_delete_mover_->set_options( options_ );
	add_or_delete_mover_->set_submotif_library( submotif_library_ );

	resample_mover_ = ResampleMoverOP( new ResampleMover( stepwise_modeler_->clone_modeler() ) );
	resample_mover_->set_native_pose( get_native_pose() );
	resample_mover_->set_options( options_ );

	vary_loop_length_mover_ = VaryLoopLengthMoverOP( new VaryLoopLengthMover );

	stepwise_move_selector_ = StepWiseMoveSelectorOP( new StepWiseMoveSelector( options_ ) );
	stepwise_move_selector_->set_choose_random( !options_->enumerate() );
	if ( options_->new_move_selector() || options_->test_all_moves() ) stepwise_move_selector_->set_force_unique_moves( true );
	stepwise_move_selector_->set_submotif_library( submotif_library_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::initialize_pose_if_empty( pose::Pose & pose ){
	if ( pose.size() > 0 ) return;
	runtime_assert( options_->from_scratch_frequency() > 0.0 || options_->mapfile_activated() );
	add_or_delete_mover_->apply( pose );
	runtime_assert( pose.size() > 0 );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::set_minimize_single_res( bool const minimize_single_res ){
	add_or_delete_mover_->set_minimize_single_res( minimize_single_res );
	resample_mover_->set_minimize_single_res( minimize_single_res );
	minimize_single_res_ = minimize_single_res;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::preminimize_pose( pose::Pose & pose ) {
	fix_up_residue_type_variants( pose );
	stepwise_modeler_->set_moving_res_and_reset( 0 );
	stepwise_modeler_->set_working_prepack_res( get_all_residues( pose ) );
	stepwise_modeler_->set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );
	stepwise_modeler_->apply( pose );
	utility::vector1< pose::PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ) {
		preminimize_pose( *( other_pose_list[ n ] ) );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseMasterMover::do_test_move( StepWiseMove const & move,
	pose::Pose & pose ) {
	// preminimization stuff.
	if ( options_->preminimize() ) preminimize_pose( pose );
	if ( options_->test_all_moves() ) return test_all_moves( pose );
	if ( move.move_type() == NO_MOVE ) {
		if ( options_->preminimize() ) {
			std::string out_pdb( "preminimize.pdb" );
			std::string const tag = tag_from_pose( pose );
			if ( tag != "empty_tag" && tag.substr( tag.size()-4 ) == ".pdb" ) out_pdb = tag.substr( 0, tag.size()-4 ) +  ".preminimize.pdb";
			TR << TR.Green << "Created: " <<  out_pdb << TR.Reset << std::endl;
			pose.dump_pdb( out_pdb );
			return true; // did something.
		}
		return false;
	}

	do_the_move( move, pose );

	return true;
}

std::string
name_from_move( StepWiseMove const & stepwise_move, Pose const & start_pose ) {
	std::stringstream ss;
	ss << "seq_rebuild_temp_"
		<< start_pose.pdb_info()->chain( stepwise_move.moving_res() )
		<< ':' << start_pose.pdb_info()->number( stepwise_move.moving_res() )
		<< ".out";
	return ss.str();
}

/////////////////////////////////////////////////////////
void
StepWiseMasterMover::resample_full_model( pose::Pose const & start_pose, pose::Pose & output_pose, bool const checkpointing_breadcrumbs ) {
	utility::vector1< Size > residues_to_rebuild;
	for ( Size ii = 1; ii <= start_pose.size(); ++ii ) {
		residues_to_rebuild.push_back( ii );
	}
	resample_full_model( start_pose, output_pose, checkpointing_breadcrumbs, residues_to_rebuild );
}

void
StepWiseMasterMover::resample_full_model( pose::Pose const & start_pose, pose::Pose & output_pose, bool const checkpointing_breadcrumbs, utility::vector1< Size > const & residues_to_rebuild ) {
	// Generates a set of 'bad residues' based on DOFs and other criteria. Then, only rebuilds the ones also present in
	// residues_to_rebuild.

	using namespace options;
	output_pose = start_pose;

	initialize();

	// AMW: This was the old, move selector dependent method.
	//stepwise_move_selector_->figure_out_all_possible_moves( start_pose );
	//for ( StepWiseMove const & stepwise_move : stepwise_move_selector_->swa_moves() ) {
	// TR << "Possible Move: " << stepwise_move << "." << std::endl;
	//}

	utility::vector1< StepWiseMove > stepwise_moves;//, internal_moves, terminal_moves;
	//stepwise_move_selector_->get_resample_terminal_move_elements( output_pose, terminal_moves );
	//stepwise_move_selector_->get_resample_internal_local_move_elements( output_pose, internal_moves );
	//for ( StepWiseMove const & stepwise_move : terminal_moves ) {
	// TR << "Applying terminal move: " << stepwise_move << "." << std::endl;
	//}
	//for ( StepWiseMove const & stepwise_move : internal_moves ) {
	// TR << "Applying internal move: " << stepwise_move << "." << std::endl;
	//}

	// Pose analysis.
	for ( Size ii = 1; ii <= start_pose.size(); ++ii ) {
		if ( !residues_to_rebuild.contains( ii ) ) continue;
		if ( start_pose.residue( ii ).is_virtual_residue() ) continue;

		if ( !start_pose.residue( ii ).is_polymer() ) {
			// Ligand -- jump_dock
			TR << "Found ligand -- " << start_pose.residue( ii ).name() << " -- at " << ii << std::endl;
			TR << start_pose.fold_tree() << std::endl;
			MoveElement move_element;
			move_element.push_back( ii );
			utility::vector1< Attachment > attachments;

			bool connected_by_jump = false;

			// Somehow parent might not be. Might need to look for child?
			Size attached_res = Size( start_pose.fold_tree().get_parent_residue( ii, connected_by_jump ) );

			TR << "Found an attached res POSSIBLY connected by a jump: " << attached_res << std::endl;
			if ( !connected_by_jump ) {
				for ( Size jj = 1; jj <= start_pose.size(); ++jj ) {
					bool connected_by_jump2 = false;
					Size const possibly_ii = Size( start_pose.fold_tree().get_parent_residue( jj, connected_by_jump2 ) );
					if ( ii == possibly_ii ) {
						connected_by_jump = connected_by_jump2;
						attached_res = jj;
						break;
					}
				}
			}
			TR << "Found an attached res connected by a jump: " << attached_res << std::endl;
			runtime_assert( connected_by_jump );
			attachments.emplace_back( attached_res, JUMP_DOCK );
			stepwise_moves.emplace_back( move_element, attachments, RESAMPLE );

			continue;
		}

		//if ( start_pose.residue( ii - 1 ).is_virtual_residue() ) continue;
		//if ( start_pose.residue( ii + 1 ).is_virtual_residue() ) continue;
		MoveElement move_element;
		move_element.push_back( ii );
		utility::vector1< Attachment > attachments;
		TR.Trace << "Triplet of concerned residues:" << std::endl;
		if ( ii > 1 ) TR.Trace << ii - 1 << " " << start_pose.residue( ii - 1 ).name() << std::endl;
		TR.Trace << ii << " " << start_pose.residue( ii ).name() << std::endl;
		if ( ii < start_pose.size() ) TR.Trace << ii + 1 << " " << start_pose.residue( ii + 1 ).name() << std::endl;
		if ( !start_pose.residue( ii ).has_lower_connect() ) {
			// Doesn't have to be a lower terminus variant -- what matters is specifically lacking lower.
			if ( start_pose.residue( ii ).is_carbohydrate() ) {
				attachments.emplace_back( ii + 1, BOND_TO_NEXT );
				stepwise_moves.emplace_back( move_element, attachments, RESAMPLE );
			}
		} else if ( !start_pose.residue( ii ).has_upper_connect() ) {
			// Doesn't have to be a upper terminus variant -- what matters is specifically lacking upper.
			if ( start_pose.residue( ii ).is_carbohydrate() ) {
				attachments.emplace_back( ii - 1, BOND_TO_PREVIOUS );
				stepwise_moves.emplace_back( move_element, attachments, RESAMPLE );
			}
		} else {
			attachments.emplace_back( ii - 1, BOND_TO_PREVIOUS );
			attachments.emplace_back( ii + 1, BOND_TO_NEXT );
			stepwise_moves.emplace_back( move_element, attachments, RESAMPLE_INTERNAL_LOCAL );
			TR << "Adding Move: " << stepwise_moves[ stepwise_moves.size() ] << "." << std::endl;
		}
	}


	// TODO: make this happen at the very end. Right now, do it first so we can see the cool effects.
	// AMW: Don't support this yet.
	// bool has_mg = false;
	// for ( Size ii = 1; ii <= output_pose.size(); ++ii ) {
	//  if ( output_pose.residue_type( ii ).name() == "MG" ) has_mg = true;
	// }
	// if ( has_mg ) {
	// MgMonteCarlo mg;
	// mg.set_output_pdb( true );
	// mg.apply( output_pose );
	//}

	// do moves in serial
	for ( StepWiseMove const & stepwise_move : stepwise_moves ) {
		TR.Debug << "Applying Move: " << stepwise_move << "." << std::endl;
		//continue;
		apply( output_pose, stepwise_move, true /* figure_out_all_possible_moves */ );
		scorefxn_->show( output_pose );
		if ( checkpointing_breadcrumbs ) {
			std::ofstream out( name_from_move( stepwise_move, start_pose ) );
			out << "DONE!\n";
			out.close();
		}
	}
}

/////////////////////////////////////////////////////////
// Called by build_full_model() in stepwise/monte_carlo/util.cc
void
StepWiseMasterMover::build_full_model( pose::Pose const & start_pose, pose::Pose & full_model_pose ) {
	using namespace options;
	full_model_pose = start_pose;

	runtime_assert( options_->skip_deletions() ); // totally inelegant, must be set outside.
	initialize();
	add_or_delete_mover_->set_choose_random( false );
	add_mover_->set_start_added_residue_in_aform( true );
	add_mover_->set_presample_added_residue( false );

	std::string move_type_string;
	while ( add_or_delete_mover_->apply( full_model_pose, move_type_string ) ) {
		TR.Debug << "Building full model: " << move_type_string << std::endl;
	}

	// check that everything was indeed built.
	utility::vector1< Size > const & working_res = const_full_model_info( full_model_pose ).working_res();
	utility::vector1< Size > const & res_list    = const_full_model_info( full_model_pose ).res_list();
	utility::vector1< Size > const & bulge_res   = const_full_model_info( full_model_pose ).rna_bulge_res();

	for ( Size n = 1; n <= working_res.size(); n++ ) {
		if ( !stepwise_addable_residue( working_res[ n ],
				const_full_model_info( start_pose ).full_model_parameters()->non_standard_residue_map() ) ) continue;
		runtime_assert( res_list.has_value( working_res[ n ] ) || bulge_res.has_value( working_res[ n ] ) );
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseMasterMover::set_options( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options )
{
	options_ = options;
	initialize(); // all movers accept options -- got to make them again.
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// used in all movers.
modeler::StepWiseModelerOP
setup_unified_stepwise_modeler( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP mc_options, core::scoring::ScoreFunctionCOP scorefxn ){
	using namespace modeler;
	using namespace modeler::rna;
	using namespace modeler::protein;
	using namespace modeler::precomputed;

	StepWiseModelerOP stepwise_modeler( new StepWiseModeler( scorefxn ) );
	protocols::stepwise::modeler::options::StepWiseModelerOptionsOP options = mc_options->setup_modeler_options();
	stepwise_modeler->set_options( options );

	// Precomputed library will probably be deprecated in favor of SubmotifLibrary soon. -- rhiju, jan 2015.
	if ( mc_options->use_precomputed_library() ) {
		PrecomputedLibraryMoverOP precomputed_library_mover( new PrecomputedLibraryMover );
		stepwise_modeler->set_precomputed_library_mover( precomputed_library_mover );
	}
	return stepwise_modeler;
}


} //mover
} //monte_carlo
} //stepwise
} //protocols
