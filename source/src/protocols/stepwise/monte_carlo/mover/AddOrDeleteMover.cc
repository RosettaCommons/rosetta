// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AddOrDeleteMover
/// @brief AddOrDeletes an RNA residue from a chain terminus.
/// @details
/// @author Rhiju Das

#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

using namespace core;
using core::Real;
using namespace core::pose::full_model_info;

//////////////////////////////////////////////////////////////////////////
// Makes a choice, based on current pose, and information in full_model_info
//  as to whether to add or delete nucleotide and chunks, and where.
//
// This may be deprecated soon, with development of StepWiseMasterMover, which can
//  make the choice of StepWiseMove (keeping track of probabilities needed
//  for detailed balance) and could then take the job of running
//  AddMover or DeleteMover.
//
//////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.AddOrDeleteMover" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {


//////////////////////////////////////////////////////////////////////////
//constructor!
AddOrDeleteMover::AddOrDeleteMover( AddMoverOP rna_add_mover,
	DeleteMoverOP rna_delete_mover,
	FromScratchMoverOP rna_from_scratch_mover ) :
	rna_add_mover_( rna_add_mover ),
	rna_delete_mover_( rna_delete_mover ),
	rna_from_scratch_mover_( rna_from_scratch_mover ),
	disallow_deletion_of_last_residue_( false ),
	swa_move_selector_( StepWiseMoveSelectorOP( new StepWiseMoveSelector ) ),
	choose_random_( true )
{}

//////////////////////////////////////////////////////////////////////////
//destructor
AddOrDeleteMover::~AddOrDeleteMover()
{}

void
AddOrDeleteMover::apply( core::pose::Pose & pose ){
	std::string move_type = "";
	apply( pose, move_type );
}

///////////////////////////////////////////////////////////////////////////////
void
AddOrDeleteMover::apply( core::pose::Pose & pose, StepWiseMove const & swa_move ){
	TR << swa_move << std::endl;
	if ( options_->filter_complex_cycles() ) runtime_assert( swa_move_selector_->just_simple_cycles( swa_move, pose, false /*verbose*/ ) );
	TR.Debug << "Starting from: " << pose.annotated_sequence() << std::endl;
	if ( swa_move.move_type() == DELETE ) {
		rna_delete_mover_->apply( pose, swa_move.move_element() );
	} else if ( swa_move.move_type() == FROM_SCRATCH ) {
		rna_from_scratch_mover_->apply( pose, swa_move.move_element() );
	} else {
		if ( swa_move.move_type() == ADD_SUBMOTIF ) {
			runtime_assert( submotif_library_ != 0 );

			// could encapsulate in full_model_info.add_submotif();
			FullModelInfo & full_model_info = nonconst_full_model_info( pose );
			full_model_info.add_other_pose(
				submotif_library_->create_new_submotif( swa_move.move_element(), swa_move.submotif_tag() , pose )
			);
		} else {
			runtime_assert( swa_move.move_type() == ADD );
		}
		rna_add_mover_->apply( pose, swa_move );
	}
	TR.Debug << "Ended with: " << pose.annotated_sequence() << std::endl;
}

//////////////////////////////////////////////////////////////////////////
// This may be deprecated soon since we are moving move selection out
// to StepWiseMasterMover.
//////////////////////////////////////////////////////////////////////////
bool
AddOrDeleteMover::apply( core::pose::Pose & pose, std::string & move_type_string /* just used by monte carlo*/ )
{
	utility::vector1< Size > const moving_res_list = core::pose::full_model_info::get_moving_res_from_full_model_info( pose );

	bool disallow_delete  = disallow_deletion_of_last_residue_ && ( moving_res_list.size() <= 1 );
	if ( options_->skip_deletions() ||  options_->rebuild_bulge_mode() ) disallow_delete = true;

	StepWiseMove swa_move;
	swa_move_selector_->set_options( options_ );
	swa_move_selector_->set_allow_delete( !disallow_delete );
	swa_move_selector_->set_choose_random( choose_random_ );
	swa_move_selector_->set_submotif_library( submotif_library_ );

	swa_move_selector_->get_add_or_delete_element( pose, swa_move );

	move_type_string = get_move_type_string( swa_move );

	if ( swa_move.move_type() == NO_MOVE ) return false;
	apply( pose, swa_move );
	return true;
}

///////////////////////////////////////////////////////////////////////////////
void
AddOrDeleteMover::set_minimize_single_res( bool const setting ){
	rna_add_mover_->set_minimize_single_res( setting );
	rna_delete_mover_->set_minimize_after_delete( !setting );
}

///////////////////////////////////////////////////////////////////////////////
std::string
AddOrDeleteMover::get_name() const {
	return "AddOrDeleteMover";
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
AddOrDeleteMover::set_options( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options ){
	options_ = options;
}


} //mover
} //monte_carlo
} //stepwise
} //protocols
