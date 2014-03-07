// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_AddOrDeleteMover
/// @brief AddOrDeletes an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/stepwise/monte_carlo/rna/RNA_AddOrDeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_AddMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_FromScratchMover.hh>
#include <protocols/stepwise/monte_carlo/SWA_MoveSelector.hh>
#include <protocols/stepwise/monte_carlo/rna/StepWiseRNA_MonteCarloOptions.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

using namespace core;
using core::Real;
using namespace core::pose::full_model_info;

//////////////////////////////////////////////////////////////////////////
// Makes a choice, based on current pose, and information in full_model_info
//  as to whether to add or delete nucleotide and chunks, and where.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.monte_carlo.RNA_AddOrDeleteMover" ) ;

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_AddOrDeleteMover::RNA_AddOrDeleteMover( RNA_AddMoverOP rna_add_mover,
																							RNA_DeleteMoverOP rna_delete_mover,
																							RNA_FromScratchMoverOP rna_from_scratch_mover ) :
		rna_add_mover_( rna_add_mover ),
		rna_delete_mover_( rna_delete_mover ),
		rna_from_scratch_mover_( rna_from_scratch_mover ),
		disallow_deletion_of_last_residue_( false ),
		swa_move_selector_( new SWA_MoveSelector )
	{}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_AddOrDeleteMover::~RNA_AddOrDeleteMover()
  {}

  void
  RNA_AddOrDeleteMover::apply( core::pose::Pose & pose ){
		std::string move_type = "";
		apply( pose, move_type );
	}

	///////////////////////////////////////////////////////////////////////////////
	void
	RNA_AddOrDeleteMover::apply( core::pose::Pose & pose, SWA_Move const & swa_move ){
		TR << swa_move << std::endl;
		TR.Debug << "Starting from: " << pose.annotated_sequence() << std::endl;
		if ( swa_move.move_type() == DELETE ) {
			rna_delete_mover_->apply( pose, swa_move.move_element() );
		} else if ( swa_move.move_type() == FROM_SCRATCH ) {
			rna_from_scratch_mover_->apply( pose, swa_move.move_element() );
		} else {
			runtime_assert( swa_move.move_type() == ADD );
			rna_add_mover_->apply( pose, swa_move );
		}
		TR.Debug << "Ended with: " << pose.annotated_sequence() << std::endl;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  RNA_AddOrDeleteMover::apply( core::pose::Pose & pose, std::string & move_type /* just used by monte carlo*/ )
	{
		utility::vector1< Size > const moving_res_list = core::pose::full_model_info::get_moving_res_from_full_model_info( pose );

		bool disallow_delete  = disallow_deletion_of_last_residue_ && ( moving_res_list.size() <= 1 );
		if ( options_->skip_deletions() || 	options_->rebuild_bulge_mode() ) disallow_delete = true;

		SWA_Move swa_move;
		swa_move_selector_->set_allow_delete( !disallow_delete );
		swa_move_selector_->set_allow_skip_bulge( options_->allow_skip_bulge() );
		swa_move_selector_->set_from_scratch_frequency( options_->from_scratch_frequency() );
		swa_move_selector_->set_intermolecular_frequency( options_->intermolecular_frequency() );

		utility::vector1< Size > const actual_sample_res = figure_out_actual_sample_res( pose );
		swa_move_selector_->get_random_add_or_delete_element( pose, swa_move, actual_sample_res );

		move_type = to_string( swa_move.move_type() );
		std::transform(move_type.begin(), move_type.end(), move_type.begin(), ::tolower); // this is why we love C

		if ( swa_move.move_type() == NO_MOVE ) return false;
		apply( pose, swa_move );
		return true;
	}

	///////////////////////////////////////////////////////////////////////////////
	void
	RNA_AddOrDeleteMover::set_minimize_single_res( bool const setting ){
		rna_add_mover_->set_minimize_single_res( setting );
		rna_delete_mover_->set_minimize_after_delete( !setting );
	}

	///////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size >
	RNA_AddOrDeleteMover::figure_out_actual_sample_res( pose::Pose const & pose ) const{

		utility::vector1< Size > sample_res;
		if ( options_->sample_res().size() > 0 ){ // user provided
			sample_res = options_->sample_res();
		} else { // all are allowed
			std::string const & full_sequence = const_full_model_info( pose ).full_sequence();
			for ( Size n = 1; n <= full_sequence.size(); n++ )	sample_res.push_back( n );
		}

		// get rid of bulge_res.
		utility::vector1< Size > actual_sample_res;
		for ( Size n = 1; n <= sample_res.size(); n++ ){
			if ( options_->bulge_res().has_value( sample_res[n] ) ) continue;
			actual_sample_res.push_back( sample_res[n] );
		}
		return actual_sample_res;
	}


	///////////////////////////////////////////////////////////////////////////////
	std::string
	RNA_AddOrDeleteMover::get_name() const {
		return "RNA_AddOrDeleteMover";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_AddOrDeleteMover::set_options( StepWiseRNA_MonteCarloOptionsCOP options ){
		options_ = options;
	}


} //rna
} //monte_carlo
} //stepwise
} //protocols
