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

#include <protocols/swa/monte_carlo/RNA_AddOrDeleteMover.hh>
#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/RNA_DeleteMover.hh>
#include <protocols/swa/monte_carlo/SWA_MoveSelector.hh>

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
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose full_model_info object.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.monte_carlo.RNA_AddOrDeleteMover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_AddOrDeleteMover::RNA_AddOrDeleteMover( RNA_AddMoverOP rna_add_mover,
																							RNA_DeleteMoverOP rna_delete_mover ) :
		rna_add_mover_( rna_add_mover ),
		rna_delete_mover_( rna_delete_mover ),
		disallow_deletion_of_last_residue_( false ),
		skip_deletions_( false ),
		disallow_skip_bulge_( true ),
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

  //////////////////////////////////////////////////////////////////////////
	bool
  RNA_AddOrDeleteMover::apply( core::pose::Pose & pose, std::string & move_type /* just used by monte carlo*/ )
	{
		utility::vector1< Size > const moving_res_list = core::pose::full_model_info::get_moving_res_from_full_model_info( pose );

		//always have something in play!!?? Or permit removal??!! need to check this carefully.
		bool disallow_delete  = disallow_deletion_of_last_residue_ && ( moving_res_list.size() <= 1 );
		if ( skip_deletions_ ) disallow_delete = true;

		SWA_Move swa_move;
		swa_move_selector_->set_disallow_delete( disallow_delete );
		swa_move_selector_->set_disallow_skip_bulge( disallow_skip_bulge_ );
		swa_move_selector_->set_disallow_resample( true );
		utility::vector1< Size > const actual_sample_res = figure_out_actual_sample_res( pose );
		swa_move_selector_->get_random_move_element_at_chain_terminus( pose, swa_move, actual_sample_res );

		if ( swa_move.move_type() == NO_ADD_OR_DELETE ) {
			move_type = "no move";
			return false;
		}

		TR << swa_move << std::endl;
		TR.Debug << "Starting from: " << pose.annotated_sequence() << std::endl;
		if ( swa_move.move_type() == DELETE ) {
			move_type = "delete";
			rna_delete_mover_->apply( pose, swa_move.move_element() );
		} else {
			runtime_assert( swa_move.move_type() == ADD );
			move_type = "add";
			rna_add_mover_->apply( pose, swa_move.moving_res(), swa_move.attached_res() );
		}
		TR.Debug << "Ended with: " << pose.annotated_sequence() << std::endl;

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
		if ( sample_res_.size() > 0 ){ // user provided
			sample_res = sample_res_;
		} else { // all are allowed
			std::string const & full_sequence = const_full_model_info( pose ).full_sequence();
			for ( Size n = 1; n <= full_sequence.size(); n++ )	sample_res.push_back( n );
		}

		// get rid of bulge_res.
		utility::vector1< Size > actual_sample_res;
		for ( Size n = 1; n <= sample_res.size(); n++ ){
			if ( bulge_res_.has_value( sample_res[n] ) ) continue;
			actual_sample_res.push_back( sample_res[n] );
		}
		return actual_sample_res;
	}


	///////////////////////////////////////////////////////////////////////////////
	std::string
	RNA_AddOrDeleteMover::get_name() const {
		return "RNA_AddOrDeleteMover";
	}


}
}
}
