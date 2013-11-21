// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeleteMover
/// @brief Deletes an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/swa/monte_carlo/RNA_DeleteMover.hh>
#include <protocols/swa/monte_carlo/SWA_MonteCarloUtil.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/StepWiseUtil.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>


using namespace core;
using namespace core::pose::full_model_info;
using core::Real;
using utility::make_tag_with_dashes;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose full_model_info object.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_delete_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_DeleteMover::RNA_DeleteMover( core::pose::PoseOP native_pose, core::Real constraint_x0, core::Real constraint_tol ):
		minimize_after_delete_( true ),
		native_pose_ ( native_pose ),
		constraint_x0_( constraint_x0 ),
		constraint_tol_( constraint_tol )
  {}
	
	RNA_DeleteMover::RNA_DeleteMover( ):
	minimize_after_delete_( true )
	{}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_DeleteMover::~RNA_DeleteMover()
  {}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_DeleteMover::apply( core::pose::Pose &  )
	{
		std::cout << "not defined yet" << std::endl;
	}


	//////////////////////////////////////////////////////////////////////
  void
  RNA_DeleteMover::apply( core::pose::Pose & pose, Size const res_to_delete_in_full_model_numbering ) const
	{
		apply( pose, utility::tools::make_vector1( res_to_delete_in_full_model_numbering ) );
	}


	//////////////////////////////////////////////////////////////////////
  void
  RNA_DeleteMover::apply( core::pose::Pose & pose, utility::vector1< Size > const & residues_to_delete_in_full_model_numbering ) const
	{
		using namespace core::pose;

		FullModelInfo & full_model_info = nonconst_full_model_info( pose );
		utility::vector1< Size > const residues_to_delete = full_model_info.full_to_sub( residues_to_delete_in_full_model_numbering );

		PoseOP sliced_out_pose_op = new Pose;
		slice_out_pose( pose, *sliced_out_pose_op, residues_to_delete );
		if ( sliced_out_pose_op->total_residue() > 1 ) full_model_info.add_other_pose( sliced_out_pose_op );

		fix_up_residue_type_variants( *sliced_out_pose_op ); // now make this include chain terminus!
		fix_up_residue_type_variants( pose ); // now make this include chain terminus!
		
		if ( native_pose_ ) {
			clear_constraints_recursively( pose );
			superimpose_recursively_and_add_constraints( pose, *native_pose_, constraint_x0_, constraint_tol_ );
		}

		if ( minimize_after_delete_ ) minimize_after_delete( pose );

	}


	//////////////////////////////////////////////////////////////////////
	// following should be deprecated by fix_up_residue_type_variants
  void
  RNA_DeleteMover::remove_cutpoint_variants_at_res_to_delete( core::pose::Pose & pose, Size const & res_to_delete ) const {

		using namespace core::chemical;
		using namespace core::pose;

		if ( pose.residue_type( res_to_delete ).has_variant_type( CUTPOINT_UPPER ) ){

			runtime_assert( res_to_delete > 1 );
			runtime_assert( pose.residue_type( res_to_delete - 1 ).has_variant_type( CUTPOINT_LOWER ) );

			remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, res_to_delete - 1 );
			remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, res_to_delete ); // this is actually gratuitous as we are about to delete.

		}

		if ( pose.residue_type( res_to_delete ).has_variant_type( CUTPOINT_LOWER ) ){

			runtime_assert( res_to_delete < pose.total_residue() );
			runtime_assert( pose.residue_type( res_to_delete+1 ).has_variant_type( CUTPOINT_UPPER ) );

			remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, res_to_delete ); // this is actually gratuitous as we are about to delete.
			remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, res_to_delete + 1 );

		}

	}


	//////////////////////////////////////////////////////////////////////
	void
	RNA_DeleteMover::wipe_out_moving_residues( pose::Pose & pose ) {

		// don't do any minimizing -- just get rid of everything...
		bool const minimize_after_delete_save( minimize_after_delete_ );
		minimize_after_delete_ = false;

		utility::vector1< SWA_Move > swa_moves;
		get_delete_move_elements( pose, swa_moves);

		if ( swa_moves.size() > 0 ){ // recursively delete all residues.
			apply( pose, swa_moves[1].move_element() );
			wipe_out_moving_residues( pose );
		}

		minimize_after_delete_ = minimize_after_delete_save;

	}

	////////////////////////////////////////////////////////////////////
	void
	RNA_DeleteMover::minimize_after_delete( pose::Pose & pose ) const{

		using namespace core::pose::full_model_info;

		stepwise_rna_modeler_->set_moving_res_and_reset( 0 );
		stepwise_rna_modeler_->set_skip_sampling( true );
		stepwise_rna_modeler_->set_minimize_res( get_moving_res_from_full_model_info( pose ) );
		stepwise_rna_modeler_->apply( pose );

	}

	///////////////////////////////////////////////////////////////////
	void
	RNA_DeleteMover::set_stepwise_rna_modeler( protocols::swa::rna::StepWiseRNA_ModelerOP stepwise_rna_modeler ){
		stepwise_rna_modeler_ = stepwise_rna_modeler;
	}

	//////////////////////////////////////////////////////////////////////
	std::string
	RNA_DeleteMover::get_name() const {
		return "RNA_DeleteMover";
	}

}
}
}
