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
#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloUtil.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>


using namespace core;
using namespace core::pose::full_model_info;
using core::Real;

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
	RNA_DeleteMover::RNA_DeleteMover():
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
  RNA_DeleteMover::apply( core::pose::Pose & pose, Size const res_to_delete, MovingResidueCase const moving_residue_case ) const
	{

		remove_cutpoint_variants_at_res_to_delete( pose, res_to_delete );

		pose.delete_polymer_residue( res_to_delete );

		if ( moving_residue_case == CHAIN_TERMINUS_5PRIME )	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_delete );

		// important book-keeping.
		reorder_full_model_info_after_delete( pose, res_to_delete );

		if ( minimize_after_delete_ ) minimize_after_delete( pose );
	}


	//////////////////////////////////////////////////////////////////////
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

		utility::vector1< Size >  possible_res;
		utility::vector1< MovingResidueCase > moving_residue_cases;
		utility::vector1< AddOrDeleteChoice > add_or_delete_choices;

		// don't do any minimizing -- just get rid of everything...
		bool const minimize_after_delete_save( minimize_after_delete_ );
		minimize_after_delete_ = false;

		get_potential_delete_residues( pose,
																	 possible_res,
																	 moving_residue_cases,
																	 add_or_delete_choices );

		if ( possible_res.size() > 0 ){ // recursively delete all residues.
			apply( pose, possible_res[1], moving_residue_cases[1] );
			wipe_out_moving_residues( pose );
		}

		minimize_after_delete_ = minimize_after_delete_save;

	}

	////////////////////////////////////////////////////////////////////
	void
	RNA_DeleteMover::set_minimize_scorefxn( core::scoring::ScoreFunctionOP minimize_scorefxn ){
		minimize_scorefxn_ = minimize_scorefxn;
	}

	////////////////////////////////////////////////////////////////////
	void
	RNA_DeleteMover::minimize_after_delete( pose::Pose & pose ) const{

		using namespace core::pose::full_model_info;

		runtime_assert( minimize_scorefxn_ != 0 );

		TR << "Minimizing after delete " << std::endl;

		TR << "Initial: " << ( *minimize_scorefxn_) ( pose ) << std::endl;

		// following is a bit of a hack.
		// however, I wanted to make sure that I use the same minimizer as AddMover.

		// need a 'sampling' residue to send into StepWiseRNA_Modeler. It actually wont be sampled
		// because of skip_sampling.
		utility::vector1< Size > possible_res;
		get_potential_resample_residues( pose, possible_res );
		if ( possible_res.size() == 0 ) return;
		Size const some_residue = possible_res[ 1 ];

		swa::rna::StepWiseRNA_Modeler stepwise_rna_modeler( some_residue, minimize_scorefxn_ );
		stepwise_rna_modeler.set_skip_sampling( true );
		stepwise_rna_modeler.set_minimize_res( nonconst_full_model_info_from_pose( pose ).moving_res_list() );

		stepwise_rna_modeler.apply( pose );

		TR << "Final: " << ( *minimize_scorefxn_) ( pose ) << std::endl;

	}

	//////////////////////////////////////////////////////////////////////
	std::string
	RNA_DeleteMover::get_name() const {
		return "RNA_DeleteMover";
	}

}
}
}
