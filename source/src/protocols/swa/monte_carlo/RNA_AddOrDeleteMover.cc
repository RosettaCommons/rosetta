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
#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloUtil.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
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

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_add_or_delete_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_AddOrDeleteMover::RNA_AddOrDeleteMover( RNA_AddMoverOP rna_add_mover,
																							RNA_DeleteMoverOP rna_delete_mover ) :
		rna_add_mover_( rna_add_mover ),
		rna_delete_mover_( rna_delete_mover ),
		allow_deletion_of_last_residue_( false )
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
  void
  RNA_AddOrDeleteMover::apply( core::pose::Pose & pose, std::string & move_type )
	{

		// should stuff into AddOrDeleteMover
		Size res_at_terminus;
		MovingResidueCase moving_residue_case;
		AddOrDeleteChoice add_or_delete_choice;

		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		utility::vector1< Size > const & moving_res_list = full_model_info.moving_res_list();

		//always have something in play!!?? Or permit removal??!! need to check this carefully.
		bool disallow_delete  = allow_deletion_of_last_residue_ && ( moving_res_list.size() <= 1 );

		get_random_residue_at_chain_terminus( pose, res_at_terminus, moving_residue_case, add_or_delete_choice, disallow_delete  );

		//		TR.Debug << "ADD/DELETE move ==> res: " << res_at_terminus << "  case: " << moving_residue_case << "  add/delete: " << add_or_delete_choice << std::endl;
		//		TR.Debug << std::endl;
		TR << "Move: add_or_delete " <<  add_or_delete_choice << " with moving residue case " << moving_residue_case <<  " at " << res_at_terminus << " starting from: " << pose.annotated_sequence() << std::endl;

		if ( add_or_delete_choice == DELETE ) {
			move_type = "delete";
			//std::cout << "Before delete: " << (*scorefxn)( pose ) << std::endl;
			rna_delete_mover_->apply( pose, res_at_terminus, moving_residue_case );
			TR.Debug << std::cout << pose.annotated_sequence() << std::endl;
			//std::cout << "After delete: " << (*scorefxn)( pose ) << std::endl << std::endl;
		} else {
			runtime_assert( add_or_delete_choice == ADD );
			// try to add a residue that is supposed to be sampled.
			move_type = "add";
			//std::cout << "Before adding onto " << res_at_terminus << " : " << (*scorefxn)( pose ) << std::endl;
			rna_add_mover_->apply( pose, res_at_terminus, moving_residue_case );
			TR.Debug << pose.annotated_sequence() << std::endl;
			// std::cout << "After add: " << (*scorefxn)( pose ) << std::endl << std::endl;
			//pose.dump_pdb( "after_add.pdb" );
		}

		TR << "Move: " <<  add_or_delete_choice << " at " << res_at_terminus << " resulting in : " << pose.annotated_sequence() << std::endl;

	}


	std::string
	RNA_AddOrDeleteMover::get_name() const {
		return "RNA_AddOrDeleteMover";
	}

}
}
}
