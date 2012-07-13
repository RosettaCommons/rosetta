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
#include <protocols/swa/monte_carlo/SubToFullInfo.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose sub_to_full_info object.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_delete_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	//  RNA_DeleteMover::RNA_DeleteMover()  {}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_DeleteMover::~RNA_DeleteMover()
  {}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_DeleteMover::apply( core::pose::Pose & pose )
	{
		std::cout << "not defined yet" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////
  void
  RNA_DeleteMover::apply( core::pose::Pose & pose, Size const res_to_delete, MovingResidueCase const moving_residue_case )
	{

		pose.delete_polymer_residue( res_to_delete );
		if ( moving_residue_case == CHAIN_TERMINUS_5PRIME )	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_delete );

		// important book-keeping.
		reorder_sub_to_full_info_after_delete( pose, res_to_delete );
	}

	std::string
	RNA_DeleteMover::get_name() const {
		return "RNA_DeleteMover";
	}

}
}
}
