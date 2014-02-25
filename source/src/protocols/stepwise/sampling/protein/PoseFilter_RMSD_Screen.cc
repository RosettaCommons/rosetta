// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseFilter_RMSD_Screen
/// @brief Subclass of PoseFilter
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/PoseFilter_RMSD_Screen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <utility/exit.hh>
#include <ObjexxFCL/string.functions.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	PoseFilter_RMSD_Screen::PoseFilter_RMSD_Screen( utility::vector1< Size > const calc_rms_res,
																									core::pose::PoseCOP native_pose,
																									Real const rmsd_cutoff,
																									bool const force_align ):
		calc_rms_res_( calc_rms_res ),
		native_pose_( native_pose ),
		rmsd_cutoff_( rmsd_cutoff ),
		force_align_( force_align ),
		cluster_by_all_atom_rmsd_( false ) // currently disabled -- need to fix for RNA.
  {
		initialize_corresponding_atom_id_map( *native_pose );
  }

  //////////////////////////////////////////////////////////////////////////
	bool
	PoseFilter_RMSD_Screen::passes_filter( core::pose::Pose & pose ){

		using namespace core::scoring;

		Real rmsd( 0.0 );

		if ( calc_rms_res_.size() == 0 ) {
			rmsd = rms_at_corresponding_atoms( pose, *native_pose_, corresponding_atom_id_map_ );
		} else if ( force_align_ ) {
			rmsd = rms_at_corresponding_atoms( pose, *native_pose_, corresponding_atom_id_map_, calc_rms_res_ );
		} else {
			// assumes prealignment of poses!!!
			rmsd = rms_at_corresponding_atoms_no_super( pose, *native_pose_,
																									corresponding_atom_id_map_, calc_rms_res_ );
		}

		if ( rmsd < rmsd_cutoff_ ) {
			//std::cout << "YES DID PASS RMSD FILTER: " << rmsd << " " << rmsd_cutoff_ << std::endl;
			return true;
		}

		//std::cout << "DID NOT PASS RMSD FILTER: " << rmsd << " " << rmsd_cutoff_ << std::endl;
		return false;
	}



	/////////////////////////////////////////////////////////////////////
	void
	PoseFilter_RMSD_Screen::initialize_corresponding_atom_id_map( core::pose::Pose const & pose ){
		using namespace core::scoring;

		// presumably can figure this out based on whether protein or RNA. anyway.
		if( cluster_by_all_atom_rmsd_ ) {
			setup_matching_heavy_atoms( pose, pose, corresponding_atom_id_map_ );
		} else {
			setup_matching_protein_backbone_heavy_atoms( pose, pose, corresponding_atom_id_map_ );
		}

	}


} //protein
} //sampling
} //stepwise
} //protocols
