// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Aroop Sircar

#ifndef INCLUDED_protocols_antibody_legacy_AntibodyClass_hh
#define INCLUDED_protocols_antibody_legacy_AntibodyClass_hh

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>

// C++ Headers

// Utility Headers

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody_legacy {

/// antibody definition
class Antibody {

public:
	/// default constructor
	Antibody();

	/// constructor with arguments
	Antibody( core::pose::Pose& pose_in );
	Antibody( core::pose::Pose& pose_in, bool camelid );
	Antibody( core::pose::Pose& pose_in, std::string cdr_name );

	void set_defaults();

	void set_Fv( core::pose::Pose& pose_in );

	void set_Fv( core::pose::Pose& pose_in, bool camelid );

	void populate_all_cdrs();

	void all_cdr_fold_tree();

	/// align current Fv to native.Fv
	void align_to_native( Antibody & native );

	core::Size cdrl_[4][3];
	core::Size cdrh_[4][3];
	core::Size lfr_[8][3];
	core::Size hfr_[7][3];

	// Start coordinates of active loop
	core::Size current_start;
	// End coordinates of active loop
	core::Size current_end;

	// Pose containing antibody variable region, Fv
	core::pose::Pose Fv;

	bool kinked_;
	bool extended_;
	utility::vector1< char > Fv_sequence_;

	loops::Loops all_cdr_loops;

	core::kinematics::MoveMap ab_movemap;

private:

	core::Size cdr_h3_cut_;
	bool camelid_;

	void detect_CDR_H3_stem_type();
	void detect_camelid_CDR_H3_stem_type();
	void detect_regular_CDR_H3_stem_type();
	void update_sequence();
};


} //namespace antibody
} //namespace protocols


#endif //INCLUDED_protocols_loops_AntibodyClass_HH
