/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/util.cc
/// @brief util functions for denovo design of structures
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/util.hh>

//Project Headers

//Protocol Headers
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

//Core Headers
#include <core/pose/Pose.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>

// Boost/ObjexxFCL Headers

//C++ Headers

static thread_local basic::Tracer TR("protocols.denovo_design.components.util");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Tells whether the two given poses are identical based on # resides and dihedrals
bool same_pose( core::pose::Pose const & pose1, core::pose::Pose const & pose2 )
{
	if (pose1.total_residue() != pose2.total_residue()) {
		return false;
	}

	for (core::Size i = 1; i <= pose1.total_residue(); ++i) {
		if ( pose1.residue(i).name() != pose2.residue(i).name() )
			return false;
		if ( pose1.residue(i).is_protein() && !pose2.residue(i).is_protein() )
			return false;
		if ( !pose1.residue(i).is_protein() && pose2.residue(i).is_protein() )
			return false;
		if ( !pose1.residue(i).is_protein() && !pose2.residue(i).is_protein() )
			continue;

		if (std::abs(pose1.phi(i) - pose2.phi(i)) > 0.0001) {
			return false;
		}
		if (std::abs(pose1.psi(i) - pose2.psi(i)) > 0.0001) {
			return false;
		}
		if (std::abs(pose1.omega(i) - pose2.omega(i)) > 0.0001) {
			return false;
		}
	}
	return true;
}

/// @brief creates a poly-ala pose where every non-gly, non-cyd, protein residue except those in the given set are converted to alanine
void construct_poly_ala_pose( core::pose::Pose & pose, bool const keep_disulf )
{
	utility::vector1< core::Size > positions;
	for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
		if ( pose.residue(i).is_protein() ) {
			positions.push_back(i);
		}
		if ( !keep_disulf && ( pose.residue(i).name() == "CYD" ) ) {
			protocols::simple_moves::MutateResidue mut( i, "ALA" );
			mut.apply( pose );
		}
	}
	protocols::toolbox::pose_manipulation::construct_poly_ala_pose(
			pose, positions,
			false, // bool keep_pro,
			true, // bool keep_gly,
			keep_disulf ); // bool keep_disulfide_cys
}

} // namespace denovo_design
} // namespace protocols
