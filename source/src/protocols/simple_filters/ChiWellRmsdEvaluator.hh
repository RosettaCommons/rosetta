// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_simple_filters_ChiWellRmsdEvaluator_hh
#define INCLUDED_protocols_simple_filters_ChiWellRmsdEvaluator_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <list>

#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

class ChiWellRmsdEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	ChiWellRmsdEvaluator( core::pose::PoseCOP, core::Size nchi_max, core::Real sasa_threshold, std::string column_tag );
	ChiWellRmsdEvaluator( core::pose::PoseCOP, core::Size nchi_max, core::Real sasa_threshold, utility::vector1< core::Size> const& selection, std::string column_tag );

	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;

private:
	core::pose::PoseCOP rmsd_pose_;
	core::scoring::ResidueSelection selection_;
	core::Size nchi_max_;
	core::Real sasa_threshold_;
	std::string tag_;
};

}
}

#endif
