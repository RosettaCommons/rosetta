// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/chemically_conjugated_docking/Gp_quantification_metrics.hh
/// @brief contains helper quantification metrics for the original publication of the UBQ_Gp code
/// @author Steven Lewis and Hope Anderson

#ifndef INCLUDED_protocols_chemically_conjugated_docking_Gp_quantification_metrics_HH
#define INCLUDED_protocols_chemically_conjugated_docking_Gp_quantification_metrics_HH

// Unit Headers

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

//JD headers
#include <protocols/jd2/Job.fwd.hh>

namespace protocols {
namespace chemically_conjugated_docking {


/// @details do ubiquitin-ras pair distance measurement and reporting
void ubq_ras_distance(core::pose::Pose & pose);


void ubq_ras_rotation_angle(
	core::pose::Pose const & pose,
	protocols::jd2::JobOP job_me,
	core::Size const GTPase_target,
	bool const ubiquitin);

/// @details This function is specific to the original system for which this code
///was written - if you are not trying to duplicate the initial results you
///should remove it!
void create_extra_output(
	core::pose::Pose & pose,
	bool const ubiquitin /*do pair-matching or not*/,
	core::Size const GTPase_target);

}//chemically_conjugated_docking
}//protocols

#endif //INCLUDED_protocols_chemically_conjugated_docking_Gp_quantification_metrics_HH
