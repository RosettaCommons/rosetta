// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_designability_score_hh
#define INCLUDED_protocols_sic_dock_designability_score_hh

#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

namespace protocols{
namespace sic_dock{

void
get_xform_stats(
	core::kinematics::Stub const & sir,
	core::kinematics::Stub const & sjr,
	core::Real& dx, core::Real& dy, core::Real& dz,
	core::Real& ex, core::Real& ey, core::Real& ez
);

struct
XfoxmScore
{
	char *hh,*he,*hl,*ee,*el,*ll;

	XfoxmScore(
		std::string datadir
	);

	void
	fillarray(
		char *a, std::string fname
	);

	void
	makebinary(
		char *a, std::string fname
	);

	float
	score(
		core::kinematics::Stub const & s1,
		core::kinematics::Stub const & s2,
		char ss1, char ss2
	) const;

	float
	score(
		core::pose::Pose const & pose,
		core::Size rsd1,
		core::Size rsd2
	) const;

	float
	score(
		core::pose::Pose & pose,
		bool compute_ss = true		
	) const;

	float
	score(
		core::pose::Pose const & pose
	) const;
};

} // end namespace sic_dock
} // end namespace protocols

#endif // INCLUDED_protocols_sic_dock_designability_score_hh
