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

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/types.hh>

namespace protocols {
namespace sic_dock {

void
get_xform_stats(
	Xform const & sir,
	Xform const & sjr,
	platform::Real& dx, platform::Real& dy, platform::Real& dz,
	platform::Real& ex, platform::Real& ey, platform::Real& ez
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
		Xform const & s1,
		Xform const & s2,
		char ss1, char ss2
	) const;

	float
	score(
		core::pose::Pose const & pose,
		platform::Size rsd1,
		platform::Size rsd2
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
