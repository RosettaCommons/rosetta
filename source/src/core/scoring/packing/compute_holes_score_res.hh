// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/compute_holes_score_res.hh
/// @brief  Packing Score
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_packing_compute_holes_score_res_hh
#define INCLUDED_core_scoring_packing_compute_holes_score_res_hh

#include <core/scoring/packing/PoseBalls.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.hh>

#include <numeric/xyzVector.hh>

#include <core/kinematics/Jump.hh>
#include <core/scoring/packing/HolesParamsRes.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packing {


Real
compute_holes_score_res(
	pose::Pose const & pose,
	HolesParamsRes const & params
);

Real
compute_holes_deriv_res(
	pose::Pose const & pose,
	HolesParamsRes const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & derivs
);

Real
compute_holes_score_res(
	pose::Pose const & pose,
	PoseBalls const & pb,
	HolesParamsRes const & params
);

Real
compute_holes_deriv_res(
	pose::Pose const & pose,
	PoseBalls const & pb,
	HolesParamsRes const & params,
	core::id::AtomID_Map< numeric::xyzVector<core::Real> > & derivs
);


}
}
}

#endif // INCLUDED_core_scoring_packing_compute_holes_score_HH
