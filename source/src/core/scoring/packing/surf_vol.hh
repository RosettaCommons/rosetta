// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/surf_vol.hh
/// @brief  Packing Score
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_packing_surf_vol_hh
#define INCLUDED_core_scoring_packing_surf_vol_hh

//Unit headers

//Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>

//numeric headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packing {

	struct SurfVol {
		Real tot_surf, tot_vol;
		core::id::AtomID_Map< core::Real> surf;
		core::id::AtomID_Map< core::Real> vol;
	};

	struct SurfVolDeriv {
		Real tot_surf, tot_vol;
		core::id::AtomID_Map< core::Real > surf;
		core::id::AtomID_Map< core::Real > vol;
		core::id::AtomID_Map< numeric::xyzVector<core::Real> > dsurf;
		core::id::AtomID_Map< numeric::xyzVector<core::Real> > dvol;
	};

	Real
	get_surf_tot(
		pose::Pose const & pose,
		core::Real const   probe_radius
	);

	SurfVol
	get_surf_vol(
		pose::Pose const & pose,
		core::Real const   probe_radius = 1.4
	);

	SurfVol
	get_surf_vol(
		pose::Pose const & pose,
		core::id::AtomID_Mask const & whichatoms,
		core::Real const   probe_radius = 1.4
	);

	SurfVolDeriv
	get_surf_vol_deriv(
		pose::Pose const & pose,
		core::Real const   probe_radius = 1.4
	);

	SurfVolDeriv
	get_surf_vol_deriv(
		pose::Pose const & pose,
		core::id::AtomID_Mask const & whichatoms,
		core::Real const   probe_radius = 1.4
	);


}
}
}

#endif // INCLUDED_core_scoring_packing_surf_vol_HH
