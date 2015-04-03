// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packing/compute_holes_score.hh
/// @brief  Packing Score
/// @author Will Sheffler


#ifndef INCLUDED_core_scoring_packing_compute_holes_score_hh
#define INCLUDED_core_scoring_packing_compute_holes_score_hh


// you cannot #include yourself #include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packing/PoseBalls.hh>
//#include <core/scoring/ScoringManager.hh>
#include <iomanip>
#include <iostream>
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <utility/exit.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packing {

	/// the result class holding the three scores and the per-atom scores
	class HolesResult : public utility::pointer::ReferenceCount {
	public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~HolesResult();
		HolesResult() : score(0.0),decoy_score(0.0),resl_score(0.0),dec15_score(0.0) {}
		core::Real score, decoy_score, resl_score, dec15_score, natom;
		core::id::AtomID_Map< core::Real > atom_scores;
	};

	/// for the standard scores
	HolesResult
	compute_rosettaholes_score(
		pose::Pose const & pose
	);

	Real
	compute_dec15_score(
		pose::Pose const & pose
	);

	/// if you have custom parameters, or want per-atom scores for a specific score
	/// result goes into the "score" field
	HolesResult
	compute_holes_score(
		pose::Pose  const & pose,
		HolesParams const & params
	);

	/// computes the cartesian space derivative WRT the given params
	HolesResult
	compute_holes_deriv(
		pose::Pose  const & pose,
		HolesParams const & params,
		core::id::AtomID_Map< numeric::xyzVector<core::Real> > & deriv
	);

//////

	HolesResult
	compute_holes_deriv(
		pose::Pose  const & pose,
		PoseBalls         & pb,
		HolesParams const & params,
		core::id::AtomID_Map< numeric::xyzVector<core::Real> > & deriv
	);

	HolesResult
	compute_rosettaholes_score(
		pose::Pose const & pose,
		PoseBalls & pb,
		HolesParams const & resl_params,
		HolesParams const & dec_params,
		HolesParams const & dec15_params,
		bool use_cached_surfs = false,
		std::string cmd = ""
	);

   HolesResult
   compute_rosettaholes_score(
	   pose::Pose const & pose,
		PoseBalls & pb
   );

	HolesResult
	compute_rosettaholes_score(
		pose::Pose  const & pose,
		HolesParams const & resl_params,
		HolesParams const & dec_params,
		HolesParams const & dec15_params
   );

   HolesResult
   compute_holes_score(
	   pose::Pose  const & pose,
	   PoseBalls         & pb,
	   HolesParams const & params,
		 bool use_cached_surfs = false,
		 std::string cmd = ""
	);

   HolesResult
   compute_holes_score(
	   pose::Pose  const & pose,
		 std::string const & cmd
   );


}
}
}

#endif // INCLUDED_core_scoring_packing_compute_holes_score_HH
