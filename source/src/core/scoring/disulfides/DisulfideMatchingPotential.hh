// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideMatchingPotential.hh
/// @brief  Centroid Disulfide Energy Potentials
/// @author rvernon@u.washington.edu
/// @date   02/10/10

#ifndef INCLUDED_core_scoring_disulfides_DisulfideMatchingPotential_hh
#define INCLUDED_core_scoring_disulfides_DisulfideMatchingPotential_hh

//Unit headers
#include <core/scoring/disulfides/DisulfideMatchingPotential.fwd.hh>
#include <core/scoring/disulfides/DisulfideMatchingDatabase.hh>
#include <utility/pointer/ReferenceCount.hh>

//Project headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>
// AUTO-REMOVED #include <core/scoring/func/Func.hh>

//Utility Headers
// AUTO-REMOVED #include <numeric/interpolation/Histogram.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>

namespace core {
namespace scoring {
namespace disulfides {

/**
 * @details This class scores centroid disulfide bonds
 * It is intended to be a singleton with a single instance held by ScoringManager.
 *
 * The energy functions are derived from those present in Rosetta++
 */
class DisulfideMatchingPotential : public utility::pointer::ReferenceCount {
public:
	DisulfideMatchingPotential();
	virtual ~DisulfideMatchingPotential();

	/**
	 * @brief Calculates scoring terms for the disulfide bond specified
	 */
	void
	score_disulfide(core::conformation::Residue const & res1,
			core::conformation::Residue const & res2,
									core::Energy & match_t,
									core::Energy & match_r,
									core::Energy & match_rt

			) const;

	// Not used by scoring machinery, exists so that other apps can compute the score directly
	Energy compute_matching_energy( pose::Pose const & pose ) const;


private:

	/**
	 * @brief calculates RT object for residue pair
	 */
	core::kinematics::RT
	disulfide_RT(
			core::conformation::Residue const& res1,
			core::conformation::Residue const& res2) const;

private:

	mutable disulfides::DisulfideMatchingDatabase matching_database_; //HACK!

}; // DisulfideMatchingPotential


// RT Helper makes an RT object out of two Epos arrays
// I should probably move this somewhere more appropriate...
// -rvernon@u.washington.edu
class RT_helper {
public:

static
core::kinematics::RT
RT_from_epos( ObjexxFCL::FArray2A_float Epos1, ObjexxFCL::FArray2A_float Epos2);

private:

static
numeric::xyzMatrix_double
get_ncac ( ObjexxFCL::FArray2A_float pos );

static
void
get_ncac(
	ObjexxFCL::FArray2A_float pos,
	numeric::xyzMatrix_double & p
	);

static
void
get_coordinate_system(
	numeric::xyzMatrix_double const & p, //FArray2A_double p, // input
	numeric::xyzMatrix_double & m //FArray2A_double m // output
	);

}; // RT_helper

} // disulfides
} // scoring
} // core

#endif
