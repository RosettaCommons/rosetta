// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/TorsionFeature.hh
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_features_TorsionFeature_HH
#define INCLUDED_protocols_comparative_modeling_features_TorsionFeature_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/comparative_modeling/features/ResidueFeature.hh>
#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>
#include <protocols/comparative_modeling/features/TorsionFeature.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace comparative_modeling {
namespace features {

enum TorsionBin {
	A = 1,
	B,
	E,
	G,
	O,
	X,
	INVALID,
	n_torsion_bins = INVALID
};

class TorsionFeature : public ResidueFeature {
public:
	TorsionFeature();
	TorsionFeature( TorsionFeature const & other );
	TorsionFeature( TorsionBin bin );

	std::string type() const;

	utility::vector1< ResidueFeatureOP >
	values_from_pose( core::pose::Pose & pose ) const;

	ResidueFeatureOP clone() const;

	TorsionBin torsion_bin() const;


	static TorsionBin
	torsion2big_bin(
		core::Real const phi,
		core::Real const psi,
		core::Real const omega
	);

private:
	TorsionBin torsion_bin_;
}; // TorsionFeature

} // features
} // comparative_modeling
} // protocols

#endif
