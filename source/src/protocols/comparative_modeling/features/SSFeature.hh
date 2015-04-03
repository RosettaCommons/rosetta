// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_features_SSFeature_HH
#define INCLUDED_protocols_comparative_modeling_features_SSFeature_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/comparative_modeling/features/ResidueFeature.hh>
#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>
#include <protocols/comparative_modeling/features/SSFeature.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace comparative_modeling {
namespace features {

enum SSType {
	H_SS = 1,
	E_SS,
	L_SS,
	INVALID_SS,
	n_ss_types = INVALID_SS
};

class SSFeature : public ResidueFeature {
public:
	SSFeature();
	SSFeature( SSFeature const & other );
	SSFeature( SSType bin );

	utility::vector1< ResidueFeatureOP >
	values_from_pose( core::pose::Pose & pose ) const;

	std::string type() const;

	ResidueFeatureOP clone() const;

	SSType ss_type() const;

	static SSType char2ss_type(
		char const ss
	);

private:
	SSType ss_type_;
}; // SSFeature

} // features
} // comparative_modeling
} // protocols

#endif
