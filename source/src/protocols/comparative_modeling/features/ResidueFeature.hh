// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/features/ResidueFeature.hh
/// @brief abstract base class for per-residue features used in comparative
/// modeling.
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_features_ResidueFeature_HH
#define INCLUDED_protocols_comparative_modeling_features_ResidueFeature_HH

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/id/SequenceMapping.fwd.hh>

#include <protocols/comparative_modeling/features/ResidueFeature.fwd.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <string>
#endif

namespace protocols {
namespace comparative_modeling {
namespace features {

class ResidueFeature : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~ResidueFeature();
	virtual utility::vector1< ResidueFeatureOP >
	values_from_pose( core::pose::Pose & pose ) const = 0;

	virtual ResidueFeatureOP clone() const = 0;
	virtual std::string type() const = 0;

	core::Size resnum() const;

	void resnum( core::Size const resnum );

	void remap( core::id::SequenceMapping const & mapping );

private:
	Size resnum_;
};

} // protocols
} // comparative_modeling
} // features

#endif
