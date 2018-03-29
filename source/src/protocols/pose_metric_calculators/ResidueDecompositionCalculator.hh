// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Colin A. Smith


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_ResidueDecompositionCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_ResidueDecompositionCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <set>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_metric_calculators {

class ResidueDecompositionCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

	core::pose::metrics::PoseMetricCalculatorOP clone() const = 0;

protected:

	ResidueDecompositionCalculator();

	ResidueDecompositionCalculator( ResidueDecompositionCalculator const & calculator );

	virtual std::string print( std::string const & key ) const;

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;

	virtual void recompute( core::pose::Pose const & this_pose ) = 0;

	void
	residue_decomposition_to_set_numbers(
		core::pose::Pose const & this_pose
	);

	void
	residue_set_numbers_to_decomposition();

	utility::vector1<std::set<core::Size> > const &
	residue_decomposition() const {
		return residue_decomposition_;
	}

	utility::vector1<core::Size> const &
	residue_set_numbers() const {
		return residue_set_numbers_;
	}

	utility::vector1<std::set<core::Size> > residue_decomposition_;
	utility::vector1<core::Size> residue_set_numbers_;
	utility::vector1<std::string> set_names_;
#ifdef    SERIALIZATION
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace pose_metric_calculators
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_ResidueDecompositionCalculator )
#endif // SERIALIZATION


#endif
