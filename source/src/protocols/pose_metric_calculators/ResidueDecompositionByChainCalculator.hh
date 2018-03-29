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


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_ResidueDecompositionByChainCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_ResidueDecompositionByChainCalculator_hh
#include <protocols/pose_metric_calculators/ResidueDecompositionCalculator.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace pose_metric_calculators {

class ResidueDecompositionByChainCalculator : public ResidueDecompositionCalculator {

public:

	ResidueDecompositionByChainCalculator();

	ResidueDecompositionByChainCalculator(
		ResidueDecompositionByChainCalculator const & calculator
	);

	virtual
	core::pose::metrics::PoseMetricCalculatorOP
	clone() const;

	utility::vector1<std::set<char> > const &
	chain_letters() const {
		return chain_letters_;
	}

	void
	chain_letters(utility::vector1<std::set<char> > const & chain_letters) {
		chain_letters_ = chain_letters;
		use_numbers_ = false;
	}

	utility::vector1<std::set<core::Size> > const &
	chain_numbers() const {
		return chain_numbers_;
	}

	void
	chain_numbers(utility::vector1<std::set<core::Size> > const & chain_numbers) {
		chain_numbers_ = chain_numbers;
		use_numbers_ = true;
	}

	bool
	use_numbers() const {
		return use_numbers_;
	}

	void
	use_numbers(
		bool use_numbers
	) {
		use_numbers_ = use_numbers;
	}

protected:

	virtual void recompute( core::pose::Pose const & this_pose );

	utility::vector1<std::set<char> > chain_letters_;
	utility::vector1<std::set<core::Size> > chain_numbers_;
	bool use_numbers_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef utility::pointer::shared_ptr< ResidueDecompositionByChainCalculator > ResidueDecompositionByChainCalculatorOP;
typedef utility::pointer::shared_ptr< ResidueDecompositionByChainCalculator const > ResidueDecompositionByChainCalculatorCOP;


} // namespace pose_metric_calculators
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_ResidueDecompositionByChainCalculator )
#endif // SERIALIZATION


#endif
