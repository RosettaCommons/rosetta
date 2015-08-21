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
/// @author John Karanicolas


#ifndef INCLUDED_core_pose_metrics_simple_calculators_InterfaceDefinitionCalculatorBase_HH
#define INCLUDED_core_pose_metrics_simple_calculators_InterfaceDefinitionCalculatorBase_HH

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace pose {
namespace metrics {
namespace simple_calculators {

class InterfaceDefinitionCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

	InterfaceDefinitionCalculator( core::Size const chain1_number, core::Size const chain2_number );

	InterfaceDefinitionCalculator( char const chain1_letter, char const chain2_letter );

	core::pose::metrics::PoseMetricCalculatorOP clone() const = 0;

protected:

	virtual std::string print( std::string const & key ) const = 0;

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const = 0;

	virtual void recompute( core::pose::Pose const & this_pose ) = 0;

	core::Size ch1_begin_num_, ch1_end_num_, ch2_begin_num_, ch2_end_num_;

	core::Size chain1_number_, chain2_number_;
	char chain1_letter_, chain2_letter_;
	bool got_chain_numbers_;

	virtual void verify_chain_setup( core::pose::Pose const & pose );

	virtual core::Size chain_letter_to_number( core::pose::Pose const & pose, char const chain_id );

	virtual void fill_in_chain_terminii( core::pose::Pose const & pose );

};


} // namespace simple_calculators
} // namespace metrics
} // namespace pose
} // namespace core

#endif //INCLUDED_core_pose_metrics_simple_calculators_InterfaceDefinitionCalculatorBase_HH
