// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Aliza Rubenstein aliza.rubenstein@gmail.com


#ifndef INCLUDED_protocols_mean_field_EnergiesByTaskCalculator_hh
#define INCLUDED_protocols_mean_field_EnergiesByTaskCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/vector1.hh>

#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace mean_field {

class EnergiesByTaskCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	// preferred constructor - use an existing InterfaceNeighborDefinitionCalculator
	EnergiesByTaskCalculator(
		core::pack::task::PackerTaskCOP task
	);

	EnergiesByTaskCalculator(
		EnergiesByTaskCalculator const & calculator
	);

	core::pose::metrics::PoseMetricCalculatorOP clone() const;

	core::Real total() const { return total_score_; }
	core::pack::task::PackerTaskCOP task() const { return task_; }

	void
	show(
		std::ostream & out
	) const;

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	core::pack::task::PackerTaskCOP task_;
	core::Real total_score_;
};

typedef utility::pointer::shared_ptr< EnergiesByTaskCalculator > EnergiesByTaskCalculatorOP;
typedef utility::pointer::shared_ptr< EnergiesByTaskCalculator const > EnergiesByTaskCalculatorCOP;


} // namespace mean_field
} // namespace protocols

#endif
