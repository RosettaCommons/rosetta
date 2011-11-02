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


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_InterfaceDeltaEnergeticsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_InterfaceDeltaEnergeticsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/vector1.hh>


// AUTO-REMOVED #include <core/id/AtomID_Map.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class InterfaceDeltaEnergeticsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	// preferred constructor - use an existing InterfaceNeighborDefinitionCalculator
	InterfaceDeltaEnergeticsCalculator( std::string const & NameOfInterfaceNeighborDefinitionCalculator );

	// less preferred constructor - create a new InterfaceNeighborDefinitionCalculator
	InterfaceDeltaEnergeticsCalculator( core::Size const chain1_number, core::Size const chain2_number );

	// less preferred constructor - create a new InterfaceNeighborDefinitionCalculator
	InterfaceDeltaEnergeticsCalculator( char const chain1_letter, char const chain2_letter );

	core::pose::metrics::PoseMetricCalculatorOP clone() const
	{ return new InterfaceDeltaEnergeticsCalculator( name_of_InterfaceNeighborDefinitionCalculator_ ); }

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	std::string name_of_InterfaceNeighborDefinitionCalculator_;

	core::scoring::EnergyMap delta_energies_unweighted_;
	core::scoring::EnergyMap weights_;
	core::Real weighted_total_;

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
