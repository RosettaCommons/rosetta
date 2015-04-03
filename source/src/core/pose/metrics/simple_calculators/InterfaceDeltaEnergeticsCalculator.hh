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
/// @author Roland A Pache


#ifndef INCLUDED_core_pose_metrics_simple_calculators_InterfaceDeltaEnergeticsCalculator_HH
#define INCLUDED_core_pose_metrics_simple_calculators_InterfaceDeltaEnergeticsCalculator_HH

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/vector1.hh>
#include <core/scoring/ScoreType.hh>


namespace core{
namespace pose {
namespace metrics {
namespace simple_calculators {

class InterfaceDeltaEnergeticsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	// preferred constructor - use an existing InterfaceNeighborDefinitionCalculator
	InterfaceDeltaEnergeticsCalculator( std::string const & NameOfInterfaceNeighborDefinitionCalculator );
    
    // preferred alternative constructor - use an existing InterfaceNeighborDefinitionCalculator and define a set of score types to ignore
	InterfaceDeltaEnergeticsCalculator( std::string const & NameOfInterfaceNeighborDefinitionCalculator, utility::vector1<core::scoring::ScoreType> const & score_types_to_ignore );

	// less preferred constructor - create a new InterfaceNeighborDefinitionCalculator
	InterfaceDeltaEnergeticsCalculator( core::Size const chain1_number, core::Size const chain2_number );

	// less preferred constructor - create a new InterfaceNeighborDefinitionCalculator
	InterfaceDeltaEnergeticsCalculator( char const chain1_letter, char const chain2_letter );

	core::pose::metrics::PoseMetricCalculatorOP clone() const
	{ return core::pose::metrics::PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( name_of_InterfaceNeighborDefinitionCalculator_, score_types_to_ignore_ ) ); }

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	std::string name_of_InterfaceNeighborDefinitionCalculator_;
    utility::vector1<core::scoring::ScoreType> score_types_to_ignore_;

	core::scoring::EnergyMap delta_energies_unweighted_;
	core::scoring::EnergyMap weights_;
	core::Real weighted_total_;

};


} // namespace simple_calculators
} // namespace metrics
} // namespace pose
} // namespace core

#endif //INCLUDED_core_pose_metrics_simple_calculators_InterfaceDeltaEnergeticsCalculator_HH
