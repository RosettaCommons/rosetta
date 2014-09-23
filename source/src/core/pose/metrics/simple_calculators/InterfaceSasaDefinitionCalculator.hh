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


#ifndef INCLUDED_core_pose_metrics_simple_calculators_InterfaceSasaDefinitionCalculator_HH
#define INCLUDED_core_pose_metrics_simple_calculators_InterfaceSasaDefinitionCalculator_HH

#include <core/pose/metrics/simple_calculators/InterfaceDefinitionCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

#include <set>

#include <utility/vector1.hh>


namespace core{
namespace pose {
namespace metrics {
namespace simple_calculators {

class InterfaceSasaDefinitionCalculator : public InterfaceDefinitionCalculator {

public:

	InterfaceSasaDefinitionCalculator( core::Size const chain1_number, core::Size const chain2_number ) :
		InterfaceDefinitionCalculator(chain1_number, chain2_number) {};

	InterfaceSasaDefinitionCalculator( char const chain1_letter, char const chain2_letter ) :
		InterfaceDefinitionCalculator(chain1_letter, chain2_letter) {};

	core::pose::metrics::PoseMetricCalculatorOP clone() const
	{ return core::pose::metrics::PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator( chain1_number_, chain2_number_ ) ); }

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	core::Real delta_sasa_;
	core::Real fraction_chain1_delta_sasa_;
	core::Real fraction_chain2_delta_sasa_;
	core::id::AtomID_Map< core::Real > atom_delta_sasa_;
	utility::vector1< core::Real > residue_delta_sasa_;
	core::id::AtomID_Map< bool > interface_atoms_;
	std::set< core::Size > interface_residues_;
	std::set< core::Size > chain1_interface_residues_;
	std::set< core::Size > chain2_interface_residues_;
	core::Size num_interface_residues_;
	core::Size num_chain1_interface_residues_;
	core::Size num_chain2_interface_residues_;

};


} // namespace simple_calculators
} // namespace metrics
} // namespace pose
} // namespace core

#endif //INCLUDED_core_pose_metrics_simple_calculators_InterfaceSasaDefinitionCalculator_HH
