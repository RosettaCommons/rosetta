// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// How many salt bridge interactions are there?
///
/// @details
/// Not much detailed here. Iterate through the oxygens of acidic residues and compare that to
/// the distance of the polar hydrogens in basic residues. Default distance is 3.2A.
/// Wait, you want to know how to use this? Well, within your protocol, you need to do the following:
/// First, create the calculator. To do this, see below:
/// core::pose::metrics::PoseMetricCalculatorOP sb_calculator = new protocols::pose_metric_calculators::SaltBridgeCalculator();
/// Then you must register this so that the pose understands it. See below:
/// core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sb_metric", sb_calculator );
/// To actually get the metric, you have to print it. For example:
/// core::pose::Pose pose;
/// pose.print_metric("sb_metric", "salt_bridge")
/// Where sb_metric is the name that it is registered under and "salt_bridge" is the key, seen below.
///
///
///
/// @author
/// Steven Combs
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_SaltBridgeCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_SaltBridgeCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
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

class SaltBridgeCalculator : public core::pose::metrics::StructureDependentCalculator{
public:
	//default constructor where distance_cutoff is = to 3.2
	SaltBridgeCalculator();

	//constructor where you define what the distance cutoff is for the Hydrogen and Acceptor atoms
	SaltBridgeCalculator(core::Real dist_cutoff);

	core::pose::metrics::PoseMetricCalculatorOP clone() const {
		return core::pose::metrics::PoseMetricCalculatorOP( new SaltBridgeCalculator( distance_cutoff_) ); };

private:
	core::Real distance_cutoff_; //distance cutoff between the Hydrogen and Acceptor atoms. Default is 3.2
	core::Size salt_bridge_total_;

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


}
}


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_SaltBridgeCalculator )
#endif // SERIALIZATION


#endif
