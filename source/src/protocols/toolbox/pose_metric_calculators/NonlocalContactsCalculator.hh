// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NonlocalContactsCalculator.hh
/// @brief
/// @author Florian Richter


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_NonlocalContactsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_NonlocalContactsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/graph/Graph.fwd.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>

#include <set>


// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>


namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class NonlocalContactsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	NonlocalContactsCalculator(
		core::Size min_sequence_separation = basic::options::option[basic::options::OptionKeys::pose_metrics::min_sequence_separation],
		core::Real contact_cutoffE = basic::options::option[basic::options::OptionKeys::pose_metrics::contact_cutoffE]
	);


	NonlocalContactsCalculator(
		std::set< core::Size > const & special_region,
		core::Size min_sequence_separation = basic::options::option[basic::options::OptionKeys::pose_metrics::min_sequence_separation],
		core::Real contact_cutoffE = basic::options::option[basic::options::OptionKeys::pose_metrics::contact_cutoffE]
	);

	NonlocalContactsCalculator(
		std::set< core::Size > const & special_region1,
		std::set< core::Size > const & special_region2,
		core::Size min_sequence_separation = basic::options::option[basic::options::OptionKeys::pose_metrics::min_sequence_separation],
		core::Real contact_cutoffE = basic::options::option[basic::options::OptionKeys::pose_metrics::contact_cutoffE]
	);

	~NonlocalContactsCalculator();


	core::pose::metrics::PoseMetricCalculatorOP clone() const {
		return core::pose::metrics::PoseMetricCalculatorOP( new NonlocalContactsCalculator( special_region1_, special_region2_, min_seq_separation_, cutoffE_) ); };

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );


private:

	core::Size total_nlcontacts_;
	core::Size special_region1_nlcontacts_;
	core::Size special_region2_nlcontacts_;
	core::Size special_region1_intra_nlcontacts_;
	core::Size special_region1_to_other_nlcontacts_;
	core::Size region1_region2_nlcontacts_;

	utility::vector1< core::Size > residue_nlcontacts_;
	utility::vector1< core::Real > residue_nlscore_;

	core::graph::GraphOP nlcontacts_graph_;

	//how far two residues need to be apart in sequence to count as nonlocal
	core::Size min_seq_separation_;

	//minimum energy between two residues to count as interacting
	core::Real cutoffE_;

	std::set< core::Size > special_region1_;
	std::set< core::Size > special_region2_;

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
