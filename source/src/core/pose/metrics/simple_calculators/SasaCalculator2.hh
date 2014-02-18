// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/metrics/simple_calculators/SasaCalculator2
/// @brief Calculator for SASA.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com).  Based on SasaCalculatorLegacy


#ifndef INCLUDED_core_pose_metrics_simple_calculators_SasaCalculator2_HH
#define INCLUDED_core_pose_metrics_simple_calculators_SasaCalculator2_HH

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/scoring/sasa/SasaCalc.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>

#include <basic/options/option.hh>

#include <utility/vector1.hh>


// option key includes

#include <basic/options/keys/sasa.OptionKeys.gen.hh>


namespace core{
namespace pose {
namespace metrics {
namespace simple_calculators {
	using utility::vector1;
	
class SasaCalculator2 : public core::pose::metrics::StructureDependentCalculator {

public:

	SasaCalculator2( core::Real probe_r = basic::options::option[basic::options::OptionKeys::sasa::probe_radius]() );

	core::pose::metrics::PoseMetricCalculatorOP clone() const { return new core::pose::metrics::simple_calculators::SasaCalculator2(); };

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:
	
	core::scoring::sasa::SasaCalcOP sasa_calc_;
	
	//These need to be specified and returned from sasacalc or else clang gives error: taking the address of a temporary object of type...
	
	core::Real total_sasa_;
	core::Real total_hsasa_;
	core::id::AtomID_Map<core::Real > atom_sasa_;
	vector1< core::Real > residue_sasa_;
	vector1< core::Real > residue_hsasa_;
	vector1< core::Real > residue_rel_hsasa_;
	
};


} // namespace simple_calculators
} // namespace metrics
} // namespace pose
} // namespace core

#endif //INCLUDED_core_pose_metrics_simple_calculators_SasaCalculator2_HH


