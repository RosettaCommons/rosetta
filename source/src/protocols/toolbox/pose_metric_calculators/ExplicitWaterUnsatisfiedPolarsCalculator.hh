// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin ExplicitWaterUnsatisfiedPolarsCalculator
///
/// @brief This Calculator tries to solvate all polar groups by
///         docking explicit TP5 water molecules, then counts
///         unsatisfied hydrogen bonds using the same criteria in 
///         BuriedUnsatisfiedHydrogenBondCalculator
///       
///
///
///
///
///
///
///
///
///
///
///
///
///
/// @author Chris King - dr.chris.king@gmail.com
/// 
///
/// @last_modified 4.8.2011
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_ExplicitWaterUnsatisfiedPolarsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_ExplicitWaterUnsatisfiedPolarsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class ExplicitWaterUnsatisfiedPolarsCalculator : public core::pose::metrics::StructureDependentCalculator{
public:
	//default constructor shell_cutoff is = to 4.0
	ExplicitWaterUnsatisfiedPolarsCalculator( core::scoring::ScoreFunctionOP scorefxn );

	//constructor where you define what the distance cutoff is for the Hydrogen and Acceptor atoms
	ExplicitWaterUnsatisfiedPolarsCalculator( core::scoring::ScoreFunctionOP scorefxn, core::Real shell_cutoff );

	core::pose::metrics::PoseMetricCalculatorOP clone() const {
	    return new ExplicitWaterUnsatisfiedPolarsCalculator( scorefxn_, shell_cutoff_ ); };

private:
	core::Real shell_cutoff_; //water shell approx cutoff (for initial position in docking)
	core::Size all_unsat_polars_;
	core::scoring::ScoreFunctionOP scorefxn_;
protected:
	  virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	  virtual std::string print( std::string const & key ) const;
	  virtual void recompute( core::pose::Pose const & this_pose );

};




}
}
}






#endif
