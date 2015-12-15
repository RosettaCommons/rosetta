// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/SemiExplicitWaterUnsatisfiedPolarsCalculator.hh
/// @brief
/// @author Chris King, templated from BuriedUnsatisfiedPolarsCalulator by Florian Richter


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_SemiExplicitWaterUnsatisfiedPolarsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_SemiExplicitWaterUnsatisfiedPolarsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/conformation/Residue.fwd.hh>


#include <basic/options/option.hh>

#include <utility/vector1.hh>

#include <set>


// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class SemiExplicitWaterUnsatisfiedPolarsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	SemiExplicitWaterUnsatisfiedPolarsCalculator(
		std::string hbond_calc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::Real semiexpl_water_cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::semiex_water_burial_cutoff]
	);


	SemiExplicitWaterUnsatisfiedPolarsCalculator(
		std::string hbond_calc,
		core::scoring::ScoreFunctionOP scorefxn,
		std::set< core::Size > const & special_region,
		core::Real semiexpl_water_cutoff = basic::options::option[basic::options::OptionKeys::pose_metrics::semiex_water_burial_cutoff]
	);


	core::Real
	semiexpl_water_hbgeom_score(
		core::pose::Pose pose,
		core::scoring::ScoreFunctionOP scorefxn,
		core::Size seqpos,
		core::Size atomno,
		core::conformation::Residue new_rsd,
		core::Size new_atomno
	);

	core::pose::metrics::PoseMetricCalculatorOP clone() const {
		return core::pose::metrics::PoseMetricCalculatorOP( new SemiExplicitWaterUnsatisfiedPolarsCalculator( name_of_hbond_calc_, scorefxn_, semiexpl_water_cutoff_ ) ); };

	std::string const & name_of_hbond_calc() const { return name_of_hbond_calc_; }

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );


private:

	void assert_calculators();

	static
	core::Size satisfaction_cutoff( std::string atom_type );


	core::scoring::hbonds::HBondDatabaseCOP hb_database_;
	core::Size all_unsat_polars_;
	core::Size special_region_unsat_polars_;
	core::id::AtomID_Map< bool > atom_unsat_;
	utility::vector1< core::Size > residue_unsat_polars_;
	utility::vector1< core::Real > residue_semiexpl_score_;
	core::id::AtomID_Map< core::Real > atom_semiexpl_score_;
	core::Real semiexpl_water_cutoff_;

	//holds the atom hbonds calculators necessary for this calculator
	std::string name_of_hbond_calc_;
	core::scoring::ScoreFunctionOP scorefxn_;

	std::set< core::Size > special_region_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SemiExplicitWaterUnsatisfiedPolarsCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_SemiExplicitWaterUnsatisfiedPolarsCalculator )
#endif // SERIALIZATION


#endif
