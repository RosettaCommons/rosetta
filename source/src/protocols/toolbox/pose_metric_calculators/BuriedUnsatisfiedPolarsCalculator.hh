// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/BuriedUnsatisfiedPolarsCalculator.hh
/// @brief calculates buried polar atoms that are not satisfied by hydrogen bonds
/// @details significantly updated in 2017 to have more generous definition of h-bonds
///  (previously, legit h-bonds were excluded because of sfxn exceptions); users can now choose
///  between legacy SASA and VSASA for burial; poses with more than 3 chains now supported; the
///  way unsats are counted and reported is now different (before, all unsats were counted as equal,
///  which is misleading); users can choose different reporting behaviours; legacy options can be
///  restored by setting legacy=true, but this is only recommended for benchmarking purposes
/// @author Florian Richter
/// @author Scott Boyken (sboyken@gmail.com), major updates and refactoring

#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_BuriedUnsatisfiedPolarsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_BuriedUnsatisfiedPolarsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>

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

class BuriedUnsatisfiedPolarsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	///@brief default constructor ( keeping this for historical reasons and compatibility with other code that calls this )
	BuriedUnsatisfiedPolarsCalculator(
		std::string const & sasa_calc,
		std::string const & hbond_calc,
		core::Real const burial_cutoff = -1.0, // must be 3rd to ensure legacy compatible // -1.0 dummy so we check if user defined
		bool const generous = true,
		bool const legacy = false,
		bool const vsasa = true
	);

	///@brief default constructor with special_region ( keeping this for historical reasons and compatibility with other code that calls this )
	BuriedUnsatisfiedPolarsCalculator(
		std::string const & sasa_calc,
		std::string const & hbond_calc,
		std::set< core::Size > const & special_region, // must be 3rd
		core::Real const burial_cutoff = -1.0, // must be 3rd to ensure legacy compatible // -1.0 dummy so we check if user defined
		bool const generous = true,
		bool const legacy = false,
		bool const vsasa = true
	);

	///@brief constructor where all options must be defined, for calling from code
	BuriedUnsatisfiedPolarsCalculator(
		std::string const & sasa_calc,
		std::string const & hbond_calc,
		std::set< core::Size > const & special_region,
		core::Real const burial_cutoff,
		core::Real const probe_r,
		core::Real const residue_surface_cutoff,
		bool const generous,
		bool const legacy,
		bool const vsasa,
		bool const use_sc_neighbors,
		bool const skip_surface_res
	);


	///@brief copy constructor implicit in clone call
	core::pose::metrics::PoseMetricCalculatorOP clone() const {
		return core::pose::metrics::PoseMetricCalculatorOP( new BuriedUnsatisfiedPolarsCalculator( name_of_sasa_calc_, name_of_hbond_calc_, special_region_, burial_cutoff_, probe_radius_, residue_surface_cutoff_, generous_hbonds_, legacy_counting_, vsasa_, use_sc_neighbors_, skip_surface_res_ ) ); };

	std::string const & name_of_hbond_calc() const { return name_of_hbond_calc_; }
	std::string const & name_of_sasa_calc() const { return name_of_sasa_calc_; }

	void set_special_region( std::set< core::Size > const & special_region ) {
		special_region_ = special_region;
		this->notify_energy_change();
	}


protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );


private:

	void assert_calculators();

	///@brief ONLY USED FOR LEGACY BEHAVIOR (legacy=true), cutoff to determine if unsat
	static core::Size satisfaction_cutoff( std::string atom_type );

	core::Size all_bur_unsat_polars_;
	core::Size bb_heavy_unsats_;
	core::Size all_heavy_unsats_;
	core::Size countable_nonheavy_unsats_;
	core::id::AtomID_Map< bool > atom_bur_unsat_;
	utility::vector1< core::Size > residue_bur_unsat_polars_;
	//holds the sasa and atom hbonds calculators necessary for this calculator
	std::string name_of_hbond_calc_, name_of_sasa_calc_;
	std::set< core::Size > special_region_;
	core::Real burial_cutoff_;
	core::Real probe_radius_;
	core::Real residue_surface_cutoff_;
	bool generous_hbonds_; // if false, will only count the h-bonds that the scorefxn does
	bool legacy_counting_; // revert to legacy options
	bool vsasa_;
	bool use_sc_neighbors_;
	bool skip_surface_res_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	BuriedUnsatisfiedPolarsCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_BuriedUnsatisfiedPolarsCalculator )
#endif // SERIALIZATION


#endif
