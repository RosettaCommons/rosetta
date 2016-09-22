// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/RotamerBoltzmannWeight2.hh
/// @brief Next-generation RotamerBoltzmannWeight filter
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_hh
#define INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_hh

// Unit headers
#include <protocols/simple_filters/RotamerBoltzmannWeight2.fwd.hh>
#include <protocols/filters/Filter.hh>

// Protocol headers
#include <protocols/toolbox/EnergyLandscapeEvaluator.fwd.hh>
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/SingletonBase.hh>

// C++ headers
#include <set>

namespace protocols {
namespace simple_filters {

template< class T >
class IdManager : public utility::SingletonBase< IdManager< T > > {
public:
	IdManager();

	bool id_exists( T const & id ) const;
	void register_id( T const & id );
	void unregister_id( T const & id );
	T const & register_new_id();

private:
	std::set< T > used_ids_;
};

///@brief Next-generation RotamerBoltzmannWeight filter
class RotamerBoltzmannWeight2 : public protocols::filters::Filter {
public:
	enum ScoreType {
		MEAN_PROBABILITY,
		MAX_PROBABILITY,
		MODIFIED_DDG
	};
	typedef std::map< core::Size, core::Real > ResidueDDGs;

public:
	RotamerBoltzmannWeight2();

	/// @brief copy constructor -- required because each instance should have unique calculator id
	RotamerBoltzmannWeight2( RotamerBoltzmannWeight2 const & rval );

	// destructor (important for properly forward-declaring smart-pointer members)
	~RotamerBoltzmannWeight2() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public:
	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::filters::FilterOP
	clone() const override;

public:
	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	void
	set_energy_landscape_evaluator( protocols::toolbox::EnergyLandscapeEvaluatorCOP evaluator );

	void
	set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

	void
	set_score_type( ScoreType const scoretype );

public:
	std::string const &
	calculator_id() const;

private:
	core::Real
	compute_score(
		core::pose::Pose const & pose,
		toolbox::pose_metric_calculators::RotamerProbabilities const & probabilities ) const;

	void
	register_calculator() const;

	void
	unregister_calculator();

	static std::string
	new_calculator_id();

	core::Real
	compute_modified_ddg(
		core::pose::Pose const & pose,
		toolbox::pose_metric_calculators::RotamerProbabilities const & probs ) const;

	core::Real
	compute_ddg( core::pose::Pose const & pose ) const;

	core::Real
	compute_modified_residue_energy(
		toolbox::pose_metric_calculators::RotamerProbability const & probability ) const;

private:
	ScoreType score_type_;
	std::string calculator_id_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	protocols::toolbox::pose_metric_calculators::RotamerBoltzCalculatorOP calculator_;

private:
	// constants
	static core::Real const default_temperature_;
	static core::Real const default_lambda_;

};

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_hh
