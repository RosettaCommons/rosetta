// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SimpleMetricSelector.hh
/// @brief  Allows selecting residues based on the result of a PerResidueRealSimpleMetric
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_core_select_residue_selector_SimpleMetricSelector_HH
#define INCLUDED_core_select_residue_selector_SimpleMetricSelector_HH

// Unit headers
#include <core/select/residue_selector/SimpleMetricSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/simple_metrics/PerResidueRealMetric.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/numbers.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

// forward declaration for testing
class SimpleMetricSelectorTests;

namespace core {
namespace select {
namespace residue_selector {



/// @brief The SimpleMetricSelector allows conditional application of a residue selector
class SimpleMetricSelector : public ResidueSelector {
public:
	friend class ::SimpleMetricSelectorTests;

	/// @brief Copy constructor
	///
	SimpleMetricSelector( SimpleMetricSelector const &src );

	// derived from base class
	~SimpleMetricSelector() override;

	SimpleMetricSelector & operator=( SimpleMetricSelector const & ot );

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	ResidueSelectorOP clone() const override;

	SimpleMetricSelector(
		core::simple_metrics::PerResidueRealMetricCOP metric = nullptr,
		Real lower_bound = utility::get_undefined_real(),
		Real upper_bound = utility::get_undefined_real(),
		bool outside_bounds = false
	);

	ResidueSubset apply( core::pose::Pose const & pose ) const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	set_metric( core::simple_metrics::PerResidueRealMetricCOP const & metric );

private:
	core::simple_metrics::PerResidueRealMetricCOP metric_;
	Real lower_bound_;
	Real upper_bound_;
	bool outside_bounds_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_SimpleMetricSelector )
#endif // SERIALIZATION


#endif
