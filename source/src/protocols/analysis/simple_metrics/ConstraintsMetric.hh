// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/simple_metrics/ConstraintsMetric.hh
/// @brief A simple metric that writes out all the constraints in a pose or sub-region of a pose,
/// in a format matching the (non-enzdes) constraints file format.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetric_HH
#define INCLUDED_protocols_analysis_simple_metrics_ConstraintsMetric_HH

#include <protocols/analysis/simple_metrics/ConstraintsMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace analysis {
namespace simple_metrics {

/// @brief A simple metric that writes out all the constraints in a pose or sub-region of a pose,
/// in a format matching the (non-enzdes) constraints file format.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class ConstraintsMetric : public core::simple_metrics::StringMetric {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ConstraintsMetric();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ConstraintsMetric() override;

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///@brief Calculate the metric.
	std::string
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Name of the class
	std::string
	name() const override;

	///@brief Name of the class for creator.
	static
	std::string
	name_static();

	///@brief Name of the metric
	std::string
	metric() const override;

	/// @brief This simple metric is unpublished.  It returns Vikram K. Mulligan as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

public:

	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

public: //Setters

	/// @brief Set the residue selector.  Use nullptr to clear the residue selector.
	/// @details Does not clone the input residue selector, but stores the const owning pointer directly.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

public: //Getters

	/// @brief Get the residue selector.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; }

private: //Calculators

	/// @brief Of the residues that a particular constraint acts on, do any lie in the selected residues?
	bool
	selected_residue_included_in_cst(
		core::scoring::constraints::Constraint const & cst,
		core::select::residue_selector::ResidueSubset const & residue_selection
	) const;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

private: //Data

	/// @brief A residue selector.  If provided, those constraints that involve at least one residue
	/// in the selection are written out.  If not provided, all constraints in the pose are written out.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

};

} //simple_metrics
} //analysis
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_analysis_simple_metrics_ConstraintsMetric )
#endif // SERIALIZATION


#endif //protocols_analysis_simple_metrics_ConstraintsMetric_HH
