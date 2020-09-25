// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/drug_design/RDKitMetric.hh
/// @brief A SimpleMetric which measures properties calcualted by RDKit on a ligand.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_RDKitMetric_HH
#define INCLUDED_protocols_drug_design_RDKitMetric_HH

#include <protocols/drug_design/RDKitMetric.fwd.hh>
#include <core/simple_metrics/RealMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace drug_design {

/// @brief A SimpleMetric which measures properties calcualted by RDKit on a ligand.
class RDKitMetric : public core::simple_metrics::RealMetric {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	RDKitMetric();

	RDKitMetric( core::select::residue_selector::ResidueSelectorCOP residue, std::string const & metric_name );

	void residue_selector( core::select::residue_selector::ResidueSelectorCOP residue ) { residue_ = residue; }
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_; }

	void rdkit_metric( std::string const & setting );
	std::string const & rdkit_metric() const { return rdkit_metric_; }


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in RealMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( core::pose::Pose & pose, prefix="", suffix="" ) override;

	/// @brief Calculate the metric.
	core::Real
	calculate( core::pose::Pose const & pose ) const override;

public:

	/// @brief Name of the class
	std::string
	name() const override;

	/// @brief Name of the class for creator.
	static
	std::string
	name_static();

	/// @brief Name of the metric
	std::string
	metric() const override;

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

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

private:
	core::select::residue_selector::ResidueSelectorCOP residue_;
	std::string rdkit_metric_;
};

} //drug_design
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_drug_design_RDKitMetric )
#endif // SERIALIZATION

#endif //protocols_drug_design_RDKitMetric_HH





