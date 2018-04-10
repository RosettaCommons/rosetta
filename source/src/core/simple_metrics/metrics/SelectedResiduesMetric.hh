// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/metrics/SelectedResiduesMetric.hh
/// @brief Output residue-selected residues to a score file as rosetta resnums or pdbnums.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_metrics_SelectedResiduesMetric_HH
#define INCLUDED_core_simple_metrics_metrics_SelectedResiduesMetric_HH

#include <core/simple_metrics/metrics/SelectedResiduesMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace simple_metrics {
namespace metrics {

///@brief Output residue-selected residues to a score file as rosetta resnums or pdbnums.
class SelectedResiduesMetric : public core::simple_metrics::StringMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SelectedResiduesMetric();

	SelectedResiduesMetric( select::residue_selector::ResidueSelectorCOP selector );


	/// @brief Copy constructor (not needed unless you need deep copies)
	SelectedResiduesMetric( SelectedResiduesMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SelectedResiduesMetric() override;

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
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Output a comma-separated list of residues either in PDB (default) or Rosetta numbers.
	std::string
	calculate( core::pose::Pose const & pose ) const override;

public:

	///@brief Set the required residue selector.
	void
	set_residue_selector( select::residue_selector::ResidueSelectorCOP selector);

	///@brief Set the mode in which we will output the selected residues.
	void
	set_output_in_rosetta_num( bool output_in_rosetta_num);

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

private:

	select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	bool rosetta_nums_ = false;

};

} //core
} //simple_metrics
} //metrics



#endif //core_simple_metrics_metrics_SelectedResiduesMetric_HH





