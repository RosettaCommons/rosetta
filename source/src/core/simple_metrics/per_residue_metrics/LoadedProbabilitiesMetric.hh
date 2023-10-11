// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/LoadedProbabilitiesMetric.hh
/// @brief A class to load a probabilities weights file into a PerResidueProbabilitiesMetric
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifndef INCLUDED_core_simple_metrics_per_residue_metrics_LoadedProbabilitiesMetric_HH
#define INCLUDED_core_simple_metrics_per_residue_metrics_LoadedProbabilitiesMetric_HH

#include <core/simple_metrics/per_residue_metrics/LoadedProbabilitiesMetric.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

///@brief A class to load a probabilities weights file into a PerResidueProbabilitiesMetric
class LoadedProbabilitiesMetric : public core::simple_metrics::PerResidueProbabilitiesMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	LoadedProbabilitiesMetric();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~LoadedProbabilitiesMetric() override;


public:

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///@brief Calculate the metric.
	/// This map is Rosetta Resnum->(AA, value)
	///
	///@details
	/// Return by value as this function can not STORE the result, it only calculates.
	/// Store the result in the pose by using the apply method, which calls this method and stores the result
	/// in the pose as ExtraScoreValues.
	std::map< core::Size, std::map< core::chemical::AA, core::Real >>
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

	/// @brief This simple metric is unpublished.  It returns Moritz Ertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

	///@brief set the filename of the probabilities weights file to be loaded
	void set_filename( std::string const & filename );

	///@brief load amino acid probabilities from a given vectors of lines in format (POSNUM, AA_TYPE, WEIGHT)
	void set_probabilities_from_lines( utility::vector1< std::string > const & lines );

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

	// name of the weights file to load probabilities from
	std::string filename_;
	// probabilities loaded from a file
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> probabilities_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //per_residue_metrics
} //simple_metrics
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_LoadedProbabilitiesMetric )
#endif // SERIALIZATION

#endif //core_simple_metrics_per_residue_metrics_LoadedProbabilitiesMetric_HH
