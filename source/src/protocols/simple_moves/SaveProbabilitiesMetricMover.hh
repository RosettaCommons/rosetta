// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SaveProbabilitiesMetricMover.hh
/// @brief A class to save a PerResidueProbabilitiesMetric to a weight file.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifndef INCLUDED_protocols_simple_moves_SaveProbabilitiesMetricMover_HH
#define INCLUDED_protocols_simple_moves_SaveProbabilitiesMetricMover_HH

// Unit headers
#include <protocols/simple_moves/SaveProbabilitiesMetricMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

namespace protocols {
namespace simple_moves {

///@brief A class to save a PerResidueProbabilitiesMetric to a weight file.
class SaveProbabilitiesMetricMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SaveProbabilitiesMetricMover();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SaveProbabilitiesMetricMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	///@brief Set the PerResidueProbabilitiesMetric that will be used to calculate the pseudo-perplexity.
	void
	set_metric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric );

	///@brief Set the name of the output file.
	void
	set_filename( std::string const & filename );

	///@brief Set the type of the output file.
	void
	set_filetype( std::string const & filetype );

	///@brief Set a boolean to attempt to find cached data matching the name/custom_type of the passed in simple_metric.
	/// Optionally pass any set prefix/suffix.
	///
	/// This will allow the filter to re-use previously calculated data.
	///
	void
	set_use_cached_data( bool use_cache, std::string const & prefix="", std::string const & suffix="");

	///@brief If use_cache is set to false, do we fail if no data is found in the pose?
	/// Default True
	void
	set_fail_on_missing_cache(bool fail);

	///@brief Save a set of amino acid probabilities to a weight file with formatting: PoseNum ResidueType Weight
	static void save_aa_probabilities_to_file( std::string const & weights_file, std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & prob_set);

	///@brief Get the sequence of the selection present in the probability map
	static std::string get_selection_sequence(
		std::string const &pose_sequence,
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> const &position_map
	);

	///@brief convert the map of probabilities to one of logits for writing out to a PSSM
	static std::map< core::Size, std::map<core::chemical::AA, core::Real>> convert_probabilities_to_logits(
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> const &probabilities_map);

public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//SaveProbabilitiesMetricMover & operator=( SaveProbabilitiesMetricMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: //Function overrides needed for the citation manager:

	/// @brief This mover is unpublished.  It returns moritzertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

private: // methods

private: // data

	// metric to be used
	core::simple_metrics::PerResidueProbabilitiesMetricCOP metric_ = nullptr;

	// filename or output file
	std::string filename_;
	// filetype, either weights or pssm file
	std::string filetype_ = "weights";

	//Accessing cached data from the pose
	bool use_cache_ = false;
	std::string cache_prefix_;
	std::string cache_suffix_;
	bool fail_on_missing_cache_ = true;
};

std::ostream &
operator<<( std::ostream & os, SaveProbabilitiesMetricMover const & mover );

} //simple_moves
} //protocols

#endif //protocols_simple_moves_SaveProbabilitiesMetricMover_HH
