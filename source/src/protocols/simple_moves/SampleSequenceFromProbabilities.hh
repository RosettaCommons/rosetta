// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SampleSequenceFromProbabilities.hh
/// @brief A class to sample sequences from a PerResidueProbabilitiesMetric and thread them onto the pose.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com), University of Leipzig

#ifndef INCLUDED_protocols_simple_moves_SampleSequenceFromProbabilities_HH
#define INCLUDED_protocols_simple_moves_SampleSequenceFromProbabilities_HH

// Unit headers
#include <protocols/simple_moves/SampleSequenceFromProbabilities.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/simple_metrics/PerResidueProbabilitiesMetric.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/TaskFactory.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

namespace protocols {
namespace simple_moves {

///@brief A class to sample sequences from a PerResidueProbabilitiesMetric and thread them onto the pose.
class SampleSequenceFromProbabilities : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SampleSequenceFromProbabilities();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SampleSequenceFromProbabilities() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///@brief Set the PerResidueProbabilitiesMetric that will be used to calculate the pseudo-perplexity.
	void
	set_metric( core::simple_metrics::PerResidueProbabilitiesMetricCOP metric );

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

	///@brief Set the positional temperature option
	void set_pos_temp(core::Real pos_temp);

	///@brief Set the amino acid temperature option
	void set_aa_temp(core::Real aa_temp);

	///@brief Set the probability cutoff option
	void set_prob_cutoff(core::Real prob_cutoff);

	///@brief Set the delta probability cutoff option
	void set_delta_prob_cutoff(core::Real delta_prob_cutoff);

	///@brief Set the maximum number of mutations allowed option
	void set_max_mutations(core::Size max_mutations);

	///@brief Set a bool to define whether we repack
	void set_packing(bool packing);

	///@brief Set the taskfactory option
	void set_task_factory(core::pack::task::TaskFactoryOP & task_factory);

public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	//SampleSequenceFromProbabilities & operator=( SampleSequenceFromProbabilities const & src );

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

	/// @brief This mover is unpublished. It returns Moritz Ertelt as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

private: // methods

	///@brief Main function to sample mutations from a probability distribution
	/*!@details Creates a discrete distribution based on the delta_prob (curr_AA - max_prob_AA) of each position
	and draws from that distribution taking into a account the pos_temp_ which controls the randomness. The idea behind
	this is that positions with the larges disagreement between current AA and predicted AAs will be mutated first.
	Also checks the probability_cutoff_, delta_probability_cutoff and present TaskOperations and sets probabilities not
	meeting the requirements to zero. Returns as many positions with their proposed AA as are set through
	the max_mutation_ option set by the user. The randomness of choice of AA is controlled by the aa_temp_ option.
	*/
	///@param[in] pose Uses the pose to get the current residues and check for TaksOps.
	///@param[in] values The values of a PerResidueProbabilitiesMetric.
	///@returns A map containing positions and their proposed AA drawn from the input probabilities.
	std::map<core::Size, core::chemical::AA> sample_mutations(
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> & values,
		core::pose::Pose const & pose
	) const;

	///@brief Helper function to get the ranked positions based on maximum difference in probabilities.
	std::vector<core::Size> sample_positions(
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> & values,
		core::pose::Pose const & pose
	) const;

	///@brief Helper function to sample an amino acid from a transformed probability distribution
	core::chemical::AA sample_amino_acid(
		std::map<core::chemical::AA, core::Real> const & aa_probs
	) const;

	///@brief Helper function to return a string of three amino acid letters from the proposed mutations for the SimpleThreadingMover
	static std::string construct_modified_sequence(core::pose::Pose& pose, std::map< core::Size, core::chemical::AA > & mutations) ;

	///@brief Helper function to check whether the AA is set to designable in the Packer
	static bool is_aa_allowed_by_task( core::pack::task::ResidueLevelTask const &rlt, core::chemical::AA aa) ;

	///@brief Helper function to calculate differences between current and other AAs, as well as disabling unwanted AAs
	utility::vector1<std::pair<core::Size, core::Real>> calculate_position_diffs(
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> & probabilities,
		core::pose::Pose const& pose
	) const;

	///@brief Helper function to return positions ranked based on their delta_probability value to the current AA
	std::vector<core::Size> get_ranked_positions(
		utility::vector1<std::pair<core::Size, core::Real>>& position_diffs,
		core::pose::Pose const& pose
	) const;

private: // data

	// metric to be used
	core::simple_metrics::PerResidueProbabilitiesMetricCOP metric_ = nullptr;

	// sampling options
	core::Real pos_temp_ = 1.0;
	core::Real aa_temp_ = 1.0;
	core::Real prob_cutoff_ = 0.0001;
	core::Real delta_prob_cutoff_ = 0.0;
	core::Size max_mutations_ = 10000;
	bool packing_ = true;

	//Accessing cached data from the pose
	bool use_cache_ = false;
	std::string cache_prefix_;
	std::string cache_suffix_;
	bool fail_on_missing_cache_ = true;

	// TaskOperation support
	core::pack::task::TaskFactoryOP task_factory_;

};

std::ostream &
operator<<( std::ostream & os, SampleSequenceFromProbabilities const & mover );

} //simple_moves
} //protocols

#endif //protocols_simple_moves_SampleSequenceFromProbabilities_HH
