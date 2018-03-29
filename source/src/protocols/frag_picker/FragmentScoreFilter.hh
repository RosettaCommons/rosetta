// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author Cody Krivacic (krivacic@berkeley.edu)


#ifndef INCLUDED_protocols_frag_picker_FragmentScoreFilter_hh
#define INCLUDED_protocols_frag_picker_FragmentScoreFilter_hh

// Unit headers
#include <protocols/frag_picker/FragmentScoreFilter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/Tag.fwd.hh> //transcluded from Filter.hh
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Filter.hh
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentCrmsd.hh>

namespace protocols {
namespace frag_picker {

void run_command( std::string const & command );

void convert_binary_checkpoint( std::string const & check_filename );

///@brief --brief--
class FragmentScoreFilter : public protocols::filters::Filter {

public:
	FragmentScoreFilter();

	// destructor (important for properly forward-declaring smart-pointer members)
	~FragmentScoreFilter() override;

	/// @brief returns true if the structure passes the filter, false otherwise
	bool
	apply( core::pose::Pose const & pose ) const override;

	void rescore_fragment_crmsd( core::pose::PoseOP const pose, utility::vector1<frag_picker::Candidate> candidates,  frag_picker::scores::FragmentCrmsd* fc )const;

	void setup_fragment_picker( core::pose::PoseOP const pose, frag_picker::FragmentPickerOP picker ) const;

	core::Real
	compute( core::pose::Pose const pose ) const;

	core::Real
	get_result( utility::vector1<core::Real> ) const;

	/// @brief required for reporting score values
	core::Real
	report_sm( core::pose::Pose const & pose ) const override;

	/// @brief allows printing data to a stream
	void
	report( std::ostream & os, core::pose::Pose const & pose ) const override;

public:
	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief parse XML tag (to use this Filter in Rosetta Scripts)
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

private:
	/// @brief Threshold number for passing or failing filter
	core::Real threshold_;

	/// @brief Choose whether numbers that are higher or numbers that are lower than the filter pass
	std::string direction_;

	/// @brief Choose which score type to filter based on
	std::string score_type_;

	/// @brief Start evaluating fragments at this residue
	core::Size start_res_;
	/// @brief Stop evaluating fragments after this residue
	core::Size end_res_;

	/// @brief Size of fragments to evaluate
	core::Size fragment_size_ = 9;

	/// @brief Should the final score be the minimum, maximum, or average of the positions queried?
	std::string compute_ = "maximum";

	/// @brief How to choose the best fragment at each position (default FragmentCrmsd; TotalScore can also be used)
	std::string sort_by_ = "FragmentCrmsd";

	/// @brief Folder to place sequence profile files
	std::string outputs_folder_;

	/// @brief Basename of sequence profile files
	std::string outputs_name_ = "pose";

	/// @brief Location of CSBLAST program
	std::string csblast_;

	/// @brief Location of BLAST PGP program
	std::string blast_pgp_;

	/// @brief Location of sequence database for BLAST PGP
	std::string placeholder_seqs_;

	/// @brief Location of SPARKS-X directory
	std::string sparks_x_;

	/// @brief Location of SPARKS-X query script
	std::string sparks_x_query_;

	/// @brief Path to run_psipred_single script
	std::string psipred_;

	/// @brief Path to vall database
	std::string vall_path_;

	/// @brief Path to scoring config file
	std::string frags_scoring_config_;

	/// @brief How many fragments per position?
	core::Size n_frags_;

	/// @brief How many candidates per position?
	core::Size n_candidates_;

};
}
}


#endif //INCLUDED_--path_underscore--_--class--_hh
