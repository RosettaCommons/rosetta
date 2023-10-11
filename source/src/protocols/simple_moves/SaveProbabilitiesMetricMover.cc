// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SaveProbabilitiesMetricMover.cc
/// @brief A class to save a PerResidueProbabilitiesMetric to a weight file.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

// Unit headers
#include <protocols/simple_moves/SaveProbabilitiesMetricMover.hh>
#include <protocols/simple_moves/SaveProbabilitiesMetricMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/simple_metrics/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/io/ozstream.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.simple_moves.SaveProbabilitiesMetricMover" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SaveProbabilitiesMetricMover::SaveProbabilitiesMetricMover():
	protocols::moves::Mover( SaveProbabilitiesMetricMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SaveProbabilitiesMetricMover::~SaveProbabilitiesMetricMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
SaveProbabilitiesMetricMover::apply( core::pose::Pose& pose ){
	using namespace core::chemical;

	// check that it's a PerResidueProbabilitiesMetric
	if ( metric_ == nullptr ) {
		utility_exit_with_message("SaveProbabilitiesMetricMover: This mover requires a PerResidueProbabilitiesMetric!");
	}
	// get values from the metric
	TR << "Calculating/Fetching probabilities from metric..." << std::endl;
	std::map< core::Size, std::map< AA, core::Real >> values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );
	TR << "Done!" << std::endl;
	TR << "Saving probabilities to file " << filename_ << " formatted as " << filetype_ << std::endl;

	// if filetype is both create a name with .pssm and one with .weights, otherwise leave filename alone
	std::string weights_filename = filename_;
	std::string pssm_filename = filename_;

	if ( filetype_ == "both" ) {
		core::Size lastdot = filename_.find_last_of(".");
		if ( lastdot != std::string::npos ) {
			std::string base_filename = filename_.substr(0, lastdot);  // Remove extension if any
			weights_filename = base_filename + ".weights";
			pssm_filename = base_filename + ".pssm";
		} else {
			weights_filename += ".weights";
			pssm_filename += ".pssm";
		}
	}

	if ( filetype_ == "weights" || filetype_ == "both" ) {
		save_aa_probabilities_to_file(weights_filename, values);
	}

	if ( filetype_ == "pssm" || filetype_ == "both" ) {
		std::string const selection_sequence = get_selection_sequence(pose.sequence(), values);
		std::map<core::Size, std::map<core::chemical::AA, core::Real>> logit_map = convert_probabilities_to_logits(values);
		core::simple_metrics::PerResidueProbabilitiesMetric::output_sequence_profile(selection_sequence, logit_map, pssm_filename);
	}
	TR << "Done! Saved as: " << filename_ << std::endl;
}
void
SaveProbabilitiesMetricMover::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}
void
SaveProbabilitiesMetricMover::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}
void
SaveProbabilitiesMetricMover::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}
void
SaveProbabilitiesMetricMover::set_filename( std::string const & filename ){
	filename_ = filename;
}
void
SaveProbabilitiesMetricMover::set_filetype( std::string const & filetype ){
	if ( filetype != "weights" && filetype != "pssm" && filetype != "both" ) {
		utility_exit_with_message( "Allowed filetypes are pssm, weights or both, unrecognized type " + filetype );
	}
	filetype_ = filetype;
}

void
SaveProbabilitiesMetricMover::save_aa_probabilities_to_file( std::string const & weights_file, std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & prob_set) {

	utility::io::ozstream out_stream(weights_file);
	// Write a header line to the file for clarity
	out_stream << "#POSNUM RESIDUETYPE WEIGHT\n";

	// Iterate through the probabilities and write each one to the file
	for ( auto const & residue_entry : prob_set ) {
		core::Size resi = residue_entry.first;
		auto const & aa_probs = residue_entry.second;

		for ( auto const & aa_entry : aa_probs ) {
			auto aa = aa_entry.first;
			core::Real weight = aa_entry.second;

			core::Real rounded_weight = std::round(weight * 1000000.0) / 1000000.0; // rounding to 6 decimal places

			if ( rounded_weight < 0.0 || rounded_weight > 1.0 ) {
				utility_exit_with_message("Probability value out of range [0,1]: " + std::to_string(rounded_weight));
			}

			// Convert the AA enum to a string
			std::string resn = core::chemical::name_from_aa( aa );

			// Write the data to the file
			out_stream << resi << " " << resn << " " << rounded_weight << "\n";
		}
	}
	out_stream.close();
}

///@brief Get the sequence of the selection present in the probability map
std::string
SaveProbabilitiesMetricMover::get_selection_sequence(
	std::string const & pose_sequence,
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & position_map
) {
	std::string selection_sequence;
	// maps in C++ have sorted keys
	for ( auto const & pair : position_map ) {
		// pair.first is the position (of type core::Size)
		selection_sequence += pose_sequence[pair.first - 1];
	}
	return selection_sequence;
}

///@brief convert the map of probabilities to one of logits for writing out to a PSSM
std::map<core::Size, std::map<core::chemical::AA, core::Real>>
SaveProbabilitiesMetricMover::convert_probabilities_to_logits(
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> const & probabilities_map)
{
	const core::Real epsilon = 1e-9;  // Small constant to avoid log(0) and log(1)

	std::map<core::Size, std::map<core::chemical::AA, core::Real>> logits_map;

	for ( const auto& pos_and_probs : probabilities_map ) {
		core::Size position = pos_and_probs.first;

		for ( const auto& aa_and_prob : pos_and_probs.second ) {
			core::chemical::AA aa = aa_and_prob.first;
			core::Real prob = aa_and_prob.second;

			// Add epsilon to avoid log(0) and log(1)
			prob = std::max(epsilon, std::min(1.0 - epsilon, prob));

			core::Real logit = std::log(prob / (1.0 - prob));

			logits_map[position][aa] = logit;
		}
	}
	return logits_map;
}


////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
SaveProbabilitiesMetricMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data
) {

	// PerResidueProbabilitiesMetric we will sample from
	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("The SaveProbabilitiesMetricMover only works with PerResidueProbabilitiesMetrics!");
	}

	core::simple_metrics::PerResidueProbabilitiesMetricCOP res_metric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::PerResidueProbabilitiesMetric const>( metric );

	set_metric( res_metric );

	set_filename( tag->getOption< std::string >( "filename") );

	if ( tag->hasOption( "filetype" ) ) {
		set_filetype( tag->getOption< std::string >( "filetype" ) );
	}

	// options for using cached data
	bool use_cache = tag->getOption< bool >("use_cached_data", false);
	std::string prefix;
	std::string suffix;

	if ( tag->hasOption("cache_prefix") ) {
		prefix = tag->getOption< std::string >("cache_prefix");
	}
	if ( tag->hasOption("cache_suffix") ) {
		suffix = tag->getOption< std::string >("cache_suffix");
	}
	set_use_cached_data(use_cache, prefix, suffix);
	set_fail_on_missing_cache(tag->getOption< bool>("fail_on_missing_cache", fail_on_missing_cache_));

}

void SaveProbabilitiesMetricMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to calculate the pseudo-perplexity from." );

	attlist + XMLSchemaAttribute::required_attribute( "filename", xs_string, "The name of the output weights file storing the probabilities." );

	attlist + XMLSchemaAttribute::attribute_w_default( "filetype",  xs_string, "The output filetype, either psi-blast pssm (pssm), weights (weights) file or both (both) (default is weights)", "weights");

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A class to save the probabilities from a PerResidueProbabilitiesMetric in a weights file.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SaveProbabilitiesMetricMover::fresh_instance() const
{
	return utility::pointer::make_shared< SaveProbabilitiesMetricMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SaveProbabilitiesMetricMover::clone() const
{
	return utility::pointer::make_shared< SaveProbabilitiesMetricMover >( *this );
}

std::string SaveProbabilitiesMetricMover::get_name() const {
	return mover_name();
}

std::string SaveProbabilitiesMetricMover::mover_name() {
	return "SaveProbabilitiesMetricMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
SaveProbabilitiesMetricMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< SaveProbabilitiesMetricMover >();
}

std::string
SaveProbabilitiesMetricMoverCreator::keyname() const
{
	return SaveProbabilitiesMetricMover::mover_name();
}

void SaveProbabilitiesMetricMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SaveProbabilitiesMetricMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns moritzertelt as its author.
void
SaveProbabilitiesMetricMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"SaveProbabilitiesMetricMover", basic::citation_manager::CitedModuleType::Mover,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		"Wrote the SaveProbabilitiesMetricMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, SaveProbabilitiesMetricMover const & mover )
{
	mover.show(os);
	return os;
}


} //simple_moves
} //protocols
