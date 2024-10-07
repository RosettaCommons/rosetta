// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SampleSequenceFromProbabilities.cc
/// @brief A class to sample sequences from a PerResidueProbabilitiesMetric and thread them onto the pose.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com), University of Leipzig

// Unit headers
#include <protocols/simple_moves/SampleSequenceFromProbabilities.hh>
#include <protocols/simple_moves/SampleSequenceFromProbabilitiesCreator.hh>

// protocol headers
#include <protocols/simple_moves/SimpleThreadingMover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/simple_metrics/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

// Std headers
#include <random>
#include <unordered_set>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.simple_moves.SampleSequenceFromProbabilities" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SampleSequenceFromProbabilities::SampleSequenceFromProbabilities():
	protocols::moves::Mover( SampleSequenceFromProbabilities::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SampleSequenceFromProbabilities::~SampleSequenceFromProbabilities()= default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
SampleSequenceFromProbabilities::apply( core::pose::Pose& pose ){
	using namespace core::chemical;

	// check that it's a PerResidueProbabilitiesMetric
	if ( metric_ == nullptr ) {
		utility_exit_with_message("SampleSequenceFromProbabilities: This mover requires a PerResidueProbabilitiesMetric!");
	}
	// get values from the metric
	TR << "Calculating/Fetching probabilities from metric..." << std::endl;
	std::map< core::Size, std::map< AA, core::Real >> values = metric_->cached_calculate( pose, use_cache_, cache_prefix_, cache_suffix_, fail_on_missing_cache_ );
	TR << " Done!" << std::endl;

	TR << "Sampling mutations...." << std::endl;
	// sample mutations from probabilities
	std::map< core::Size, AA > mutations = sample_mutations( values, pose );
	// Construct the modified sequence string
	std::string modified_sequence = construct_modified_sequence(pose, mutations );
	TR << " Done!" << std::endl;

	// Check if the sequence is unmutated (only hyphens)
	bool is_unmutated = true;
	for ( std::size_t i = 0; i < modified_sequence.length(); i += 4 ) {
		if ( i >= modified_sequence.length() || modified_sequence[i] != '-' ) {
			is_unmutated = false;
			break;
		}
	}

	// Call SimpleThreadingMover with sampled sequence if mutated
	if ( !is_unmutated ) {
		TR << "Threading sequence onto pose..." << std::endl;
		SimpleThreadingMoverOP threader(utility::pointer::make_shared<SimpleThreadingMover>(modified_sequence, 1));
		threader->set_pack_neighbors(true);
		// disable packing if specified by user
		if ( !packing_ ) {
			threader->set_pack_rounds(0);
		}
		threader->set_sequence_mode("threeletter");
		threader->apply(pose);
		TR.Info << " Done!" << std::endl;
	} else {
		TR.Warning << "No mutations match the thresholds set. Skipping SimpleThreadingMover." << std::endl;
	}
}

void
SampleSequenceFromProbabilities::set_use_cached_data(bool use_cache, std::string const & prefix, std::string const & suffix){
	use_cache_ = use_cache;
	cache_prefix_ = prefix;
	cache_suffix_ = suffix;
	if ( use_cache_ ) {
		TR << "Attempting to use cached data with set prefix/suffix:" <<prefix <<" "<<suffix << std::endl;
	}
}

void
SampleSequenceFromProbabilities::set_fail_on_missing_cache(bool fail){
	fail_on_missing_cache_ = fail;
}

void
SampleSequenceFromProbabilities::set_metric(core::simple_metrics::PerResidueProbabilitiesMetricCOP metric){
	metric_ = metric;
}
void SampleSequenceFromProbabilities::set_pos_temp(core::Real pos_temp) {
	if ( pos_temp == 0.0 ) {
		pos_temp += 0.000001; // to avoid 0 division
	}
	pos_temp_ = pos_temp;
}
void SampleSequenceFromProbabilities::set_aa_temp(core::Real aa_temp) {
	if ( aa_temp == 0.0 ) {
		aa_temp += 0.000001; // to avoid 0 division
	}
	aa_temp_ = aa_temp;
}
void SampleSequenceFromProbabilities::set_prob_cutoff(core::Real prob_cutoff) {
	prob_cutoff_ = prob_cutoff;
}
void SampleSequenceFromProbabilities::set_delta_prob_cutoff(core::Real delta_prob_cutoff) {
	delta_prob_cutoff_ = delta_prob_cutoff;
}
void SampleSequenceFromProbabilities::set_max_mutations(core::Size max_mutations) {
	max_mutations_ = max_mutations;
}
void SampleSequenceFromProbabilities::set_packing(bool packing) {
	packing_ = packing;
}
void SampleSequenceFromProbabilities::set_task_factory(core::pack::task::TaskFactoryOP & task_factory) {
	task_factory_ = task_factory;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
SampleSequenceFromProbabilities::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
SampleSequenceFromProbabilities::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data
) {
	//SimpleMetric::parse_base_tag( tag );

	// PerResidueProbabilitiesMetric we will sample from
	core::simple_metrics::SimpleMetricCOP metric = core::simple_metrics::get_metric_from_datamap_and_subtags(tag, data, "metric");

	if ( metric->simple_metric_type() != "PerResidueProbabilitiesMetric" ) {
		utility_exit_with_message("The SampleSequenceFromProbabilities mover only works with PerResidueProbabilitiesMetrics!");
	}

	core::simple_metrics::PerResidueProbabilitiesMetricCOP res_metric = utility::pointer::dynamic_pointer_cast< core::simple_metrics::PerResidueProbabilitiesMetric const>( metric );

	set_metric( res_metric );

	// sampling options
	set_pos_temp( tag->getOption< core::Real >( "pos_temp", pos_temp_ ) );
	set_aa_temp( tag->getOption< core::Real >( "aa_temp", aa_temp_ ) );
	set_prob_cutoff( tag->getOption< core::Real >( "prob_cutoff", prob_cutoff_ ) );
	set_delta_prob_cutoff( tag->getOption< core::Real>( "delta_prob_cutoff", delta_prob_cutoff_ ) );
	set_max_mutations( tag->getOption< core::Size >( "max_mutations", max_mutations_ ) );
	set_packing( tag->getOption< bool >( "packing", true ) );

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

	/// Parse TaskOps
	core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	if ( new_task_factory != nullptr ) {
		set_task_factory( new_task_factory );
	}

}
void SampleSequenceFromProbabilities::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::required_attribute( "metric", xs_string, "A PerResidueProbabilitiesMetric to calculate the pseudo-perplexity from." );
	attlist + XMLSchemaAttribute::attribute_w_default(
		"pos_temp", xsct_real,
		"Sampling temperature for choosing positions. T greater 1 is more deterministic, T smaller 1 increases diversity. Defaults to 1", "1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"aa_temp", xsct_real,
		"Sampling temperature for choosing amino acids. T greater 1 is more deterministic, T smaller 1 increases diversity. Defaults to 1", "1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"prob_cutoff", xsct_real,
		"Probability cutoff for amino acids to consider for sampling. Defaults to 0.0001", "0.0001");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"delta_prob_cutoff", xsct_real,
		"Delta probability (current-proposed) cutoff for amino acids to consider for sampling. Defaults to 0 (at least as likely as the current)", "0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_mutations", xsct_non_negative_integer,
		"The maximal amount of mutations allowed. Defaults to basically unlimited (10000)", "10000");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"packing", xsct_rosetta_bool,
		"Whether to do any packing at all, if set to false results in basic rotamers. You might want to set this to false if you e.g. relax after this mover anyways. (Default true)", "true");

	attlist + XMLSchemaAttribute::attribute_w_default( "use_cached_data",  xsct_rosetta_bool, "Use any data stored in the datacache that matches the set metrics name (and any prefix/suffix.)  Data is stored during a SimpleMetric's apply function, which is called during RunSimpleMetrics", "false");
	attlist + XMLSchemaAttribute("cache_prefix", xs_string, "Any prefix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute("cache_suffix", xs_string, "Any suffix used during apply (RunSimpleMetrics), that we will match on if use_cache is true");
	attlist + XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If use_cached_data is True and cache is not found, should we fail?", "true");

	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A class to sample sequences from a PerResidueProbabilitiesMetric and thread them onto the pose.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SampleSequenceFromProbabilities::fresh_instance() const
{
	return utility::pointer::make_shared< SampleSequenceFromProbabilities >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SampleSequenceFromProbabilities::clone() const
{
	return utility::pointer::make_shared< SampleSequenceFromProbabilities >( *this );
}

std::string SampleSequenceFromProbabilities::get_name() const {
	return mover_name();
}

std::string SampleSequenceFromProbabilities::mover_name() {
	return "SampleSequenceFromProbabilities";
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
SampleSequenceFromProbabilitiesCreator::create_mover() const
{
	return utility::pointer::make_shared< SampleSequenceFromProbabilities >();
}

std::string
SampleSequenceFromProbabilitiesCreator::keyname() const
{
	return SampleSequenceFromProbabilities::mover_name();
}

void SampleSequenceFromProbabilitiesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SampleSequenceFromProbabilities::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Moritz Ertelt as its author.
void
SampleSequenceFromProbabilities::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"SampleSequenceFromProbabilities", basic::citation_manager::CitedModuleType::Mover,
		"Moritz Ertelt",
		"University of Leipzig",
		"moritz.ertelt@gmail.com",
		"Wrote the SampleSequenceFromProbabilities."
		)
	);
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

std::vector<core::Size>
SampleSequenceFromProbabilities::sample_positions(
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> & probabilities,
	core::pose::Pose const & pose
) const {
	using namespace core::chemical;
	utility::vector1<std::pair<core::Size, core::Real>> position_diffs = calculate_position_diffs( probabilities, pose );
	return get_ranked_positions( position_diffs, pose );
}

core::chemical::AA
SampleSequenceFromProbabilities::sample_amino_acid(
	std::map<core::chemical::AA, core::Real> const & aa_probs
) const {
	std::vector<core::Real> accumulated_weights;
	std::vector<core::chemical::AA> aas;

	core::Real total_weight = 0.0;
	for ( auto const & aa_prob : aa_probs ) {
		core::Real weight = std::pow(aa_prob.second, 1.0 / aa_temp_);
		total_weight += weight;
		accumulated_weights.push_back(total_weight);
		TR.Debug << aa_prob.first << " " << aa_prob.second << " weight " << weight << std::endl;
		aas.push_back(aa_prob.first);
	}
	// Generate a random number between 0 and the total accumulated weight.
	core::Real random_weight = total_weight * numeric::random::uniform();

	// Use binary search to find the appropriate index.
	// Note: I use this instead of std::discrete_distribution as it handles the cases where all but one weight equal 0.
	auto it = std::lower_bound(accumulated_weights.begin(), accumulated_weights.end(), random_weight);
	core::Size sampled_idx = std::distance(accumulated_weights.begin(), it);

	return aas[sampled_idx];
}

std::map<core::Size, core::chemical::AA>
SampleSequenceFromProbabilities::sample_mutations(
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> & values,
	core::pose::Pose const & pose
) const {
	std::map<core::Size, core::chemical::AA> sampled_mutations;
	for ( core::Size const &position : sample_positions(values, pose) ) {
		auto const & aa_probs = values.at(position);
		core::chemical::AA sampled_aa = sample_amino_acid(aa_probs);
		TR.Debug << "THE SAMPLED AA AT POSITION " << position << " IS " << sampled_aa << std::endl;
		sampled_mutations[position] = sampled_aa;
	}
	return sampled_mutations;
}

std::string
SampleSequenceFromProbabilities::construct_modified_sequence(core::pose::Pose& pose, std::map< core::Size, core::chemical::AA > & mutations) {
	using namespace core::chemical;
	std::vector<std::string> modified_sequence;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		AA aa = pose.residue(i).aa();
		if ( mutations.find(i) != mutations.end() ) {
			aa = mutations[i];
		}
		std::string aa_name = name_from_aa(aa);
		if ( aa == pose.residue(i).aa() ) {
			aa_name = "-"; // if current aa just leave it, formatted for SimpleThreadingMover
		}
		modified_sequence.push_back(aa_name);
	}

	std::string modified_sequence_str;
	for ( core::Size i = 0; i < modified_sequence.size(); ++i ) {
		modified_sequence_str += modified_sequence[i];
		// don't add a comma to the last position
		if ( i < modified_sequence.size() - 1 ) {
			modified_sequence_str += ",";
		}
	}
	return modified_sequence_str;
}

bool
SampleSequenceFromProbabilities::is_aa_allowed_by_task( core::pack::task::ResidueLevelTask const &rlt, core::chemical::AA aa) {
	if ( !rlt.being_packed() ) return false;

	if ( rlt.being_designed() ) {
		auto const &allowed_residues = rlt.allowed_residue_types();
		// returns true if the AA is in the list of allowed residues
		return std::any_of(allowed_residues.begin(), allowed_residues.end(),
			[&aa]( std::shared_ptr<const core::chemical::ResidueType> const &res_type_ptr) { return res_type_ptr->aa() == aa; });
	}
	return false;
}

utility::vector1<std::pair<core::Size, core::Real>>
SampleSequenceFromProbabilities::calculate_position_diffs(
	std::map<core::Size, std::map<core::chemical::AA, core::Real>> & probabilities,
	core::pose::Pose const& pose
) const {
	using namespace core::chemical;

	// setup task factory
	core::pack::task::PackerTaskCOP packer_task;
	if ( task_factory_ != nullptr ) {
		packer_task = task_factory_->create_task_and_apply_taskoperations(pose);
	}
	// Calculate the maximum difference for each position and disable unwanted AAs by setting them to zero
	utility::vector1<std::pair<core::Size, core::Real>> position_diffs; // Contains pairs of position and their max_diff
	for ( auto & pos_and_probs : probabilities ) {
		core::Size position = pos_and_probs.first;
		std::map<core::chemical::AA, core::Real> & aa_probs = pos_and_probs.second;

		// get the probability of the amino acid currently present in the pose at that position
		AA current_aa = pose.residue(position).aa();

		core::Real current_prob;
		if ( aa_probs.find(current_aa) != aa_probs.end() ) {
			current_prob = aa_probs.at( current_aa );
		} else {
			current_prob = 0.0;
			TR.Error << "The amino acid type at position " << position << " was not found in the probabilities map!" << std::endl;
		}

		core::Real max_prob = 0.0;
		for ( auto & aa_and_prob : aa_probs ) {
			core::chemical::AA aa = aa_and_prob.first;
			core::Real prob = aa_and_prob.second;

			// checking task operations and setting disallowed AAs to zero
			// Check against task operations
			if ( task_factory_ != nullptr && !is_aa_allowed_by_task( packer_task->residue_task( position ), aa ) ) {
				aa_and_prob.second = 0.0;
				continue;
			}
			// setting probs below user cutoff to zero
			if ( prob < prob_cutoff_ || prob - current_prob < delta_prob_cutoff_ ) {
				aa_and_prob.second = 0.0;
				continue;
			}
			// else check if its max
			if ( prob > max_prob ) {
				max_prob = prob;
			}
		}

		core::Real max_diff = max_prob - current_prob;

		// Check if the maximum difference (delta) for this position is above the threshold
		if ( max_diff >= delta_prob_cutoff_ && max_prob >= prob_cutoff_ ) {
			position_diffs.emplace_back(position, max_diff);
		}
	}

	return position_diffs;
}

std::vector<core::Size>
SampleSequenceFromProbabilities::get_ranked_positions(
	utility::vector1<std::pair<core::Size, core::Real>>& position_diffs,
	core::pose::Pose const& pose
) const {
	// Calculate weights using the position temperature to control randomness
	// std::vector<core::Real> weights;
	utility::vector1< core::Real > weights;
	for ( auto const & pos_and_diff : position_diffs ) {
		core::Real diff = pos_and_diff.second;
		weights.push_back(std::pow(diff, 1.0 / pos_temp_));
	}
	std::vector<core::Size> ranked_positions;
	std::unordered_set<core::Size> selected_positions;
	while ( ranked_positions.size() < max_mutations_ && ranked_positions.size() < pose.total_residue() && !weights.empty() ) {
		// I create this each time since we change the weights to avoid re-sampling the same position
		// std::discrete_distribution breaks the mac build so I use the WeightedSampler instead
		// core::Size sampled_idx = std::discrete_distribution<>(weights.begin(), weights.end())(numeric::random::rg());
		numeric::random::WeightedSampler sampler;
		sampler.weights( weights );
		core::Size sampled_idx = sampler.random_sample( numeric::random::rg());

		// Ensure we don't re-sample the same position
		if ( selected_positions.find(position_diffs[sampled_idx].first) == selected_positions.end() ) {
			ranked_positions.push_back(position_diffs[sampled_idx].first);
			selected_positions.insert(position_diffs[sampled_idx].first);
		}
		// Remove the sampled position from the list
		std::swap(position_diffs[sampled_idx], position_diffs.back());
		position_diffs.pop_back();

		std::swap(weights[sampled_idx], weights.back());
		weights.pop_back();
	}
	return ranked_positions;
}

std::ostream &
operator<<( std::ostream & os, SampleSequenceFromProbabilities const & mover )
{
	mover.show(os);
	return os;
}


} //simple_moves
} //protocols
