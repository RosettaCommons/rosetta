// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/ResidueProbDesignOperation.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/task_operations/ResidueProbDesignOperation.hh>
#include <protocols/task_operations/ResidueProbDesignOperationCreator.hh>

#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>
#include <utility/io/util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


static basic::Tracer TR( "protocols.TaskOperations.ResidueProbDesignOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

using namespace core::chemical;
using core::pack::task::PackerTask;
using core::pack::task::operation::TaskOperationOP;
using core::Size;
using core::Real;
using utility::vector1;
typedef std::map< core::chemical::AA, Real > AAProbabilities; //Map of an amino acid and it's probability.
typedef std::map< Size, AAProbabilities > PerResidueAAProbSet; //Amino acid probabilities for a particular residue number.

core::pack::task::operation::TaskOperationOP
ResidueProbDesignOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ResidueProbDesignOperation );
}

void ResidueProbDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueProbDesignOperation::provide_xml_schema( xsd );
}

std::string ResidueProbDesignOperationCreator::keyname() const
{
	return ResidueProbDesignOperation::keyname();
}


ResidueProbDesignOperation::ResidueProbDesignOperation() : core::pack::task::operation::TaskOperation()
{
	set_defaults();
}

void
ResidueProbDesignOperation::set_defaults() {
	set_include_native_restype(true);
	set_keep_task_allowed_aas(false);
	set_sample_zero_probs_at(0.0);
	set_picking_rounds(1);
}

/*ResidueProbDesignOperation::ResidueProbDesignOperation(PerResidueAAProbSet aa_probs) {
ResidueProbDesignOperation();
prob_set_ = aa_probs;
}*/

ResidueProbDesignOperation::~ResidueProbDesignOperation()= default;

ResidueProbDesignOperation::ResidueProbDesignOperation(ResidueProbDesignOperation const & rhs): core::pack::task::operation::TaskOperation(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

core::pack::task::operation::TaskOperationOP
ResidueProbDesignOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ResidueProbDesignOperation( *this ) );
}

/*
ResidueProbDesignOperation & ResidueProbDesignOperation::operator =(ResidueProbDesignOperation const & rhs){
if ( this == &rhs ) return *this;
core::pack::task::operation::TaskOperation::operator=( rhs );
init_for_equal_operator_and_copy_constructor( *this, rhs );
return *this;
}
*/

void
ResidueProbDesignOperation::init_for_equal_operator_and_copy_constructor(ResidueProbDesignOperation& lhs, ResidueProbDesignOperation const & rhs){
	lhs.include_native_restype_ = rhs.include_native_restype_;
	lhs.keep_task_allowed_aa_ = rhs.keep_task_allowed_aa_;
	lhs.overall_prob_set_ = rhs.overall_prob_set_;
	lhs.picking_rounds_ = rhs.picking_rounds_;
	lhs.prob_set_ = rhs.prob_set_;
	lhs.zero_probs_overwrite_ = rhs.zero_probs_overwrite_;
	lhs.no_probability_ = rhs.no_probability_;
}

////////////////////////////////////////////////////////////////////////////
// Probability Settings
//

void
ResidueProbDesignOperation::set_aa_probabilities(const Size resnum, AAProbabilities aa_probs) {
	prob_set_[resnum] = aa_probs;
}

void
ResidueProbDesignOperation::set_aa_probability_set(PerResidueAAProbSet per_res_aa_probs) {
	prob_set_ = per_res_aa_probs;
}

void
ResidueProbDesignOperation::set_overall_aa_probabilities(AAProbabilities aa_probs) {
	overall_prob_set_ = aa_probs;
}

void
ResidueProbDesignOperation::set_aa_probabilities_from_file( const std::string& weights_file) {
	auto lines = utility::io::get_lines_from_file_data(weights_file);

	PerResidueAAProbSet prob_set;

	// FIXME: assign all weights
	for ( std::string const& line : lines ) {
		auto columns = utility::split(line);
		if ( line[0] == '#' ) continue;
		if ( columns.size() < 3 ) {
			utility_exit_with_message("Weights has to be specified in the following format: POSNUM RESIDUETYPE WEIGHT");
		}

		core::Size resi = 0;
		resi = utility::from_string(columns[1], resi);
		auto resn = columns[2];
		core::Real weight = 0.0;
		weight = utility::from_string(columns[3], weight);
		TR << resi << " " << resn << " " << weight << std::endl;
		auto aa = core::chemical::aa_from_one_or_three(resn);

		if ( prob_set.find(resi) == prob_set.end() ) {
			prob_set.insert(std::make_pair(resi, AAProbabilities()));
		}
		prob_set[resi][aa] = weight;
	}

	set_aa_probability_set(prob_set);
}

////////////////////////////////////////////////////////////////////////////
// Reset Probabilities
//

void
ResidueProbDesignOperation::clear_prob_set() {
	prob_set_.clear();
}

void
ResidueProbDesignOperation::clear_overall_prob_set() {
	overall_prob_set_.clear();
}


////////////////////////////////////////////////////////////////////////////
// General Options
//

void
ResidueProbDesignOperation::set_include_native_restype(bool include) {
	include_native_restype_ = include;
}

void
ResidueProbDesignOperation::set_keep_task_allowed_aas(bool setting) {
	keep_task_allowed_aa_ = setting;
}

void
ResidueProbDesignOperation::set_picking_rounds(core::Size picking_rounds) {
	picking_rounds_ = picking_rounds;
}

void
ResidueProbDesignOperation::set_sample_zero_probs_at(const core::Real probability) {
	zero_probs_overwrite_ = probability;
}

bool
ResidueProbDesignOperation::include_native_restype() const {
	return include_native_restype_;
}
bool
ResidueProbDesignOperation::keep_task_allowed_aas() const {
	return keep_task_allowed_aa_;
}

void
ResidueProbDesignOperation::set_no_probability(bool no_probability){
	no_probability_ = no_probability;
}

Size
ResidueProbDesignOperation::picking_rounds() const {
	return picking_rounds_;
}
Real
ResidueProbDesignOperation::sample_zero_probs_at() const {
	return zero_probs_overwrite_;
}


bool
ResidueProbDesignOperation::resnum_exists_in_set(core::Size const resnum) const {
	auto iter(prob_set_.find(resnum));

	return iter != prob_set_.end();
}

bool
ResidueProbDesignOperation::aa_exists(AA amino_acid, std::map< AA, Real > const & aa_prob_set) const {
	auto iter(aa_prob_set.find(amino_acid));

	return iter != aa_prob_set.end();
}

vector1< Real >
ResidueProbDesignOperation::get_weights(AAProbabilities & aa_prob_set) const{

	vector1<Real> aa_weights(20, zero_probs_overwrite_);

	//Non-cannonicals will have to come later
	for ( core::Size i = 1; i <= 20; ++i ) {
		auto amino = static_cast< core::chemical::AA >(i);
		if ( ! aa_exists(amino, aa_prob_set) || aa_prob_set[amino] == 0.0 ) { continue; }
		aa_weights[i] = aa_prob_set[amino];
	}

	return aa_weights;

}

void
ResidueProbDesignOperation::apply(core::pose::Pose const & pose, core::pack::task::PackerTask& task) const {

	using namespace core::chemical;

	vector1<bool > design_positions = task.designing_residues();

	//Copy class variables to keep function const.
	PerResidueAAProbSet prob_set(prob_set_);
	AAProbabilities overall_prob_set(overall_prob_set_);


	numeric::random::WeightedSampler sampler;
	vector1<Real> aa_weights;
	//Setup the distribution once if we have a single distribution to choose aa from.
	if ( ! overall_prob_set_.empty() ) {
		TR << "Overall AA Probability set found.  Using." << std::endl;
		aa_weights = get_weights(overall_prob_set);
		sampler.weights(aa_weights);
	} else {
		// Fill potentially missing weights
		for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
			AAProbabilities probs;
			if ( prob_set.find(resi) != prob_set.end() ) {
				probs = prob_set[resi];
			}
			for ( core::Size aai = core::chemical::first_l_aa;
					aai <= core::chemical::num_canonical_aas; ++aai ) {
				auto aa = static_cast<core::chemical::AA>(aai);
				if ( probs.find(aa) == probs.end() ) {
					probs[aa] = 1.0;
				}
			}
			prob_set[resi] = probs;
		}
	}

	for ( core::Size i = 1; i <=pose.size(); ++i ) {

		//If Residue is not designable or ResidueProb not set + there is no overall probability set to use, keep going.
		if ( ( ! design_positions[i] ) || ( ! resnum_exists_in_set(i) && overall_prob_set.empty() ) ) {
			continue;
		}

		if ( overall_prob_set.empty() ) {
			aa_weights = get_weights( prob_set[i] );
			sampler.weights(aa_weights);
		}
		vector1<bool> allowed_aminos(20, false);


		//Get the Amino Acid from the sampler.  Add to task depending on options.
		if ( include_native_restype_ ) {
			allowed_aminos[task.nonconst_residue_task(i).get_original_residue()] = true;
		}

		///Add all residues with non-zero probability for the position
		if ( no_probability_ ) {
			for ( core::Size x = 1; x <= allowed_aminos.size(); ++x ) {
				if ( aa_weights[x] > 0 ) {
					allowed_aminos[x] = true;
				}
			}

			task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aminos);
		} else {
			for ( Size picks = 1; picks <= picking_rounds_; ++picks ) {

				core::Size aa_num = sampler.random_sample(numeric::random::rg());
				allowed_aminos[aa_num] = true;
			}

			//Add to already allowed aminos
			if ( keep_task_allowed_aa_ ) {
				for ( core::Size aa_num = 1; aa_num <= 20; ++aa_num ) {
					auto amino = static_cast<core::chemical::AA>(aa_num);
					if ( allowed_aminos[aa_num] ) {
						task.nonconst_residue_task(i).allow_aa(amino);
					}
				}
			} else {
				//Replace the current aminos
				task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aminos);
			}
		}
	}
}

void
ResidueProbDesignOperation::parse_tag( TagCOP tag , DataMap & ) {
	set_aa_probabilities_from_file(tag->getOption<std::string>("weights_file"));
	set_keep_task_allowed_aas(tag->getOption<bool>("keep_task_allowed_aas"));
	set_include_native_restype(tag->getOption<bool>("include_native_restype"));
	set_sample_zero_probs_at(tag->getOption<core::Real>("sample_zero_probs"));
	set_picking_rounds(tag->getOption<core::Size>("picking_rounds"));
	set_no_probability(tag->getOption<bool>("no_probabilities", no_probability_));

}

void
ResidueProbDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::required_attribute(
		"weights_file", xs_string,
		"Path to a file containting weights for each residue type at each position in the following format:\nPOSNUM RESIDUETYPE WEIGHT\nWith weight as a real number between 0.0 and 1.0. Unspecified residues and/or residue types will automatically default to a weight of 1.0.")
		+ XMLSchemaAttribute::attribute_w_default(
		"include_native_restype", xsct_rosetta_bool,
		"The native residue type is always an allowed residue type.", "true")
		+ XMLSchemaAttribute::attribute_w_default(
		"keep_task_allowed_aas", xsct_rosetta_bool,
		"If set to true, the sampled residue types will not replace all other previously allowed residue types.", "false")
		+ XMLSchemaAttribute::attribute_w_default(
		"sample_zero_probs", xsct_real,
		"Overwrite the sampling probability for all residue types with a weight of zero with the given weight. For example, if you have a probability of zero, you can sample this instead at like .3 or .5 or whatever you want. The idea is that you don't need to go and change your input data to add some variability in the data - useful if you have a very low sampling rate of the input data.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default( "no_probabilities", xsct_rosetta_bool, "Should we sample ALL AA that does not have prob of 0 at 1.0 instead?  This basically works to add ALL AA seen at a given position in a particular cluster to the set of design residues.  Used to increase variability of designs or for testing purposes", "false")
		+ XMLSchemaAttribute::attribute_w_default(
		"picking_rounds", xsct_positive_integer,
		"Allowed residue types can be sampled multiple times. Especially of interest in combination with the option 'keep_task_allowed_aas'.",
		"1");

	task_op_schema_w_attributes(xsd, keyname(), attributes, "Randomly adds allowed residue types to the PackerTask. Based on the user specified weights [0.0-1.0] for each residue type per position this operation will add the residue type to the list of allowed design residues in a non-deterministic manner.");
}

} //task_operation
} //protocols
