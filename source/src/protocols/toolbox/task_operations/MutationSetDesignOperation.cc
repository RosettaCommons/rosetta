// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/MutationSetDesignOperation.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/toolbox/task_operations/MutationSetDesignOperation.hh>

// Numeric Includes
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

static basic::Tracer TR("protocols.toolbox.task_operations.MutationSetDesignOperation");

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;
using utility::vector1;

/* AMW: no Creator
void MutationSetDesignOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
MutationSetDesignOperation::provide_xml_schema( xsd );
}

std::string MutationSetDesignOperationCreator::keyname() const
{
return MutationSetDesignOperation::keyname();
}
*/

MutationSetDesignOperation::~MutationSetDesignOperation()= default;

MutationSetDesignOperation::MutationSetDesignOperation():
	TaskOperation()
{
	set_defaults();
}

MutationSetDesignOperation::MutationSetDesignOperation(utility::vector1<MutationSet> mutation_sets):
	TaskOperation()
{
	set_defaults();
	set_mutation_sets(mutation_sets);
}

MutationSetDesignOperation::MutationSetDesignOperation(
	utility::vector1<MutationSet> mutation_sets,
	utility::vector1<core::Real> weights ) :
	TaskOperation()
{
	set_defaults();
	set_mutation_sets(mutation_sets, weights);
}


MutationSetDesignOperation::MutationSetDesignOperation(MutationSetDesignOperation const & /*src*/) = default;

TaskOperationOP
MutationSetDesignOperation::clone() const{
	return TaskOperationOP(new MutationSetDesignOperation( *this));
}

void
MutationSetDesignOperation::set_defaults(){
	add_to_allowed_aas_ = false;
	include_native_aa_ = false;
	picking_rounds_ = 1;
	reset_sample_index();

	//sample_number_ = 0;
}

void
MutationSetDesignOperation::set_mutation_sets(utility::vector1<MutationSet> mutation_sets) {
	mutation_sets_ = mutation_sets;
	weights_.clear();
	weights_.resize(mutation_sets_.size(), 1.0); //Set all weights to one.
}

void
MutationSetDesignOperation::set_mutation_sets(
	utility::vector1<MutationSet> mutation_sets,
	utility::vector1<core::Real> mutation_set_weights)
{
	mutation_sets_ = mutation_sets;
	weights_ = mutation_set_weights;
}

void
MutationSetDesignOperation::add_mutation_set(MutationSet mutation_set) {
	mutation_sets_.push_back(mutation_set);
	weights_.push_back(1.0);
}

void
MutationSetDesignOperation::add_mutation_set(MutationSet mutation_set, core::Real weight){
	mutation_sets_.push_back(mutation_set);
	weights_.push_back(weight);
}

void
MutationSetDesignOperation::clear_mutation_sets(){
	mutation_sets_.clear();
	weights_.clear();
}

void
MutationSetDesignOperation::set_picking_rounds(core::Size picking_rounds){
	picking_rounds_ = picking_rounds;
}

void
MutationSetDesignOperation::add_to_allowed_aas(const bool& setting){
	add_to_allowed_aas_ = setting;
}

void
MutationSetDesignOperation::include_native_aa(const bool& setting){
	include_native_aa_ = setting;
}

core::Size
MutationSetDesignOperation::get_sample_index() const{
	return sample_index_;
}

core::Size
MutationSetDesignOperation::get_total_mutation_sets() const {
	return mutation_sets_.size();
}

void
MutationSetDesignOperation::reset_sample_index(){
	sample_index_ = 0;
}

void
MutationSetDesignOperation::set_sample_index(core::Size sample_index) {
	sample_index_ = sample_index;
}

void
MutationSetDesignOperation::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{
	if ( mutation_sets_.size() != weights_.size() ) {
		utility_exit_with_message("Mutation sets and weights must be equivalent in length!");
	}

	//Setup what we have now for each residue
	vector1< vector1< bool > > pose_allowed_aminos;
	vector1< bool > sampled_resnums(pose.size(), false);

	for ( core::Size i = 1;  i <= pose.size(); ++i ) {
		vector1<bool> allowed_aminos(20, false);
		if ( include_native_aa_ ) {
			allowed_aminos[ task.residue_task(i).get_original_residue() ] = true;
		}
		pose_allowed_aminos.push_back(allowed_aminos);
	}

	core::Size picking_rounds = picking_rounds_;
	if ( sample_index_ != 0 ) { picking_rounds = 1;}

	for ( core::Size round = 1; round <= picking_rounds; ++round ) {
		MutationSet current_set;
		if ( sample_index_ != 0 ) {
			if ( sample_index_ > mutation_sets_.size() ) {
				utility_exit_with_message("Sample index > length of mutation_sets!");
			}

			current_set = mutation_sets_[sample_index_];
		} else {

			numeric::random::WeightedSampler sampler;
			sampler.weights(weights_);
			current_set = mutation_sets_[sampler.random_sample(numeric::random::rg())];

		}

		//Iterate through the map, setting positions.
		//std::string mutant;
		utility::vector1<core::Size> positions;

		for ( std::map< core::Size, core::chemical::AA>::const_iterator ele = current_set.begin(), ele_end = current_set.end(); ele != ele_end; ++ele ) {
			core::Size resnum = ele->first;
			core::chemical::AA allowed_mutation = ele->second;
			//mutant = mutant + oneletter_code_from_aa(allowed_mutation);
			pose_allowed_aminos[ resnum ][ core::Size(allowed_mutation)] = true;
			sampled_resnums[ resnum ] = true;
			positions.push_back(resnum);

		}
		//TR << "positions: "<< utility::to_string(positions) << std::endl;
		//TR << "mutant: " << mutant << std::endl;


	} //For picking rounds

	//Add to already allowed aminos
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! sampled_resnums[ i ] ) { continue;}

		if ( add_to_allowed_aas_ ) {
			for ( core::Size aa_num = 1; aa_num <= 20; ++aa_num ) {
				auto amino = static_cast<core::chemical::AA>(aa_num);
				if ( pose_allowed_aminos[ i ][ aa_num ] ) {
					task.nonconst_residue_task(i).allow_aa(amino);
				}
			}
		} else {
			//Replace the current aminos
			task.nonconst_residue_task(i).restrict_absent_canonical_aas(pose_allowed_aminos[ i ]);
		}
	}
}

// No parse_tag
void MutationSetDesignOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}


} //design
} //antibody
} //protocols
