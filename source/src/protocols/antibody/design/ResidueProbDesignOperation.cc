// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/ResidueProbDesignOperation.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/ResidueProbDesignOperation.hh>
#include <protocols/antibody/design/ResidueProbDesignOperationCreator.hh>

#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>



static basic::Tracer TR( "protocols.toolbox.TaskOperations.ResidueProbDesignOperation" );
static numeric::random::RandomGenerator RG(136548);

namespace protocols {
namespace antibody {
namespace design {
	
	using namespace core::chemical;
	using core::pack::task::PackerTask;
	using core::pack::task::operation::TaskOperationOP;
	using core::Size;
	using core::Real;

	typedef std::map< core::chemical::AA, Real > AAProbabilities; //Map of an amino acid and it's probability.
	typedef std::map< Size, AAProbabilities > PerResidueAAProbSet; //Amino acid probabilities for a particular residue number.
	
ResidueProbDesignOperation::ResidueProbDesignOperation() : core::pack::task::operation::TaskOperation()	
{
	init();
}

void
ResidueProbDesignOperation::init() {
	set_include_native_restype(true);
	set_keep_task_allowed_aas(false);
	set_sample_zero_probs_at(0.0);
	set_picking_rounds(1);
}

/*ResidueProbDesignOperation::ResidueProbDesignOperation(PerResidueAAProbSet aa_probs) {
	ResidueProbDesignOperation();
	prob_set_ = aa_probs;
}*/

ResidueProbDesignOperation::~ResidueProbDesignOperation(){}

ResidueProbDesignOperation::ResidueProbDesignOperation(ResidueProbDesignOperation const & rhs): core::pack::task::operation::TaskOperation(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

core::pack::task::operation::TaskOperationOP
ResidueProbDesignOperation::clone() const
{
	return new ResidueProbDesignOperation( *this );
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
}

////////////////////////////////////////////////////////////////////////////
// Probability Settings
//
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


////////////////////////////////////////////////////////////////////////////
// Reset Probabilities
//
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
	std::map< Size, std::map< AA, Real > >::const_iterator iter(prob_set_.find(resnum));
	
	return iter != prob_set_.end();
}

bool
ResidueProbDesignOperation::aa_exists(AA amino_acid, std::map< AA, Real > const & aa_prob_set) const {
	std::map< core::chemical::AA, Real >::const_iterator iter(aa_prob_set.find(amino_acid));
	
	return iter != aa_prob_set.end();
}

vector1< Real >
ResidueProbDesignOperation::get_weights(AAProbabilities & aa_prob_set) const{
	
	vector1<Real> aa_weights(20, zero_probs_overwrite_);
	
	//Non-cannonicals will have to come later
	for (core::Size i = 1; i <= 20; ++i){
		core::chemical::AA amino = static_cast< core::chemical::AA >(i);
		if (! aa_exists(amino, aa_prob_set) || aa_prob_set[amino] == 0.0 ){ continue; }
		aa_weights[i] = aa_prob_set[amino];
	}
	
	return aa_weights;
	
}

void
ResidueProbDesignOperation::apply(core::pose::Pose const & pose, core::pack::task::PackerTask& task) const {
	
	using namespace core::chemical;
	
	vector1<bool > design_positions = task.designing_residues();
	
	//Copy class variables to keep function const.  
	std::map< Size, std::map< AA, Real > > prob_set(prob_set_);
	std::map< AA, Real > overall_prob_set(overall_prob_set_);
	

	numeric::random::WeightedSampler sampler;
	
	//Setup the distribution once if we have a single distribution to choose aa from.
	if (! overall_prob_set_.empty() ){
		TR << "Overall AA Probability set found.  Using." << std::endl;
		vector1<Real> aa_weights = get_weights(overall_prob_set);
		sampler.weights(aa_weights);
	}
	
	for (core::Size i = 1; i <=pose.total_residue(); ++i){
		
		//If Residue is not designable or ResidueProb not set + there is no overall probability set to use, keep going.
		if (( ! design_positions[i] ) || ( ! resnum_exists_in_set(i) && overall_prob_set.empty() )) {
			continue;
		}
			
		if (overall_prob_set.empty()){
			vector1< Real > aa_weights = get_weights( prob_set[i] );
			sampler.weights(aa_weights);
		}
		vector1<bool> allowed_aminos(20, false);
		
		
		//Get the Amino Acid from the sampler.  Add to task depending on options.
		if (include_native_restype_){
			allowed_aminos[task.nonconst_residue_task(i).get_original_residue()] = true;	
		}
		
		for (Size picks = 1; picks <= picking_rounds_; ++picks){
			
			core::Size aa_num = sampler.random_sample(RG);
			allowed_aminos[aa_num] = true;
		}
		
		if (keep_task_allowed_aa_){
			for (core::Size aa_num = 1; aa_num <= 20; ++aa_num){
				core::chemical::AA amino = static_cast<core::chemical::AA>(aa_num);
				if (allowed_aminos[aa_num]){
					task.nonconst_residue_task(i).allow_aa(amino);
				}
			}
		}
		else{
			task.nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aminos);
		}
	}
}


} //task_operations
} //toolbox
} //protocols
