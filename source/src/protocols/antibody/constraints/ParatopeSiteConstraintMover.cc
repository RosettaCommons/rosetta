// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/ParatopeSiteConstraintMover.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/AntibodyInfo.hh>

#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/SiteConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>

static thread_local basic::Tracer TR("antibody.constraints.ParatopeSiteConstraintMover");

namespace protocols {
namespace antibody {
namespace constraints {
	using utility::vector1;
ParatopeSiteConstraintMover::ParatopeSiteConstraintMover() :
	protocols::moves::Mover(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

ParatopeSiteConstraintMover::ParatopeSiteConstraintMover(AntibodyInfoCOP ab_info) :
	protocols::moves::Mover(),
	current_func_(/* NULL */)
{
	ab_info_ = ab_info;
	set_defaults();
}

ParatopeSiteConstraintMover::~ParatopeSiteConstraintMover(){}

void
ParatopeSiteConstraintMover::set_defaults() {
	cdrs_to_apply_.clear();
	cdrs_to_apply_.resize(6, true);
	
	antigen_chains_.clear();
	paratope_residues_.clear();
	interface_distance_ = 10;
	//cst_set_ = new core::scoring::constraints::ConstraintSet();
}

void
ParatopeSiteConstraintMover::constrain_to_paratope_cdrs(const vector1<CDRNameEnum>& paratope_cdrs){
	cdrs_to_apply_.clear();
	cdrs_to_apply_.resize(6, false);
	for (core::Size i = 1; i <= paratope_cdrs.size(); ++i){
		cdrs_to_apply_[core::Size(paratope_cdrs[i])] = true;
	}
	
}

void
ParatopeSiteConstraintMover::constrain_to_paratope_cdrs(const vector1<bool>& paratope_cdrs) {
	cdrs_to_apply_ = paratope_cdrs;
}

void
ParatopeSiteConstraintMover::constrain_to_antigen_chains(const vector1<core::Size>& antigen_chains){
	antigen_chains_ = antigen_chains;
}

void
ParatopeSiteConstraintMover::constrain_to_paratope_residues(const vector1<bool>& paratope_residues) {
	paratope_residues_ = paratope_residues;
}

void
ParatopeSiteConstraintMover::set_constraint_func(core::scoring::func::FuncOP constraint_func){
	current_func_ = constraint_func;
}

void
ParatopeSiteConstraintMover::set_interface_distance(core::Real interface_distance){
	interface_distance_ = interface_distance;
}

void
ParatopeSiteConstraintMover::remove(core::pose::Pose& pose, bool reset_paratope_residues){
	using namespace core::scoring::constraints;
	
	if (! reset_paratope_residues){
		assert(paratope_residues_.size() == pose.total_residue());
	}
	else{
		this->setup_paratope_residues_from_cdrs(pose);
	}
	
	if (antigen_chains_.size() == 0) antigen_chains_ = ab_info_->get_antigen_chain_ids(pose);
	
	vector1<ConstraintOP> csts_to_be_removed;
	for (core::Size i = 1; i <= pose.total_residue(); ++i){
		if (! paratope_residues_[i]) continue;
		
		for (core::Size x = 1; x <= antigen_chains_.size(); ++x){
			SiteConstraintOP res_constraint = setup_constraints(pose, i, utility::to_string(core::pose::get_chain_from_chain_id(x, pose)));
			csts_to_be_removed.push_back(res_constraint);
		}
	}
	pose.remove_constraints(csts_to_be_removed, true);
}

void
ParatopeSiteConstraintMover::setup_paratope_residues_from_cdrs(core::pose::Pose const & pose){
	
	paratope_residues_.clear();
	paratope_residues_.resize(pose.total_residue(), false);
	for (core::Size i =1; i <= core::Size(CDRNameEnum_total); ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		
		if (cdrs_to_apply_[cdr]){
			for (core::Size x = ab_info_->get_CDR_start(cdr, pose); x <= ab_info_->get_CDR_end(cdr, pose); ++x){
				paratope_residues_[x] = true;
			}
		}
	}
}

void
ParatopeSiteConstraintMover::apply(core::pose::Pose& pose){

	using namespace core::scoring::constraints;
	
	if (! ab_info_){
		ab_info_ = AntibodyInfoCOP( AntibodyInfoOP( new AntibodyInfo(pose) ) );
	}
	//Check if antigen is present
	if (! ab_info_->antigen_present())
	{
		TR <<"Antigen not present!  Could not apply constraints" << std::endl;
		set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
		return;
	}
	//If the antibody is camelid, remove trying to setup constraints to light chain
	if (ab_info_->is_camelid())
	{
		cdrs_to_apply_[l1] = false;
		cdrs_to_apply_[l2] = false;
		cdrs_to_apply_[l3] = false;
	}
	
	//Check any settings, set defaults from our antibody info.
	if (antigen_chains_.size() == 0) antigen_chains_ = ab_info_->get_antigen_chain_ids(pose);
	if (paratope_residues_.size() == 0) setup_paratope_residues_from_cdrs(pose);
	
	//Check any set function
	if (! current_func_){
		current_func_ = core::scoring::func::FuncOP( new core::scoring::func::FlatHarmonicFunc(0, 1, interface_distance_) );
	}
	
	assert(paratope_residues_.size() == pose.total_residue());
	
	//Ready to go!
	ConstraintCOPs current_csts = pose.constraint_set()->get_all_constraints();

	//pose.constraint_set()->show(TR);
	for (core::Size i = 1; i <= pose.total_residue(); ++i){
		if (! paratope_residues_[i]) continue;
		
		for (core::Size x = 1; x <= antigen_chains_.size(); ++x){
			SiteConstraintOP res_constraint = setup_constraints(pose, i, utility::to_string(core::pose::get_chain_from_chain_id(antigen_chains_[x], pose)));
			
			//Use find  - this may take some time.  Also, the func comparisons seem to be by value - I really don't think they would ever be the same since we clone.
			//Looks like clone of atom pair constraint DOES NOT clone the func, so we should be ok.
			if (std::find(current_csts.begin(), current_csts.end(), res_constraint) == current_csts.end()){
				pose.add_constraint(res_constraint);
				//TR<< "Added constraint: " << i << " to chain "<<utility::to_string(core::pose::get_chain_from_chain_id(antigen_chains_[x], pose)) << std::endl;
			}
		}
	}
}

core::scoring::constraints::SiteConstraintOP
ParatopeSiteConstraintMover::setup_constraints(core::pose::Pose const & pose, core::Size resnum, const std::string chain){
	using namespace core::scoring::constraints;
	
	//core::scoring::constraints::ConstraintSetOP atom_constraints = new ConstraintSet();
	//AmbiguousConstraintOP res_constraint = new AmbiguousConstraint();
	

	SiteConstraintOP atom_constraint( new SiteConstraint() );
	atom_constraint->setup_csts(
				resnum, 
				"CA", 
				chain, 
				pose, 
				current_func_);
		
	return atom_constraint;
	
}


} //constraints
} //antibody
} //protocols
