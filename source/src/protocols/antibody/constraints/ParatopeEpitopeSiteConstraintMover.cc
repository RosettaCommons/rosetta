// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/ParatopeEpitopeSiteConstraintMover.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>

#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/SiteConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR("antibody.constraints.ParatopeEpitopeSiteConstraintMover");

namespace protocols {
namespace antibody {
namespace constraints {
	using utility::vector1;

ParatopeEpitopeSiteConstraintMover::ParatopeEpitopeSiteConstraintMover() : 
	protocols::moves::Mover(),
	ab_info_(/* NULL */),
	current_func_(/* NULL */)
{
	set_defaults();
}

ParatopeEpitopeSiteConstraintMover::ParatopeEpitopeSiteConstraintMover(AntibodyInfoCOP ab_info) :
	protocols::moves::Mover(),
	current_func_(/* NULL */)
{
	ab_info_ = ab_info;
	set_defaults();
}

ParatopeEpitopeSiteConstraintMover::ParatopeEpitopeSiteConstraintMover(AntibodyInfoCOP ab_info, vector1<CDRNameEnum> paratope_cdrs):
	protocols::moves::Mover(),
	paratope_cdrs_(paratope_cdrs),
	current_func_(/* NULL */)
{
	ab_info_ = ab_info;
	constrain_to_paratope_cdrs(paratope_cdrs);
	interface_distance_ = 10.0;
}

ParatopeEpitopeSiteConstraintMover::ParatopeEpitopeSiteConstraintMover(AntibodyInfoCOP ab_info, vector1<CDRNameEnum> paratope_cdrs, vector1<bool> epitope_residues):
	protocols::moves::Mover(),
	epitope_residues_(epitope_residues),
	current_func_(/* NULL */)
{
	ab_info_ = ab_info;
	constrain_to_paratope_cdrs(paratope_cdrs);
	interface_distance_ = 10.0;
}
		
ParatopeEpitopeSiteConstraintMover::~ParatopeEpitopeSiteConstraintMover(){}

void
ParatopeEpitopeSiteConstraintMover::set_defaults(){
	
	paratope_cdrs_.clear();
	paratope_cdrs_.resize(6, true);
	paratope_residues_.clear();
	epitope_residues_.clear();
	current_func_ = NULL;
	interface_distance_ = 10.0;
}

void
ParatopeEpitopeSiteConstraintMover::set_constraint_func(core::scoring::func::FuncOP constraint_func){
	current_func_ = constraint_func;
}

void
ParatopeEpitopeSiteConstraintMover::set_interface_distance(const core::Size distance){
	interface_distance_ = distance;
}

void
ParatopeEpitopeSiteConstraintMover::constrain_to_epitope_residues(vector1<bool> const &epitope_residues) {
	epitope_residues_ = epitope_residues;
}

void
ParatopeEpitopeSiteConstraintMover::constrain_to_epitope_residues(const vector1<design::PDBNumbering >& epitope_residues, const core::pose::Pose& pose) {
	
	epitope_residues_ = protocols::antibody::design::get_resnum_from_pdb_numbering(pose, epitope_residues);
}

void
ParatopeEpitopeSiteConstraintMover::constrain_to_paratope_cdrs(vector1<CDRNameEnum> const & paratope_cdrs) {
	paratope_cdrs_.clear();
	paratope_cdrs_.resize(6, false);
	for (core::Size i = 1; i <= paratope_cdrs.size(); ++i){
		paratope_cdrs_[core::Size(paratope_cdrs[i])] = true;
	}
}

void
ParatopeEpitopeSiteConstraintMover::constrain_to_paratope_cdrs(const vector1<bool>& paratope_cdrs) {
	
	if (paratope_cdrs.size() != 6){
		utility_exit_with_message("Passed paratope cdrs does not equal the total number of cdrs!");
	}
	paratope_cdrs_ = paratope_cdrs;
	

}

void
ParatopeEpitopeSiteConstraintMover::constrain_to_paratope_residues(vector1<bool> const & paratope_residues) {
	paratope_residues_ = paratope_residues;
}

vector1<bool>
ParatopeEpitopeSiteConstraintMover::paratope_residues_from_cdrs(core::pose::Pose const & pose, const vector1<bool>& paratope_cdrs) const {
	
	vector1<bool> residues(pose.total_residue(), false);
	for (core::Size i = 1; i <= paratope_cdrs.size(); ++i){
		if (paratope_cdrs[i]){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			core::Size cdr_start = ab_info_->get_CDR_start(cdr, pose);
			core::Size cdr_end = ab_info_->get_CDR_end(cdr, pose);
			for (core::Size x = cdr_start; x <= cdr_end; ++x){
				residues[x] = true;
			}
		}
	}
	return residues;	
}

void
ParatopeEpitopeSiteConstraintMover::apply(core::pose::Pose& pose) {
	using namespace core::scoring::constraints;
	
	if (! ab_info_){
		ab_info_ = AntibodyInfoCOP( AntibodyInfoOP( new AntibodyInfo(pose) ) );
	}
	
	//Check if antigen is present
	if (! ab_info_->antigen_present()){
		TR <<"Antigen not present!  Could not apply constraints" << std::endl;
		set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
		return;
	}
	
	//If we have a camelid, only add constraints to H:
	if ( ab_info_->is_camelid()){
		paratope_cdrs_[l1] = false;
		paratope_cdrs_[l2] = false;
		paratope_cdrs_[l3] = false;
	}
	
	if (paratope_residues_.size() == 0){
		paratope_residues_ = this->paratope_residues_from_cdrs(pose, paratope_cdrs_);
	}
	
	
	//If no constraint is set.  Use the default.
	if (!current_func_){
		current_func_ = core::scoring::func::FuncOP( new core::scoring::func::FlatHarmonicFunc(0, 1, interface_distance_) );
	}
	
	//Setup antigen paratope residues if none are set.
	if (epitope_residues_.size() == 0){
		epitope_residues_ = protocols::antibody::select_epitope_residues(ab_info_, pose, interface_distance_);
	}
	
	assert(paratope_residues_.size() == pose.total_residue());
	assert(epitope_residues_.size() == pose.total_residue());
	
	TR << "Currently added constraints: "<<std::endl;
	//pose.constraint_set()->show(TR);
	
	//Setup constraint from paratope to epitope and from epitope to paratope.
	ConstraintCOPs current_csts = pose.constraint_set()->get_all_constraints();
	for (core::Size i = 1; i <= pose.total_residue(); ++i){
		
		if ((paratope_residues_[ i ] == true) && (epitope_residues_[ i ] == true)){
			utility_exit_with_message("Cannot be both paratope and epitope residue ");
		}
		
		if (paratope_residues_[i]){
			core::scoring::constraints::SiteConstraintOP constraint = setup_constraints(pose, i, epitope_residues_);
			if (std::find(current_csts.begin(), current_csts.end(), constraint) == current_csts.end()){
				pose.add_constraint(constraint);
				//TR << "Adding paratope-> epitope constraint: "<<i <<std::endl;
			}
		}
		else if  (epitope_residues_[i]){
			core::scoring::constraints::SiteConstraintOP constraint = setup_constraints(pose, i, paratope_residues_);
			if (std::find(current_csts.begin(), current_csts.end(), constraint) == current_csts.end()){
				pose.add_constraint(constraint);
				//TR << "Adding epitope-> paratope constraint: "<<i <<std::endl;
			}
			
		}

	}
	
} // apply

void
ParatopeEpitopeSiteConstraintMover::remove(core::pose::Pose & pose){
	
	using namespace core::scoring::constraints;
	
	vector1<ConstraintOP> csts_to_be_removed;
	for (core::Size i = 1; i <= pose.total_residue(); ++i){
		
		assert(paratope_residues_[i] != true && epitope_residues_[i] != true);
		
		if (paratope_residues_[i]){
			core::scoring::constraints::SiteConstraintOP constraint = setup_constraints(pose, i, epitope_residues_);
			//constraint_map_[i].push_back(L_constraint);
			csts_to_be_removed.push_back(constraint);
		}
		else if (epitope_residues_[i]){
			core::scoring::constraints::SiteConstraintOP constraint = setup_constraints(pose, i, paratope_residues_);
			//constraint_map_[i].push_back(H_constraint);
			csts_to_be_removed.push_back(constraint);
		}
	}
	pose.remove_constraints(csts_to_be_removed, true);
}

core::scoring::constraints::SiteConstraintOP
ParatopeEpitopeSiteConstraintMover::setup_constraints(core::pose::Pose const & pose, core::Size residue, vector1<bool> const & residues) const {
	using namespace core::scoring::constraints;
	
	//core::scoring::constraints::ConstraintSetOP atom_constraints = new ConstraintSet();
	
	SiteConstraintOP atom_constraint( new SiteConstraint() );
	atom_constraint->setup_csts(
				residue, 
				"CA", 
				residues, 
				pose, 
				current_func_);
		

	return atom_constraint;
}

} //constraints
} //antibody
} //protocols
