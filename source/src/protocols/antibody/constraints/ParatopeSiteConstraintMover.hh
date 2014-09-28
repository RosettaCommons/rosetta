// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/constraints/ParatopeSiteConstraintMover.hh
/// @brief Adds and removes ambiguos site constraints for the Antibody Paratope to antigen, defined for simplicity as the CDRs.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_constraints_PARATOPESITECONSTRAINTMOVER_HH
#define INCLUDED_protocols_antibody_constraints_PARATOPESITECONSTRAINTMOVER_HH

#include <core/pose/Pose.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/constraints/AmbiguousConstraint.fwd.hh>
#include <core/scoring/constraints/SiteConstraint.fwd.hh>

#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace antibody {
namespace constraints {
	
	
///@brief Adds and removes ambiguous site constraints for the Antibody Paratope to antigen,
/// defined for simplicity as the CDRs, however a set of paratope residues can be given to the class to use instead.
///
class ParatopeSiteConstraintMover : public protocols::moves::Mover {

public:
	
	ParatopeSiteConstraintMover();
	ParatopeSiteConstraintMover(AntibodyInfoCOP ab_info);
	
	~ParatopeSiteConstraintMover();
	
	///@brief Add AmbiguousSiteConstraints to specified paratope residues, or each CDR set.  Default is all of them.
	/// If there are already exits AmbiguousSiteConstraints from the antigen to the residue, add them. 
	/// If they are missing, such as after a graft, add them.
	void
	apply(core::pose::Pose & pose);
	
	///@brief Remove constraints from each paratope residue and antigen chain set.  If reset_paratope_residues is true, then it will update the
	/// set of paratope residues it has.  This is used after pose length changes such as CDR insertion or deletion.  
	void
	remove(core::pose::Pose & pose, bool reset_paratope_residues = false);
	
	
public:
	
	///@brief Optionally use only these CDRs as the paratope.  Useful for constraining to Light or Heavy chain CDRs
	void
	constrain_to_paratope_cdrs(utility::vector1<CDRNameEnum> const & paratope_cdrs);
	
	void
	constrain_to_paratope_cdrs(utility::vector1<bool> const & paratope_cdrs);
	
	///@brief Optionally constrain to a set of pre-determined paratope residues
	void
	constrain_to_paratope_residues(utility::vector1<bool> const & paratope_residues);
	
	///@brief Optionally constrain to a set of antigen chains instead of all of them
	void
	constrain_to_antigen_chains(utility::vector1<core::Size> const & antigen_chains);
	
	
	///@brief Optionally set the Func that will be used for the constraint.  Default is the Flat_Harmonic at 0, 1, 5
	void
	set_constraint_func(core::scoring::func::FuncOP constraint_func);
	
	///@brief Set the interface distance that will be only be used as the tolerance for the linear harmonic constraint 
	/// if no constraint func is set. 
	void
	set_interface_distance(core::Real interface_distance);
	
public:
	void
	set_defaults();
	
	std::string 
	get_name() const {
		return "ParatopeSiteConstraintMover";
	}
	
	//parse_my_tag
	

	
private:
	
	void
	setup_paratope_residues_from_cdrs(core::pose::Pose const & pose);
	

	
	core::scoring::constraints::SiteConstraintOP
	setup_constraints(core::pose::Pose const & pose, core::Size resnum, std::string chain);
	
private:
	
	AntibodyInfoCOP ab_info_;
	utility::vector1<bool> cdrs_to_apply_;
	utility::vector1<bool> paratope_residues_;
	utility::vector1<core::Size> antigen_chains_;
	core::Real interface_distance_;
	
	//std::map< core::Size, vector1<core::scoring::constraints::AmbiguousConstraintCOP > > constraint_map_;
	core::scoring::func::FuncOP current_func_;
	
	//core::scoring::constraints::ConstraintSetOP cst_set_; //Cannot attach to conformation as pose only gives out references - not OPs
};	



}
}
}


#endif	//#ifndef INCLUDED_protocols/antibody_design_PARATOPESITECONSTRAINTMOVER_HH

