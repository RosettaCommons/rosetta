// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/ResidueTorsionRestraints.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/ResidueTorsionRestraints.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

//Auto Headers
namespace protocols {
namespace ligand_docking {


static basic::Tracer TR( "protocols.ligand_dock.ResidueTorsionRestraints", basic::t_debug );


ResidueTorsionRestraints::ResidueTorsionRestraints(
	core::pose::Pose & pose,
	core::Size resid,
	core::Real stddev_degrees
):
	resid_(resid),
	stddev_degrees_(stddev_degrees),
	my_constraints_(),
	old_chi_(),
	old_constraints_()
{
	// my_constraints_ is empty to start with, so all existing constraints will be kept initially.
	enable( pose );
}

bool ResidueTorsionRestraints::operator==(const ResidueTorsionRestraints &other){
	return resid_==other.resid_;
}

/// @details Adds constraints to all rotatable torsions except proton chis.
/// Conserves all existing constraints except ones previously added by this object.
void ResidueTorsionRestraints::setup_constraints(core::pose::Pose & pose)
{
	using namespace core::scoring::constraints;
	using core::chemical::ResidueType;
	using core::chemical::AtomIndices;
	using core::conformation::Residue;
	using core::id::AtomID;

	// input stddev is in degrees, but dihedral constraints deal in radians
	core::Real const stddev_radians = numeric::conversions::radians( stddev_degrees_ );

	utility::vector1< core::scoring::constraints::ConstraintCOP > dont_care;
	ConstraintSetOP new_constraints = without_my_constraints( pose.constraint_set(), dont_care );
	// Over the lifetime of this object, its original Pose may be cloned, copied over,
	// reverted to some previous version (think MC trials), etc.
	// Thus, *any* of the constraints we've *ever* created may be in the Pose passed to enable/disable().
	// So we need to keep record of all of them.  Fortunately they're stored as strings, which are small.
	//my_constraints_.clear();

	Residue const & rsd = pose.residue(resid_);
	ResidueType const & rsd_type = pose.residue_type(resid_);
	for ( core::Size j = 1, j_end = rsd_type.nchi(); j <= j_end; ++j ) {
		core::Real const curr_chi_degrees = rsd.chi(j);
		core::Real const curr_chi_radians = numeric::conversions::radians( curr_chi_degrees );
		core::scoring::func::FuncOP restr_func( new core::scoring::func::CircularHarmonicFunc( curr_chi_radians, stddev_radians ) );
		AtomIndices chi_idx = rsd_type.chi_atoms(j); // 1-based
		DihedralConstraintCOP constraint( DihedralConstraintOP( new DihedralConstraint(
			AtomID(chi_idx[1], resid_),
			AtomID(chi_idx[2], resid_),
			AtomID(chi_idx[3], resid_),
			AtomID(chi_idx[4], resid_),
			restr_func
			) ) );
		TR << "Constraint: " << curr_chi_degrees << " deg, " << constraint->atom(1) << " "
			<< constraint->atom(2) << " " << constraint->atom(3) << " " << constraint->atom(4) << std::endl;
		// Is this still necessary (no) or advisable (maybe)?  Constraint will be removed before packing...
		if ( rsd.atom_type(chi_idx[1]).is_hydrogen() || rsd.atom_type(chi_idx[4]).is_hydrogen() ) {
			TR << "Constraint involves hydrogen atom; skipping it for PROTON_CHI." << std::endl;
		} else {
			new_constraints->add_constraint( constraint );
			my_constraints_.insert( constraint->to_string() );
		}
	}

	pose.constraint_set( new_constraints );
}


core::scoring::constraints::ConstraintSetOP
ResidueTorsionRestraints::without_my_constraints(
	core::scoring::constraints::ConstraintSetCOP old_constraints,
	utility::vector1< core::scoring::constraints::ConstraintCOP > & removed_constraints
)
{
	using namespace core::scoring::constraints;
	ConstraintSetOP new_constraints( new ConstraintSet() );

	// Cycle through all existing constraints, and keep all except the ones we added:
	if ( old_constraints.get() != nullptr ) {
		utility::vector1< ConstraintCOP > old_constr = old_constraints->get_all_constraints();
		for ( Size i = 1; i <= old_constr.size(); ++i ) {
			if ( my_constraints_.find( old_constr[i]->to_string() ) == my_constraints_.end() ) {
				new_constraints->add_constraint( old_constr[i] );
				TR << "Keeping old constraint " << old_constr[i]->to_string() << std::endl;
			} else {
				removed_constraints.push_back( old_constr[i] );
				TR << "Removing old constraint " << old_constr[i]->to_string() << std::endl;
			}
		}
	}
	return new_constraints;
}


void ResidueTorsionRestraints::enable( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;
	TR.Trace << "enable(), # constraints before: " << pose.constraint_set()->get_all_constraints().size() << std::endl;

	// Has our target residue changed conformation significantly?
	utility::vector1< core::Real > new_chi = pose.residue(resid_).chi(); // need a copy, not a reference
	bool resid_has_changed = (old_chi_.size() != new_chi.size());
	if ( !resid_has_changed ) {
		for ( core::Size i = 1; i <= old_chi_.size(); ++i ) {
			if ( std::abs( old_chi_[i] - new_chi[i] ) > 1e-1 ) resid_has_changed = true;
		}
	}

	TR.Trace << "Residue " << resid_ << " has changed conformation? " << resid_has_changed << std::endl;
	if ( !resid_has_changed && !old_constraints_.empty() ) {
		// No:  restore constraints exactly as they were
		ConstraintSetOP new_constraints = pose.constraint_set()->clone(); // deep copy, constraints and their funcs are cloned too
		for ( core::scoring::constraints::ConstraintCOP constraint : old_constraints_ ) {
			new_constraints->add_constraint( constraint );
		}
		pose.constraint_set( new_constraints );
	} else {
		// Yes:  generate new constraints for our DOFs
		setup_constraints( pose );
	}
	old_constraints_.clear();
	TR.Trace << "enable(), # constraints after: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
}


void ResidueTorsionRestraints::disable( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;
	// Memorize initial conformation of our target residue
	old_chi_ = pose.residue(resid_).chi(); // need a copy, not a reference
	// Remove constraints
	TR.Trace << "disable(), # constraints before: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
	old_constraints_.clear(); // Guard against double-disable
	ConstraintSetOP new_constraints = without_my_constraints( pose.constraint_set(), old_constraints_ );
	pose.constraint_set( new_constraints );
	TR.Trace << "disable(), # constraints after: " << pose.constraint_set()->get_all_constraints().size() << std::endl;
}


} // namespace ligand_docking
} // namespace protocols
