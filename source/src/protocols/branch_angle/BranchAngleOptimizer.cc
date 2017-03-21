// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/branch_angle/BranchAgnleOptimizer.cc
/// @brief implementation of BranchAgnleOptimizer methods
/// @author Colin A. Smith (colin.smith@ucsf.edu)

#include <protocols/branch_angle/BranchAngleOptimizer.hh>

// Protocols Headers
#include <protocols/branch_angle/BranchCoef1.hh>
#include <protocols/branch_angle/BranchCoef2.hh>
#include <protocols/branch_angle/BranchParam1.hh>
#include <protocols/branch_angle/BranchParam2.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID.hh>
#include <basic/database/open.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/conversions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/xyz.functions.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>

using namespace core;
using namespace core::scoring::mm;

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.branch_angle.BranchAngleOptimizer" );

namespace protocols {
namespace branch_angle {

BranchAngleOptimizer::BranchAngleOptimizer(
	core::scoring::mm::MMBondAngleLibrary const & mm_bondangle_library
):
	mm_bondangle_library_( mm_bondangle_library ),
	bond_angle_residue_type_param_set_(/* NULL */),
	tolerance_(0.0000001),
	initialized_(false)
{}

BranchAngleOptimizer::BranchAngleOptimizer():
	mm_bondangle_library_( scoring::ScoringManager::get_instance()->get_MMBondAngleLibrary() ),
	bond_angle_residue_type_param_set_(/* NULL */),
	tolerance_(0.0000001),
	initialized_(false)
{}

BranchAngleOptimizer::BranchAngleOptimizer( BranchAngleOptimizer const & ) = default;

BranchAngleOptimizer::~BranchAngleOptimizer() = default;

core::scoring::mm::MMBondAngleResidueTypeParamSetOP
BranchAngleOptimizer::bond_angle_residue_type_param_set()
{
	return bond_angle_residue_type_param_set_;
}

core::scoring::mm::MMBondAngleResidueTypeParamSetCOP
BranchAngleOptimizer::bond_angle_residue_type_param_set() const
{
	return bond_angle_residue_type_param_set_;
}

void
BranchAngleOptimizer::bond_angle_residue_type_param_set(
	core::scoring::mm::MMBondAngleResidueTypeParamSetOP param_set
)
{
	bond_angle_residue_type_param_set_ = param_set;
}

void
BranchAngleOptimizer::bond_angle_residue_type_param_set(
	core::scoring::mm::MMBondAngleResidueTypeParamSetCOP param_set
)
{
	bond_angle_residue_type_param_set_ = core::scoring::mm::MMBondAngleResidueTypeParamSetOP( new core::scoring::mm::MMBondAngleResidueTypeParamSet( *param_set ) );
}

/// @details
/// Add a branch point to the object, if optimization parameters were found, returns an index
/// to them, otherwise returns 0. If other failure conditions are met, a warning is output
/// and the function returns 0.
Size
BranchAngleOptimizer::optimize_angles(
	pose::Pose & pose,
	id::AtomID main_atomid1,
	id::AtomID center_atomid,
	id::AtomID main_atomid2,
	bool optimize_for_minimum // = false
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());

	//TR << "Optimizing Backbone: " << main_atomid1 << center_atomid << main_atomid2 << std::endl;

	kinematics::tree::AtomCOP main_atom1( pose.atom_tree().atom(main_atomid1).get_self_ptr() );
	kinematics::tree::AtomCOP const center_atom( pose.atom_tree().atom(center_atomid).get_self_ptr() );
	kinematics::tree::AtomCOP main_atom2( pose.atom_tree().atom(main_atomid2).get_self_ptr() );

	Size coef_index(0);

	// make sure that the center atom has a complete set of connections
	runtime_assert(!pose.residue(center_atomid.rsd()).has_incomplete_connection(center_atomid.atomno()));

	// count the number of atoms bonded to the central atom
	Size num_neighbors(pose.residue(center_atomid.rsd()).n_bonded_neighbor_all_res(center_atomid.atomno()));

	Real m2_bond_angle(numeric::angle_radians(pose.xyz(main_atomid1), pose.xyz(center_atomid),
		pose.xyz(main_atomid2)));

	if ( num_neighbors == 3 ) {

		// lookup the branching atom
		id::AtomID branch_atomid1;
		branching_atomid1(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1);

		//TR << "Optimizing 1 Branching Atom: " << branch_atomid1 << std::endl;

		kinematics::tree::AtomCOP const branch_atom1( pose.atom_tree().atom(branch_atomid1).get_self_ptr() );

		// switch main atoms if necessary for having a working stub
		if ( branch_atom1->input_stub_atom2() == main_atom2 ) {
			id::AtomID temp_atomid(main_atomid1);
			kinematics::tree::AtomCOP temp_atom(main_atom1);
			main_atomid1 = main_atomid2;
			main_atom1 = main_atom2;
			main_atomid2 = temp_atomid;
			main_atom2 = temp_atom;
		}

		// make sure both branching atom builds off main_atom1 & center_atom
		if ( !(branch_atom1->input_stub_atom1() == center_atom && branch_atom1->input_stub_atom2() == main_atom1 &&
				branch_atom1->input_stub_atom3() == main_atom2) ) {

			TR.Warning << "Branching atom (" << branch_atomid1 << ") for segment" << main_atomid1 << center_atomid
				<< main_atomid2 << "does not use segment atoms for input stub, ignoring" << std::endl;
			return 0;
		}

		// lookup the branching angle coefficients, returning 0 if they can't be found
		BranchParam1 const params(param1(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1));

		auto coef_map_iter(coef_map1_.find(params));
		if ( coef_map_iter == coef_map1_.end() ) {
			undefined_coef1_.insert(params);
			return 0;
		} else {
			coef_index = coef_map_iter->second;
		}

		if ( optimize_for_minimum ) m2_bond_angle = coef1_[coef_index].overall_theta0();

		// calculate the new torsion offsets and bond angles
		Real b1_torsion_offset;
		Real b1_bond_angle;
		coef1_[coef_index].evaluate(m2_bond_angle, b1_torsion_offset, b1_bond_angle);

		// set the corresponding degrees of freedom
		id::DOF_ID const branch_atom1_PHI(branch_atomid1, id::PHI);
		id::DOF_ID const branch_atom1_THETA(branch_atomid1, id::THETA);
		pose.set_dof(branch_atom1_PHI, b1_torsion_offset);
		pose.set_dof(branch_atom1_THETA, pi - b1_bond_angle);

	} else if ( num_neighbors == 4 ) {

		// lookup the branching atoms
		id::AtomID branch_atomid1;
		id::AtomID branch_atomid2;
		branching_atomids2(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1, branch_atomid2);

		//TR << "Optimizing 2 Branching Atoms: " << branch_atomid1 << branch_atomid2 << std::endl;

		kinematics::tree::AtomCOP branch_atom1( pose.atom_tree().atom(branch_atomid1).get_self_ptr() );
		kinematics::tree::AtomCOP branch_atom2( pose.atom_tree().atom(branch_atomid2).get_self_ptr() );

		// switch main atoms (and branching atom chirality) if necessary for having a working stub
		if ( branch_atom1->input_stub_atom2() == main_atom2 && branch_atom2->input_stub_atom2() == main_atom2 ) {
			id::AtomID temp_atomid(main_atomid1);
			kinematics::tree::AtomCOP temp_atom(main_atom1);
			main_atomid1 = main_atomid2;
			main_atom1 = main_atom2;
			main_atomid2 = temp_atomid;
			main_atom2 = temp_atom;
			temp_atomid = branch_atomid1;
			temp_atom = branch_atom1;
			branch_atomid1 = branch_atomid2;
			branch_atom1 = branch_atom2;
			branch_atomid2 = temp_atomid;
			branch_atom2 = temp_atom;
		}

		// make sure both branching atom builds off main_atom1 & center_atom
		if ( !(branch_atom1->input_stub_atom1() == center_atom && branch_atom1->input_stub_atom2() == main_atom1 &&
				branch_atom2->input_stub_atom1() == center_atom && branch_atom2->input_stub_atom2() == main_atom1) ) {

			TR.Warning << "Branching atoms (" << branch_atomid1 << branch_atomid2 << ") for segment" << main_atomid1
				<< center_atomid << main_atomid2 << "do not use segment atoms for input stub, ignoring" << std::endl;
			return 0;
		}

		// lookup the branching angle coefficients, returning 0 if they can't be found
		BranchParam2 const params(param2(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1, branch_atomid2));

		auto coef_map_iter(coef_map2_.find(params));
		if ( coef_map_iter == coef_map2_.end() ) {
			undefined_coef2_.insert(params);
			return 0;
		} else {
			coef_index = coef_map_iter->second;
		}

		if ( optimize_for_minimum ) m2_bond_angle = coef2_[coef_index].overall_theta0();

		// calculate the new torsion offsets and bond angles
		Real b1_torsion_offset;
		Real b1_bond_angle;
		Real b2_torsion_offset;
		Real b2_bond_angle;
		coef2_[coef_index].evaluate(m2_bond_angle, b1_torsion_offset, b1_bond_angle, b2_torsion_offset, b2_bond_angle);

		// set the corresponding degrees of freedom
		id::DOF_ID const branch_atom1_PHI(branch_atom1->id(), id::PHI);
		id::DOF_ID const branch_atom1_THETA(branch_atom1->id(), id::THETA);
		id::DOF_ID const branch_atom2_PHI(branch_atom2->id(), id::PHI);
		id::DOF_ID const branch_atom2_THETA(branch_atom2->id(), id::THETA);

		// use slightly different behavior depending on how the torsion offsets are ordered
		if ( branch_atom1->input_stub_atom3() == main_atom2 && branch_atom2->input_stub_atom3() == branch_atom1 ) {
			pose.set_dof(branch_atom1_PHI, b1_torsion_offset);
			pose.set_dof(branch_atom2_PHI, b2_torsion_offset - b1_torsion_offset);
			// make sure there aren't any children after branch_atom2
			// needs a const compatible next_child method
			//runtime_assert(!center_atom->next_child(branch_atom2));
		} else if ( branch_atom1->input_stub_atom3() == branch_atom2 && branch_atom2->input_stub_atom3() == main_atom2 ) {
			pose.set_dof(branch_atom1_PHI, b1_torsion_offset - b2_torsion_offset);
			pose.set_dof(branch_atom2_PHI, b2_torsion_offset);
			// make sure there aren't any children after branch_atom1
			// needs a const compatible next_child method
			//runtime_assert(!center_atom->next_child(branch_atom1));
		} else {
			TR.Warning << "Branching atoms (" << branch_atomid1 << branch_atomid2 << ") for segment" << main_atomid1
				<< center_atomid << main_atomid2 << "do not have main atom 2 as their direct descendent, ignoring"
				<< std::endl;
			return 0;
		}

		pose.set_dof(branch_atom1_THETA, pi - b1_bond_angle);
		pose.set_dof(branch_atom2_THETA, pi - b2_bond_angle);

	} else {

		TR.Warning << "Segment" << main_atomid1 << center_atomid << main_atomid2
			<< "does not have 1 or 2 branching atoms, ignoring" << std::endl;
	}

	return coef_index;
}

/// @details
/// if overall parameters are not available, single bond parameters will be returned
core::Size
BranchAngleOptimizer::overall_params(
	core::pose::Pose const & pose,
	core::id::AtomID main_atomid1,
	core::id::AtomID center_atomid,
	core::id::AtomID main_atomid2,
	core::Real & Ktheta,
	core::Real & theta0,
	core::Real & energy0
)
{
	//kinematics::tree::AtomCOP main_atom1(& pose.atom_tree().atom(main_atomid1));
	//kinematics::tree::AtomCOP const center_atom(& pose.atom_tree().atom(center_atomid));
	kinematics::tree::AtomCOP main_atom2( pose.atom_tree().atom(main_atomid2).get_self_ptr() );

	//TR << "Getting Overall Parameters: " << main_atomid1 << center_atomid << main_atomid2 << std::endl;

	Size coef_index(0);

	// make sure that the center atom has a complete set of connections
	runtime_assert(!pose.residue(center_atomid.rsd()).has_incomplete_connection(center_atomid.atomno()));

	// count the number of atoms bonded to the central atom
	Size num_neighbors(pose.residue(center_atomid.rsd()).n_bonded_neighbor_all_res(center_atomid.atomno()));

	if ( num_neighbors == 3 ) {

		// lookup the branching atom
		id::AtomID branch_atomid1;
		branching_atomid1(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1);

		//TR << "1 Branching Atom: " << branch_atomid1 << std::endl;

		kinematics::tree::AtomCOP const branch_atom1( pose.atom_tree().atom(branch_atomid1).get_self_ptr() );

		// switch main atoms if necessary for having a working stub
		if ( branch_atom1->input_stub_atom2() == main_atom2 ) {
			id::AtomID temp_atomid(main_atomid1);
			main_atomid1 = main_atomid2;
			main_atomid2 = temp_atomid;
		}

		// lookup the branching angle coefficients
		BranchParam1 const params(param1(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1));

		auto coef_map_iter(coef_map1_.find(params));
		if ( coef_map_iter == coef_map1_.end() ) {
			undefined_coef1_.insert(params);
			Ktheta = params.m1_m2_Ktheta();
			theta0 = params.m1_m2_theta0();
			energy0 = 0;
		} else {
			coef_index = coef_map_iter->second;
			Ktheta = coef1_[coef_index].overall_Ktheta();
			theta0 = coef1_[coef_index].overall_theta0();
			energy0 = coef1_[coef_index].overall_energy0();
		}

	} else if ( num_neighbors == 4 ) {

		// lookup the branching atoms
		id::AtomID branch_atomid1;
		id::AtomID branch_atomid2;
		branching_atomids2(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1, branch_atomid2);

		//TR << "2 Branching Atoms: " << branch_atomid1 << branch_atomid2 << std::endl;

		kinematics::tree::AtomCOP branch_atom1( pose.atom_tree().atom(branch_atomid1).get_self_ptr() );
		kinematics::tree::AtomCOP branch_atom2( pose.atom_tree().atom(branch_atomid2).get_self_ptr() );

		// switch main atoms (and branching atom chirality) if necessary for having a working stub
		if ( branch_atom1->input_stub_atom2() == main_atom2 && branch_atom2->input_stub_atom2() == main_atom2 ) {
			id::AtomID temp_atomid(main_atomid1);
			main_atomid1 = main_atomid2;
			main_atomid2 = temp_atomid;
			temp_atomid = branch_atomid1;
			branch_atomid1 = branch_atomid2;
			branch_atomid2 = temp_atomid;
		}

		// lookup the branching angle coefficients
		BranchParam2 const params(param2(pose, main_atomid1, center_atomid, main_atomid2, branch_atomid1, branch_atomid2));

		auto coef_map_iter(coef_map2_.find(params));
		if ( coef_map_iter == coef_map2_.end() ) {
			undefined_coef2_.insert(params);
			Ktheta = params.m1_m2_Ktheta();
			theta0 = params.m1_m2_theta0();
			energy0 = 0;
		} else {
			coef_index = coef_map_iter->second;
			Ktheta = coef2_[coef_index].overall_Ktheta();
			theta0 = coef2_[coef_index].overall_theta0();
			energy0 = coef2_[coef_index].overall_energy0();
		}

	} else {

		TR.Warning << "Segment" << main_atomid1 << center_atomid << main_atomid2
			<< "does not have 1 or 2 branching atoms, ignoring" << std::endl;
	}

	return coef_index;
}

/// @details
/// use the MMBondAngleLibrary referenced by this object
BranchParam1
BranchAngleOptimizer::param1(
	pose::Pose const & pose,
	id::AtomID const & main_atomid1,
	id::AtomID const & center_atomid,
	id::AtomID const & main_atomid2,
	id::AtomID const & branch_atomid1
) const
{
	if ( bond_angle_residue_type_param_set_ ) {

		Real m1_m2_Ktheta, m1_m2_theta0, m1_b1_Ktheta, m1_b1_theta0, m2_b1_Ktheta, m2_b1_theta0;

		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid1, center_atomid, main_atomid2, m1_m2_Ktheta, m1_m2_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid1, center_atomid, branch_atomid1, m1_b1_Ktheta, m1_b1_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid2, center_atomid, branch_atomid1, m2_b1_Ktheta, m2_b1_theta0
		);

		return BranchParam1(m1_m2_Ktheta, m1_m2_theta0, m1_b1_Ktheta, m1_b1_theta0, m2_b1_Ktheta, m2_b1_theta0,
			tolerance_);
	}

	int const main_mmtype1(pose.residue_type(main_atomid1.rsd()).atom(main_atomid1.atomno()).mm_atom_type_index());
	int const center_mmtype(pose.residue_type(center_atomid.rsd()).atom(center_atomid.atomno()).mm_atom_type_index());
	int const main_mmtype2(pose.residue_type(main_atomid2.rsd()).atom(main_atomid2.atomno()).mm_atom_type_index());
	int const branch_mmtype1(pose.residue_type(branch_atomid1.rsd()).atom(branch_atomid1.atomno()).mm_atom_type_index());

	mm_bondangle_library_citer_pair m1_m2(mm_bondangle_library_.lookup(main_mmtype1, center_mmtype, main_mmtype2));
	mm_bondangle_library_citer_pair m1_b1(mm_bondangle_library_.lookup(main_mmtype1, center_mmtype, branch_mmtype1));
	mm_bondangle_library_citer_pair m2_b1(mm_bondangle_library_.lookup(main_mmtype2, center_mmtype, branch_mmtype1));

	BranchParam1 const params((m1_m2.first->second).key1(), (m1_m2.first->second).key2(),
		(m1_b1.first->second).key1(), (m1_b1.first->second).key2(),
		(m2_b1.first->second).key1(), (m2_b1.first->second).key2(),
		tolerance_);

	return params;
}

/// @details
/// use the MMBondAngleLibrary referenced by this object
BranchParam2
BranchAngleOptimizer::param2(
	pose::Pose const & pose,
	id::AtomID const & main_atomid1,
	id::AtomID const & center_atomid,
	id::AtomID const & main_atomid2,
	id::AtomID const & branch_atomid1,
	id::AtomID const & branch_atomid2
) const
{
	if ( bond_angle_residue_type_param_set_ ) {

		Real m1_m2_Ktheta, m1_m2_theta0, m1_b1_Ktheta, m1_b1_theta0, m2_b1_Ktheta, m2_b1_theta0,
			m1_b2_Ktheta, m1_b2_theta0, m2_b2_Ktheta, m2_b2_theta0, b1_b2_Ktheta, b1_b2_theta0;

		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid1, center_atomid, main_atomid2, m1_m2_Ktheta, m1_m2_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid1, center_atomid, branch_atomid1, m1_b1_Ktheta, m1_b1_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid2, center_atomid, branch_atomid1, m2_b1_Ktheta, m2_b1_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid1, center_atomid, branch_atomid2, m1_b2_Ktheta, m1_b2_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), main_atomid2, center_atomid, branch_atomid2, m2_b2_Ktheta, m2_b2_theta0
		);
		bond_angle_residue_type_param_set_->lookup(
			pose.conformation(), branch_atomid1, center_atomid, branch_atomid2, b1_b2_Ktheta, b1_b2_theta0
		);

		return BranchParam2(m1_m2_Ktheta, m1_m2_theta0, m1_b1_Ktheta, m1_b1_theta0, m2_b1_Ktheta, m2_b1_theta0,
			m1_b2_Ktheta, m1_b2_theta0, m2_b2_Ktheta, m2_b2_theta0, b1_b2_Ktheta, b1_b2_theta0,
			tolerance_);
	}

	int const main_mmtype1(pose.residue_type(main_atomid1.rsd()).atom(main_atomid1.atomno()).mm_atom_type_index());
	int const center_mmtype(pose.residue_type(center_atomid.rsd()).atom(center_atomid.atomno()).mm_atom_type_index());
	int const main_mmtype2(pose.residue_type(main_atomid2.rsd()).atom(main_atomid2.atomno()).mm_atom_type_index());
	int const branch_mmtype1(pose.residue_type(branch_atomid1.rsd()).atom(branch_atomid1.atomno()).mm_atom_type_index());
	int const branch_mmtype2(pose.residue_type(branch_atomid2.rsd()).atom(branch_atomid2.atomno()).mm_atom_type_index());

	mm_bondangle_library_citer_pair m1_m2(mm_bondangle_library_.lookup(main_mmtype1, center_mmtype, main_mmtype2));
	mm_bondangle_library_citer_pair m1_b1(mm_bondangle_library_.lookup(main_mmtype1, center_mmtype, branch_mmtype1));
	mm_bondangle_library_citer_pair m2_b1(mm_bondangle_library_.lookup(main_mmtype2, center_mmtype, branch_mmtype1));
	mm_bondangle_library_citer_pair m1_b2(mm_bondangle_library_.lookup(main_mmtype1, center_mmtype, branch_mmtype2));
	mm_bondangle_library_citer_pair m2_b2(mm_bondangle_library_.lookup(main_mmtype2, center_mmtype, branch_mmtype2));
	mm_bondangle_library_citer_pair b1_b2(mm_bondangle_library_.lookup(branch_mmtype1, center_mmtype, branch_mmtype2));

	BranchParam2 const params((m1_m2.first->second).key1(), (m1_m2.first->second).key2(),
		(m1_b1.first->second).key1(), (m1_b1.first->second).key2(),
		(m2_b1.first->second).key1(), (m2_b1.first->second).key2(),
		(m1_b2.first->second).key1(), (m1_b2.first->second).key2(),
		(m2_b2.first->second).key1(), (m2_b2.first->second).key2(),
		(b1_b2.first->second).key1(), (b1_b2.first->second).key2(),
		tolerance_);

	return params;
}

/// @details
/// get number of single branching atom coefficients
core::Size
BranchAngleOptimizer::num_coef1() const
{
	return coef1_.size();
}

/// @details
/// get number of double branching atom coefficients
core::Size
BranchAngleOptimizer::num_coef2() const
{
	return coef2_.size();
}

/// @details
/// get number of undefined single branching atom coefficients
core::Size
BranchAngleOptimizer::num_undefined_coef1() const
{
	return undefined_coef1_.size();
}

/// @details
/// get number of undefined double branching atom coefficients
core::Size
BranchAngleOptimizer::num_undefined_coef2() const
{
	return undefined_coef2_.size();
}

/// @details
/// read known parameters from the database
void
BranchAngleOptimizer::read_database()
{
	read_coef1_default();
	read_coef2_default();
	read_coef1_user();
	read_coef2_user();
	initialized_ = true;
}

/// @details
/// write undefined parameters to the database
void
BranchAngleOptimizer::write_database() const
{
	if ( undefined_coef1_.size() ) {
		write_undefined_coef1();
		TR << "Wrote " << undefined_coef1_.size() << " undefined single branch parameters to the database" << std::endl;
	}
	if ( undefined_coef2_.size() ) {
		write_undefined_coef2();
		TR << "Wrote " << undefined_coef2_.size() << " undefined double branch parameters to the database" << std::endl;
	}
}

/// @details
/// reads records of the format:
///
/// m1_m2_Ktheta(kcal/radians^2) m1_m2_theta0(degrees)
/// m1_b1_Ktheta(kcal/radians^2) m1_b1_theta0(degrees)
/// m2_b1_Ktheta(kcal/radians^2) m2_b1_theta0(degrees)
///
/// overall_Ktheta(kcal/radians^2) overall_theta0(degrees) overall_energy0(kcal/mol)
///
/// b1_torsion_offset_A(radians) b1_torsion_offset_B(unitless) b1_torsion_offset_C(radians^-1)
/// b1_bond_angle_A(radians) b1_bond_angle_B(unitless) b1_bond_angle_C(radians^-1)
void
BranchAngleOptimizer::read_coef1(
	std::istream & in
)
{
	Real m1_m2_Ktheta;
	Real m1_m2_theta0;
	Real m1_b1_Ktheta;
	Real m1_b1_theta0;
	Real m2_b1_Ktheta;
	Real m2_b1_theta0;

	Real overall_Ktheta;
	Real overall_theta0;
	Real overall_energy0;
	Real b1_torsion_offset_A;
	Real b1_torsion_offset_B;
	Real b1_torsion_offset_C;
	Real b1_bond_angle_A;
	Real b1_bond_angle_B;
	Real b1_bond_angle_C;

	while ( in >> m1_m2_Ktheta >> m1_m2_theta0
			>> m1_b1_Ktheta >> m1_b1_theta0
			>> m2_b1_Ktheta >> m2_b1_theta0
			>> overall_Ktheta >> overall_theta0 >> overall_energy0
			>> b1_torsion_offset_A >> b1_torsion_offset_B >> b1_torsion_offset_C
			>> b1_bond_angle_A >> b1_bond_angle_B >> b1_bond_angle_C ) {

		numeric::conversions::to_radians(m1_m2_theta0);
		numeric::conversions::to_radians(m1_b1_theta0);
		numeric::conversions::to_radians(m2_b1_theta0);
		numeric::conversions::to_radians(overall_theta0);

		BranchParam1 const params(m1_m2_Ktheta, m1_m2_theta0,
			m1_b1_Ktheta, m1_b1_theta0,
			m2_b1_Ktheta, m2_b1_theta0,
			tolerance_);

		BranchCoef1 const coefs(overall_Ktheta, overall_theta0, overall_energy0,
			b1_torsion_offset_A, b1_torsion_offset_B, b1_torsion_offset_C,
			b1_bond_angle_A, b1_bond_angle_B, b1_bond_angle_C);

		auto iter(coef_map1_.find(params));

		if ( iter == coef_map1_.end() ) {
			// if the parameter key doesn't exist, add a new coefficient record and map entry
			coef1_.push_back(coefs);
			coef_map1_[params] = coef1_.size();
		} else {
			// if the parameter key already exists, replace the existing coefficients
			coef1_[iter->second] = coefs;
		}
	}
}

/// @details
/// read from an uncompressed or gzip-compressed file, returns false on failure
bool
BranchAngleOptimizer::read_coef1(
	std::string const & filename
)
{
	utility::io::izstream infile(filename);

	if ( infile ) {
		read_coef1(infile);
		infile.close();
		return true;
	}

	return false;
}

/// @details
/// read from sampling/branch_angle/branch_angle_1.txt
void
BranchAngleOptimizer::read_coef1_default()
{
	utility::io::izstream infile;
	basic::database::open(infile, "sampling/branch_angle/branch_angle_1.txt");
	read_coef1(infile);
	infile.close();
}

/// @details
/// read from branch_angle/branch_angle_1_user.txt
bool
BranchAngleOptimizer::read_coef1_user()
{
	return read_coef1(basic::database::full_name("branch_angle/branch_angle_1_user.txt", false));
}

/// @details
/// reads records of the format:
///
/// m1_m2_Ktheta(kcal/radians^2) m1_m2_theta0(degrees)
/// m1_b1_Ktheta(kcal/radians^2) m1_b1_theta0(degrees)
/// m2_b1_Ktheta(kcal/radians^2) m2_b1_theta0(degrees)
/// m1_b2_Ktheta(kcal/radians^2) m1_b2_theta0(degrees)
/// m2_b2_Ktheta(kcal/radians^2) m2_b2_theta0(degrees)
/// b1_b2_Ktheta(kcal/radians^2) b1_b2_theta0(degrees)
///
/// overall_Ktheta(kcal/radians^2) overall_theta0(degrees) overall_energy0(kcal/mol)
///
/// b1_torsion_offset_A(radians) b1_torsion_offset_B(unitless) b1_torsion_offset_C(radians^-1)
/// b1_bond_angle_A(radians) b1_bond_angle_B(unitless) b1_bond_angle_C(radians^-1)
/// b2_torsion_offset_A(radians) b2_torsion_offset_B(unitless) b2_torsion_offset_C(radians^-1)
/// b2_bond_angle_A(radians) b2_bond_angle_B(unitless) b2_bond_angle_C(radians^-1)
void
BranchAngleOptimizer::read_coef2(
	std::istream & in
)
{
	Real m1_m2_Ktheta;
	Real m1_m2_theta0;
	Real m1_b1_Ktheta;
	Real m1_b1_theta0;
	Real m2_b1_Ktheta;
	Real m2_b1_theta0;
	Real m1_b2_Ktheta;
	Real m1_b2_theta0;
	Real m2_b2_Ktheta;
	Real m2_b2_theta0;
	Real b1_b2_Ktheta;
	Real b1_b2_theta0;

	Real overall_Ktheta;
	Real overall_theta0;
	Real overall_energy0;
	Real b1_torsion_offset_A;
	Real b1_torsion_offset_B;
	Real b1_torsion_offset_C;
	Real b1_bond_angle_A;
	Real b1_bond_angle_B;
	Real b1_bond_angle_C;
	Real b2_torsion_offset_A;
	Real b2_torsion_offset_B;
	Real b2_torsion_offset_C;
	Real b2_bond_angle_A;
	Real b2_bond_angle_B;
	Real b2_bond_angle_C;

	while ( in >> m1_m2_Ktheta >> m1_m2_theta0
			>> m1_b1_Ktheta >> m1_b1_theta0
			>> m2_b1_Ktheta >> m2_b1_theta0
			>> m1_b2_Ktheta >> m1_b2_theta0
			>> m2_b2_Ktheta >> m2_b2_theta0
			>> b1_b2_Ktheta >> b1_b2_theta0
			>> overall_Ktheta >> overall_theta0 >> overall_energy0
			>> b1_torsion_offset_A >> b1_torsion_offset_B >> b1_torsion_offset_C
			>> b1_bond_angle_A >> b1_bond_angle_B >> b1_bond_angle_C
			>> b2_torsion_offset_A >> b2_torsion_offset_B >> b2_torsion_offset_C
			>> b2_bond_angle_A >> b2_bond_angle_B >> b2_bond_angle_C ) {

		numeric::conversions::to_radians(m1_m2_theta0);
		numeric::conversions::to_radians(m1_b1_theta0);
		numeric::conversions::to_radians(m2_b1_theta0);
		numeric::conversions::to_radians(m1_b2_theta0);
		numeric::conversions::to_radians(m2_b2_theta0);
		numeric::conversions::to_radians(b1_b2_theta0);
		numeric::conversions::to_radians(overall_theta0);

		BranchParam2 const params(m1_m2_Ktheta, m1_m2_theta0,
			m1_b1_Ktheta, m1_b1_theta0,
			m2_b1_Ktheta, m2_b1_theta0,
			m1_b2_Ktheta, m1_b2_theta0,
			m2_b2_Ktheta, m2_b2_theta0,
			b1_b2_Ktheta, b1_b2_theta0,
			tolerance_);

		BranchCoef2 const coefs(overall_Ktheta, overall_theta0, overall_energy0,
			b1_torsion_offset_A, b1_torsion_offset_B, b1_torsion_offset_C,
			b1_bond_angle_A, b1_bond_angle_B, b1_bond_angle_C,
			b2_torsion_offset_A, b2_torsion_offset_B, b2_torsion_offset_C,
			b2_bond_angle_A, b2_bond_angle_B, b2_bond_angle_C);

		auto iter(coef_map2_.find(params));

		if ( iter == coef_map2_.end() ) {
			// if the parameter key doesn't exist, add a new coefficient record and map entry
			coef2_.push_back(coefs);
			coef_map2_[params] = coef2_.size();
		} else {
			// if the parameter key already exists, replace the existing coefficients
			coef2_[iter->second] = coefs;
		}
	}
}

/// @details
/// read from an uncompressed or gzip-compressed file, returns false on failure
bool
BranchAngleOptimizer::read_coef2(
	std::string const & filename
)
{
	utility::io::izstream infile(filename);

	if ( infile ) {
		read_coef2(infile);
		infile.close();
		return true;
	}

	return false;
}

/// @details
/// read from sampling/branch_angle/branch_angle_2.txt, fails hard
void
BranchAngleOptimizer::read_coef2_default()
{
	utility::io::izstream infile;
	basic::database::open(infile, "sampling/branch_angle/branch_angle_2.txt");
	read_coef2(infile);
	infile.close();
}

/// @details
/// read from branch_angle/branch_angle_2_user.txt
bool
BranchAngleOptimizer::read_coef2_user()
{
	return read_coef2(basic::database::full_name("branch_angle/branch_angle_2_user.txt", false));
}

/// @details
/// reads records of the format:
///
/// m1_m2_Ktheta(kcal/radians^2) m1_m2_theta0(degrees)
/// m1_b1_Ktheta(kcal/radians^2) m1_b1_theta0(degrees)
/// m2_b1_Ktheta(kcal/radians^2) m2_b1_theta0(degrees)
void
BranchAngleOptimizer::read_undefined_coef1(
	std::istream & in
)
{
	Real m1_m2_Ktheta;
	Real m1_m2_theta0;
	Real m1_b1_Ktheta;
	Real m1_b1_theta0;
	Real m2_b1_Ktheta;
	Real m2_b1_theta0;

	while ( in >> m1_m2_Ktheta >> m1_m2_theta0
			>> m1_b1_Ktheta >> m1_b1_theta0
			>> m2_b1_Ktheta >> m2_b1_theta0 ) {

		numeric::conversions::to_radians(m1_m2_theta0);
		numeric::conversions::to_radians(m1_b1_theta0);
		numeric::conversions::to_radians(m2_b1_theta0);

		BranchParam1 const params(m1_m2_Ktheta, m1_m2_theta0,
			m1_b1_Ktheta, m1_b1_theta0,
			m2_b1_Ktheta, m2_b1_theta0,
			tolerance_);

		undefined_coef1_.insert(params);
	}
}

/// @details
/// read from an uncompressed or gzip-compressed file, returns false on failure
bool
BranchAngleOptimizer::read_undefined_coef1(
	std::string const & filename
)
{
	utility::io::izstream infile(filename);

	if ( infile ) {
		read_undefined_coef1(infile);
		infile.close();
		return true;
	}

	return false;
}

/// @details
/// read from branch_angle/branch_angle_1_undefined.txt
bool
BranchAngleOptimizer::read_undefined_coef1()
{
	return read_undefined_coef1(basic::database::full_name("branch_angle/branch_angle_1_undefined.txt", false));
}

/// @details
/// for every set of parameters, dump out three lines in the format:
///
/// m1_m2_Ktheta(kcal/radians^2) m1_m2_theta0(degrees)
/// m1_b1_Ktheta(kcal/radians^2) m1_b1_theta0(degrees)
/// m2_b1_Ktheta(kcal/radians^2) m2_b1_theta0(degrees)
///
/// each set of parameters will be followed by a blank line
void
BranchAngleOptimizer::write_undefined_coef1(
	std::ostream & out
) const
{
	std::streamsize oldprecision = out.precision();
	out << std::setprecision(16);
	for ( auto const & params : undefined_coef1_ ) {
		out << params.m1_m2_Ktheta() << " " << numeric::conversions::degrees(params.m1_m2_theta0()) << std::endl;
		out << params.m1_b1_Ktheta() << " " << numeric::conversions::degrees(params.m1_b1_theta0()) << std::endl;
		out << params.m2_b1_Ktheta() << " " << numeric::conversions::degrees(params.m2_b1_theta0()) << std::endl;
		out << std::endl;
	}
	out << std::setprecision(oldprecision);
}

/// @details
/// write to an uncompressed or gzip-compressed file, returns false on failure
bool
BranchAngleOptimizer::write_undefined_coef1(
	std::string const & filename
) const
{
	utility::io::ozstream outfile(filename);

	if ( outfile ) {
		write_undefined_coef1(outfile);
		outfile.close();
		return true;
	}

	return false;
}

/// @details
/// overwrite branch_angle/branch_angle_1_undefined.txt if undefined parameters exist
bool
BranchAngleOptimizer::write_undefined_coef1() const
{
	if ( undefined_coef1_.size() ) {
		return write_undefined_coef1(basic::database::full_name("branch_angle/branch_angle_1_undefined.txt", false));
	}

	return true;
}

/// @details
/// reads records of the format:
///
/// m1_m2_Ktheta(kcal/radians^2) m1_m2_theta0(degrees)
/// m1_b1_Ktheta(kcal/radians^2) m1_b1_theta0(degrees)
/// m2_b1_Ktheta(kcal/radians^2) m2_b1_theta0(degrees)
/// m1_b2_Ktheta(kcal/radians^2) m1_b2_theta0(degrees)
/// m2_b2_Ktheta(kcal/radians^2) m2_b2_theta0(degrees)
/// b1_b2_Ktheta(kcal/radians^2) b1_b2_theta0(degrees)
void
BranchAngleOptimizer::read_undefined_coef2(
	std::istream & in
)
{
	Real m1_m2_Ktheta;
	Real m1_m2_theta0;
	Real m1_b1_Ktheta;
	Real m1_b1_theta0;
	Real m2_b1_Ktheta;
	Real m2_b1_theta0;
	Real m1_b2_Ktheta;
	Real m1_b2_theta0;
	Real m2_b2_Ktheta;
	Real m2_b2_theta0;
	Real b1_b2_Ktheta;
	Real b1_b2_theta0;

	while ( in >> m1_m2_Ktheta >> m1_m2_theta0
			>> m1_b1_Ktheta >> m1_b1_theta0
			>> m2_b1_Ktheta >> m2_b1_theta0
			>> m1_b2_Ktheta >> m1_b2_theta0
			>> m2_b2_Ktheta >> m2_b2_theta0
			>> b1_b2_Ktheta >> b1_b2_theta0 ) {

		numeric::conversions::to_radians(m1_m2_theta0);
		numeric::conversions::to_radians(m1_b1_theta0);
		numeric::conversions::to_radians(m2_b1_theta0);
		numeric::conversions::to_radians(m1_b2_theta0);
		numeric::conversions::to_radians(m2_b2_theta0);
		numeric::conversions::to_radians(b1_b2_theta0);

		BranchParam2 const params(m1_m2_Ktheta, m1_m2_theta0,
			m1_b1_Ktheta, m1_b1_theta0,
			m2_b1_Ktheta, m2_b1_theta0,
			m1_b2_Ktheta, m1_b2_theta0,
			m2_b2_Ktheta, m2_b2_theta0,
			b1_b2_Ktheta, b1_b2_theta0,
			tolerance_);

		undefined_coef2_.insert(params);
	}
}

/// @details
/// read from an uncompressed or gzip-compressed file, returns false on failure
bool
BranchAngleOptimizer::read_undefined_coef2(
	std::string const & filename
)
{
	utility::io::izstream infile(filename);

	if ( infile ) {
		read_undefined_coef2(infile);
		infile.close();
		return true;
	}

	return false;
}

/// @details
/// read from branch_angle/branch_angle_2_undefined.txt
bool
BranchAngleOptimizer::read_undefined_coef2()
{
	return read_undefined_coef2(basic::database::full_name("branch_angle/branch_angle_2_undefined.txt", false));
}

/// @details
/// for every set of parameters, dump out three lines in the format:
///
/// m1_m2_Ktheta(kcal/radians^2) m1_m2_theta0(degrees)
/// m1_b1_Ktheta(kcal/radians^2) m1_b1_theta0(degrees)
/// m2_b1_Ktheta(kcal/radians^2) m2_b1_theta0(degrees)
/// m1_b2_Ktheta(kcal/radians^2) m1_b2_theta0(degrees)
/// m2_b2_Ktheta(kcal/radians^2) m2_b2_theta0(degrees)
/// b1_b2_Ktheta(kcal/radians^2) b1_b2_theta0(degrees)
///
/// each set of parameters will be followed by a blank line
void
BranchAngleOptimizer::write_undefined_coef2(
	std::ostream & out
) const
{
	std::streamsize oldprecision = out.precision();
	out << std::setprecision(16);
	for ( auto const & params : undefined_coef2_ ) {
		out << params.m1_m2_Ktheta() << " " << numeric::conversions::degrees(params.m1_m2_theta0()) << std::endl;
		out << params.m1_b1_Ktheta() << " " << numeric::conversions::degrees(params.m1_b1_theta0()) << std::endl;
		out << params.m2_b1_Ktheta() << " " << numeric::conversions::degrees(params.m2_b1_theta0()) << std::endl;
		out << params.m1_b2_Ktheta() << " " << numeric::conversions::degrees(params.m1_b2_theta0()) << std::endl;
		out << params.m2_b2_Ktheta() << " " << numeric::conversions::degrees(params.m2_b2_theta0()) << std::endl;
		out << params.b1_b2_Ktheta() << " " << numeric::conversions::degrees(params.b1_b2_theta0()) << std::endl;
		out << std::endl;
	}
	out << std::setprecision(oldprecision);
}

/// @details
/// write to an uncompressed or gzip-compressed file, returns false on failure
bool
BranchAngleOptimizer::write_undefined_coef2(
	std::string const & filename
) const
{
	utility::io::ozstream outfile(filename);

	if ( outfile ) {
		write_undefined_coef2(outfile);
		outfile.close();
		return true;
	}

	return false;
}

/// @details
/// overwrite branch_angle/branch_angle_2_undefined.txt if undefined parameters exist
bool
BranchAngleOptimizer::write_undefined_coef2() const
{
	if ( undefined_coef2_.size() ) {
		return write_undefined_coef2(basic::database::full_name("branch_angle/branch_angle_2_undefined.txt", false));
	}

	return true;
}

/// @details
/// get ordered branching atom around an atom with 3 neighbors
void
branching_atomid1(
	pose::Pose const & pose,
	id::AtomID main_atomid1,
	id::AtomID center_atomid,
	id::AtomID main_atomid2,
	id::AtomID & branch_atomid1
)
{
	utility::vector1<id::AtomID> neighbors(pose.conformation().bonded_neighbor_all_res(center_atomid));

	runtime_assert(neighbors.size() == 3);

	bool found_main_atomid1(false);
	bool found_main_atomid2(false);

	for ( Size i = 1; i <= 3; ++i ) {
		id::AtomID & atomid(neighbors[i]);
		if ( atomid == main_atomid1 ) {
			found_main_atomid1 = true;
		} else if ( atomid == main_atomid2 ) {
			found_main_atomid2 = true;
		} else {
			branch_atomid1 = atomid;
		}
	}

	runtime_assert(found_main_atomid1 && found_main_atomid2);
}

/// @details
/// get ordered branching atoms such that the torsion offset from main_atomid2 to branch_atomid1
/// is less than that to branch_atomid2, given that both are in the range [0, 2*pi). In other
/// words the atoms will be returned clockwise from main_atomid2, when looking from main_atomid1
/// to center_atomid.
void
branching_atomids2(
	pose::Pose const & pose,
	id::AtomID main_atomid1,
	id::AtomID center_atomid,
	id::AtomID main_atomid2,
	id::AtomID & branch_atomid1,
	id::AtomID & branch_atomid2
)
{
	Real static const pi_2(numeric::NumericTraits<Real>::pi_2());

	utility::vector1<id::AtomID> neighbors(pose.conformation().bonded_neighbor_all_res(center_atomid));

	runtime_assert(neighbors.size() == 4);

	bool found_main_atomid1(false);
	bool found_main_atomid2(false);
	bool found_branch_atomid1(false);

	for ( Size i = 1; i <= 4; ++i ) {
		id::AtomID & atomid(neighbors[i]);
		if ( atomid == main_atomid1 ) {
			found_main_atomid1 = true;
		} else if ( atomid == main_atomid2 ) {
			found_main_atomid2 = true;
		} else if ( !found_branch_atomid1 ) {
			branch_atomid1 = atomid;
			found_branch_atomid1 = true;
		} else {
			branch_atomid2 = atomid;
		}
	}

	runtime_assert(found_main_atomid1 && found_main_atomid2);

	// get dihedral offsets of the two branching atoms adjusted to [0, 2*pi)
	Real dihedral1(numeric::dihedral_radians(pose.xyz(main_atomid2), pose.xyz(main_atomid1), pose.xyz(center_atomid),
		pose.xyz(branch_atomid1)));
	if ( dihedral1 < 0 ) dihedral1 += pi_2;
	Real dihedral2(numeric::dihedral_radians(pose.xyz(main_atomid2), pose.xyz(main_atomid1), pose.xyz(center_atomid),
		pose.xyz(branch_atomid2)));
	if ( dihedral2 < 0 ) dihedral2 += pi_2;

	// switch the atoms if their dihedral offsets are in the wrong order
	if ( dihedral2 < dihedral1 ) {
		id::AtomID temp_atomid(branch_atomid1);
		branch_atomid1 = branch_atomid2;
		branch_atomid2 = temp_atomid;
	}
}

/// @details
/// get ordered branching atoms such that the torsion offset from main_atom2 to branch_atom1
/// is less than that to branch_atom2, given that both are in the range [0, 2*pi). In other
/// words the atoms will be returned clockwise from main_atom2, when looking from main_atom1
/// to the central atom.
void
get_branching_atoms2(
	kinematics::tree::AtomCOP const main_atom2,
	kinematics::tree::AtomCOP & branch_atom1,
	kinematics::tree::AtomCOP & branch_atom2
)
{
	Real static const pi_2(numeric::NumericTraits<Real>::pi_2());

	kinematics::tree::AtomCOP const parent(main_atom2->parent());

	// check to see if parent exsits and that the correct number of bonded atoms are present
	runtime_assert(parent != nullptr);
	runtime_assert(parent->get_nonjump_atom(3) != nullptr);
	runtime_assert(!parent->get_nonjump_atom(4));

	branch_atom1 = nullptr;
	branch_atom2 = nullptr;

	// iterate through the bonded atoms and find the two siblings
	for ( Size i = 1; i <= 3; ++i ) {
		kinematics::tree::AtomCOP const sibling(parent->get_nonjump_atom(i));
		TR << sibling->id() << std::endl;
		if ( sibling != main_atom2 ) {
			if ( !branch_atom1 ) {
				branch_atom1 = sibling;
			} else {
				branch_atom2 = sibling;
			}
		}
	}

	runtime_assert(branch_atom1 != nullptr);
	runtime_assert(branch_atom2 != nullptr);

	// get dihedral offsets of the two branching atoms adjusted to [0, 2*pi)
	Real dihedral1(parent->dihedral_between_bonded_children(*main_atom2, *branch_atom1));
	if ( dihedral1 < 0 ) dihedral1 += pi_2;
	Real dihedral2(parent->dihedral_between_bonded_children(*main_atom2, *branch_atom2));
	if ( dihedral2 < 0 ) dihedral2 += pi_2;

	// switch the atoms if their dihedral offsets are in the wrong order
	if ( dihedral2 < dihedral1 ) {
		kinematics::tree::AtomCOP const temp_atom(branch_atom1);
		branch_atom1 = branch_atom2;
		branch_atom2 = temp_atom;
	}
}

} // branch_angle
} // protocols
