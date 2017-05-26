// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/simple_moves/ConstrainToIdealMover.hh
/// @brief Adds to your pose constraints suitable for idealization of bond lengths and angles
/// @author Steven Lewis (smlewi@gmail.com); Rhiju Das

#ifndef INCLUDED_protocols_simple_moves_ConstrainToIdealMover_hh
#define INCLUDED_protocols_simple_moves_ConstrainToIdealMover_hh

// Unit Headers
#include <protocols/simple_moves/ConstrainToIdealMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh> //need to move AtomLevelDomainMap to toolbox XRW2

#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <utility/vector1.hh>


// Utility Headers
// #include <core/types.hh>

namespace protocols {
namespace simple_moves {


//////////////////WARNING WARNING WARNING////////////////////////////////////////////////////
//This code has been refactored from RNA_minimizer, but has not yet been updated to work with
//protein poses - I wouldn't suggest trying to use it yet!
//TODO for proteins
//make_pose_from_sequence from annotated sequence (?)
//test with cen and fa poses
//check RNA assumptions in i_want_this_atom_to_move

/// @details Idealization is generally performed by adding constraints to the Pose that keep bond lengths and angles within the appropriate ranges; then freeing those DOFs and performing minimization or sampling.  This protocols::moves::Mover creates bond and angle constraints compatible with idealization; it does not modify the input Pose other than by adding constraints.  If your pose is already ideal, this is likely to unidealize it.  Also, it will add a LOT of degrees of freedom to your minimization, which may lead to significant slowdowns!

class ConstrainToIdealMover : public protocols::moves::Mover {

public:

	ConstrainToIdealMover();
	// ConstrainToIdealMover(ConstrainToIdealMover const & rhs);
	ConstrainToIdealMover & operator = ( ConstrainToIdealMover const & rhs );
	~ConstrainToIdealMover() override;

	using Mover::apply;
	void apply( core::pose::Pose & pose ) override;
	void apply( core::pose::Pose & pose, core::kinematics::MoveMap & mm );

	std::string get_name() const override;

	moves::MoverOP fresh_instance() const override;
	moves::MoverOP clone() const override;

	/// @brief setter for AtomLevelDomainMap; makes a shallow copy
	void set_atom_level_domain_map( protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map );

	/// @brief getter for AtomLevelDomainMap
	protocols::toolbox::AtomLevelDomainMapCOP get_atom_level_domain_map() const;

	void set_score_type( core::scoring::ScoreType const setting );
	void set_just_rna_backbone( bool const setting ){ just_rna_backbone_ = setting; }
	void set_just_polar_hydrogens( bool const setting ){ just_polar_hydrogens_ = setting; }

	void set_fix_lengths( bool const & setting ){ fix_lengths_ = setting; }
	void set_fix_angles( bool const & setting ){ fix_angles_ = setting; }
	void set_fix_torsions( bool const & setting ){ fix_torsions_ = setting; }

private:

	/// @most of the work happens here.  This modifies the movemap and adds constraints to the Pose that will gently idealize angles and bond lengths
	void
	vary_bond_geometry(
		core::pose::Pose & pose,
		core::kinematics::MoveMap & mm,
		core::pose::Pose const & pose_reference ) const;

	/// @brief this function generates a "reference pose" for the input pose that has the same chemistry (as best as possible), but has ideal angles and bond lengths.  These are then used to generate the constraints later.  Virtual in case you want to make the reference in a different fashion.
	virtual
	void
	create_pose_reference(
		core::pose::Pose const & pose,
		core::pose::Pose & pose_reference );

	/// @brief maps to other version of function; should this type of atom be moved during idealization?
	virtual
	bool
	i_want_this_atom_to_move( core::pose::Pose const & pose, core::id::AtomID const & atom_id ) const;

	/// @brief returns whether or not this atom should move during idealiation categorically; mostly boils down to "don't move the sidechains".  Virtual in case you want to deny using a different metric.
	virtual
	bool
	i_want_this_atom_to_move( core::conformation::Residue const & residue2, core::Size const & k ) const;

	bool
	check_if_really_connected(
		core::pose::Pose const & pose,
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2 ) const;

	bool
	check_in_bonded_list(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list ) const;

	bool
	check_in_bond_angle_list(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::id::AtomID const & atom_id3,
		utility::vector1< std::pair< core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list ) const;

	void
	add_bond_angle_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::id::AtomID const & atom_id3,
		utility::vector1< std::pair < core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list,
		core::pose::Pose const & pose,
		core::pose::Pose const & pose_reference,
		core::scoring::constraints::ConstraintSetOP & cst_set ) const;

	void
	add_bond_dihedral_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::id::AtomID const & atom_id3,
		core::id::AtomID const & atom_id4,
		core::pose::Pose const & pose,
		core::pose::Pose const & pose_reference,
		core::scoring::constraints::ConstraintSetOP & cst_set ) const;

	void
	add_bond_length_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list,
		core::pose::Pose const & pose,
		core::pose::Pose const & pose_reference,
		core::scoring::constraints::ConstraintSetOP & cst_set ) const;

private:
	/// @brief atom_level_domain_map has info on which atoms should move; complementary to move_map (which instead focuses on DOFs).
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map_;

	core::Real const bond_length_sd_;
	core::Real const bond_length_sd_polar_hydrogen_;
	core::Real const bond_angle_sd_;
	core::Real const bond_angle_sd_polar_hydrogen_;
	core::Real const bond_torsion_sd_;
	core::Real const bond_torsion_sd_polar_hydrogen_;

	core::scoring::ScoreType score_type_;
	bool just_rna_backbone_;
	bool just_polar_hydrogens_;

	bool const legacy_dof_allow_move_;
	bool const verbose_;

	bool fix_lengths_;
	bool fix_angles_;
	bool fix_torsions_;
	bool disallow_vary_geometry_proton_chi_;

}; //end ConstrainToIdealMover

void
setup_vary_rna_bond_geometry( core::kinematics::MoveMap & mm, core::pose::Pose & pose,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::scoring::ScoreType score_type = core::scoring::rna_bond_geometry );

void
setup_vary_polar_hydrogen_geometry( core::kinematics::MoveMap & mm,
	core::pose::Pose & pose, toolbox::AtomLevelDomainMapCOP atom_level_domain_map );

} //namespace moves
} //namespace protocols

#endif // INCLUDED_protocols_simple_moves_ConstrainToIdealMover_HH
