// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/toolbox/AllowInsert.fwd.hh> //need to move AllowInsert to toolbox XRW2

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
//remap which score term these appear in (so it's not rna_bond_geometry)

///@details Idealization is generally performed by adding constraints to the Pose that keep bond lengths and angles within the appropriate ranges; then freeing those DOFs and performing minimization or sampling.  This protocols::moves::Mover creates bond and angle constraints compatible with idealization; it does not modify the input Pose other than by adding constraints.  If your pose is already ideal, this is likely to unidealize it.  Also, it will add a LOT of degrees of freedom to your minimization, which may lead to significant slowdowns!
class ConstrainToIdealMover : public protocols::moves::Mover {

public:

	ConstrainToIdealMover();
	ConstrainToIdealMover(ConstrainToIdealMover const & rhs);
	ConstrainToIdealMover & operator=(ConstrainToIdealMover const & rhs);
	virtual ~ConstrainToIdealMover();

	using Mover::apply;
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	virtual moves::MoverOP fresh_instance() const;
	virtual moves::MoverOP clone() const;

	///@brief supply a MoveMap to be further modified by this mover; otherwise it will modify a blank one
	void set_movemap( core::kinematics::MoveMapOP movemap );

	///@brief get the MoveMap last created during apply()
	core::kinematics::MoveMapOP get_movemap() const;

	///@brief setter for AllowInsert; makes a shallow copy
	void set_AllowInsert(protocols::toolbox::AllowInsertOP allow_insert);

	///@brief getter for AllowInsert
	protocols::toolbox::AllowInsertOP get_AllowInsert() const;

	///@brief parse XML (specifically in the context of the parser/scripting scheme)
	// virtual void parse_my_tag(
	// 	TagCOP,
	// 	basic::datacache::DataMap &,
	// 	Filters_map const &,
	// 	protocols::moves::Movers_map const &,
	// 	Pose const & );

private:

	///@most of the work happens here.  This modifies the movemap (member variable) and adds constraints to the Pose that will gently idealize angles and bond lengths
	void
	vary_bond_geometry(
		//core::kinematics::MoveMap & mm, //now a member variable
		core::pose::Pose & pose,
		core::pose::Pose const & pose_reference );

	///@brief this function generates a "reference pose" for the input pose that has the same chemistry (as best as possible), but has ideal angles and bond lengths.  These are then used to generate the constraints later.  Virtual in case you want to make the reference in a different fashion.
	virtual
	void
	create_pose_reference(
		core::pose::Pose const & pose,
		core::pose::Pose & pose_reference );

	///@brief maps to other version of function; should this type of atom be moved during idealization?
	virtual
	bool
	i_want_this_atom_to_move( core::pose::Pose const & pose, core::id::AtomID const & atom_id);

	///@brief returns whether or not this atom should move during idealiation categorically; mostly boils down to "don't move the sidechains".  Virtual in case you want to deny using a different metric.
	virtual
	bool
	i_want_this_atom_to_move( core::conformation::Residue const & residue2, core::Size const & k );

	bool
	check_if_really_connected(
		core::pose::Pose const & pose,
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2);

	bool
	check_in_bonded_list(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list );

	bool
	check_in_bond_angle_list(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::id::AtomID const & atom_id3,
		utility::vector1< std::pair< core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list );

	void
	add_bond_angle_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::id::AtomID const & atom_id3,
		utility::vector1< std::pair < core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list,
		core::pose::Pose const & pose,
		core::pose::Pose const & pose_reference,
		core::scoring::constraints::ConstraintSetOP & cst_set );

	void
	add_bond_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list,
		core::pose::Pose const & pose,
		core::pose::Pose const & pose_reference,
		core::scoring::constraints::ConstraintSetOP & cst_set );

private:
	///@brief sort of like a MoveMap?  Allowed to maintain use w/Rhiju's RNA code.
	protocols::toolbox::AllowInsertOP allow_insert_;

	///@brief MoveMap used to allow bond and angle torsions
	core::kinematics::MoveMapOP mm_;

	///@brief does the mm_ variable contain a movemap filled by apply(), or a movemap supplied to be modified?
	///Necessary to know if a self-generated movemap needs to be tossed out at the beginning of apply()
	bool supplied_movemap_;

	//Not stored locally (not reusable, no reason (?) not to save the memory by leaving it on the stack)
	//core::pose::PoseOP pose_reference_;

};//end ConstrainToIdealMover

}//namespace moves
}//namespace protocols

#endif // INCLUDED_protocols_simple_moves_ConstrainToIdealMover_HH
