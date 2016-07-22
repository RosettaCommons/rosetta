// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/ErraserMinimizerMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_ErraserMinimizerMover_HH
#define INCLUDED_protocols_farna_ErraserMinimizerMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/farna/ErraserMinimizerMover.fwd.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/id/AtomID.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <set>

namespace protocols {
namespace farna {

class ErraserMinimizerMover : public moves::Mover {

	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef utility::tag::Tag Tag;
	typedef core::id::AtomID AtomID;
	typedef moves::Movers_map Movers_map;
	typedef filters::Filters_map Filters_map;

public:

	ErraserMinimizerMover();

	virtual ~ErraserMinimizerMover() {}

	std::string get_name() const { return "ErraserMinimizerMover"; }

	void pyrimidine_flip_trial( Pose & pose );
	void setup_fold_tree( Pose & pose );

	void initialize_from_options();

	void apply( Pose & pose );

	virtual void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & /*data*/,
		Filters_map const & /*filters*/,
		Movers_map const & /*movers*/,
		Pose const & /*pose*/
	);

	bool
	check_in_bonded_list(
		AtomID const & atom_id1,
		AtomID const & atom_id2
	);

	bool
	check_in_bond_angle_list(
		AtomID const & central_atom,
		AtomID const & side_one,
		AtomID const & side_two
	);

	void
	apply_ideal_coordinates( Pose const & pose );

	void
	add_bond_constraint(
		AtomID const & atom_id1,
		AtomID const & atom_id2,
		Pose const & pose,
		core::scoring::constraints::ConstraintSetOP & cst_set
	);

	void
	add_bond_angle_constraint(
		AtomID const & atom_id1,
		AtomID const & atom_id2,
		AtomID const & atom_id3,
		Pose const & pose,
		core::scoring::constraints::ConstraintSetOP & cst_set
	);

	bool
	check_if_connected_in_atom_tree(
		Pose const & pose,
		AtomID const & atom_id1,
		AtomID const & atom_id2
	);

	// Virts and sidechain atoms that aren't the first base atom should not move
	bool
	i_want_this_atom_to_move(
		core::chemical::ResidueType const & residue_type,
		Size const & atomno
	);

	bool
	i_want_this_atom_to_move(
		Pose const & pose,
		AtomID const & atom_id
	) {
		return i_want_this_atom_to_move( pose.residue_type( atom_id.rsd() ), atom_id.atomno() );
	}

	int
	add_virtual_res( Pose & pose );

	bool
	does_atom_exist_in_reference( Pose const & pose, AtomID const & atom_id );

	void
	create_pose_reference( Pose const & pose );

	void vary_bond_geometry(
		core::kinematics::MoveMap & mm,
		Pose & pose,
		ObjexxFCL::FArray1D< bool > & allow_insert // Operationally: not fixed, cutpoint, virt
	);

	void
	add_fixed_res_constraints(
		Pose & pose,
		Size const fixed_res_num,
		Size const my_anchor
	);


private:

	utility::vector1< std::pair< AtomID, AtomID > > bonded_atom_list_;
	utility::vector1< std::pair< AtomID, std::pair< AtomID, AtomID > > > bond_angle_list_;
	Pose pose_reference_;

	bool vary_bond_geometry_;
	bool constrain_phosphate_;
	bool ready_set_only_;
	bool skip_minimize_;
	bool attempt_pyrimidine_flip_;
	std::set< core::Size > fixed_res_list_;
	utility::vector1< core::Size > cutpoint_list_;
	std::string output_pdb_name_;

	ScoreFunctionOP scorefxn_;
	ScoreFunctionOP edens_scorefxn_;
};

} //farna
} //protocols

#endif
