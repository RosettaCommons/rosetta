// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RepeatPropagationMover.cc
/// @brief Can repeat both a symmetric and non symmetric protein. Modifies the pose.
/// @author TJ Brunette

// Unit headers

#ifndef INCLUDED_protocols_simple_moves_RepeatPropagationMover_hh
#define INCLUDED_protocols_simple_moves_RepeatPropagationMover_hh

#include <protocols/simple_moves/RepeatPropagationMover.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
// Utility Headers
#include <core/types.hh>
#include <core/id/SequenceMapping.hh>

// C++ Headers
#include <string>
#include <utility/vector1.hh>

//unit header
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>


static basic::Tracer TR( "protocols.simple_moves.RepeatPropagationMover" );
namespace protocols {
namespace simple_moves {
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace protocols::moves;

class RepeatPropagationMover : public protocols::moves::Mover
{
public:
	RepeatPropagationMover();
	void apply( Pose & pose ) override;
	std::string get_name() const override;
	moves::MoverOP clone() const override { return moves::MoverOP( new RepeatPropagationMover( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
private:
	void initialize_repeat_pose( Pose & pose, Pose & repeat_pose);
	void duplicate_residues_by_type(Pose & pose, Pose & repeat_pose);
	void copy_phi_psi_omega(Pose & pose, Pose & repeat_pose);
	void add_caps(Pose & pose, Pose & repeat_pose);
	std::vector<Real> get_center_of_mass(Real* coordinates, int number_of_atoms);
	void repeat_ligand(Pose & pose, Pose & repeat_pose);
	utility::vector1<Size> initial_constrained_residues(const Pose & pose);
	void fix_ligand_residues(Pose & pose, Pose & repeat_pose);
	core::id::SequenceMapping setup_repeat_seqmap(Size repeat_number,Size ligand_in_pose, Size ligand_in_repeat);
	void repeat_ligand_constraints(Pose & pose, Pose & repeat_pose);
	void determine_overlap(const Pose & pose, Pose & parent_pose,Size overlap_max_length,Size overlap_range, std::string overlap_location_pose,Size & start_overlap_parent, Size & end_overlap_parent, Size & start_overlap_pose, Size & end_overlap_pose);
	void generate_overlap(Pose & pose, Pose & parent_pose, std::string overlap_location_pose,Size start_overlap_parent, Size end_overlap_parent, Size start_overlap_pose, Size end_overlap_pose);
	void detect_last_repeat_residue(const Pose & pose);
	void trim_back_repeat_to_repair_scar(Pose & pose);
	Size first_res_;
	Size last_res_;
	bool deal_with_ideal_loop_closure_scar_;
	Size orig_number_repeats_;
	Size orig_pdb_length_;
	Size numb_repeats_;
	bool repeat_without_replacing_pose_;
	bool maintain_cap_;
	bool maintain_ligand_;
	Size nTerm_cap_size_;
	Size cTerm_cap_size_;
};
} //protocols
} //moves
#endif
