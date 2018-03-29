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

#include <protocols/pose_creation/RepeatPropagationMover.fwd.hh>
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


namespace protocols {
namespace simple_moves {

class RepeatPropagationMover : public protocols::moves::Mover
{
public:
	RepeatPropagationMover();
	RepeatPropagationMover(core::Size numb_repeats);
	void apply( core::pose::Pose & pose ) override;
	moves::MoverOP clone() const override { return moves::MoverOP( new RepeatPropagationMover( *this ) ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	std::string get_name() const override;
	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
private:
	void initialize_repeat_pose( core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void duplicate_residues_by_type(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void copy_phi_psi_omega(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void add_cap_seq_and_structure(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void add_cap_seq(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	std::vector<core::Real> get_center_of_mass(core::Real* coordinates, int number_of_atoms);
	void repeat_ligand(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	utility::vector1<core::Size> initial_constrained_residues(const core::pose::Pose & pose);
	void fix_ligand_residues(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	core::id::SequenceMapping setup_repeat_seqmap(core::Size repeat_number,core::Size ligand_in_pose, core::Size ligand_in_repeat);
	void repeat_ligand_constraints(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void determine_overlap(const core::pose::Pose & pose, core::pose::Pose & parent_pose,core::Size overlap_max_length,core::Size overlap_range, std::string overlap_location_pose,core::Size & start_overlap_parent, core::Size & end_overlap_parent, core::Size & start_overlap_pose, core::Size & end_overlap_pose);
	void generate_overlap(core::pose::Pose & pose, core::pose::Pose & parent_pose, std::string overlap_location_pose,core::Size start_overlap_parent, core::Size end_overlap_parent, core::Size start_overlap_pose, core::Size end_overlap_pose);
	void extract_repeat_info_from_pose(const core::pose::Pose & pose);
	void trim_back_repeat_to_repair_scar(core::pose::Pose & pose);
	core::Size first_res_;
	core::Size last_res_;
	core::Size numb_repeats_;
	bool ideal_repeat_;
	bool repeat_without_replacing_pose_;
	bool maintain_cap_;
	bool maintain_cap_sequence_only_;
	bool maintain_ligand_;
	core::Size nTerm_cap_size_;
	core::Size cTerm_cap_size_;
	bool extract_repeat_info_from_pose_;
	core::Size extract_repeat_template_repeat_;
	core::Size start_pose_numb_repeats_;
	core::Size start_pose_length_;
	core::Size start_pose_duplicate_residues_;

	bool deal_with_length_change_scar_;

};
} //protocols
} //moves
#endif
