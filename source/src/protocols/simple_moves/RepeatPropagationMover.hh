// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RepeatPropagationMover.cc
/// @brief Can repeat both a symmetric and non symmetric protein. Modifies the pose.
/// @author TJ Brunette

// Unit headers

#ifndef INCLUDED_protocols_simple_moves_RepeatPropagationMover_hh
#define INCLUDED_protocols_simple_moves_RepeatPropagationMover_hh

#include <protocols/simple_moves/RepeatPropagationMover.fwd.hh>
// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <utility/vector1.hh>

//unit header
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>


static basic::Tracer TR( "protocols.simple_moves.RepeatPropagationMover" );
namespace protocols {
namespace simple_moves {
using namespace core;
using namespace std;
using namespace protocols::moves;

class RepeatPropagationMover : public protocols::moves::Mover
{
public:
	RepeatPropagationMover();
	virtual void apply( Pose & pose );
	std::string get_name() const;
	moves::MoverOP clone() const { return moves::MoverOP( new RepeatPropagationMover( *this ) ); }
	virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	void initialize_repeat_pose( core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void duplicate_residues_by_type(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	void copy_phi_psi_omega(core::pose::Pose & pose, core::pose::Pose & repeat_pose);
	Size first_res_;
	Size last_res_;
	Size numb_repeats_;
	bool repeat_without_replacing_pose_;
	bool maintain_cap_seq_and_structure_;
	bool maintain_cap_sequence_alone_;
	Size nTerm_cap_size_;
	Size cTerm_cap_size_;
};
} //protocols
} //moves
#endif
