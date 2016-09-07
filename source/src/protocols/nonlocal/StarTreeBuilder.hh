// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/StarTreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH

// Unit headers
#include <protocols/nonlocal/StarTreeBuilder.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Package headers
#include <protocols/nonlocal/TreeBuilder.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

class StarTreeBuilder : public TreeBuilder {
protected:
	static const std::string PREFIX_INITIAL;
	static const std::string PREFIX_FINAL;

public:
	StarTreeBuilder();

	/// @brief Constructs a star fold tree by placing a virtual residue at
	/// <chunks> center of mass and adding jumps from it to a stochastically
	/// chosen anchor residue in each chunk. Cutpoints are added on chunk
	/// boundaries.
	///
	/// Important: chunks must be sorted in increasing order of start position.
	/// The simplest way to achieve this is a call to Loops::sequential_order().
	///
	/// Additionally, every residue in pose must belong to one and only one chunk.
	void set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose) override;

	/// @brief Removes the virtual residue added to <pose> in calls to set_up()
	void tear_down(core::pose::Pose* pose) override;

protected:
	/// @brief Stochastically selects an anchor position
	core::Size choose_anchor_position(const protocols::loops::Loop& chunk) const;

	/// @brief When native is available, computes rmsd of jump residues, storing
	/// the results as comments in silent file output. If specified, prefix string
	/// will precede the result in the silent file output.
	void do_compute_jump_rmsd(core::pose::Pose* model, const std::string& prefix = "") const;

private:
	/// @brief Index of the virtual residue we added to the pose in set_up()
	int virtual_res_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_
