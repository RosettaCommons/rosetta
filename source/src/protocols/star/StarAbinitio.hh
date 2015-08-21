// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/StarAbinitio.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_STAR_STARABINITIO_HH
#define INCLUDED_PROTOCOLS_STAR_STARABINITIO_HH

// Unit header
#include <protocols/star/StarAbinitio.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/SaneMinMover.hh>

namespace protocols {
namespace star {

class StarAbinitio : public protocols::moves::Mover {
public:
	StarAbinitio();
	void apply(core::pose::Pose& pose);

	/// @detail Uses the copy constructor to create a new instance
	protocols::moves::MoverOP clone() const;

	/// @detail Uses the no-argument constructor to create a new instance
	protocols::moves::MoverOP fresh_instance() const;

	/// @detail Returns the name of this mover
	std::string get_name() const;

private:
	/// @detail Sets up kinematics to keep the orientation of the aligned regions
	/// fixed with respect to one another. A virtual residue is placed at the
	/// aligned regions' center of mass. A jump from the virtual residue to the
	/// midpoint of each region is added. Interior cutpoints (i.e. those between
	/// consecutive aligned regions) are retrieved from the input parameter.
	void setup_kinematics(const protocols::loops::Loops& aligned,
		const utility::vector1<unsigned>& interior_cuts,
		core::pose::Pose & pose) const;

	/// @detail Removes virtual residue, cutpoint variants and restores simple kinematics
	void tear_down_kinematics(core::pose::Pose & pose) const;

	core::fragment::FragSetOP fragments_lg_;
	core::fragment::FragSetOP fragments_sm_;
	core::fragment::SecondaryStructureOP pred_ss_;
	protocols::simple_moves::SaneMinMoverOP minimizer_;
};

}  // namespace star
}  // namespace protocols

#endif  // PROTOCOLS_ABINITIO_STAR_STAR_ABINITIO_HH_
