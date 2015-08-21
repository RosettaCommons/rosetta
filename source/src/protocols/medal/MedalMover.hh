// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_MEDAL_MEDALMOVER_HH
#define INCLUDED_PROTOCOLS_MEDAL_MEDALMOVER_HH

// Unit header
#include <protocols/medal/MedalMover.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>


namespace protocols {
namespace medal {

/// @class Alternating rigid body and fragment insertion moves with
/// a linearly increasing chainbreak
class MedalMover : public protocols::moves::Mover {
protected:
	typedef utility::vector1<double> Probabilities;

public:
	MedalMover();
	void apply(core::pose::Pose& pose);

	// -- jd2 -- //
	std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;

protected:
	/// @brief Closes chainbreaks in <pose>
	void do_loop_closure(core::pose::Pose* pose) const;

	/// @brief Computes per-residue sampling probabilities
	void compute_per_residue_probabilities(
		const unsigned num_residues,
		const core::sequence::SequenceAlignment& alignment,
		const protocols::loops::Loops& chunks,
		const core::kinematics::FoldTree& tree,
		const core::fragment::FragSet& fragments,
		Probabilities* probs) const;

	/// @brief Partitions the structure into non-overlapping chunks
	void decompose_structure(const core::pose::Pose& pose,
		protocols::loops::LoopsOP & chunks) const;

	/// @brief Configures a mover for biased fragment insertion
	protocols::moves::MoverOP create_fragment_mover(
		core::pose::PoseOP native,
		core::scoring::ScoreFunctionOP score,
		core::fragment::FragSetOP fragments,
		const Probabilities& probs,
		const std::string& policy,
		unsigned library_size) const;

	/// @brief Configures a mover for alternating biased fragment insertion
	/// and rigid body moves
	protocols::moves::MoverOP create_fragment_and_rigid_mover(
		const core::pose::Pose& pose,
		core::pose::PoseOP native,
		core::scoring::ScoreFunctionOP score,
		core::fragment::FragSetOP fragments,
		const Probabilities& probs,
		const std::string& policy,
		unsigned library_size) const;

	/// @brief Configures a mover for performing alternating small and shear moves
	protocols::moves::MoverOP create_small_mover(core::pose::PoseOP native,
		core::scoring::ScoreFunctionOP score) const;

	/// @brief Configures a basic score functions which callers can then specialize
	core::scoring::ScoreFunctionOP score_function() const;

	core::fragment::FragSetOP fragments_sm_;
	core::fragment::FragSetOP fragments_lg_;
};

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_MOVER_HH_
