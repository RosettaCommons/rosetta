// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/BiasedFragmentMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_BIASEDFRAGMENTMOVER_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_BIASEDFRAGMENTMOVER_HH

// Unit header
#include <protocols/nonlocal/BiasedFragmentMover.fwd.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>

namespace protocols {
namespace nonlocal {

/// @class A kinematically-aware mover for performing fragment insertion.
/// User-specified, per-residue sampling probabilities allow fine grained
/// control over the simulation.
class BiasedFragmentMover : public protocols::moves::Mover {
	typedef boost::unordered_map<unsigned, core::fragment::Frame> FrameMap;
	typedef core::fragment::FragSetCOP FragSetCOP;
	typedef utility::vector1<double> Probabilities;

public:
	/// @brief Creates a new BiasedFragmentMover that selects uniformly from the
	/// available fragments at the selected insertion position.
	BiasedFragmentMover(const PolicyOP& policy, const Probabilities& probs);

	~BiasedFragmentMover() override = default;

	/// @brief Inserts a single fragment into pose.
	///
	/// Insertion position is chosen in a biased manner using the per-residue
	/// probabilities provided in the constructor. The decision on which fragment
	/// to insert from the fragment library is delegated to the policy specified
	/// in the constructor.
	///
	/// Respects the underlying kinematics of the system.
	void apply(core::pose::Pose& pose) override;

	/// @brief Returns the name of this mover
	std::string get_name() const override;

private:
	/// @brief Creates a position-indexable list of Frames
	void initialize_library();

	/// @brief Generates cdf from pdf
	void initialize_probabilities();

	/// @brief Verifies that the probability of selecting invalid positions is 0
	void verify_probabilities_or_die(const core::kinematics::FoldTree& tree) const;

	/// @brief Returns a randomly chosen position according to the input probabilities
	unsigned random_position() const;


	// -- Members --

	/// @brief Avoid creating a useless MoveMap for each call to Frame::apply().
	core::kinematics::MoveMap movable_;

	/// @brief Fragment library
	FragSetCOP fragments_;

	/// @brief Position-indexable Frame lookup
	FrameMap frames_;

	/// @brief Guidance for selecting the fragment to be inserted at a given position
	PolicyOP policy_;

	/// @brief PDF of residue sampling probabilities. Must remain in sync with cdf_.
	Probabilities pdf_;

	/// @brief CDF of residue sampling probabilities. Must remain in sync with pdf_.
	Probabilities cdf_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_BIASED_FRAGMENT_MOVER_HH_
