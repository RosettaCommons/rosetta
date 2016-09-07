// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/loophash/FastGapMover.hh
/// @brief protocols for closing gaps
/// @author Ken Jung <kenjung@uw.ed>


#ifndef INCLUDED_protocols_loophash_FastGapMover_hh
#define INCLUDED_protocols_loophash_FastGapMover_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loophash/FastGapMover.fwd.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LocalInserter.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace loophash {


/// @brief Mover class for closing gaps.
/// This Mover checks for any gaps using residue residue distances
/// Then eats back at the chain surrounding it until loophash finds
/// a fragment that fits in the space without changing the rest of
/// pose too much.

class FastGapMover : public moves::Mover {

public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	FastGapMover();

	/// @brief clone has to be overridden only if clone invocation is expected.
	moves::MoverOP clone() const override {
		return moves::MoverOP( new FastGapMover( *this ) );
	}

	moves::MoverOP fresh_instance() const override {
		return moves::MoverOP( new FastGapMover );
	}

	void
	apply( Pose & pose ) override;

	std::string get_name() const override;

	FastGapMover &
	min_rms( core::Real const setting )
	{
		min_rms_ = setting;
		return *this;
	}

	FastGapMover &
	max_rms( core::Real const setting )
	{
		max_rms_ = setting;
		return *this;
	}

	FastGapMover &
	non_ideal( bool const setting )
	{
		non_ideal_ = setting;
		return *this;
	}

private:
	// methods
	void find_next_gap( Pose & pose, Size & idx, Real & gap_distance );


private:
	// data
	LocalInserter_SimpleMinOP simple_inserter_;
	LoopHashLibraryOP lhlibrary_;
	LoopHashSamplerOP lhsampler_;

	core::Size min_loop_size_;
	core::Size max_loop_size_;
	core::Real max_rms_;
	core::Real min_rms_;
	bool non_ideal_;

};

} // loophash
} // protocols

#endif
