// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_SugarRotamer.hh
/// @brief Generate sugar pucker rotamers for RNA.
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_SugarRotamer_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_SugarRotamer_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_SugarRotamer.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSized.hh>

// Project headers
#include <core/types.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

class RNA_SugarRotamer : public RotamerSized {
public:
	RNA_SugarRotamer(
		core::Size const rsd_id,
		core::Size const pucker_state
	);

	/// @brief Initialization
	void init();

	/// @brief Reset to the first (or random if random()) rotamer.
	void reset();

	/// @brief Move to next rotamer
	void operator++();

	/// @brief Check if reach the end of rotamer list
	bool not_end() const;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose ) { apply( pose, id_ ); }

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose & pose, core::Size const i );

	/// @brief Get the total number of rotamers in sampler
	 core::Size size() const {
		runtime_assert( is_init() );
		return pucker_states_.size();
	}

	/// @brief Set the residue id being sampled
	void set_rsd_id( core::Size const setting ) { rsd_id_ = setting; }

	/// @brief Get the current pucker state.
	core::Size pucker() const {
		runtime_assert( is_init() );
		return pucker_states_[id_];
	}

	/// @brief Set the pucker_state (WHATEVER / NORTH / SOUTH)
	void set_pucker_state( core::Size const setting ) {
		set_and_reinit( pucker_state_, setting );
	}

	/// @brief Set if the sampler will skip pucker applying when input pose has
	//  same pucker assginment as sampler.
	void set_skip_same_pucker( bool const setting ) {
		skip_same_pucker_ = setting;
	}

	/// @brief Set if using RNA_IdealCoord to sample puckers
	void set_idealize_coord( bool const setting ) {
		idealize_coord_ = setting;
	}

	/// @brief Name of the class
	std::string get_name() const { return "RNA_SugarRotamer"; }

private:
	utility::vector1<core::Size> pucker_states_;

	core::Size id_;

	core::Size rsd_id_, pucker_state_;

	bool skip_same_pucker_, idealize_coord_;
};

}
}
}

#endif
