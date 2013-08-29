// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSizedAny.hh
/// @brief Aggregate multiple samplers for sampling from any one of them.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSizedAny_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSizedAny_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerSizedAny.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSized.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerSizedAny : public RotamerSized {
public:
	RotamerSizedAny();

	virtual ~RotamerSizedAny();

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const i );

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const {
		runtime_assert( is_init() );
		return size_;
	}

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_rotamer( RotamerSizedOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		size_list_.clear();
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerSizedAny"; }
private:
	/// @brief Convert an id number to the sampler state pair.
	std::pair<core::Size, core::Size> id2state( core::Size const id ) const;

	core::Size size_, id_;

	std::pair<core::Size, core::Size> curr_state_;

	utility::vector1<core::Size> size_list_;

	utility::vector1<RotamerSizedOP> rotamer_list_;
};

}
}

#endif

