// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerAny.hh
/// @brief Aggregate multiple samplers for sampling from any one of them.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerAny_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerAny_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerAny.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerBase.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerAny : public RotamerBase {
public:
	RotamerAny();

	virtual ~RotamerAny();

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

	/// @brief Set the random sampling state
	virtual void set_random( bool const setting );

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_rotamer( RotamerBaseOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_rotamer(
		RotamerBaseOP const & rotamer,
		core::Real const weight
	) {
		rotamer_list_.push_back( rotamer );
		weights_.push_back( weight );
		set_init( false );
	}

	/// @brief Set the weights of each rotamer sampler
	virtual void set_weights( utility::vector1<core::Real> const & weights ) {
		weights_ = weights;
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerAny"; }

private:
	core::Size curr_rotamer_;
	bool is_weighted_, is_empty_, has_empty_;
	utility::vector1<RotamerBaseOP> rotamer_list_;
	utility::vector1<core::Real> weights_, cdf_;
};

}
}

#endif

