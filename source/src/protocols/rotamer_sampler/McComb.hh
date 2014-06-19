// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/McComb.hh
/// @brief Ensemble of Markov chain samplers for sampling combinatorially.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_McComb_HH
#define INCLUDED_protocols_rotamer_sampler_McComb_HH

// Unit headers
#include <protocols/rotamer_sampler/McComb.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/McRotamer.hh>

namespace protocols {
namespace rotamer_sampler {

class McComb : public McRotamer {
public:
	McComb();

	virtual ~McComb(){}

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Update the DOFs
	virtual void update();

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Get uniform sampling state
	virtual void set_uniform_sampling( bool const setting );

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_rotamer( McRotamerOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "McComb"; }

private:
	bool is_empty_;
	utility::vector1<McRotamerOP> rotamer_list_;
};

} //rotamer_sampler
} //protocols

#endif

