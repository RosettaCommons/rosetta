// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerSamplerComb.hh
/// @brief Aggregate of multiple rotamer samplers for sampling combinatorially.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerSamplerComb_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerSamplerComb_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerSamplerComb.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSamplerBase.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerSamplerComb : public RotamerSamplerBase {
public:
	RotamerSamplerComb();

	virtual ~RotamerSamplerComb();

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
	virtual void add_external_loop_rotamer( RotamerSamplerBaseOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerSamplerComb"; }

	/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
	virtual RotamerSamplerType type() const { return COMB; }

private:
	bool is_empty_;
	utility::vector1<RotamerSamplerBaseOP> rotamer_list_;
};

} //rotamer_sampler
} //protocols

#endif

