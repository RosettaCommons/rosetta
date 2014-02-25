// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/McOneTorsion.hh
/// @brief Markov chain sampler for one torsion.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_McOneTorsion_HH
#define INCLUDED_protocols_rotamer_sampler_McOneTorsion_HH

// Unit headers
#include <protocols/rotamer_sampler/McOneTorsion.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/McRotamer.hh>

// Project headers
#include <core/id/TorsionID.hh>

namespace protocols {
namespace rotamer_sampler {

class McOneTorsion : public McRotamer {
public:

	McOneTorsion(
		core::id::TorsionID const & tor_id,
		core::Real const start_torsion
	);

	/// @brief Initialization
	void init() {
		set_init( true );
		reset();
	}

	/// @brief Reset to current angle
	void reset() { active_angle_ = stored_angle_; }

	/// @brief Update the active angle based on stored
	/// (do not update stored_angle_)
	void operator++();

	/// @brief Update the stored angle to match active angle
	void update() { stored_angle_ = active_angle_; }

	/// @brief Apply the active torsion to pose
	void apply( core::pose::Pose & pose );

	/// @brief Get the stored angle
	core::Real stored_angle() const { return stored_angle_; }

	/// @brief Get the active angle
	core::Real active_angle() const { return active_angle_; }

	/// @brief Set the stored angle
	void set_angle( core::Real const setting	) { stored_angle_ = setting; }

	/// @brief Set the angle range
	void set_angle_range( core::Real const min, core::Real const max ) {
		angle_max_ = max;
		angle_min_ = min;
	}

	/// @brief Set the standard deviation of Gaussian sampler
	void set_gaussian_stdev( core::Real const setting ) {
		stdev_ = setting;
		// Hack: If stdev < 0, do uniform_sampling instead
		set_uniform_sampling( stdev_ < 0 );
	}

	/// @brief Set the TorsionID of the sampler
	void set_torsion_id( core::id::TorsionID const & setting ) {
		torsion_id_ = setting;
	}

	/// @brief Name of the class
	std::string get_name() const { return "McOneTorsion"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return MC_ONE_TORSION; }

private:
	bool check_angle_in_range( core::Real angle ) const;
	void regularize_angle( core::Real & angle );
	core::Real stored_angle_, active_angle_, angle_min_, angle_max_, stdev_;
	core::id::TorsionID torsion_id_;
};

}
}

#endif
