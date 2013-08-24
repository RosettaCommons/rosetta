// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerOneTorsion.hh
/// @brief Generate rotamer for one torsion angle.
/// @detailed
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerOneTorsion_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerOneTorsion_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerOneTorsion.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSized.hh>

// Project headers
#include <core/types.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerOneTorsion : public RotamerSized {
public:
	typedef utility::vector1<core::Real> TorsionList;

	RotamerOneTorsion();

	RotamerOneTorsion(
		core::id::TorsionID const tor_id,
		TorsionList const torsions
	);

	virtual ~RotamerOneTorsion();

	/// @brief Initialization
	virtual void init() {
		set_init( true );
		reset();
	}

	/// @brief Reset to the first (or random if is_random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose ) {	apply( pose, id_ ); }

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const i );

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const {
		runtime_assert( is_init() );
		return torsions_.size();
	}

	/// @brief Get the value of current torsion
	virtual core::Real value() const {
		runtime_assert( is_init() );
		return torsions_[id_];
	}

	/// @brief Get the value of i-th torsion
	virtual core::Real value( core::Size const i ) const {
		runtime_assert( is_init() );
		return torsions_[i];
	}

	/// @brief Set the allowed torsions in sampler
	virtual void
	set_torsions(
		TorsionList const & setting
	) {
		set_and_reinit<TorsionList>( torsions_, setting );
	}

	/// @brief Set the residue id being sampled
	virtual void set_rsd_id( core::Size const setting ) {
		torsion_id_ =
			core::id::TorsionID( setting, torsion_id_.type(), torsion_id_.torsion() );
	}

	/// @brief Set the TorsionID of the sampler
	virtual void set_torsion_id( core::id::TorsionID const setting ) {
		torsion_id_ = setting;
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerOneTorsion"; }
private:
	core::Size id_;

	TorsionList torsions_;

	core::id::TorsionID torsion_id_;
};

}
}

#endif
