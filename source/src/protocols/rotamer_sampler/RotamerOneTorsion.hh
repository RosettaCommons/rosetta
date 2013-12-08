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
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerOneTorsion_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerOneTorsion_HH

// Unit headers
#include <protocols/rotamer_sampler/RotamerOneTorsion.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerOneValue.hh>

// Project headers
#include <core/id/TorsionID.hh>

namespace protocols {
namespace rotamer_sampler {

class RotamerOneTorsion : public RotamerOneValue {
public:
	using RotamerBase::TorsionList;

	RotamerOneTorsion();

	RotamerOneTorsion(
		core::id::TorsionID const & tor_id,
		TorsionList const & torsions
	);

	virtual ~RotamerOneTorsion();

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose ) {	apply( pose, id_ ); }

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, Size const i );

	/// @brief Set the allowed torsions in sampler
	virtual void set_torsions( TorsionList const & setting	) {
		set_values( setting );
	}

	/// @brief Set the residue id being sampled
	virtual void set_rsd_id( core::Size const setting ) {
		torsion_id_ =
			core::id::TorsionID( setting, torsion_id_.type(),
			torsion_id_.torsion() );
	}

	/// @brief Set the TorsionID of the sampler
	virtual void set_torsion_id( core::id::TorsionID const & setting ) {
		torsion_id_ = setting;
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "RotamerOneTorsion"; }

private:

	core::id::TorsionID torsion_id_;

};

}
}

#endif
