// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/BBTorsionSRFD.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)

#ifndef INCLUDED_core_fragment_BBTorsionSRFD_HH
#define INCLUDED_core_fragment_BBTorsionSRFD_HH

// Unit Headers
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <core/fragment/SecstructSRFD.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.hh>

// C/C++ headers
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

class BBTorsionSRFD : public SecstructSRFD {
	typedef SecstructSRFD Parent;

public:
	/// @brief constructor
	BBTorsionSRFD( Size const nbb_in = 3, char secstruct = 'X', char sequence = 'X')
	: SecstructSRFD(secstruct, sequence), torsions_(nbb_in), coords_(3), has_coords_(false) {}

	/// @brief copy assignment
	BBTorsionSRFD & operator =( BBTorsionSRFD const & rval );

	/// @brief clone this object
	virtual SingleResidueFragDataOP clone() const {
		return SingleResidueFragDataOP( new BBTorsionSRFD( *this ) );
	}

	/// @brief create a new instance of this object
	virtual SingleResidueFragDataOP create() const {
		return SingleResidueFragDataOP( new BBTorsionSRFD() );
	}

	/// @brief number of backbone torsions described by this fragment
	inline Size nbb() const {
		return torsions_.size();
	}

	/// @brief set value for specific torsion in this piece of fragment.
	void set_torsion( Size const tor, Real const setting ) {
		torsions_[tor] = setting;
	}

	/// @brief get the value for a specific torsion in this fragment
	inline Real torsion( Size const torsion_number ) const {
		return torsions_[ torsion_number ];
	}

	/// @brief Returns true if this instance contains cartesian coordinates,
	/// false otherwise. Coordinates are available if the <write_ca_coords>
	/// option is enabled in the new fragment picker and rosetta++ fragments
	/// are used.
	bool has_coordinates() const {
		return has_coords_;
	}

	/// @brief Returns the x coordinate of this residue's CA
	Real x() const {
		if ( !has_coords_ ) {
			utility_exit_with_message("Cartesian coordinates uninitialized!");
		}

		return coords_[1];
	}

	/// @brief Returns the y coordinate of this residue's CA
	Real y() const {
		if ( !has_coords_ ) {
			utility_exit_with_message("Cartesian coordinates uninitialized!");
		}

		return coords_[2];
	}

	/// @brief Returns the z coordinate of this residue's CA
	Real z() const {
		if ( !has_coords_ ) {
			utility_exit_with_message("Cartesian coordinates uninitialized!");
		}

		return coords_[3];
	}

	/// @brief Convenience method for setting this residue's
	/// CA coordinates all at once
	void set_coordinates(Real x, Real y, Real z) {
		coords_[1] = x;
		coords_[2] = y;
		coords_[3] = z;
		has_coords_ = true;
	}

	/// @brief insert all backbone torsions into pose at position seq_pos
	virtual bool apply( pose::Pose&, Size seq_pos ) const;

	/// @brief insert all backbone torsions into pose at position seq_pos
	/// @param[in] movemap This MoveMap will be *ignored* at the BBTorsionSRFD level,
	///  but will be passed to any superclass <tt>apply()</tt>.
	/// @param[in,out] pose The pose to modify.
	/// @param[in] seqpos Sequence position to modify.
	/// @return True if <tt>apply()</tt> successful, False otherwise.
	/// @warning MoveMap settings at the BBTorsionSRFD level are *ignored*.
	///  For speed, does not check to see whether or not all backbone torsions
	///  are moveable in MoveMap -- use <tt>is_applicable()</tt> for this
	///  purpose prior to calling <tt>apply()</tt>.
	virtual bool apply( kinematics::MoveMap const & movemap, pose::Pose & pose, Size const seqpos ) const;

	virtual bool steal( pose::Pose const&, Size seq_pos );
	virtual bool is_compatible( SingleResidueFragData const& ) const;

	/// @brief check if all backbone torsions at the sequence position moveable
	///  in the MoveMap
	/// @return True if all backbone torsions moveable and <tt>is_applicable()</tt>
	///  succeeded for superclass, otherwise False.
	virtual bool is_applicable( kinematics::MoveMap const&, Size seq_pos ) const;

	virtual void show( std::ostream &out ) const;

	virtual void read_data( std::istream &in );

	virtual std::string type() const {
		return _static_type_name();
	}

	static std::string _static_type_name() {
		return "BBTorsion";
	}

private:
	utility::vector1<Real> torsions_;

	/// @brief Cartesian coordinates for CA
	utility::vector1<Real> coords_;

	/// @brief Indicates whether this object contains cartesian coordinates
	bool has_coords_;
};

} //fragment
} //core

#endif
