// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Florian Richter (floric@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/BBTorsionAndAnglesSRFD.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

static thread_local basic::Tracer tr( "core.fragment.BBTorsionAndAnglesSRFD" );

bool BBTorsionAndAnglesSRFD::apply( pose::Pose& pose, Size seqpos ) const {
	bool const success ( Parent::apply( pose, seqpos ) );

	// only move forward with changes if prior ops successful
	if ( success ) {
		for ( Size i = 1; i <= angles_.size(); ++i ) {
			pose.conformation().set_bond_angle(
				core::id::AtomID ((angles_[i].first)[0], seqpos),
				core::id::AtomID ((angles_[i].first)[1], seqpos ),
				core::id::AtomID ((angles_[i].first)[2], seqpos ),
				angles_[i].second
			);
		}
	}
	return success; //can something go wrong ?
}

/// @brief insert backbone torsions and angles into pose at position seqpos
///  if all bb torsions are moveable in MoveMap
/// @return True if *all* torsions and angles are inserted and superclass apply()
///  is successful, otherwise False.
/// @remarks This is currently all or nothing -- all torsions for seqpos
///  must be moveable because it's not entirely clear what the behavior
///  of partial angle insertion is.  In addition, DOF_IDs are not made
///  explicitly available within this class, meaning there is no way to
///  look them up within the MoveMap; the implementation in this class
///  must be changed if this is desired.
bool BBTorsionAndAnglesSRFD::apply( kinematics::MoveMap const & movemap, pose::Pose & pose, Size const seqpos ) const {
	// parent apply() successful?
	if ( !Parent::apply( movemap, pose, seqpos ) ) {
		return false; // punt
	}

	// first check to see if all torsions moveable
	bool success2 = true; // true only if all torsions moveable
	for ( Size j = 1, je = nbb(); j <= je; ++j ) {
		success2 = success2 && movemap.get( id::TorsionID( seqpos, id::BB, j ) );
	}

	// insert all angles only if all torsions moveable
	if ( success2 ) {
		for ( Size i = 1; i <= angles_.size(); ++i ) {
			pose.conformation().set_bond_angle(
				core::id::AtomID ((angles_[i].first)[0], seqpos),
				core::id::AtomID ((angles_[i].first)[1], seqpos ),
				core::id::AtomID ((angles_[i].first)[2], seqpos ),
				angles_[i].second
			);
		}
	}
	return success2;
}

bool BBTorsionAndAnglesSRFD::steal( pose::Pose const& pose, Size seqpos ) {
	runtime_assert( angles_.size() > 0 );
	bool success ( Parent::steal( pose, seqpos ) );

	for ( Size i=1; i<= angles_.size(); ++i ) {
		angles_[i].second = pose.conformation().bond_angle(
			core::id::AtomID ((angles_[i].first)[0], seqpos),
			core::id::AtomID ((angles_[i].first)[1], seqpos),
			core::id::AtomID ((angles_[i].first)[2], seqpos) );
	}
	return success; //can something go wrong ?
}

bool BBTorsionAndAnglesSRFD::is_compatible( SingleResidueFragData const& aSRFD) const {
	if ( dynamic_cast< BBTorsionAndAnglesSRFD const * > ( & aSRFD ) ) {
		BBTorsionAndAnglesSRFD const & bbtaasrfd = static_cast< BBTorsionAndAnglesSRFD const & > ( aSRFD );
		return ( (bbtaasrfd.nbb() == nbb()) && (bbtaasrfd.nangles() == nangles()) );
	}
	return false; //wrong SRFD-type (cast not successfull)
}

/// @details there is no information in the MoveMap as to which angles can be moved, so this function
/// @details will return its value solely based on the parent
bool BBTorsionAndAnglesSRFD::is_applicable( kinematics::MoveMap const& move_map, Size seqpos) const {
	if ( ! Parent::is_applicable( move_map, seqpos ) ) return false;
	return true;
}

void BBTorsionAndAnglesSRFD::show( std::ostream &out ) const {
	Parent::show( out );
}

void BBTorsionAndAnglesSRFD::read( std::istream &in ) {
	Parent::read_data( in );
}

} //fragment
} //core
