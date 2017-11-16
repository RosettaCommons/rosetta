// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/BBTorsionSRFD.cc
/// @author Oliver Lange (olange@u.washington.edu)

// Unit Headers
#include <core/fragment/BBTorsionSRFD.hh>

// Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {

static basic::Tracer tr( "core.fragment" );

/// @brief copy assignment
BBTorsionSRFD & BBTorsionSRFD::operator =( BBTorsionSRFD const & rval ) {
	if ( this != &rval ) {
		Parent::operator =( rval );

		// copy local data
		torsions_ = rval.torsions_;
		coords_ = rval.coords_;
		has_coords_ = rval.has_coords_;
	}
	return *this;
}

/// @brief insert all backbone torsions into pose at position seq_pos
bool BBTorsionSRFD::apply( pose::Pose& pose, Size seqpos ) const {
	bool const success ( Parent::apply( pose, seqpos ) );

	// only move forward with changes if prior ops successful
	if ( success ) {
		for ( Size j = 1, je = nbb(); j <= je; ++j ) {
			pose.set_torsion( id::TorsionID( seqpos, id::BB, j ), torsions_[j] );
		}
	}
	return success; //can something go wrong ? check has_torsion() ?
}

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
bool BBTorsionSRFD::apply( kinematics::MoveMap const & movemap, pose::Pose & pose, Size const seqpos ) const {
	// parent apply() successful?
	bool const success = Parent::apply( movemap, pose, seqpos );

	// only move forward with changes if prior ops successful
	if ( success ) {
		for ( Size j = 1, je = nbb(); j <= je; ++j ) {
			pose.set_torsion( id::TorsionID( seqpos, id::BB, j ), torsions_[ j ] );
		}
	}
	return success;
}

bool BBTorsionSRFD::steal( pose::Pose const& pose, Size seqpos ) {
	runtime_assert( nbb() > 0 );
	bool success ( Parent::steal( pose, seqpos ) );
	for ( Size j=1; j<= nbb() && success; ++j ) {
		torsions_[j] = pose.torsion( id::TorsionID( seqpos, id::BB, j ));
	}
	return success; //can something go wrong ? check has_torsion() ?
}

bool BBTorsionSRFD::is_compatible( SingleResidueFragData const& aSRFD) const {
	if ( dynamic_cast< BBTorsionSRFD const* > ( & aSRFD ) ) {
		BBTorsionSRFD const & bbtsrfd = static_cast< BBTorsionSRFD const & > ( aSRFD );
		return bbtsrfd.nbb() == nbb();
	};
	return false; //wrong SRFD-type (cast not successfull)
}

/// @brief check if all backbone torsions at the sequence position moveable
///  in the MoveMap
/// @return True if all backbone torsions moveable and <tt>is_applicable()</tt>
///  succeeded for superclass, otherwise False.
bool BBTorsionSRFD::is_applicable( kinematics::MoveMap const& move_map, Size seqpos) const {
	if ( ! Parent::is_applicable( move_map, seqpos ) ) return false;

	for ( Size j=1; j<= nbb(); ++j ) {
		// catch a user-error that is otherwise difficult to find:
		if ( j == 3 ) { //omega
			if ( !( move_map.get( id::TorsionID( seqpos, id::BB, j )) ) ) {
				tr.Warning << "MoveMap allows phi/psi motion but not omega motion --> "
					<< "Fragment cannot be applied --> is this intended ?"
					<< std::endl;
			}
		}
		if ( !( move_map.get( id::TorsionID( seqpos, id::BB, j )) ) ) return false;
	}
	return true;
}

void BBTorsionSRFD::show( std::ostream &out ) const {
	using namespace ObjexxFCL::format;
	Parent::show( out );
	runtime_assert( nbb() == torsions_.size() );

	// print torsions
	for ( Size j=1; j<= nbb() ; ++j ) {
		out << F( 10, 3, torsions_[ j ] );
	}

	// print cartesian coordinates
	if ( has_coordinates() ) {
		out << F(10, 3, x())
			<< F(10, 3, y())
			<< F(10, 3, z());
	}
}

void BBTorsionSRFD::read_data( std::istream &in ) {
	Parent::read_data( in );
	Real t;
	torsions_.clear();
	while ( in >> t ) {
		torsions_.push_back( t );
	}
	if ( in.eof() && in.fail() && !in.bad() ) {
		if ( nbb() ) {
			in.clear( std::ios_base::eofbit );
		}
	}
}

} //fragment
} //core
