// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/OriginalScaffoldBuildPoint.cc
/// @brief  Class implementations for the launch point geometry on the Scaffold.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

// Package headers
#include <protocols/match/output/PoseInserter.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>

// Numeric headers
#include <numeric/HomogeneousTransform.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

/// @details Initialize the coordinates to bogus yet numically palpable values.
/// If someone should accidentally try to construct a coordinate frame from
/// CA pointing toward CB reading from uninitialized variables,
/// they wouldn't hit a divide-by-zero error that might occur if CA and C
/// were both at the origin.
ProteinBackboneBuildPoint::ProteinBackboneBuildPoint() :
	parent(),
	phi_( 0.0 ),
	psi_( 0.0 ),
	N_pos_( 0.0 ),
	CA_pos_( 0.0, 0.0, 1.5 ),
	C_pos_( 0.0, 0.75, 2.25 ),
	O_pos_( 0.0, 1.4, 1.8 ),
	H_pos_( 0.0, 0.0, 0.0 ),
	HA_pos_( 0.0, 0.0, 0.0 )
{}

ProteinBackboneBuildPoint::ProteinBackboneBuildPoint( Size index ) :
	parent( index ),
	phi_( 0.0 ),
	psi_( 0.0 ),
	N_pos_( 0.0 ),
	CA_pos_( 0.0, 0.0, 1.5 ),
	C_pos_( 0.0, 0.75, 2.25 ),
	O_pos_( 0.0, 1.4, 1.8 ),
	H_pos_( 0.0, 0.0, 0.0 ),
	HA_pos_( 0.0, 0.0, 0.0 )
{}

ProteinBackboneBuildPoint::~ProteinBackboneBuildPoint() {}

/// Setters
void ProteinBackboneBuildPoint::phi( Real setting ) { phi_ = setting; }
void ProteinBackboneBuildPoint::psi( Real setting ) { psi_ = setting; }

void ProteinBackboneBuildPoint::N_pos( Vector const & setting ) { N_pos_  = setting; }
void ProteinBackboneBuildPoint::CA_pos( Vector const & setting ){ CA_pos_ = setting; }
void ProteinBackboneBuildPoint::C_pos( Vector const & setting ) { C_pos_  = setting; }
void ProteinBackboneBuildPoint::O_pos( Vector const & setting ) { O_pos_  = setting; }
void ProteinBackboneBuildPoint::H_pos( Vector const & setting ) { H_pos_  = setting; }
void ProteinBackboneBuildPoint::HA_pos( Vector const & setting ){ HA_pos_ = setting; }


OriginalBackboneBuildPoint::OriginalBackboneBuildPoint( Size index ):
	parent( index ),
	original_resid_( 0 )
{}

OriginalBackboneBuildPoint::~OriginalBackboneBuildPoint() {}

OriginalBackboneBuildPoint::OriginalBackboneBuildPoint(
	core::conformation::Residue const & res
) :
	parent(),
	original_resid_( res.seqpos() )
{
	initialize_from_residue( res );
}


OriginalBackboneBuildPoint::OriginalBackboneBuildPoint(
	core::conformation::Residue const & res,
	Size index
) :
	parent( index ),
	original_resid_( res.seqpos() )
{
	initialize_from_residue( res );
}

void
OriginalBackboneBuildPoint::initialize_from_residue(
	core::conformation::Residue const & res
)
{
	typedef numeric::HomogeneousTransform< core::Real > HTReal;

	input_conformation_ = core::conformation::ResidueOP( new core::conformation::Residue( res ) );


	debug_assert(   res.has( "N" ));
	debug_assert(  res.has( "CA" ));
	debug_assert(   res.has( "C" ));
	debug_assert(   res.has( "O" ));

	debug_assert( res.is_protein() );
	phi(    res.mainchain_torsion( 1 ));
	psi(    res.mainchain_torsion( 2 ));
	N_pos(  res.xyz( "N" ));
	CA_pos( res.xyz( "CA" ));
	C_pos(  res.xyz( "C" ));
	O_pos(  res.xyz( "O" ));

	if ( res.is_lower_terminus() || res.aa() == core::chemical::aa_pro ) {
		/// Build H
		core::Real res_phi = res.is_lower_terminus() ? 180 : phi();
		HTReal Nframe( C_pos(), CA_pos(), N_pos() );
		HTReal phi_rot;
		phi_rot.set_zaxis_rotation_deg( res_phi );
		HTReal twistedNframe = Nframe * phi_rot;
		twistedNframe.walk_along_z( 1.0 ); // N--H bond length
		H_pos( twistedNframe.point() );
	} else {
		H_pos(  res.xyz( "H" ));
	}

	if ( res.is_lower_terminus() ) {
		phi( core::pack::dunbrack::SingleResidueDunbrackLibrary::NEUTRAL_PHI );
	}

	if ( res.is_upper_terminus() ) {
		psi( core::pack::dunbrack::SingleResidueDunbrackLibrary::NEUTRAL_PSI );
	}

	if ( res.aa() == core::chemical::aa_gly ) {
		HA_pos( res.xyz( "2HA" ));
	} else {
		HA_pos( res.xyz( "HA" ));
	}
}

bool OriginalBackboneBuildPoint::compatible( ScaffoldBuildPoint const & other, bool ) const
{
	return other.compatible( *this, false );
}

/// @details It's ok for a single backbone build point to be used by multiple geometric constraints
/// at least from the perspective of the build point.  If two different side chain conformations
/// are needed for a single build point, then clearly that's impossible; however, if one of the
/// geometric constraints requires a particular sidechain conformation, and the other is using
/// the backbone of the original build point, then we're fine.  Therefore, this method always
/// returns true.
bool OriginalBackboneBuildPoint::compatible( OriginalBackboneBuildPoint const &, bool ) const
{
	return true;
}

OriginalBackboneBuildPoint::Size
OriginalBackboneBuildPoint::original_insertion_point() const
{
	return original_resid_;
}

void
OriginalBackboneBuildPoint::insert(
	Size seqpos_to_insert_at,
	Hit const & hit,
	UpstreamBuilderCOP builder,
	core::pose::Pose & pose
) const
{
	output::PoseInserter inserter( pose, seqpos_to_insert_at );
	builder->recover_hit( hit, *this, inserter );
}

}
}
}

