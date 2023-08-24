// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/scoring/MEnvAtomParams.cc
/// @brief A container for storing memrbane environemnt parameters and derivatives
/// @author Rebecca Alford (rfalford12@gmail.com)

#include <protocols/membrane/scoring/MEnvAtomParams.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.membrane.scoring.MEnvAtomParams" );


namespace protocols {
namespace membrane {
namespace scoring {

MEnvAtomParams::MEnvAtomParams() :
	utility::VirtualBase(),
	atom_type_name_( "" ),
	dGfreeW_( 0.0 ),
	dGfreeB_( 0.0 ),
	hyd_( 0.0 ),
	f1_( 0.0 ),
	f2_( 0.0 )
{}

MEnvAtomParams::MEnvAtomParams(
	std::string const & atom_type_name,
	core::Real const dGfreeW,
	core::Real const dGfreeB,
	core::Real const hyd,
	core::Vector const f1,
	core::Vector const f2
) : utility::VirtualBase(),
	atom_type_name_( atom_type_name ),
	dGfreeW_( dGfreeW ),
	dGfreeB_( dGfreeB ),
	hyd_( hyd ),
	f1_( f1 ),
	f2_( f2 )
{}

MEnvAtomParams::~MEnvAtomParams() {}

MEnvAtomParams::MEnvAtomParams( MEnvAtomParams const & src ) :
	utility::VirtualBase( src ),
	atom_type_name_( src.atom_type_name_ ),
	dGfreeW_( src.dGfreeW_ ),
	dGfreeB_( src.dGfreeB_ ),
	hyd_( src.hyd_ ),
	f1_( src.f1_ ),
	f2_( src.f2_ )
{}

MEnvAtomParamsOP
MEnvAtomParams::clone() const {
	return MEnvAtomParamsOP( new MEnvAtomParams( *this ) );
}

} //protocols
} //membrane
} //scoring






