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
	utility::pointer::ReferenceCount(),
	atom_type_name_( "" ),
	dGfreeW_( 0.0 ),
	dGfreeB_( 0.0 ),
	hyd_( 0.0 ),
	hyd_deriv_( 0.0 ),
	memb_coord_( 0.0 )
{}

MEnvAtomParams::MEnvAtomParams(
	std::string const & atom_type_name,
	core::Real const dGfreeW,
	core::Real const dGfreeB,
	core::Real const hyd,
	core::Real const hyd_deriv,
	core::Vector const & memb_coord
) : utility::pointer::ReferenceCount(),
	atom_type_name_( atom_type_name ),
	dGfreeW_( dGfreeW ),
	dGfreeB_( dGfreeB ),
	hyd_( hyd ),
	hyd_deriv_( hyd_deriv ),
	memb_coord_( memb_coord )
{}

MEnvAtomParams::~MEnvAtomParams() {}

MEnvAtomParams::MEnvAtomParams( MEnvAtomParams const & src ) :
	utility::pointer::ReferenceCount( src ),
	atom_type_name_( src.atom_type_name_ ),
	dGfreeW_( src.dGfreeW_ ),
	dGfreeB_( src.dGfreeB_ ),
	hyd_( src.hyd_ ),
	hyd_deriv_( src.hyd_deriv_ ),
	memb_coord_( src.memb_coord_ )
{}

MEnvAtomParamsOP
MEnvAtomParams::clone() const {
	return MEnvAtomParamsOP( new MEnvAtomParams( *this ) );
}

} //protocols
} //membrane
} //scoring






