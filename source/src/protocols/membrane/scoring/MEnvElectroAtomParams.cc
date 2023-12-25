// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/scoring/MEnvElectroAtomParams.cc
/// @brief A container for storing implicit membrane environment electrostatic parameters and derivatives
/// @author rituparna Samanta (rsamant2@jhu.edu)

#include <protocols/membrane/scoring/MEnvElectroAtomParams.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.membrane.scoring.MEnvElectroAtomParams" );


namespace protocols {
namespace membrane {
namespace scoring {

MEnvElectroAtomParams::MEnvElectroAtomParams() :
	utility::VirtualBase(),
	atom_type_name_(),
	charge_(),
	elec_field_(),
	elec_field_gradient_(),
	f1_(),
	f2_()
{}

MEnvElectroAtomParams::MEnvElectroAtomParams(
	std::string const & atom_type_name,
	core::Real const charge,
	core::Real const elec_field,
	core::Real const elec_field_gradient,
	core::Vector const f1,
	core::Vector const f2
) : utility::VirtualBase(),
	atom_type_name_( atom_type_name ),
	charge_( charge ),
	elec_field_( elec_field ),
	elec_field_gradient_( elec_field_gradient ),
	f1_( f1 ),
	f2_( f2 )
{}

MEnvElectroAtomParams::~MEnvElectroAtomParams() {}

MEnvElectroAtomParamsOP
MEnvElectroAtomParams::clone() const {
	return MEnvElectroAtomParamsOP( new MEnvElectroAtomParams( *this ) );
}

} //protocols
} //membrane
} //scoring






