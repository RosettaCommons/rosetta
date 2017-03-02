// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MainChainTorsionClasses
/// @brief a few functions used by several StepWiseProteinAnsatz classes
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/protein/MainChainTorsionSet.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/io/silent/BinarySilentStruct.hh>


#include <string>

//Auto Headers
#include <utility/fixedsizearray1.hh>


using core::Real;
using core::Size;
using core::pose::Pose;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

MainChainTorsionSet::MainChainTorsionSet( core::Real const phi, core::Real const psi, core::Real const omega )
{
	mainchain_dihedral_values_[ 1 ] = phi;
	mainchain_dihedral_values_[ 2 ] = psi;
	mainchain_dihedral_values_[ 3 ] = omega;
}

MainChainTorsionSet::MainChainTorsionSet( core::Real const phi, core::Real const psi )
{
	mainchain_dihedral_values_[ 1 ] = phi;
	mainchain_dihedral_values_[ 2 ] = psi;
	mainchain_dihedral_values_[ 3 ] = 180.0;
}
	
MainChainTorsionSet::MainChainTorsionSet( utility::fixedsizearray1< core::Real, 3 > const & mainchain_dihedral_values, core::Real const omega )
{
	mainchain_dihedral_values_[ 1 ] = mainchain_dihedral_values[ 1 ];
	mainchain_dihedral_values_[ 2 ] = mainchain_dihedral_values[ 2 ];
	mainchain_dihedral_values_[ 3 ] = mainchain_dihedral_values[ 3 ];
	mainchain_dihedral_values_[ 4 ] = omega;
}
	
MainChainTorsionSet::MainChainTorsionSet( utility::fixedsizearray1< core::Real, 3 > const & mainchain_dihedral_values )
{
	mainchain_dihedral_values_[ 1 ] = mainchain_dihedral_values[ 1 ];
	mainchain_dihedral_values_[ 2 ] = mainchain_dihedral_values[ 2 ];
	mainchain_dihedral_values_[ 3 ] = mainchain_dihedral_values[ 3 ];
	mainchain_dihedral_values_[ 4 ] = 180.0;
}


MainChainTorsionSet::~MainChainTorsionSet(){}

// ASSUME alpha
core::Real MainChainTorsionSet::phi() const{ return mainchain_dihedral_values_[ 1 ]; }
core::Real MainChainTorsionSet::psi() const{ return mainchain_dihedral_values_[ 2 ]; }
core::Real MainChainTorsionSet::omega() const{ return mainchain_dihedral_values_[ 3 ]; }

MainChainTorsionSet &
MainChainTorsionSet::operator=( MainChainTorsionSet const & src )
{
	if ( this != & src ) {
		mainchain_dihedral_values_ =  src.mainchain_dihedral_values();
	}
	return (*this );
}

} //protein
} //modeler
} //stepwise
} //protocols

