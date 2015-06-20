// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideAtomIndices.cc
/// @brief  Disulfide Energy class implementation
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/disulfides/DisulfideAtomIndices.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace disulfides {


DisulfideAtomIndices::DisulfideAtomIndices( conformation::Residue const & res ) :
	derivative_atom_types_( res.natoms(), NO_DERIVATIVES_FOR_ATOM )
{
	core::chemical::ResidueType rt = res.type();
	
	std::string disulf_atom_name = rt.get_disulfide_atom_name();
	disulf_atom_index_ = rt.has( disulf_atom_name ) ? res.atom_index( disulf_atom_name ) : 0;
	
	if ( disulf_atom_name == "CEN" ) {
	debug_assert(rt.has("CEN") );//disulfides form to SG or CEN only
		disulf_atom_index_ = res.atom_index( "CEN" );
		derivative_atom_types_[ disulf_atom_index_ ] = CYS_CEN;
	} else {
		c_beta_index_ = rt.atom_base( disulf_atom_index_ ); //rt.has("SD") ? res.atom_index( "CG" ): res.atom_index( "CB" ) );
		c_alpha_index_ = rt.atom_base( c_beta_index_ ); //rt.has("SD") ? res.atom_index( "CB" ) : res.atom_index( "CA" ) );
		derivative_atom_types_[ c_alpha_index_ ] = CYS_C_ALPHA;
		derivative_atom_types_[ c_beta_index_  ] = CYS_C_BETA;
		derivative_atom_types_[ disulf_atom_index_ ] = CYS_S_GAMMA;
	}

}

bool
DisulfideAtomIndices::atom_gets_derivatives( Size atom_index ) const
{
	return derivative_atom_types_[ atom_index ] != NO_DERIVATIVES_FOR_ATOM;
}

DisulfideDerivativeAtom
DisulfideAtomIndices::derivative_atom( Size atom_index ) const
{
	return derivative_atom_types_[ atom_index ];
}


} // namespace disulfides
} // namespace scoring
} // namespace core
