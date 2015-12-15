// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideAtomIndices.hh
/// @brief  Disulfide Atom Indices class declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_DisulfideAtomIndices_hh
#define INCLUDED_core_scoring_disulfides_DisulfideAtomIndices_hh

// Unit headers
#include <core/scoring/disulfides/DisulfideAtomIndices.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility Headers

#include <utility/vector1_bool.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace disulfides {

/// @brief This class is used by the *DisulfideEnergyContainer and the *DisulfidePotential
/// classes to rapidly index into a residue that's known to form a disulfide.  For the sake
/// of computing derivatives, there are only three atoms that need to be readily available:
/// CA, CB, and the atom which makes the disulfide bond, either SG or CEN.
/// The DisulfideEnergyContainer is responsible for keeping the indices in one of these objects
/// up-to-date with the residue it is meant to shadow.
class DisulfideAtomIndices
{

public:
	DisulfideAtomIndices( conformation::Residue const & res );

	bool atom_gets_derivatives( Size atom_index ) const;
	DisulfideDerivativeAtom derivative_atom( Size atom_index ) const;

	Size c_alpha_index() const { return c_alpha_index_; }
	Size c_beta_index()  const { return c_beta_index_;  }
	/// @brief The atom which participates in the disulfide bond
	/// SG for fullatom or CEN for centroid
	Size disulf_atom_index() const { return disulf_atom_index_; }

private:
	Size c_alpha_index_;
	Size c_beta_index_;
	Size disulf_atom_index_;

	utility::vector1< DisulfideDerivativeAtom > derivative_atom_types_;
#ifdef    SERIALIZATION
public:
	DisulfideAtomIndices();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#endif
