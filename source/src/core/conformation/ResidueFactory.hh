// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/conformation/ResidueFactory.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_conformation_ResidueFactory_hh
#define INCLUDED_core_conformation_ResidueFactory_hh

#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// #include <core/chemical/AtomTypeSet.fwd.hh>
// #include <core/chemical/MMAtomTypeSet.fwd.hh>

namespace core {
namespace conformation {

/// a collection of functions making a single residue
class ResidueFactory
{
public:
	//  static
	//  ResidueTypeOP
	//  create_residue_type( chemical::AtomTypeSetCAP atom_types, chemical::MMAtomTypeSetCAP mm_atom_types );

	/// creates residue of desired type, coords are ideal values in some default spatial orientation
	static
	ResidueOP
	create_residue( chemical::ResidueType const & rsd_type );


	/// rotamer-style creation, uses backbone of existing residue (current_rsd)
	static
	ResidueOP
	create_residue(
		chemical::ResidueType const & rsd_type,
		Residue const & current_rsd,
		Conformation const & conformation,
		bool preserve_c_beta = false
	);

};

} // namespace conformation
} // namespace core

#endif
