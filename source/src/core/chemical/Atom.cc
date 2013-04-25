// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin Atom
///
/// @brief
/// A class for defining atom parameters, known as atom_types
///
/// @details
/// This class contains the "chemical" information for atoms. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Atom.hh. The atom_type properties
/// are assigned by the class AtomSet which is initiated from the ChemicalManager. Atom type properties
/// are currently are read in from the file located chemical/atom_type_sets/fa_standard/atom_properties.txt.
/// These properties contain the the properties of LJ_RADIUS, LJ_WDEPTH, LK_DGRFREE, LK_LAMBDA, LK_VOLUME.
/// These properties are used in the scoring function fa_atr, fa_rep, fa_sol, which is located in the Etable (core/scoring/etable/Etable.hh)
/// Additional parameters are acceptor/donor, hybridzation, and orbital paramaters.
///
///
///
/// @author
/// Phil Bradley
/// Steven Combs - comments
///
///
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////

// Rosetta headers
#include <core/chemical/Atom.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <algorithm>

namespace core {
namespace chemical {

void
Atom::print(
	std::ostream & out
) const {
	out << "Name: " << name() << std::endl;
//	out << "Type: " << type() << std::endl;
	out << "MM Name: " << mm_name() << std::endl;
	out << "atom_type_index: " << atom_type_index() << std::endl;
	out << "mm_atom_type_index: " << mm_atom_type_index() << std::endl;
	out << "charge: " << charge() << std::endl;
	//out << "xyz: " << xyz() << std::endl;
	//out << "icoor: " << icoor() << std::endl;
	out << std::endl;
}

std::ostream &
operator<< (std::ostream & out, Atom const & atom ){
	atom.print( out );
	return out;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // pose
} // core
