// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin Bond
///
/// @brief
/// A class for defining Bond parameters, known as Bond_types
///
/// @details
/// This class contains the "chemical" information for Bonds. This does not contain the actual
/// xyz coordinates of the class (xyz found in core/conformation/Bond.hh. The Bond_type properties
/// are assigned by the class BondSet which is initiated from the ChemicalManager. Bond type properties
/// are currently are read in from the file located chemical/Bond_type_sets/fa_standard/Bond_properties.txt.
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
#include <core/chemical/Bond.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

namespace core {
namespace chemical {

void
Bond::print( std::ostream & out ) const {
	out << distance_ << std::endl;
}

std::ostream &
operator<< (std::ostream & out, Bond const & bond){
	bond.print( out );
	return out;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // pose
} // core
