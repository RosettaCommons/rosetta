// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////

// Unit headers
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>

// Package headers

namespace core {
namespace chemical {

/////////////////////////////////////////////////////////////
////////// PREDICATES for FILTERED GRAPHS ///////////////////
////////////////////////////////////////////////////////////

bool HeavyAtomFilter::operator()(VD const vd) const{
		return (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_heavyatom();
}

bool HeavyAtomWithPolarHydrogensFilter::operator()(VD const vd) const{

	for(HeavyAtomOutEdgeIterPair ep = boost::out_edges(vd, *graph_); ep.first != ep.second; ++ep.first){
		HeavyAtomOutEdgeIter e_iter= ep.first;
		HeavyAtomED ed = *e_iter;
		HeavyAtomVD target = boost::target(ed, *graph_);
		Atom const& a =  (*graph_)[target];
		AtomType const& at = (*atom_types_)[ a.atom_type_index() ];
		if( at.is_polar_hydrogen() ) return true;
	}
	return false;
}

}
}
///////////////////////////////////////////////////////////////


