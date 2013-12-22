// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/SimplePDB_Atom.hh
///
/// @brief
/// @author will sheffler


#ifndef INCLUDED_core_scoring_packstat_SimplePDB_Atom_hh
#define INCLUDED_core_scoring_packstat_SimplePDB_Atom_hh


// Project forward headers
#include <core/scoring/packstat/SimplePDB_Atom.fwd.hh>
#include "core/scoring/packstat/types.hh"

#include <string>

namespace core {
namespace scoring {
namespace packstat {


/// @brief
struct SimplePDB_Atom
{
	friend std::ostream & operator>> ( std::ostream & out, SimplePDB_Atom & pdb  );
	friend std::ostream & operator<< ( std::ostream & out, SimplePDB_Atom const & atom );

	SimplePDB_Atom() :
		num(-12345),
		resnum(),
		chain(),
		icode(),
		x(),
		y(),
		z(),
		occ(),
		bfac(),
		sasa(),
		radius()
	{}

	~SimplePDB_Atom() {}

	std::string ATOM, type, res, lastcol, whole_line;
	int num, resnum;
	char chain,icode;
	PackstatReal x, y, z, occ, bfac, sasa, radius;
	numeric::xyzVector<PackstatReal> xyz;

}; // SimplePDB_Atom


} // namespace packstat
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_packstat_SimplePDB_Atom_HH
