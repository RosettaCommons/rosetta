// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/ResidueSupport.hh
/// @brief support functions for class residue; functions that
/// should not be included as part of the class.
/// @author Phil Bradley


//#ifndef INCLUDED_core_chemical_ResidueSupport_HH
//#define INCLUDED_core_chemical_ResidueSupport_HH

// Package Headers
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/graph/Graph.hh>

//Auto Headers
#include <ObjexxFCL/FArray2D.hh>


// ObjexxFCL Headers
// Commented by inclean daemon #include <ObjexxFCL/FArray2D.hh>

namespace core {
namespace chemical {

ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res )
{
	using namespace graph;
	Graph g;

	g.set_num_nodes( res.natoms() );
	for ( uint ii = 1; ii <= res.natoms(); ++ii )
	{
		AtomIndices const & ii_bonded = res.bonded_neighbor( ii  );
		for ( Size jj = 1; jj <= ii_bonded.size(); ++jj)
		{
			if ( ii_bonded[ jj ] > ii )
				g.add_edge( ii, ii_bonded[ jj ] );
		}
	}
	return g.all_pairs_shortest_paths();
}

}
}

//#endif
