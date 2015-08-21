// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/conformation/AtomGraph.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_conformation_AtomGraph_hh
#define INCLUDED_core_conformation_AtomGraph_hh

#include <core/conformation/AtomGraph.fwd.hh>


#include <core/conformation/Conformation.fwd.hh>

#include <core/types.hh>

//Auto Headers
namespace core {
namespace conformation {

void
atom_graph_from_conformation(
	Conformation const & conformation,
	AtomGraphOP pg
);

/// @brief create a pointgraph which consists of 1 node for every atom, plus 1 additional node which will be added as the last node.  The index of the additional node is returned
platform::Size
annotated_atom_graph_from_conformation(
	Conformation const & conformation,
	AtomGraphOP pg,
	PointPosition const & additional_point

);

}
}

#endif /* ATOMGRAPH_HH_ */
