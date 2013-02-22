// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin ResidueType
///
/// @brief
/// A class for defining a type of residue
///
/// @details
/// This class contains the "chemical" information for residues. This does not contain the actual xyz coordinates of a
/// particular residue in a specific peptide.  (xyz coordinates are found in core/conformation/Residue.hh).  A residue
/// in Rosetta can be a ligand, DNA, amino acid, or basically anything.  A residue is read in through residue_io.cc and
/// read from parameter files, generally located in the database chemical/residue_types.  For ligands, or anything that
/// is not one of the natural 20 AAs, a parameter has to be provided to rosetta through the -extra_res_fa flag.
/// residue_io.cc sets private member data in ResidueType.  The primary data that are set are: atoms, mmatoms,
/// orbitals, and properties of the particular residue type.  These properties can be modified through patches, which
/// is controlled through PatchOperations.cc.  If the residue_type of a residue is modified, the indices of atoms and
/// mmatoms and everything associated with those indices must be redefined.  This reordering of indices is taken care
/// of with the function reorder_primary_data().
///
/// Setting of primary data and then reordering is important.  Primary data for the following are described:
///
/// Atoms: Setting of atoms includes indexing the atoms into vectors, saving their names into vectors/maps, saving the
/// associated mm_atom_type into a vector, saving bond connections into vectors, etc, etc.  Since everything is
/// allocated into vectors, it is easy to reorder those vectors.  On any given residue, the heavy atoms are put into
/// the vector first, (their indices are first,) and hydrogens are put in last.
///
/// Properties: Properties of a residue include things like DNA, PROTEIN, SC_ORBITALS, CHARGED, etc.  These properties
/// indicate the type of residue it is and what properties are associated with the residue.  They are set when read in.
/// Several lines of code must be modified to get them to work, all found here in ResidueType.cc.
///
/// Orbitals: Orbitals are indexed separately from atoms.  They function much the same way as atoms, except for some
/// key differences.  To find atoms bonded to orbitals, you must provide the atom index, not the orbital index.  (I
/// haven't figured out how to get the reverse to work because of the separate indices.)  Orbital xyz coordinates are
/// not updated when atom coordinates are.  This is to keep speed consistent with just having atoms.  To output the
/// orbitals, use the flag -output_orbitals.
///
/// @author
/// Phil Bradley
/// Steven Combs - these comments
////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_ResiduePredicates_hh
#define INCLUDED_core_chemical_ResiduePredicates_hh

// Unit headers
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <core/chemical/AtomTypeSet.fwd.hh>

// Package headers

namespace core {
namespace chemical {

/////////////////////////////////////////
/////////// Graph typedefs //////////////

typedef boost::undirected_graph<
		Atom, // struct with properties of a node
		Bond // struct with properties of an edge
		/*,ResidueType*/
> ResidueGraph;

typedef ResidueGraph::vertex_descriptor VD;
typedef ResidueGraph::edge_descriptor ED;
typedef utility::vector1< VD > VDs;

typedef boost::graph_traits<ResidueGraph>::vertex_iterator VIter;
typedef boost::graph_traits<ResidueGraph>::edge_iterator EIter;
typedef std::pair<VIter, VIter> VIterPair;

typedef boost::graph_traits<ResidueGraph>::out_edge_iterator OutEdgeIter;
typedef boost::graph_traits<ResidueGraph>::in_edge_iterator InEdgeIter;
typedef std::pair<OutEdgeIter, OutEdgeIter> OutEdgeIterPair;

typedef std::map< std::string, VD > NameVDMap;
typedef std::pair<std::string, VD> NameVDPair;
typedef std::pair<NameVDMap::iterator, bool> NameVDInserted;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
////////// PREDICATES for FILTERED GRAPHS ///////////////////
////////////////////////////////////////////////////////////

class HeavyAtomFilter{
public:
	HeavyAtomFilter(){};
	HeavyAtomFilter(ResidueGraph& graph, AtomTypeSetCAP atom_types):graph_(&graph),atom_types_(atom_types){};
	bool operator()(VD const vd) const;
private:
	ResidueGraph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
	AtomTypeSetCAP atom_types_;
};
typedef boost::filtered_graph<ResidueGraph, boost::keep_all, HeavyAtomFilter> HeavyAtomGraph;

}
}
///////////////////////////////////////////////////////////////


#endif
