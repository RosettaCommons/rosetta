// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/ResidueConnection.hh
/// @brief  Inter-residue chemical bond connection point class declaration.
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_chemical_ResidueConnection_hh
#define INCLUDED_core_chemical_ResidueConnection_hh


// Unit headers
#include <core/chemical/ResidueConnection.fwd.hh>

// Project headers
#include <core/chemical/AtomICoor.hh>

// Utility headers

// C++ headers
// Commented by inclean daemon #include <string>


namespace core {
namespace chemical {

/// @brief A simple class marking atoms at inter-residue connections.
///
/// Each residue type specifies some number of positions at which it is expecting
/// to form a chemical bond with another residue. Think of them as ports: they are
/// parts of the residue where there are chemical bonds beyond the intra-residue
/// chemical bonds are expected -- places where they can be chemically linked
/// to the outside world.  A conformation::Residue will require that its
/// ResidueConnections be fulfilled by other Residues -- the ResConnID class
/// describes how two residues are connected: e.g., the third ResConnID for
/// residue 10 would say "I connect to residue 58 at residue 58's third residue
/// connection" if residue 10 and residue 58 were disulfide bonded as the disulfide
/// connection id is "3" for two mid-protein cystine residues.  The advantages
/// of separating ResidueConnections from atoms themselves are that 1) it allows
/// multiple residue connections to stem from the same atom -- useful for single-atom
/// residues, such as coordinated metals (Zn, Mg), and 2) it allows one residue
/// to change its set of atoms without invalidating the bond information (e.g. the atom
/// index) on its partner.  For example, if a chain-break were placed between residues
/// 57 and 58, then residue 58 will get an extra C-prev virtual atom, and the index of
/// SG will change.  Residue 10, if it had recorded the SG index would have to find
/// SG's new index.  If instead, the connection point is represented simply as
/// connection point 3, and if the new residue type (the chainbreak disulfide residue)
/// has the same number of residue connections as the original residue type (it will!)
/// then nothing about residue 10 needs to be updated.

class ResidueConnection {
public:
	/// @brief default constructor
	ResidueConnection():
		atomno_( 0 ),
		icoor_(),
		index_( 0 ),
		vertex_(0)
	{}

	/// @brief constructor with atom index number
	ResidueConnection(
		int const atomno_in,
		VD const vertex
	):
		atomno_( atomno_in ),
		icoor_(),
		index_( 0 ),
		vertex_(vertex)
	{}


	//DO NOT USE unless you add a vertex constructor with this!!!!!!
/*	/// @brief constructor with atom index number and AtomICoor
	ResidueConnection(
		int const atomno_in,
		AtomICoor const & icoor_in
	):
		atomno_( atomno_in ),
		icoor_( icoor_in ),
		index_( 0 )
	{}*/



	//DO NOT USE unless you add a vertex constructor with this!!!!!!
	/// @brief constructor with atom index number, AtomICoor, and connection index
/*	ResidueConnection(
		int const atomno_in,
		AtomICoor const & icoor_in,
		int const index
	):
		atomno_( atomno_in ),
		icoor_( icoor_in ),
		index_( index )
	{}*/

	/// @brief get atom index number
	int
	atomno() const
	{
		return atomno_;
	}

	/// @brief set atom index number
	void
	atomno( Size const atomno_in )
	{
		atomno_ = atomno_in;
	}


	/// @brief get the vetex associated with this residue connection
	VD
	vertex() const
	{
		return vertex_;
	}

	/// @brief set the vertex of this residue connection
	void
	vertex(VD const vertex)
	{
		vertex_ = vertex;
	}

	/// @brief get atom's AtomICoor
	AtomICoor const &
	icoor() const
	{
		return icoor_;
	}

	/// @brief set atom's AtomICoor
	void
	icoor( AtomICoor const & ic )
	{
		icoor_ = ic;
	}

	int index() const { return index_; }
	void index( int index_in ) { index_ = index_in; }

private:

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
			ar & atomno_;
			ar & icoor_;
			ar & index_;
	}
#endif

	/// atom index number
	int atomno_;
	/// atom AtomICoor
	AtomICoor icoor_;
	/// Which residue connection # am I in my owners list of residue connections?
	int index_;
	VD vertex_; //my vertex which corresponds to the atomno
};

} // chemical
} // core



#endif
// 	int other_rsd_;
// 	std::string other_atom_name_;
// 	///
// 	int
// 	other_rsd() const
// 	{
// 		return other_rsd_;
// 	}

// 	///
// 	std::string const &
// 	other_atom_name() const
// 	{
// 		return other_atom_name_;
// 	}

