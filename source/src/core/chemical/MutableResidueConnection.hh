// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/MutableResidueConnection.hh
/// @brief  Inter-residue chemical bond connection point class declaration.
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_MutableResidueConnection_hh
#define INCLUDED_core_chemical_MutableResidueConnection_hh


// Unit headers
#include <core/chemical/MutableResidueConnection.fwd.hh>

#include <core/chemical/ResidueGraphTypes.hh>

// Project headers
#include <core/chemical/MutableICoorRecord.hh>

// Utility headers

// C++ headers

namespace core {
namespace chemical {

/// @brief A simple class marking atoms at inter-residue connections.
///
/// Each residue type specifies some number of positions at which it is expecting
/// to form a chemical bond with another residue. Think of them as ports: they are
/// parts of the residue where there are chemical bonds beyond the intra-residue
/// chemical bonds are expected -- places where they can be chemically linked
/// to the outside world.
/// The MutableResidueConnection class stores the information for available connections
/// for the MutableResidueType object.

class MutableResidueConnection {
public:
	/// @brief default constructor
	MutableResidueConnection():
		icoor_(),
		vertex_(nullptr)
	{}

	/// @brief constructor with atom index number
	MutableResidueConnection(
		VD const vertex
	):
		icoor_(),
		vertex_(vertex)
	{}

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
	MutableICoorRecord const &
	icoor() const
	{
		return icoor_;
	}

	/// @brief set atom's AtomICoor
	void
	icoor( MutableICoorRecord const & ic )
	{
		icoor_ = ic;
	}

	/// @brief Update the internal VDs based on the provide mapping
	void remap_atom_vds( std::map< VD, VD > const & old_to_new );

private:

	/// atom AtomICoor
	MutableICoorRecord icoor_;
	VD vertex_; /// Which atom on the MutableResidueType is the connection point

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // chemical
} // core


#endif  // INCLUDED_core_chemical_MutableResidueConnection_hh

