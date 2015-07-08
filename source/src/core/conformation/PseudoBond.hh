// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_conformation_PseudoBond_hh
#define INCLUDED_core_conformation_PseudoBond_hh

// Unit Headers
#include <core/conformation/PseudoBond.fwd.hh>

// Project Headers
#include <core/chemical/ResConnID.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1_bool.hh>


namespace core {
namespace conformation {

// A PseudoBond connects two residues that are not covalently attached, but are
// within four bonds of each other through some third (or third & fourth) residue.
// PseudoBonds exist in the conformation layer and not in the chemical layer
// because they are the result of specific chemical configurations and cannot
// be predicted until runtime when a PDB is read in.  A pseudobond may also connect
// a single residue back to itself.
//
// PseudoBonds keep indices to connection-ids for the two residues that they
// bridge; a psuedobond attaches between two atoms that are already designated
// as connection points for their ResidueTypes.  The Residuetypes provide mappings
// between connection point ids and the atom ids that correspond to those connections.
// For two different residue types to utilize the same pseudo-bond information, they
// must have an identical set of connection points.  This has a counter-intuitive
// result: if one wants to flip the chiral center on a residue with a three-atom
// backbone, one cannot simply twist the residue 180 degrees; the connection points
// won't correspond.  The result is that different residue types have to be defined for
// r vs s chiral centers.

class PseudoBond : public utility::pointer::ReferenceCount
{
public:

typedef chemical::ResConnID ResConnID;

public:
	PseudoBond();
	virtual ~PseudoBond();
	PseudoBond( PseudoBond const &);
	PseudoBond & operator = ( PseudoBond const & rhs );

	bool operator == ( PseudoBond const & rhs ) const;

	// lower residue
	Size lr() const;
	void lr( Size );

	// upper residue
	Size ur() const;
	void ur( Size );

	Size lr_conn_id() const;
	void lr_conn_id( Size );

	Size ur_conn_id() const;
	void ur_conn_id( Size );


	ResConnID lr_resconnid() const;
	void lr_resconnid( ResConnID );

	ResConnID ur_resconnid() const;
	void ur_resconnid( ResConnID );

	Size nbonds() const;
	void nbonds( Size );

private:
	ResConnID lr_conn_;
	ResConnID ur_conn_;

	Size nbonds_;

};


// A PBCollection stores all of the PBs between a pair of residues.
// PBs can be added to the collection, and iterated over, but cannot
// be modified.

class PseudoBondCollection : public utility::pointer::ReferenceCount
{
public:
	typedef utility::vector1< PseudoBond >::const_iterator PBIter;

public:
	PseudoBondCollection();
	virtual ~PseudoBondCollection();

	PseudoBondCollectionCOP
	clone_with_new_sequence_numbering(
		utility::vector1< int > const & old2new
	) const;

	void push_back( PseudoBond const & );

	PBIter iter_begin() const;
	PBIter iter_end() const;

	Size size() const { return pseudo_bonds_.size(); }

private:

	utility::vector1< PseudoBond > pseudo_bonds_;

};


} // conformation
} // core

#endif
