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

// Unit Headers
#include <core/conformation/PseudoBond.hh>

namespace core {
namespace conformation {

PseudoBond::PseudoBond() :
	ReferenceCount(),
	lr_conn_(),
	ur_conn_(),
	nbonds_( 0 )
{}

PseudoBond::~PseudoBond() {}

PseudoBond::PseudoBond( PseudoBond const & src ) :
	ReferenceCount(),
	lr_conn_( src.lr_conn_ ),
	ur_conn_( src.ur_conn_ ),
	nbonds_( src.nbonds_ )
{}

PseudoBond const & PseudoBond::operator = ( PseudoBond const & rhs )
{
	lr_conn_ = rhs.lr_conn_;
	ur_conn_ = rhs.ur_conn_;
	nbonds_  = rhs.nbonds_;
	return *this;
}

bool PseudoBond::operator == ( PseudoBond const & rhs ) const
{
	if ( lr_conn_ != rhs.lr_conn_ ) return false;
	if ( ur_conn_ != rhs.ur_conn_ ) return false;
	if ( nbonds_  != rhs.nbonds_ ) return false;
	return true;
}

// lower residue
Size PseudoBond::lr() const { return lr_conn_.resid(); }
void PseudoBond::lr( Size setting ) { lr_conn_.resid( setting ); }

// upper residue
Size PseudoBond::ur() const { return ur_conn_.resid(); }
void PseudoBond::ur( Size setting ) { ur_conn_.resid( setting ); }

Size PseudoBond::lr_conn_id() const { return lr_conn_.connid(); }
void PseudoBond::lr_conn_id( Size setting ) { lr_conn_.connid( setting ); }

Size PseudoBond::ur_conn_id() const { return ur_conn_.connid(); }
void PseudoBond::ur_conn_id( Size setting ) { ur_conn_.connid( setting ); }

PseudoBond::ResConnID PseudoBond::lr_resconnid() const { return lr_conn_; }
void      PseudoBond::lr_resconnid( ResConnID setting ) { lr_conn_ = setting; }

PseudoBond::ResConnID PseudoBond::ur_resconnid() const { return ur_conn_; }
void      PseudoBond::ur_resconnid( ResConnID setting ) { ur_conn_ = setting; }

Size PseudoBond::nbonds() const { return nbonds_; }
void PseudoBond::nbonds( Size setting ) { nbonds_ = setting; }


// A PBCollection stores all of the PBs between a pair of residues.
// PBs can be added to the collection, and iterated over, but cannot
// be modified.
PseudoBondCollection::PseudoBondCollection() {}
PseudoBondCollection::~PseudoBondCollection() {}

PseudoBondCollectionCOP
PseudoBondCollection::clone_with_new_sequence_numbering(
	utility::vector1< int > const & old2new
) const
{
	PseudoBondCollectionOP new_pbc( new PseudoBondCollection );
	for ( PBIter pbiter = iter_begin(), pbiter_end = iter_end();
			pbiter != pbiter_end; ++pbiter ) {
		PseudoBond new_pb( *pbiter );
		new_pb.lr( old2new[ pbiter->lr() ] );
		new_pb.ur( old2new[ pbiter->ur() ] );
		new_pbc->push_back( new_pb );
	}
	return new_pbc;
}

void PseudoBondCollection::push_back( PseudoBond const & pb )  {
	pseudo_bonds_.push_back( pb );
}

PseudoBondCollection::PBIter PseudoBondCollection::iter_begin() const {
	return pseudo_bonds_.begin();
}

PseudoBondCollection::PBIter PseudoBondCollection::iter_end() const {
	return pseudo_bonds_.end();
}


} // conformation
} // core
