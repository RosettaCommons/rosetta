// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemica/ResConnID.hh
/// @brief  Declaration of class to represent the chemical bond between two
///         ResidueConnections from one conformation::Residue to another.
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_chemical_ResConnID_hh
#define INCLUDED_core_chemical_ResConnID_hh

// Unit Headers
#include <core/chemical/ResConnID.fwd.hh>

// Project Headers
#include <core/types.hh>

namespace core {
namespace chemical {


/// @brief The ResConnID could more properly be called the ResidueConnector.  It stores the
/// data necessary to describe how one ResidueConnection on a conformation::Residue is connected
/// to the rest of the structure (Pose), by listing the other Residue's index and the ResidueConnection
/// index.
class ResConnID {

public:

	ResConnID();
	ResConnID( ResConnID const & );
	ResConnID( Size resid, Size connid );

	ResConnID const & operator = ( ResConnID const & );
	friend bool operator < ( ResConnID const & lhs, ResConnID const & rhs );
	friend bool operator == ( ResConnID const & lhs, ResConnID const & rhs );
	friend bool operator != ( ResConnID const & lhs, ResConnID const & rhs );

	Size resid() const;
	void resid( Size );

	Size connid() const;
	void connid( Size );

	bool incomplete() const;
	void mark_incomplete();

private:
	Size res_id_;
	Size conn_id_;

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) {
			ar & res_id_;
			ar & conn_id_;
	}
#endif
};

bool operator < ( ResConnID const & lhs, ResConnID const & rhs );
bool operator == ( ResConnID const & lhs, ResConnID const & rhs );
bool operator != ( ResConnID const & lhs, ResConnID const & rhs );

}
}

#endif
