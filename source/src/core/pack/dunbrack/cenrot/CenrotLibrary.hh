// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/cenrot/CenrotLibrary.hh
/// @brief  Centroid Rotamer Library classes and utility functions
/// @author


#ifndef INCLUDED_core_pack_dunbrack_cenrot_CenrotLibrary_hh
#define INCLUDED_core_pack_dunbrack_cenrot_CenrotLibrary_hh

// Unit headers
#include <core/pack/dunbrack/cenrot/CenrotLibrary.fwd.hh>

// Package headers
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.fwd.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>
#endif

// Utility headers
#include <utility/SingletonBase.hh>

// Numeric headers

// C++ headers
#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

///////////////////////////////////////////////////////////////////////////////
/// @brief Stores and handles loading of centroid rotamers for the canonical amino acids.
class CenrotLibrary : public utility::SingletonBase< CenrotLibrary >
{
public:
	friend class utility::SingletonBase< CenrotLibrary >;

	typedef chemical::AA AA;
	typedef conformation::Residue Residue;
	typedef chemical::ResidueType ResidueType;

public:
	virtual ~CenrotLibrary();

	SingleResidueCenrotLibraryCOP
	get_cenrot_library_by_aa( chemical::AA const & aa ) const;

private:

	/// @brief Add a centroid canonical AA Dunbrack library
	///
	/// @details For thread safety, should only be called during singleton construction.
	/// If you make this public, you need to add a mutex
	void
	add_cenrot_residue_library(
		AA const & aa,
		SingleResidueCenrotLibraryCOP rot_lib
	);

	/// @brief Initialize library from the appropriate database file
	/// Called during singleton construction.
	void create_centroid_rotamer_libraries_from_ASCII();

private:

	CenrotLibrary();
	CenrotLibrary( CenrotLibrary const & ); // unimplemented
	CenrotLibrary const & operator = ( CenrotLibrary const & ); // unimplemented

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static CenrotLibrary * create_singleton_instance();

private:

	utility::vector1< SingleResidueCenrotLibraryCOP > cenrot_libraries_;

};


} // cenrot
} // dunbrack
} // pack
} // core


#endif // INCLUDED_core_pack_dunbrack_cenrot_CenrotLibrary_HH
