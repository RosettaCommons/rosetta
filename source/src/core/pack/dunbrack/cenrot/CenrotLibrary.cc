// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


// Unit headers
#include <core/pack/dunbrack/cenrot/CenrotLibrary.hh>

// Package headers
#include <core/pack/dunbrack/cenrot/SingleResidueCenrotLibrary.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Numeric headers

// Utility headers
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>

// External headers

// C++ Headers

// Boost Headers
#include <boost/bind.hpp>

using basic::T;

using basic::Error;
using basic::Warning;

// Singleton instance and mutex static data members
namespace utility {

using core::pack::dunbrack::cenrot::CenrotLibrary;

#if defined MULTI_THREADED
template <> std::mutex utility::SingletonBase< CenrotLibrary >::singleton_mutex_{};
template <> std::atomic< CenrotLibrary * > utility::SingletonBase< CenrotLibrary >::instance_( 0 );
#else
template <> CenrotLibrary * utility::SingletonBase< CenrotLibrary >::instance_( 0 );
#endif

}

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

static THREAD_LOCAL basic::Tracer TR("core.pack.dunbrack.cenrot.CenrotLibrary");

CenrotLibrary *
CenrotLibrary::create_singleton_instance()
{
	CenrotLibrary * rotamer_library = new CenrotLibrary();

	rotamer_library->create_centroid_rotamer_libraries_from_ASCII();

	return rotamer_library;
}

CenrotLibrary::CenrotLibrary():
	cenrot_libraries_( chemical::num_canonical_aas, 0 )
{}

CenrotLibrary::~CenrotLibrary()
{}


void
CenrotLibrary::add_cenrot_residue_library(
	AA const & aa,
	SingleResidueCenrotLibraryCOP rot_lib
)
{
	if ( aa > chemical::num_canonical_aas ) {
		TR.Error << "ERROR: Cannot add centroid Dunbrack rotamer library of type " << aa << " not a canonical amino acid." << std::endl;
		utility_exit_with_message("Cannot add a non-canonical centroid Dunbrack library.");
	}
	if ( cenrot_libraries_[ aa ] != 0 ) {
		TR.Error << "ERROR: Cannot add cenroid Dunbrack rotamer library of type " << aa << ": library already loaded." << std::endl;
		utility_exit_with_message("Can't add centroid rsd library twice");
	}
	cenrot_libraries_[ aa ] = rot_lib;
}

SingleResidueCenrotLibraryCOP
CenrotLibrary::get_cenrot_library_by_aa( chemical::AA const & aa ) const
{
	if ( (Size) aa <= cenrot_libraries_.size() ) {
		return cenrot_libraries_[ aa ];
	}
	TR.Error << "ERROR: Cannot get centroid Dunbrack rotamer library of type " << aa << ": not a canonical amino acid." << std::endl;
	utility_exit_with_message("Cannot get non-canonical centroid Dunbrack library.");
	return 0;
}


void CenrotLibrary::create_centroid_rotamer_libraries_from_ASCII()
{
	using namespace chemical;
	using namespace core::pack::dunbrack::cenrot;

	/// Now read in the cenrot library
	clock_t starttime = clock();
	utility::io::izstream libstream(basic::database::full_name("rotamer/cenrot_dunbrack.lib"));
	//std::cout << basic::database::full_name("rotamer/centroid_rotlibs") << std::endl;
	ResidueTypeSetCAP rsd_set=ChemicalManager::get_instance()->residue_type_set( "centroid_rot" );

	std::string nextaa;
	libstream >> nextaa;

	Size count_libraries_read( 0 );
	while ( nextaa != "" ) {
		chemical::AA aan = chemical::aa_from_name( nextaa );
		SingleResidueCenrotLibraryOP newlib( new SingleResidueCenrotLibrary(aan) );
		/// read the rotlib for current aa and save the name of the next one
		nextaa = newlib->read_from_file( libstream, true );
		++count_libraries_read;

		// We put the test here so that we eat the remaining portion of the block.
		if ( aan > chemical::num_canonical_aas ) {
			TR.Warning << "Skipping cenrot library for non-canonical amino acid " << nextaa << " (" << aan << ")" << std::endl;
			continue;
		}
		add_cenrot_residue_library( aan, newlib );
	}

	libstream.close();

	clock_t stoptime = clock();
	TR << "Cenrot library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from ASCII" << std::endl;
}

} // cenrot
} // dunbrack
} // namespace scoring
} // namespace core
