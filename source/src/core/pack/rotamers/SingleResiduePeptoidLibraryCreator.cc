// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamers/SingleResiduePeptoidLibraryCreator.hh
/// @brief  Class for instantiating a particular SingleResidueRotamerLibrary
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <core/pack/rotamers/SingleResiduePeptoidLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>

// Program header
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.hh>
#include <core/pack/rotamers/RotamericSingleResiduePeptoidLibrary.tmpl.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>
#include <core/chemical/rotamers/PeptoidRotamerLibrarySpecification.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// External headers
#include <boost/lexical_cast.hpp>

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace rotamers {

static thread_local basic::Tracer TR("core.pack.rotamers.SingleResiduePeptoidLibraryCreator");

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleResiduePeptoidLibraryCreator::create( core::chemical::ResidueType const & restype) const {
	using namespace core::chemical::rotamers;
	using namespace core::pack::dunbrack;

	RotamerLibrarySpecificationCOP libspec( restype.rotamer_library_specification() );
	// If the factory system is sound, these two checks should work.
	assert( libspec );
	PeptoidRotamerLibrarySpecificationCOP peptoid_libspec = utility::pointer::dynamic_pointer_cast< PeptoidRotamerLibrarySpecification const >(libspec);
	assert( peptoid_libspec );

	// Some basic error checking against restype
	Size pose_n_rotlib_chi( restype.nchi() - restype.n_proton_chi() );
	Size n_rotlib_chi( peptoid_libspec->peptoid_rotlib_n_rotameric_bins() );
	if( n_rotlib_chi != pose_n_rotlib_chi ) {
		TR.Error << "Number of chi mismatch. Expected " << n_rotlib_chi << " rotatable heavy atom chis, found " << pose_n_rotlib_chi << std::endl;
		utility_exit_with_message("Number of chi mismatch in Peptoid rotlib loading.");
	}

	// Don't use any restype info in actually making the SingleResidueRotamerLibraryCOP, to allow for robust caching.

	using namespace utility::options;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// create izstream from path
	std::string dir_name = basic::database::full_name( "/rotamer/peptoid_rotlibs/" );
	std::string file_name = peptoid_libspec->peptoid_rotlib_path();
	if( ! file_name.size() ) {
		utility_exit_with_message("Unspecified Peptoid rotlib path with residue type " + restype.name() );
	}
	std::string full_path = dir_name + file_name;
	utility::io::izstream rotlib_in( full_path );

	// if we cannot open in regular path, try alternate paths
	utility::options::PathVectorOption & pvec = option[ in::file::extra_rot_lib_path ];
	Size pveci(1);
	while ( !rotlib_in && pveci <= pvec.size() ) {
		full_path = pvec[ pveci ].name() + file_name;
		rotlib_in.open( full_path );
		pveci++;
	}

	TR << "Reading in rot lib " << full_path << "...";

	// get an instance of RotamericSingleResiduePeptoidLibrary, but need a RotamerLibrary to do it
	// this means that when ever you read in the Peptoid libraries you will also read in the Peptoid libraries
	// this may need to be a pointer to the full type and not just a SRRLOP

	// this comes almost directally from RotmerLibrary.cc::create_rotameric_dunlib()
	SingleResidueRotamerLibraryOP peptoid_rotlib;

	switch ( n_rotlib_chi ) {
	case 1: {
		RotamericSingleResiduePeptoidLibrary< ONE, THREE > * r1 =
			new RotamericSingleResiduePeptoidLibrary< ONE, THREE >();
		r1->set_n_chi_bins( peptoid_libspec->peptoid_rotlib_n_bin_per_rot() );
		r1->read_from_file( rotlib_in );
		peptoid_rotlib = SingleResidueRotamerLibraryOP(r1);
		break;
	}
	case 2: {
		RotamericSingleResiduePeptoidLibrary< TWO, THREE > * r2 =
			new RotamericSingleResiduePeptoidLibrary< TWO, THREE >();
		r2->set_n_chi_bins( peptoid_libspec->peptoid_rotlib_n_bin_per_rot() );
		r2->read_from_file( rotlib_in );
		peptoid_rotlib = SingleResidueRotamerLibraryOP(r2);
		break;
	}
	case 3: {
		RotamericSingleResiduePeptoidLibrary< THREE, THREE > * r3 =
			new RotamericSingleResiduePeptoidLibrary< THREE, THREE >();
		r3->set_n_chi_bins( peptoid_libspec->peptoid_rotlib_n_bin_per_rot() );
		r3->read_from_file( rotlib_in );
		peptoid_rotlib = SingleResidueRotamerLibraryOP(r3);
		break;
	}
	case 4: {
		RotamericSingleResiduePeptoidLibrary< FOUR, THREE > * r4 =
			new RotamericSingleResiduePeptoidLibrary< FOUR, THREE >();
		r4->set_n_chi_bins( peptoid_libspec->peptoid_rotlib_n_bin_per_rot() );
		r4->read_from_file( rotlib_in );
		peptoid_rotlib = SingleResidueRotamerLibraryOP(r4);
		break;
	}
	default:
		TR.Error << "ERROR: too many chi angles desired for peptoid library: " << n_rotlib_chi << std::endl;
		utility_exit_with_message( "ERROR: too many chi angles desired for peptoid library." );
		break;
	}

	TR << "done!" << std::endl;
	return peptoid_rotlib;
}

std::string
SingleResiduePeptoidLibraryCreator::keyname() const {
	return core::chemical::rotamers::PeptoidRotamerLibrarySpecification::library_name();
}

} //namespace rotamers
} //namespace pack
} //namespace core

