// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////////////////////////////
/// @file protocols/drug_design/bcl/BCLReferenceSDFilesManager
///
/// @brief A singleton class for managing supplementary SDFs to ensure that they are loaded once and only once from disk.
///
/// @details
/// This class is an experimental design choice for tracking reference and scaffold molecules that are read in
/// from disk during small molecule design. Each of the BCLFragmentMutateMover movers is going to perform some very specific
/// type of perturbation/mutation to the chemical structure of your molecule. This means that in a standard drug design XML
/// file you may define 10-20 drug design movers, all of which act on a different subset of atoms or perform different
/// types of moves or have different filters. One of the ways in which you can specify which atoms are mutable
/// (i.e. subject to chemical modification/replacement) is by supplying a reference structure in SDF format in the object data label
/// with which the ::bcl::util::Implementation< ::bcl::chemistry::FragmentMutateInterface> object is constructed. The goal here
/// is to avoid loading into memory more reference structures than you need. So if you have 20 movers and they all use the same 3
/// reference fragments, you only want to load 3 reference fragments, not 20 with a bunch of redundancy.
///
/// @author Benjamin P. Brown (benjamin.p.brown17@gmail.com)
////////////////////////////////////////////////////////////////////////

// Unit headers
#include <protocols/drug_design/bcl/BCLReferenceSDFilesManager.hh>

// Project header
#include <core/types.hh>

// Utility headers
#include <utility/pointer/memory.hh>
#include <utility/thread/threadsafe_creation.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/GeneralFileManager.hh>

// Core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>
#include <fstream>

// BCL includes
#ifdef USEBCL
#include <bcl/include/chemistry/bcl_chemistry_atom_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_atom_vector.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_track_mutable_atoms.h>
#include <bcl/include/chemistry/bcl_chemistry_rotamer_library_file.h>
#include <bcl/include/io/bcl_io_ofstream.h>
#include <bcl/include/io/bcl_io_file.h>
#include <bcl/include/linal/bcl_linal_vector_3d.h>
#include <bcl/include/storage/bcl_storage_vector.h>
#include <bcl/include/util/bcl_util_format.h>
#endif

// Boost headers
#include <functional>

// Tracer
static basic::Tracer TR( "protocols.drug_design.bcl.BCLReferenceSDFilesManager" );

namespace protocols
{
namespace drug_design
{
namespace bcl
{

// BCLReferenceSDFilesManager Public methods /////////////////////////////////////////////////////////////
// Static constant data access

/// @brief return the reference fragment from a filename if it exists
#ifdef USEBCL
::bcl::util::ShPtr< ::bcl::chemistry::FragmentComplete> BCLReferenceSDFilesManager::get_fragment_from_file( std::string const & filename) const
{
	// initialize a null ptr
	::bcl::util::ShPtr< ::bcl::chemistry::FragmentComplete> frag;

	// if the file exists then return mapped fragment
	if ( reference_fragments_.Has( filename) ) {
		frag = reference_fragments_.Find( filename)->second;
	}
	return frag;
}
#endif

//! @brief add a fragment to our collection of reference fragments
void BCLReferenceSDFilesManager::add_reference_fragment(
#ifdef USEBCL
	std::string const &reference_fragment_filename)
#else
	std::string const & /* OBJECT_DATA_LABEL */ )
#endif
{
#ifdef USEBCL
	// read in reference filename as a FragmentComplete
	if ( reference_fragment_filename.size() ) {
		// read input
		::bcl::io::IFStream input;
		::bcl::io::File::MustOpenIFStream( input, reference_fragment_filename);
		::bcl::chemistry::FragmentEnsemble reference_mol;
		reference_mol.ReadMoreFromMdl( input);
		::bcl::io::File::CloseClearFStream( input);

		// note that we want each filename to be associated with just 1 molecule, so we assume only 1 molecule per SDF
		reference_fragments_.Insert( std::make_pair
			(
			reference_fragment_filename,
			::bcl::util::ShPtr< ::bcl::chemistry::FragmentComplete>( new ::bcl::chemistry::FragmentComplete( reference_mol.GetMolecules().FirstElement()) ) )
		);
	}
#else
	utility_exit_with_message("Use of extras=bcl build required.");
#endif
}

// BCLReferenceSDFilesManager Private methods ////////////////////////////////////////////////////////////

/// @brief empty constructor
BCLReferenceSDFilesManager::BCLReferenceSDFilesManager() :
	SingletonBase< BCLReferenceSDFilesManager >()
#ifdef USEBCL
	,reference_fragments_()
#endif
{}

#ifdef USEBCL
/// @brief mutex for thread safety
::bcl::sched::Mutex &BCLReferenceSDFilesManager::GetMutex()
{
	static ::bcl::sched::Mutex mutex;
	return mutex;
}
#endif


} // bcl
} // protocols
} // drug_design
