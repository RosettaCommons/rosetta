// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////////////////////////////
/// @file protocols/drug_design/bcl/BCLSampleConfsManager
///
/// @brief A singleton class for managing BCL SampleConformations objects
///
/// @details
/// This class is an experimental design choice for managing SampleConformations and related BCL class objects.
/// Specifically, it is usually useful for memory management purposes to avoid loading a rotamer library and
/// sample conformations object more than once per application. Periodically we do need multiple unique sample conformations
/// objects, but we still generally just use the same rotamer library.
///
/// @author Benjamin P. Brown (benjamin.p.brown17@gmail.com)
////////////////////////////////////////////////////////////////////////

// Unit headers
#include <protocols/drug_design/bcl/BCLSampleConfsManager.hh>

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
static basic::Tracer TR( "protocols.drug_design.bcl.BCLSampleConfsManager" );

namespace protocols
{
namespace drug_design
{
namespace bcl
{

// BCLSampleConfsManager Public methods /////////////////////////////////////////////////////////////
// Static constant data access

/// @brief initialize the rotamer library file
#ifdef USEBCL
void BCLSampleConfsManager::init_rotamer_library_file() const
{
	rotlib_ = ::bcl::chemistry::RotamerLibraryFile();
}
#endif

/// @brief return the rotamer library
#ifdef USEBCL
::bcl::chemistry::RotamerLibraryFile BCLSampleConfsManager::get_rotamer_library_file() const
{
return rotlib_;
}
#endif

/// @brief return the sample conformations object
#ifdef USEBCL
::bcl::chemistry::SampleConformations BCLSampleConfsManager::get_sample_confs() const
{
return conf_sampler_;
}
#endif

/// @brief set a sample conformations object to member
#ifdef USEBCL
void BCLSampleConfsManager::set_sample_confs( ::bcl::chemistry::SampleConformations const &sample_confs) const
{
	conf_sampler_ = sample_confs;
}
#endif

// BCLSampleConfsManager Private methods ////////////////////////////////////////////////////////////

/// @brief empty constructor
BCLSampleConfsManager::BCLSampleConfsManager() :
	SingletonBase< BCLSampleConfsManager >()
#ifdef USEBCL
	,rotlib_(),
	conf_sampler_()
#endif
{}

#ifdef USEBCL
/// @brief mutex for thread safety
::bcl::sched::Mutex &BCLSampleConfsManager::GetMutex()
{
	static ::bcl::sched::Mutex mutex;
	return mutex;
}
#endif


} // bcl
} // protocols
} // drug_design
