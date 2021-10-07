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


#ifndef INCLUDED_protocols_drug_design_bcl_BCLSampleConfsManager_hh
#define INCLUDED_protocols_drug_design_bcl_BCLSampleConfsManager_hh

// Unit headers

// Protocols headers

// Core headers

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// C++ header
#include <map>
#include <tuple>

// BCL includes
#ifdef USEBCL
#include <bcl/include/chemistry/bcl_chemistry_fragment_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_mutate_interface.h>
#include <bcl/include/chemistry/bcl_chemistry_rotamer_library_file.h>
#include <bcl/include/chemistry/bcl_chemistry_sample_conformations.h>
#include <bcl/include/math/bcl_math_mutate_interface.h>
#include <bcl/include/sched/bcl_sched_mutex.h>
#include <bcl/include/sdf/bcl_sdf.h>
#include <bcl/include/util/bcl_util_implementation.h>
#include <bcl/include/util/bcl_util_loggers.h>
#endif

namespace protocols
{
namespace drug_design
{
namespace bcl
{

/// @brief A singleton class for managing supplementary SDFs to ensure that they are loaded once and only once from disk.
/// @author Benjamin P. Brown (benjamin.p.brown17@gmail.com)
class BCLSampleConfsManager : public utility::SingletonBase< BCLSampleConfsManager > {
	friend class utility::SingletonBase< BCLSampleConfsManager >;

public:

	////////////////
	// operations //
	////////////////

#ifdef USEBCL
	/// @brief initialize the rotamer library file
	void init_rotamer_library_file() const;

	/// @brief return the reference fragment from a filename if it exists
	::bcl::chemistry::RotamerLibraryFile get_rotamer_library_file() const;

	/// @brief return the reference fragment from a filename if it exists
	::bcl::chemistry::SampleConformations get_sample_confs() const;

	/// @brief set a sample conformations object to member
	void set_sample_confs( ::bcl::chemistry::SampleConformations const &sample_confs) const;
#endif

private:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief empty constructor.
	BCLSampleConfsManager();

	/// @brief explicitly deleted copy constructor.
	BCLSampleConfsManager( BCLSampleConfsManager const & ) = delete;

	/// @brief explicitly deleted assignment operator.
	BCLSampleConfsManager operator=( BCLSampleConfsManager const & ) = delete;

	//////////////////////
	// helper functions //
	//////////////////////

#ifdef USEBCL
	/// @brief mutex for thread safety
	static ::bcl::sched::Mutex &GetMutex();
#endif

	//////////
	// data //
	//////////

#ifdef USEBCL
	//! @brief rotamer library object for small molecule conformation sampling
	mutable ::bcl::chemistry::RotamerLibraryFile rotlib_;

	//! @brief SampleConformations object for small molecule conformation sampling
	mutable ::bcl::chemistry::SampleConformations conf_sampler_;
#endif


};

} // bcl
} // protocols
} // drug_design

#endif //INCLUDED_protocols/relax_BCLSampleConfsManager_fwd_hh



