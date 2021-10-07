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


#ifndef INCLUDED_protocols_drug_design_bcl_BCLReferenceSDFilesManager_hh
#define INCLUDED_protocols_drug_design_bcl_BCLReferenceSDFilesManager_hh

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
class BCLReferenceSDFilesManager : public utility::SingletonBase< BCLReferenceSDFilesManager > {
	friend class utility::SingletonBase< BCLReferenceSDFilesManager >;

public:

	////////////////
	// operations //
	////////////////

	/// @brief return the reference fragment from a filename if it exists
#ifdef USEBCL
	::bcl::util::ShPtr< ::bcl::chemistry::FragmentComplete> get_fragment_from_file( std::string const & filename) const;
#endif

	//! @brief add a fragment to our collection of reference fragments
	void add_reference_fragment( std::string const &reference_fragment_filename);

private:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief empty constructor.
	BCLReferenceSDFilesManager();

	/// @brief explicitly deleted copy constructor.
	BCLReferenceSDFilesManager( BCLReferenceSDFilesManager const & ) = delete;

	/// @brief explicitly deleted assignment operator.
	BCLReferenceSDFilesManager operator=( BCLReferenceSDFilesManager const & ) = delete;

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

	//! @brief A map of filename to file contents where the file is an SDF and the contents are reference fragments
#ifdef USEBCL
	mutable ::bcl::storage::Map< std::string, ::bcl::util::ShPtr< ::bcl::chemistry::FragmentComplete> > reference_fragments_;
#endif


};

} // bcl
} // protocols
} // drug_design

#endif //INCLUDED_protocols/relax_BCLReferenceSDFilesManager_fwd_hh



