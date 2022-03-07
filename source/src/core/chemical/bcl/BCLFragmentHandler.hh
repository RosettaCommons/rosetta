// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
/// @file
///
/// @brief
/// A class for converting between BCl and Rosetta small molecule objects
///
/// @details
/// This class converts between BCL and Rosetta small molecule objects and sets
/// small molecule BCL conformers as rotamers in a Rosetta pose residue.
/// This class is meant to deprecate FragmentToRestype and RestypeToFragment
/// and make one umbrella class that can be extended and/or generalized.
///
/// @author
/// Sandeep Kothiwale
/// Rocco Moretti
/// Steven Combs
/// Benjamin P. Brown (benjamin.p.brown17@gmail.com)
////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_bcl_BCLFragmentHandler_hh
#define INCLUDED_core_chemical_bcl_BCLFragmentHandler_hh

// core includes
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/StoredRotamerLibrarySpecification.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/sdf/bcl_sdf.h>
#include <bcl/include/sdf/bcl_sdf_bond_info.h>
#include <bcl/include/sdf/bcl_sdf_atom_info.h>
#include <bcl/include/chemistry/bcl_chemistry_atoms_complete_standardizer.h>
#include <bcl/include/chemistry/bcl_chemistry_bond_isometry_handler.h>
#include <bcl/include/chemistry/bcl_chemistry_stereocenters_handler.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_ensemble.h>
#include <bcl/include/chemistry/bcl_chemistry_rotamer_library_file.h>
#include <bcl/include/chemistry/bcl_chemistry_sample_conformations.h>
#include <bcl/include/util/bcl_util_loggers.h>
#endif

// utility includes
#include <utility/VirtualBase.hh>
#include <utility/vector0.hh>
#include <utility/exit.hh>
#include <string>

namespace core {
namespace chemical {
namespace bcl {

class BCLFragmentHandler : public utility::VirtualBase
{

public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief default constructor
	BCLFragmentHandler();

	/// @brief construct with BCL fragment
#ifdef USEBCL
	explicit BCLFragmentHandler( ::bcl::chemistry::FragmentComplete const &fragment);
#endif

	/// @brief construct with Rosetta residue
#ifdef USEBCL
	explicit BCLFragmentHandler( MutableResidueTypeCOP restype);
#endif

	/////////////////
	// data access //
	/////////////////

	//! @brief Return the BCL molecule as a standard FragmentComplete
#ifdef USEBCL
	::bcl::chemistry::FragmentComplete const &get_bcl_fragment() const;
#endif

	//! @brief Return a Rosetta molecule as a MutableResidueType
	MutableResidueTypeCOP get_rosetta_restype() const;

	//! @brief Get how the most recently created ResidueType corresponds to the underlying fragment.
	IndexVDMapping const &get_index_to_vd() const;

	//! @brief Get mapping of restype vertex descriptors to indices of the bcl fragments
	VDIndexMapping const &get_vd_to_index() const;

	//! @brief Get the atom used as the neighbor atom when restype is generated
	core::Size get_nbr() const;


	////////////////
	// operations //
	////////////////

	//! @brief Set the BCL molecule
#ifdef USEBCL
	void set_bcl_fragment( ::bcl::chemistry::FragmentComplete const &fragment);
#endif

	//! @brief Set the Rosetta molecule
	void set_rosetta_restype( MutableResidueTypeCOP restype);

	/// @brief Which atom in the fragment to use as the neighbor atom when the a restype is generated.
	void set_nbr( core::Size nbr );

#ifdef USEBCL
	//! @brief Convert BCL member Fragment to Rosetta MutableResidueType
	//! @return mutable residue type object for Rosetta
	MutableResidueTypeOP fragment_to_restype();
#endif

#ifdef USEBCL
	//! @brief Convert BCL Fragment to Rosetta MutableResidueType
	//! @param FRAGMENT the BCL FragmentComplete to be converted
	//! @param CONFS fragment conformers to save as rotamers
	//! @return mutable residue type object for Rosetta
	MutableResidueTypeOP fragment_to_restype
	(
		::bcl::chemistry::FragmentComplete const &fragment,
		::bcl::chemistry::FragmentEnsemble const &confs = ::bcl::chemistry::FragmentEnsemble()
	);
#endif

#ifdef USEBCL
	//! @brief Convert BCL Fragment to Rosetta MutableResidueType
	//! @param FRAGMENT the BCL FragmentComplete to be converted
	//! @param RESTYPE the the restype from which to obtain additional restype info
	//! @param MAPPING the VD index mapping object
	//! @param CONFS fragment conformers to save as rotamers
	//! @return mutable residue type object for Rosetta
	MutableResidueTypeOP fragment_to_restype
	(
		::bcl::chemistry::FragmentComplete const &fragment,
		MutableResidueType const &restype,
		VDIndexMapping const &mapping,
		::bcl::chemistry::FragmentEnsemble const &confs = ::bcl::chemistry::FragmentEnsemble()
	);
#endif

#ifdef USEBCL
	//! @brief Convert ResidueType to a BCL Fragment.
	//! @param RESTYPE the restype to be converted into a FragmentComplete
	//! @param NORMALIZE if true, remove charges and reprotonate based on BCL heuristics.
	//! Note, charges due to heavy atoms, like quaternary amines, should still be present.
	::bcl::chemistry::FragmentComplete restype_to_fragment
	(
		MutableResidueType const &restype,
		bool normalize = false
	);
#endif

private:

	//////////
	// data //
	//////////

#ifdef USEBCL
	//! @brief molecule in BCL format
	::bcl::chemistry::FragmentComplete bcl_fragment_;
#endif

	//! @brief molecule in Rosetta format
	MutableResidueTypeCOP rosetta_restype_;

	//! @brief How the fragment indices map to the most recently created ResidueType
	IndexVDMapping index_to_vd_;

	//! @brief Mapping of restype vertex descriptors to indices of the bcl fragments
	//! Uses utility::get_undefined_size() for non-represented BCL indices
	VDIndexMapping vd_to_index_;

	//! @brief Which index in the fragment is used for the neighbor atom.
	//! utility::get_undefined_size() means autodetermine.
	core::Size nbr_;

};
} // bcl
} // chemical
} // core


#endif // INCLUDED_core_chemical_bcl_BCLFragmentHandler_hh
