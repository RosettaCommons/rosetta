// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/PeptoidRotamerLibrarySpecification.hh
/// @brief  The PeptoidRotamerLibrarySpecification class says to build a rotamer library from a Peptoid rotlib
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_PeptoidRotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_PeptoidRotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/PeptoidRotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Package headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Basic headers

// Utility Headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class PeptoidRotamerLibrarySpecification : public RotamerLibrarySpecification {
public:
	PeptoidRotamerLibrarySpecification();
	PeptoidRotamerLibrarySpecification( std::string const & peptoid_rotlib_path );
	PeptoidRotamerLibrarySpecification( std::istream & input );
	virtual ~PeptoidRotamerLibrarySpecification();

	/// @brief Sets the path to the Peptoid rotlib for the residue type
	void
	peptoid_rotlib_path( std::string const & path ) { peptoid_rotlib_path_ = path; }

//	/// @brief Sets the number of rotatable bonds described by the Peptoid rotlib  (not nesesarily equal to nchi)
//	void
//	peptoid_rotlib_n_rotameric_bins( core::Size nbins ) { peptoid_rotlib_n_rots_ = nbins; }

	/// @brief Sets the number of rotamers for each rotatable bond described by the Peptoid rotlib
	void
	peptoid_rotlib_n_bin_per_rot( utility::vector1<Size> const & n_bins_per_rot) { peptoid_rotlib_n_bins_per_rot_ = n_bins_per_rot; }

	/// @brief Returns the path to the Peptoid rotlib for the residue type
	std::string const &
	peptoid_rotlib_path() const { return peptoid_rotlib_path_; }

	/// @brief Returns the number of rotatable bonds described by the Peptoid rotlib  (not nesesarily equal to nchi)
	core::Size
	peptoid_rotlib_n_rotameric_bins() const { return peptoid_rotlib_n_bins_per_rot_.size(); }

	/// @brief Returns the number of rotamers for each rotatable bond described by the Peptoid rotlib for all bonds
	utility::vector1<Size> const &
	peptoid_rotlib_n_bin_per_rot() const { return peptoid_rotlib_n_bins_per_rot_; }

	virtual
	std::string
	keyname() const;

	virtual
	std::string
	cache_tag(core::chemical::ResidueType const &) const;

	static
	std::string
	library_name();

private:

	/// @brief path to the Peptoid rotlib
	std::string peptoid_rotlib_path_;
//	/// @brief the number of non-hydrogen chi angles in the Peptoid rotlib
//	core::Size peptoid_rotlib_n_rots_;
	/// @brief the number of rotamer bins for each chi angle in the Peptoid rotlib
	utility::vector1< Size > peptoid_rotlib_n_bins_per_rot_;

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
