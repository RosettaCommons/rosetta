// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/rotamers/NCAARotamerLibrarySpecification.hh
/// @brief  The NCAARotamerLibrarySpecification class says to build a rotamer library from a NCAA rotlib
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_rotamers_NCAARotamerLibrarySpecification_HH
#define INCLUDED_core_chemical_rotamers_NCAARotamerLibrarySpecification_HH

// Unit headers
#include <core/chemical/rotamers/NCAARotamerLibrarySpecification.fwd.hh>
#include <core/chemical/rotamers/RotamerLibrarySpecification.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/types.hh>

// Basic headers

// Utility Headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <istream>

namespace core {
namespace chemical {
namespace rotamers {

class NCAARotamerLibrarySpecification : public RotamerLibrarySpecification {
public:
	NCAARotamerLibrarySpecification();
	NCAARotamerLibrarySpecification( std::string const & ncaa_rotlib_path );
	NCAARotamerLibrarySpecification( std::istream & );
	virtual ~NCAARotamerLibrarySpecification();

	/// @brief Sets the path to the NCAA rotlib for the residue type
	void
	ncaa_rotlib_path( std::string const & path ) { ncaa_rotlib_path_ = path; }

	// /// @brief Sets the number of rotatable bonds described by the NCAA rotlib  (not nesesarily equal to nchi)
	// void
	// ncaa_rotlib_n_rotameric_bins( core::Size nbins ) { ncaa_rotlib_n_rots_ = nbins; }

	/// @brief Sets the number of rotamers for each rotatable bond described by the NCAA rotlib
	void
	ncaa_rotlib_n_bin_per_rot( utility::vector1<Size> const & n_bins_per_rot) { ncaa_rotlib_n_bins_per_rot_ = n_bins_per_rot; }

	void
	semirotameric_ncaa_rotlib(bool setting) { semirotameric_ncaa_rotlib_ = setting; }

	void
	nrchi_symmetric(bool setting) { nrchi_symmetric_ = setting; }

	void
	nrchi_start_angle(Real setting) { nrchi_start_angle_ = setting; }

	/// @brief Returns the path to the NCAA rotlib for the residue type
	std::string const &
	ncaa_rotlib_path() const { return ncaa_rotlib_path_; }

	/// @brief Returns the number of rotatable bonds described by the NCAA rotlib  (not nesesarily equal to nchi)
	core::Size
	ncaa_rotlib_n_rotameric_bins() const { return ncaa_rotlib_n_bins_per_rot_.size(); }

	/// @brief Returns the number of rotamers for each rotatable bond described by the NCAA rotlib for all bonds
	utility::vector1<Size> const &
	ncaa_rotlib_n_bin_per_rot() const { return ncaa_rotlib_n_bins_per_rot_; }

	bool
	semirotameric_ncaa_rotlib() const { return semirotameric_ncaa_rotlib_; }

	bool
	nrchi_symmetric() const { return nrchi_symmetric_; }

	Real
	nrchi_start_angle() const { return nrchi_start_angle_; }

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

	/// @brief path to the NCAA rotlib
	std::string ncaa_rotlib_path_;
	// /// @brief the number of non-hydrogen chi angles in the NCAA rotlib
	// core::Size ncaa_rotlib_n_rots_;
	/// @brief the number of rotamer bins for each chi angle in the NCAA rotlib
	utility::vector1< Size > ncaa_rotlib_n_bins_per_rot_;

	/// @brief whether or not the NCAA rotlib is semirotameric
	bool semirotameric_ncaa_rotlib_;

	bool nrchi_symmetric_;
	Real nrchi_start_angle_;

};


} //namespace rotamers
} //namespace chemical
} //namespace core


#endif
