// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Class to write kinemage-formatted output for ResidueType
/// @file   core/chemical/ResidueTypeKinWriter.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_ResidueTypeKinWriter_hh
#define INCLUDED_core_chemical_ResidueTypeKinWriter_hh


// Unit headers
#include <core/chemical/ResidueTypeKinWriter.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iosfwd>

namespace core {
namespace chemical {

class ResidueTypeKinWriter : public utility::pointer::ReferenceCount
{
public:
	ResidueTypeKinWriter();
	virtual ~ResidueTypeKinWriter();

	/// @brief write the header for the kinemage to center on this residue
	void write_kin_header(
		std::ostream & ostr,
		core::chemical::ResidueType const & restype,
		core::Size which_kin = 1 // multiple kineamges can be put in a single kinemage file.
	) const;

	/// @brief Write out settings for a particular ResidueType
	/// This should be similar to the kinemage files that molfile_to_params.py writes
	void
	write_restype(
		std::ostream & ostr,
		core::chemical::ResidueType const & restype,
		core::Size which_kin = 1 // multiple kineamges can be put in a single kinemage file.
	) const;

	// void master(    std::string const & setting );
	// void dominant(                 bool setting );
	// void animate(                  bool setting );
	// void group(                    bool setting );
	// void write_virtual_atoms(      bool setting );
	//
	// /// @brief Calling this function with the setting "true", turns on polar, apolar, and backbone hydrogen writing.
	// /// Calling this function with the setting "false", turns off polar, apolar, and backbone hydrogen writing.
	// void write_hydrogens(          bool setting );
	//
	// void write_apolar_hydrogens(   bool setting );
	// void write_polar_hydrogens(    bool setting );
	// void write_backbone_hydrogens( bool setting );

private:
	// std::string master_;
	// bool dominant_;
	// bool animate_;
	// bool group_; // false for subgroup
	// bool write_apolar_hydrogens_;
	// bool write_polar_hydrogens_;
	// bool write_backbone_hydrogens_;
	// bool write_virtual_atoms_;
};


} // chemical
} // core


#endif
