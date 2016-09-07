// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Class to write kinemage-formatted output for Residue and Conformation
/// @file   core/conformation/ResidueKinWriter.hh
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_conformation_ResidueKinWriter_hh
#define INCLUDED_core_conformation_ResidueKinWriter_hh


// Unit headers
#include <core/conformation/ResidueKinWriter.fwd.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iosfwd>

namespace core {
namespace conformation {

void write_kinemage_header(
	std::ostream & ostr,
	Size const kin_number,
	std::string const & title,
	Vector const & ctr
);

class ResidueKinWriter : public utility::pointer::ReferenceCount
{
public:
	typedef core::Size   Size;

public:
	ResidueKinWriter();
	~ResidueKinWriter() override;

	/// @brief write the header for the kinemage to center on this residue
	void write_kin_header(
		std::ostream & ostr,
		core::conformation::Residue const & rsd,
		core::Size atom_to_center_on = 0, // keep this 0 to center the kinemage on the neighbor atom
		core::Size which_kin = 1 // multiple kineamges can be put in a single kinemage file.
	) const;

	/// @brief Write out the coordinates for a particular residue; the kinemage tag
	/// is assumed to have been writen already.
	void
	write_rsd_coords(
		std::ostream & ostr,
		core::conformation::Residue const & rsd,
		bool is_instance = false
	) const;

	void master(    std::string const & setting );
	void dominant(                 bool setting );
	void animate(                  bool setting );
	void group(                    bool setting );
	void write_virtual_atoms(      bool setting );

	/// @brief Calling this function with the setting "true", turns on polar, apolar, and backbone hydrogen writing.
	/// Calling this function with the setting "false", turns off polar, apolar, and backbone hydrogen writing.
	void write_hydrogens(          bool setting );

	void write_apolar_hydrogens(   bool setting );
	void write_polar_hydrogens(    bool setting );
	void write_backbone_hydrogens( bool setting );

private:
	std::string master_;
	bool dominant_;
	bool animate_;
	bool group_; // false for subgroup
	bool write_apolar_hydrogens_;
	bool write_polar_hydrogens_;
	bool write_backbone_hydrogens_;
	bool write_virtual_atoms_;
};

class ConformationKinWriter : public utility::pointer::ReferenceCount
{
public:

	ConformationKinWriter():
		write_virtual_atoms_(false)
	{}
	~ConformationKinWriter() override;

	/// @brief Write out the coordinates for an entire conformation; this includes
	/// inter-residue bonds that would be missed by the ResidueKinWriter.
	void
	write_coords(
		std::ostream & ostr,
		core::conformation::Conformation const & conf,
		bool is_instance = false
	) const;

	void write_virtual_atoms( bool setting );
	void master( std::string const & setting );

private:

	std::string master_;
	bool write_virtual_atoms_;

};


} // conformation
} // core


#endif
