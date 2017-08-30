// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/ResidueInformation.hh
/// @brief  Class definition for ResidueInformation
/// @author Sergey Lyskov


#ifndef INCLUDED_core_io_ResidueInformation_HH
#define INCLUDED_core_io_ResidueInformation_HH

// Unit headers
#include <core/io/ResidueInformation.fwd.hh>
#include <core/io/AtomInformation.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <string>


namespace core {
namespace io {

// Intermediate format for easy construction of
// core::conformation::Residue objects from a StructFileRep.
/// @details Subset of data from "ATOM" lines that is shared by all atoms in a
/// residue.\n
/// This class is used for temporary storage of information during Pose building
/// and "unbuilding".
class ResidueInformation {
public:
	/// @brief default constructor to initialize all values
	ResidueInformation();
	/// @brief Initialize the resName/chainID/etc. data from an atom,
	/// but do not add the atom to the atoms_ array.
	ResidueInformation(AtomInformation const & ai);

	bool operator==( ResidueInformation const & that ) const;
	bool operator!=( ResidueInformation const & that ) const;

	std::string const & resName() const;
	char chainID()  const;
	int  resSeq()   const;
	char iCode()    const;
	int  terCount() const;
	std::string const & segmentID() const;

	std::string const & rosetta_resName() const;

	void resName(  std::string const & setting );
	void chainID(  char setting );
	void resSeq(   int setting );
	void iCode(    char setting );
	void terCount( int setting );
	void set_xyz( std::string const &atomname, Vector const &vect );
	void set_temp( std::string const &atomname, core::Real const &val );
	void segmentID( std::string const & setting );

	utility::vector1< AtomInformation > const & atoms() const;
	void append_atom( AtomInformation const & atoms );
	void rename_atom( std::string const & orig_name, std::string const & new_name );
	void append_atoms( ResidueInformation const & source );

	/// @brief A convenience accessor for the coordinates.
	std::map< std::string, Vector > const & xyz() const;
	std::map< std::string, double > const & temps() const; //< map of names to B-factors;  redundant but used a lot in reader

	/// @brief Returns a short, printable designation for this residue.
	std::string resid() const;


private:
	/// For now, all member names have the same names as fields in PDB standard.
	std::string resName_;
	char chainID_;
	int resSeq_;
	char iCode_;
	int terCount_; //< number of TER or END cards encountered prior to this
	utility::vector1< AtomInformation > atoms_;
	std::map< std::string, Vector > xyz_; //< map of names to coords;  redundant but used a lot in reader
	std::map< std::string, double > temps_; //< map of names to B-factors;  redundant but used a lot in reader

	std::string segmentID_;

	// Additional annotations

	/// @brief Some Residues have a different three letter code in Rosetta versus the PDB. Store that here.
	std::string rosetta_resName_;

};  // class ResidueInformation

} // namespace io
} // namespace core

#endif  // INCLUDED_core_io_ResidueInformation_HH
