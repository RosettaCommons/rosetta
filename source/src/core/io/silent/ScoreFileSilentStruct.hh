// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/ScoreFileSilentStruct.hh
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_ScoreFileSilentStruct_hh
#define INCLUDED_core_io_silent_ScoreFileSilentStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/io/silent/SilentStruct.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>


// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {

class ScoreFileSilentStruct : public SilentStruct {

public:

	/// @brief Constructors.
	ScoreFileSilentStruct( SilentFileOptions const & opts );

	ScoreFileSilentStruct(
		SilentFileOptions const & opts,
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	);

	SilentStructOP clone() const override {
		return utility::pointer::make_shared< ScoreFileSilentStruct >( *this );
	}

	/// @brief Re-dimension the storage capacity of this ScoreFileSilentStruct to the given number of
	/// residues.
	//void resize(Size const nres_in );

	// destructor
	~ScoreFileSilentStruct() override {}

	/// @brief Test if this ScoreFileSilentStruct is equal to the given ScoreFileSilentStruct in
	/// terms of conformation. Doesn't check energies.
	ScoreFileSilentStruct & operator= (
		ScoreFileSilentStruct const & src
	);

	/// @brief Tells this ScoreFileSilentStruct object to initialize itself from the given set of lines.
	/// Only initializes energies.
	bool init_from_lines(
		utility::vector1< std::string > const & lines,
		SilentFileData & container
	) override;

	/// @brief Configure the conformation of the given Pose with the
	/// conformational data within this ScoreFileSilentStruct. Calls pose.clear() and
	/// rebuilds Pose from scratch using the / user-specified residue types.
	void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set,
		bool const metapatches = true
	) const override;

	/// @brief opposite of fill_pose
	void fill_struct(
		core::pose::Pose const & pose,
		std::string tag
	) override;

	/// @brief Prints the conformation information within this ScoreFileSilentStruct
	/// to the given std::ostream.
	void print_conformation( std::ostream & output ) const override;

	/// @brief Prints the header information within this ScoreFileSilentStruct
	/// to the given std::ostream.
	void print_header( std::ostream & out ) const override;

	/// @brief returns the positions of the CA atoms in this ScoreFileSilentStruct.
	/// Useful for RMS calculations.
	ObjexxFCL::FArray2D< Real > get_CA_xyz() const override;

	/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from the torsions in this
	/// ScoreFileSilentStruct and the C-alpha atoms from this ScoreFileSilentStruct.
	Real get_debug_rmsd() override;

}; // class ScoreFileSilentStruct

} // namespace silent
} // namespace io
} // namespace core

#endif
