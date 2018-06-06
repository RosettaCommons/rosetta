// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/PDBSilentStruct.hh
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author James Thompson

#ifndef INCLUDED_core_import_pose_PDBSilentStruct_hh
#define INCLUDED_core_import_pose_PDBSilentStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileOptions.fwd.hh>
#include <core/io/StructFileRep.hh>


#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <iosfwd>

#include <utility/vector1.fwd.hh>


namespace core {
namespace import_pose {

class PDBSilentStruct : public core::io::silent::SilentStruct {

public:

	/// @brief Constructors.
	PDBSilentStruct( core::io::silent::SilentFileOptions const & opts );

	PDBSilentStruct(
		core::io::silent::SilentFileOptions const & opts,
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	);

	core::io::silent::SilentStructOP clone() const override {
		return core::io::silent::SilentStructOP( new PDBSilentStruct( *this ) );
	}

	// destructor
	~PDBSilentStruct() override = default;

	using core::io::silent::SilentStruct::print_header;

	virtual void print_header( std::ostream & out );

	/// @brief Test if this PDBSilentStruct is equal to the given PDBSilentStruct in terms of conformation.
	/// Doesn't check energies.
	PDBSilentStruct & operator= (
		PDBSilentStruct const & src
	);

	/// @brief Tells this PDBSilentStruct object to initialize itself from the given set of lines. Lines should
	/// be of the format
	bool init_from_lines(
		utility::vector1< std::string > const & lines,
		core::io::silent::SilentFileData & container
	) override;

	/// @brief Configure the conformation of the given Pose with the conformational data within this PDBSilentStruct.
	/// Calls pose.clear() and rebuilds Pose from scratch using FA_STANDARD residue types.
	void fill_pose(
		core::pose::Pose & pose,
		bool const metapatches = true
	) const override;

	/// @brief Configure the conformation of the given Pose with the
	/// conformational data within this PDBSilentStruct. Calls pose.clear() and
	/// rebuilds Pose from scratch using the / user-specified residue types.
	void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set,
		bool const metapatches = true
	) const override;

	/// @brief opposite of fill_pose
	void fill_struct(
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	) override;

	/// @brief Prints the conformation information within this PDBSilentStruct
	/// to the given std::ostream.
	void print_conformation( std::ostream & output ) const override;

	/// @brief data getters/setters
	core::io::StructFileRep const
	struct_file_rep() const {
		return *sfr_;
	}

	void struct_file_rep( core::io::StructFileRep new_sfr ) {
		sfr_ = new_sfr.clone();
	}

	std::string get_pdb_lines() const {
		return pdb_lines_;
	}

	/// @brief returns the positions of the CA atoms in this PDBSilentStruct.
	/// Useful for RMS calculations.
	ObjexxFCL::FArray2D< Real > get_CA_xyz() const override;

	/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from the torsions in this
	/// PDBSilentStruct and the C-alpha atoms from this PDBSilentStruct.
	Real get_debug_rmsd() override;

protected:
	core::io::StructFileRepOP sfr_;
	std::string pdb_lines_; // a concatenated version of all the lines in the pdb file
}; // class PDBSilentStruct

} // namespace import_pose
} // namespace core

#endif
