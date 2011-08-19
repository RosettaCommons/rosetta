// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/BinaryRNASilentStruct.hh
///
/// @brief
/// @author Frank DiMaio
/// @author Mike Tyka
/// @author Rhiju Das

#ifndef INCLUDED_core_io_silent_BinaryRNASilentStruct_hh
#define INCLUDED_core_io_silent_BinaryRNASilentStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// AUTO-REMOVED #include <core/conformation/Residue.fwd.hh>

// AUTO-REMOVED #include <core/io/silent/RNA_SilentStruct.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

// C++ Headers
#include <string>

//Auto Headers
#include <core/io/silent/SilentStruct.hh>


namespace core {
namespace io {
namespace silent {

class BinaryRNASilentStruct : public SilentStruct {

public:

	/// @brief Constructors.
	BinaryRNASilentStruct( Size const nres_in );

	BinaryRNASilentStruct();

	BinaryRNASilentStruct(
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	);

	/// @brief Re-dimension the storage capacity of this BinaryRNASilentStruct to the given number of residues.
	void resize(
		Size const nres_in
	);

	virtual SilentStructOP clone() const {
		return new BinaryRNASilentStruct( *this );
	};

	// destructor
	~BinaryRNASilentStruct() {}

	/// @brief Test if this BinaryRNASilentStruct is equal to the given BinaryRNASilentStruct in terms of conformation.
	/// Doesn't check energies.
	/* BinaryRNASilentStruct & operator= (
		BinaryRNASilentStruct const & src
	); */

	/// @brief Initialize object from a set of lines
	virtual bool init_from_lines(
		utility::vector1< std::string > const & lines,
		SilentFileData & container
	);

	/// @brief Configure the conformation of the given Pose with the conformational data within this BinaryRNASilentStruct.
	/// Calls pose.clear() and rebuilds Pose from scratch using FA_STANDARD residue types.
	virtual void fill_pose(
		core::pose::Pose & pose
	) const ;

	/// @brief Configure the conformation of the given Pose with the conformational data within
	/// this BinaryRNASilentStruct. Calls pose.clear() and rebuilds Pose from scratch using the
	/// user-specified residue types.
	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set
	) const;

	/// @brief opposite of fill_pose
	virtual void fill_struct( core::pose::Pose const & pose, std::string tag );

	/// @brief print header information
	virtual void print_header( std::ostream& out ) const;

	/// @brief Prints the conformation information within this BinaryRNASilentStruct to the given std::ostream.
	virtual void print_conformation( std::ostream & output ) const;

	/// @brief returns the positions of the CA atoms in this RNA_SilentStruct.
	/// Useful for RMS calculations.
	virtual ObjexxFCL::FArray2D< Real > get_CA_xyz() const;

	// model quality-related methods.
	virtual Real CA_rmsd( RNA_SilentStruct other_pss );

	/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from the torsions in this
	/// RNA_SilentStruct and the C-alpha atoms from this RNA_SilentStruct.
	virtual Real get_debug_rmsd();


protected:
	void add_jump( kinematics::Jump jump ) {
		jumps_.push_back( jump.rt() );
	}

	kinematics::RT const & jump( Size jump_num ) const {
		return jumps_[ jump_num ];
	}

	void set_fold_tree( kinematics::FoldTree const & f ) {
		fold_tree_ = f;
	}
	kinematics::FoldTree const& fold_tree( ) const {
		return fold_tree_;
	}

	char secstruct( unsigned int seqpos ) const {
		return secstruct_[seqpos];
	}

	void secstruct( unsigned int seqpos, char ss ) {
		secstruct_[seqpos] = ss;
	}

	bool fullatom_;

	utility::vector1< char > secstruct_;
	utility::vector1< utility::vector1< numeric::xyzVector<float> > > atm_coords_;
	utility::vector1< kinematics::RT > jumps_;
	bool bJumps_use_IntraResStub_;
	kinematics::FoldTree fold_tree_;

}; // class BinaryRNASilentStruct

} // namespace silent
} // namespace io
} // namespace core

#endif
