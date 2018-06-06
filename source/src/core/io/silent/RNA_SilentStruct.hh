// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/silent/RNA_SilentStruct.hh
///
/// @brief Representation of rosetta++ RNA silent-file structures.
/// @author James Thompson
/// @author Rhiju Das

#ifndef INCLUDED_core_io_silent_RNA_SilentStruct_hh
#define INCLUDED_core_io_silent_RNA_SilentStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


#include <core/io/silent/SilentStruct.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>


#include <utility/pointer/ReferenceCount.hh>


// C++ Headers
#include <cstdlib>
#include <iostream>
#include <utility/assert.hh>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include <utility/vector1.hh>


/////////////////////////////////////////////////////////////////////////
// Following should be easy to generalize for protein, RNA, DNA.
// This may eventually be critical as we start to look at mixed systems.
// For now, just for safety and (perhaps) to avoid confusion,
// we'll go ahead and make this a separate class.
//   -- Rhiju, April 2008
/////////////////////////////////////////////////////////////////////////

namespace core {
namespace io {
namespace silent {

class RNA_SilentStruct : public SilentStruct {

public:

	/// @brief Constructors.
	RNA_SilentStruct( SilentFileOptions const & opts, Size const nres_in )
	: SilentStruct( opts )
	{
		nres( nres_in );
		fullatom_    = true;
		non_main_chain_sugar_coords_defined_ = false;
		resize( nres_in );
	}

	// RNA_SilentStruct( RNA_SilentStruct const & src );

	RNA_SilentStruct( SilentFileOptions const & opts ) :
		SilentStruct( opts )
	{
		nres( 0 );
		fullatom_    = true;
		non_main_chain_sugar_coords_defined_ = false;
		decoy_tag( "empty" );
	}

	RNA_SilentStruct(
		SilentFileOptions const & opts,
		core::pose::Pose const & pose,
		std::string tag = "empty_tag",
		bool fa = true
	);

	/// @brief Re-dimension the storage capacity of this RNA_SilentStruct to the given number of residues.
	void resize(
		Size const nres_in
	);

	virtual SilentStructOP clone() const {
		return SilentStructOP( new RNA_SilentStruct( *this ) );
	};

	// destructor
	~RNA_SilentStruct() {}

	/// @brief Test if this RNA_SilentStruct is equal to the given RNA_SilentStruct in terms of conformation.
	/// Doesn't check energies.
	RNA_SilentStruct & operator= (
		RNA_SilentStruct const & src
	);

	/// @brief Tells this RNA_SilentStruct object to initialize itself from the given set of lines. Lines should
	/// be of the format
	virtual bool init_from_lines(
		utility::vector1< std::string > const & lines,
		SilentFileData & container
	);

	/// @brief Configure the conformation of the given Pose with the conformational data within this RNA_SilentStruct.
	/// Calls pose.clear() and rebuilds Pose from scratch using FA_STANDARD residue types.
	virtual void fill_pose(
		core::pose::Pose & pose,
		bool const metapatches = true
	) const;

	/// @brief Configure the conformation of the given Pose with the conformational data within
	/// this RNA_SilentStruct. Calls pose.clear() and rebuilds Pose from scratch using the
	/// user-specified residue types.
	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set,
		bool const metapatches = true
	) const;

	//virtual void fill_pose(
	// core::pose::Pose & pose,
	// core::chemical::ResidueTypeSet const & residue_set,
	// bool const use_input_pose
	//) const;

	/// @brief opposite of fill_pose
	virtual void fill_struct( core::pose::Pose const & pose, std::string tag );

	/// @brief print header information
	virtual void print_header( std::ostream& out ) const;

	/// @brief Prints the conformation information within this RNA_SilentStruct to the given std::ostream.
	virtual void print_conformation( std::ostream & output ) const;

	/// @brief data getters/setters
	bool fullatom() const {
		return fullatom_;
	}

	void fullatom( bool fullatom ) {
		fullatom_ = fullatom;
	}

	char secstruct( Size seqpos ) const {
		return secstruct_[seqpos];
	}

	utility::vector1< Real > mainchain_torsions( Size seqpos ) const {
		return mainchain_torsions_[ seqpos ];
	}

	utility::vector1< Real > chi_torsions( Size seqpos ) const {
		return chi_torsions_[ seqpos ];
	}

	Real mainchain_torsion( Size const & seqpos, Size const & torsion_num ) const {
		return mainchain_torsions_[ seqpos ][ torsion_num ];
	}

	Real chi( Size const & seqpos, Size const & torsion_num ) const {
		return chi_torsions_[ seqpos ][ torsion_num ];
	}

	Vector coords( Size seqpos ) const {
		return coords_[seqpos];
	}

	utility::vector1< Vector > coords() const {
		return coords_;
	}

	void set_secstruct( Size const & seqpos, char const & ss ) {
		secstruct_[seqpos] = ss;
	}

	void set_mainchain_torsions( Size const & seqpos, utility::vector1< Real > & torsions ) {
		mainchain_torsions_[seqpos] = torsions;
	}

	void set_chi_torsions( Size const & seqpos, utility::vector1< Real > & torsions ) {
		chi_torsions_[seqpos] = torsions;
	}

	void set_coords( Size const & seqpos, Vector const & coords ) {
		coords_[seqpos] = coords;
	}

	void set_non_main_chain_sugar_coords( Size const & seqpos, utility::vector1< Vector > const & vecs ) {
		non_main_chain_sugar_coords_defined_ = true;
		non_main_chain_sugar_coords_[seqpos] = vecs;
	}

	void set_fold_tree( kinematics::FoldTree const & f ) {
		fold_tree_ = f;
	}

	kinematics::FoldTree const& fold_tree( ) const {
		return fold_tree_;
	}

	void add_jump( kinematics::Jump const & jump ) {
		jumps_.push_back( jump );
	}

	kinematics::Jump const & jump( Size const & jump_num ) const {
		return jumps_[ jump_num ];
	}


	/// @brief returns the positions of the CA atoms in this RNA_SilentStruct.
	/// Useful for RMS calculations.
	virtual ObjexxFCL::FArray2D< Real > get_CA_xyz() const;

	// model quality-related methods.
	virtual Real CA_rmsd( RNA_SilentStruct other_pss );

	/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from the torsions in this
	/// RNA_SilentStruct and the C-alpha atoms from this RNA_SilentStruct.
	virtual Real get_debug_rmsd();

protected:
	bool fullatom_;
	bool non_main_chain_sugar_coords_defined_;

	utility::vector1< char > secstruct_;
	utility::vector1< utility::vector1< Real > > mainchain_torsions_;
	utility::vector1< utility::vector1< Real > > chi_torsions_;
	utility::vector1< Vector > coords_;
	utility::vector1< utility::vector1< Vector > > non_main_chain_sugar_coords_;
	utility::vector1< kinematics::Jump > jumps_;

	kinematics::FoldTree fold_tree_;

}; // class RNA_SilentStruct

} // namespace silent
} // namespace io
} // namespace core

#endif
