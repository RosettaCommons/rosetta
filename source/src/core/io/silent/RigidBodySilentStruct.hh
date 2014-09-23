// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/ProteinSilentStruct.hh
///
/// @brief Representation of rosetta++ protein silent-file structures.
/// @author James Thompson, Mike Tyka

#ifndef INCLUDED_core_io_silent_RigidBodySilentStruct_hh
#define INCLUDED_core_io_silent_RigidBodySilentStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace silent {



class RigidBodySilentStruct : public SilentStruct {

public:
	RigidBodySilentStruct(
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	) {
		symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
		symminfo_->set_use_symmetry(false);
		fill_struct( pose, tag );
		write_fold_tree_ = false;
	} // RigidBodySilentStruct

	RigidBodySilentStruct()
	{
		nres( 0 );
		decoy_tag( "empty_tag" );
		symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo() );
		symminfo_->set_use_symmetry(false);
		write_fold_tree_ = false;
	}

	/// @brief Returns a new RigidBodySilentStruct with a copy of the information
	/// in this RigidBodySilentStruct.
	virtual SilentStructOP clone() const {
		return SilentStructOP( new RigidBodySilentStruct( *this ) );
	}

	// destructor
	~RigidBodySilentStruct() {}

	/// @brief Test if this RigidBodySilentStruct is equal to the given
	/// RigidBodySilentStruct in terms of conformation. Doesn't check energies.
	RigidBodySilentStruct & operator= (
		RigidBodySilentStruct const & src
	);

	/// @brief Tells this RigidBodySilentStruct object to initialize itself from
	//the given set of lines.
	virtual bool init_from_lines(
		utility::vector1< std::string > const & lines,
		SilentFileData & container
	);

	/// @brief Configure the conformation of the given Pose with the
	/// conformational data within this RigidBodySilentStruct.
	/// sets the jump stored in the RigidBodySilentStruct
	virtual void fill_pose(
		core::pose::Pose & pose
	) const;

	/// @brief Configure the conformation of the given Pose with the
	/// conformational data within this RigidBodySilentStruct.
	// invalid to use
	virtual void fill_pose(
		core::pose::Pose & pose,
		core::chemical::ResidueTypeSet const & residue_set
	) const;

	/// @brief opposite of fill_pose
	virtual void fill_struct(
		core::pose::Pose const & pose,
		std::string tag = "empty_tag"
	);

	/// @brief Prints the conformation information within this
	// RigidBodySilentStruct to the given std::ostream.
	virtual void print_conformation( std::ostream & output ) const;

	//lin Symmetry
	// @lin - move these to the .cc file so you can only include SymmetryInfo.fwd.hh!
	bool is_symmetric() const { return symminfo_->get_use_symmetry();	}

	void symmetry_info( core::conformation::symmetry::SymmetryInfo & s ) {
		symminfo_ = core::conformation::symmetry::SymmetryInfoOP( new core::conformation::symmetry::SymmetryInfo( s ) );
	}

	core::conformation::symmetry::SymmetryInfoCOP symmetry_info( ) const {
		return symminfo_;
	}

	void add_jump( kinematics::Jump const & jump ) {
		jumps_.push_back( jump.rt() );
	}

	void add_rt( kinematics::RT const & rt ) {
		jumps_.push_back( rt );
	}

	/// @brief returns the number of jumps held in this container.
	Size njumps() const {
		return jumps_.size();
	}

	// it's really odd that this function is called jump, but returns an RT.
	kinematics::RT const & jump( Size const jump_num ) const {
		return jumps_[ jump_num ];
	}

	/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from
	/// the torsions in this RigidBodySilentStruct and the C-alpha atoms from this
	/// RigidBodySilentStruct.
	virtual Real get_debug_rmsd() {
		return 1000.0; //cannot do this with this bare-bones silent-struct type
	}

	void fold_tree( kinematics::FoldTree const& );
	kinematics::FoldTree const& fold_tree() const;

	virtual ObjexxFCL::FArray2D< Real > get_CA_xyz() const;

protected:

private: // private member functions

public:
	//	virtual core::Size mem_footprint() const;

private:
	utility::vector1< kinematics::RT > jumps_;
	core::conformation::symmetry::SymmetryInfoOP symminfo_;
	kinematics::FoldTreeOP fold_tree_;
	bool write_fold_tree_;
}; // class RigidBodySilentStruct


} // namespace silent
} // namespace io
} // namespace core

#endif

// I will be removing this #include in not too long; if you need to use
// this class, you'll also have to #include the .tmpl.hh file.
// #include <core/io/silent/RigidBodySilentStruct.tmpl.hh>
