// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/jumping/ResiduePairJump.hh
/// @brief a class to create jump transform from a pair of residues
/// @details
///  This class is to create possible jump transforms between a pair of residues
///  when their sidechains are locked in certain geometry constraints. It starts
///  from the predefined constraints and takes backbone-independent rotamer
///  conformation into account and reversely generate their backbone positions.
///  Then a jump transform is measured.
/// @author Chu Wang


#ifndef INCLUDED_protocols_jumping_ResiduePairJump_hh
#define INCLUDED_protocols_jumping_ResiduePairJump_hh

// Unit Headers
#include <protocols/jumping/ResiduePairJump.fwd.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#ifdef WIN32
#include <core/pack/rotamer_set/RotamerSet.hh>
#endif

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/Jump.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <map>

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jumping {
///////////////////////////////////////////////////////////////////////////////////
/// @brief a single residue component of a ResiduePairJump class
///
/// @details a residue type with information which atoms to be used to define a jump
///and cst.
class ResiduePairJumpSingle : public utility::pointer::ReferenceCount {

public:
	/// @brief empty constructor
	ResiduePairJumpSingle();

	virtual ~ResiduePairJumpSingle();

	/// @brief constructed by residue_type
	ResiduePairJumpSingle(
		core::chemical::ResidueType const & residue_type
	);

	/// @brief have jumpAtoms been defined?
	inline
	bool jumpAtoms_defined() const { return jumpAtoms_.size() == 3; }

	/// @brief have cstAtoms been defined?
	inline
	bool cstAtoms_defined() const { return cstAtoms_.size() == 3; }

	/// @brief whether this residue has internal flexibility
	inline
	bool fixResidue() const { return fixResidue_; }

	/// @brief access residue type
	core::chemical::ResidueTypeOP residueType() const;

	/// @brief access all three jumpAtoms at once
	inline
	utility::vector1< std::string > const & jumpAtoms() const { return jumpAtoms_; }

	/// @brief access single jumpAtom one at a time
	inline
	std::string const jumpAtoms( core::Size i ) const { return jumpAtoms_[i]; }

	/// @brief access all three cstAtoms at once
	inline
	utility::vector1< std::string > const & cstAtoms() const { return cstAtoms_; }

	/// @brief access single cstAtom one at a time
	inline
	std::string const cstAtoms( core::Size i ) const { return cstAtoms_[i]; }

	/// @brief set all three jumpAtoms at once
	void set_jumpAtoms( utility::vector1< std::string > const & jump_atoms );

	/// @brief set single jumpAtom one at a time
	void set_jumpAtoms( core::Size i, std::string const & atom_name );

	/// @brief set all three cstAtoms at once
	void set_cstAtoms( utility::vector1< std::string > const & cst_atoms );

	/// @brief set all single cstAtom one at a time
	void set_cstAtoms( core::Size i, std::string const & atom_name );

private:
	/// @brief internal data -- residue type
	core::chemical::ResidueTypeOP residueType_;
	/// @brief internal data -- atom names (3) for creating the jump
	utility::vector1< std::string >  jumpAtoms_;
	/// @brief internal data -- atom names (3) for defining cst
	utility::vector1< std::string >  cstAtoms_;
	/// @brief internal data -- whether to consider internal flexibility within this residue
	bool fixResidue_;
};


enum cstType {
	disAB = 1,
	angleA,
	angleB,
	dihedralA,
	dihedralAB,
	dihedralB,
	num_cstType = dihedralB
};

enum dofType {
	rot1 = 1,
	rot2,
	cstJump,
	num_dofType = cstJump
};

class ResiduePairJump : public utility::pointer::ReferenceCount {

	typedef std::map< cstType, utility::vector1< core::Real > > cstInfoMap;
	typedef std::map< cstType, utility::vector1< core::Real > >::iterator cstInfoMapIterator;

	//typedef std::map< cstType, core::id::DOF_ID > cstTypeToDofMap;
	//typedef std::map< cstType, core::id::DOF_ID >::iterator cstTypeToDofMapIterator;

public:

	// empty constructor
	ResiduePairJump();
	virtual ~ResiduePairJump();

	// constructed by two input ResidueTypes
	ResiduePairJump(
		core::chemical::ResidueType const & residue1,
		core::chemical::ResidueType const & residue2
	);

	// add two input ResidueTypes to define this pair ( erase old data );
	void add_residue_pair(
		core::chemical::ResidueType const & residue1,
		core::chemical::ResidueType const & residue2
	);

	void add_residue_single(
		core::chemical::ResidueType const & residue
	);

	inline
	bool jumpAtoms_defined () const {
		return residues_[1]->jumpAtoms_defined() && residues_[2]->jumpAtoms_defined();
	}

	inline
	bool cstAtoms_defined () const {
		return residues_[1]->cstAtoms_defined() && residues_[2]->cstAtoms_defined();
	}

	inline
	bool cstInfoMap_defined() const {
		return cstInfoMap_.size() >= core::Size( dihedralB );
	}

	inline
	void set_cstAtoms(
		core::Size rsd,
		core::Size atom,
		std::string name
	)
	{
		residues_[rsd]->set_cstAtoms( atom, name );
	}

	inline
	utility::vector1< std::string > const & cstAtoms( core::Size i ) const {
		return residues_[i]->cstAtoms();
	}

	inline
	void set_jumpAtoms(
		core::Size rsd,
		core::Size atom,
		std::string name
	)
	{
		residues_[rsd]->set_jumpAtoms( atom, name );
	}

	inline
	utility::vector1< std::string > const & jumpAtoms( core::Size i ) const {
		return residues_[i]->jumpAtoms();
	}

	void set_cstInfo( cstType type, core::Real value );

	void set_cstInfo( cstType type, utility::vector1< core::Real >  const & values );

	core::fragment::FrameOP generate_frame();

	void init_mini_pose();

	void diversify_dof_conformers( dofType type, core::Size max_index );

	void diversify_dof_conformers();

	core::pose::PoseOP apply_dof_conformer( std::map< dofType, core::Size > const & conformer_map );

	//void setup_cstTypeToDofMap();

	void build_sidechain_rotamers();

	void build_cst_conformer_jumps();

	void diversify_cst_conformers( cstType type, core::Size max_index );

	void diversify_cst_conformers();

private:
	// residue types
	utility::vector1< ResiduePairJumpSingleOP > residues_;
	// discrete conformation for each cst parameters
	std::map< cstType, utility::vector1< core::Real > > cstInfoMap_;
	// a vector of all possible conformers derived from the combination of different cst types
	utility::vector1< std::map< cstType, core::Size > > cst_conformers_;
	// jumpSet derived from all cst_conformers
	utility::vector1< core::kinematics::Jump > cst_jumps_;
	// rotamerSet for each residues if possible
	utility::vector1< core::pack::rotamer_set::RotamerSetOP > rotsets_;
	// a vector of all possible conformers derived from the combination of sidechain rotamers and cst_jumps
	utility::vector1< std::map< dofType, core::Size > > dof_conformers_;
	// data, miniPose consisting of these two residues
	core::pose::PoseOP miniPose_;
};


} //protocols
} //jumping

#endif

