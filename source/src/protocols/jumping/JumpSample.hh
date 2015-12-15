// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file JumpSetup
/// @brief read jump-definition file   setups fold tree an chainbreak variants
/// loop code didn't work because fold-tree to complicated ( overlapping loops )
/// @details
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jumping_JumpSample_hh
#define INCLUDED_protocols_jumping_JumpSample_hh

// Unit Headers
#include <protocols/jumping/JumpSetup.fwd.hh>

// Package Headers
#include <protocols/jumping/PairingLibrary.hh>
#include <core/scoring/dssp/PairingsList.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/func/Func.hh>
#include <core/fragment/FrameList.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jumping {

class ChainbreakDistFunc : public core::scoring::func::Func {
public:
	ChainbreakDistFunc( core::Real const x0_in );
	virtual ~ChainbreakDistFunc();

	core::scoring::func::FuncOP
	clone() const;

	virtual bool operator == ( core::scoring::func::Func const & other ) const;
	virtual bool same_type_as_me( core::scoring::func::Func const & other ) const;

	virtual core::Real func( core::Real const x ) const;
	virtual core::Real dfunc( core::Real const x ) const;

private:
	core::Real d2target_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ChainbreakDistFunc();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class Interval;
typedef utility::vector1< jumping::Interval > JumpList;

//@brief JumpSample with a random choice of cut-points according to allowed cut_reg in JumpSetup
//@details this class can handle all the intricacies of actually applying a choice of jumps and cuts to a pose
class JumpSample  {
public:
	JumpSample ()
	: njump_( 0 ),
		bValidTree_( false )
	{};

	~JumpSample() {};
	JumpSample ( JumpSetup const& );

	// generate fold-tree from jumps and cuts
	JumpSample ( core::Size total_residue, core::Size njump, ObjexxFCL::FArray2D_int jumps, ObjexxFCL::FArray1D_int cuts, core::Size root = 1);

	// generate fold-tree from jumps, jump_atoms and cuts
	JumpSample ( core::Size total_residue, core::Size njump, ObjexxFCL::FArray2D_int jumps, ObjexxFCL::FArray2D<std::string> jump_atoms, ObjexxFCL::FArray1D_int cuts, core::Size root = 1);

	// this one is not safe due to FArray,  generate random fold-tree  ( as used by SheetBuilder )
	JumpSample ( core::Size total_residue, core::scoring::dssp::PairingsList const& jumps, ObjexxFCL::FArray1D_float const& cut_probability, core::Size root=1);

	// generate random fold-tree  ( as used by SheetBuilder )
	JumpSample ( core::Size total_residue, core::scoring::dssp::PairingsList const& jumps, core::fragment::SecondaryStructure const& ss_def, core::Size root=1 );


	// wrap to a copy of an existing Fold-Tree
	JumpSample ( core::kinematics::FoldTree const&f );

	// wrap an existing Fold-Tree
	JumpSample( core::kinematics::FoldTreeOP f );

	bool
	is_valid() const {
		return bValidTree_;
	}

	void
	resize( core::Size njump );

	//void generate_tree( core::Size nres, core::kinematics::FoldTree &f ) const;

	void
	set_fold_tree_in_pose( core::pose::Pose &pose ) const;

	//@brief lengthen secondary structure elements around jump-points to be at least 5 residues ( -2..+2 ) long.
	void
	safe_secstruct( core::pose::Pose &pose ) const;

	//@brief transfer native jump RT to poes, Sideeffect: changes fold-tree of native_pose.
	void
	transfer_jumps( core::pose::Pose &pose, core::pose::Pose &native_pose ) const;

	//bool
	// close_chainbreaks( core::pose::Pose &pose ) const;

	//@brief transfer native jump RT to FrameList
	void
	steal_jumps(
		core::pose::Pose &native_pose,
		core::fragment::FrameIterator const& begin,
		core::fragment::FrameIterator const& end
	) const;


	//@brief generate list of frames ( one for each jump )
	void
	generate_jump_frames( core::fragment::FrameList&, core::kinematics::MoveMap const&, bool bWithTorsion = true ) const;

	void
	steal_orientation_and_pleating( core::pose::Pose &native_pose );

	//@brief generate fragset with RTs from library
	void
	generate_jump_frags(
		PairingLibrary const& pairings,
		core::kinematics::MoveMap const&,
		bool bWithTorsion,
		core::fragment::FrameList&
	) const;

	//@brief switch on chainbreak scores for all cuts ( add chainbreak variant )
	void
	add_chainbreaks( core::pose::Pose &pose ) const;

	//@brief switch on chainbreaks-scores only if the separation in the foldtree is smaller than the threshold <max_dist>
	void
	add_chainbreaks( core::pose::Pose &pose, core::Size max_dist, core::kinematics::ShortestPathInFoldTree const& ) const;

	void
	add_chainbreaks_as_distance_constraint( core::pose::Pose &pose ) const;

	void
	add_chainbreaks_as_distance_constraint(
		core::pose::Pose &pose,
		core::Size max_dist,
		core::kinematics::ShortestPathInFoldTree const& sp
	) const;

	void
	remove_chainbreaks( core::pose::Pose &pose ) const;

	void
	set_fold_tree_and_chainbreaks( core::pose::Pose &pose ) const {
		set_fold_tree_in_pose( pose );
		add_chainbreaks( pose );
	}

	core::Size
	downstream_res_nr( core::Size jump_nr ) const {
		return  fold_tree_->downstream_jump_residue( jump_nr );
	}

	core::Size size() const
	{
		return njump_;
	}

	ObjexxFCL::FArray2D_int const& jumps() const {
		return jumps_;
	}

	ObjexxFCL::FArray1D_int const& cuts() const {
		return cuts_;
	}

	ObjexxFCL::FArray2D< std::string > const& jump_atoms() const {
		return jump_atoms_;
	}

	bool has_orientation_and_pleating() const;

	/// @brief output operator
	friend std::ostream & operator <<(std::ostream & os, JumpSample const & t);

	core::kinematics::FoldTree const& fold_tree() const {
		return *fold_tree_;
	}

	core::Size total_residue() const { return total_residue_; }

	/// @brief dump file with little pymol script that shows jumps and cuts
	void
	dump_pymol( std::string fn ) const;

	core::scoring::dssp::Pairing
	get_pairing( core::Size res1, core::Size res2 ) const;

private:
	void
	generate_tree_from_jumps_and_cuts(core::Size root=1);

	void
	generate_random_tree_from_jumps( ObjexxFCL::FArray1D_float const& prob, core::Size root=1 );

	void
	apply_to( core::kinematics::FoldTreeOP f );

	void jumps2pairings();

	void
	correct_jump_atoms_for_fragments() const;

	core::Size total_residue_;
	core::Size njump_;

	// the latter two contain redundant information
	core::scoring::dssp::PairingsList jump_pairings_;
	ObjexxFCL::FArray2D_int jumps_;
	ObjexxFCL::FArray2D< std::string > jump_atoms_;

	ObjexxFCL::FArray1D_int cuts_;
	core::kinematics::FoldTreeOP fold_tree_;
	bool bValidTree_;

};

} //jumping
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jumping_JumpSample )
#endif // SERIALIZATION


#endif
