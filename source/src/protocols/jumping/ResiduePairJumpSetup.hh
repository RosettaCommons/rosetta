// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ResiduePairJumpSetup
/// @brief setup jumps specifically to ResiduePairJump, read in definition file
/// @details derived from BaseJumpSetup
/// @author Chu Wang


#ifndef INCLUDED_protocols_jumping_ResiduePairJumpSetup_hh
#define INCLUDED_protocols_jumping_ResiduePairJumpSetup_hh

// Unit Headers
#include <protocols/jumping/ResiduePairJumpSetup.fwd.hh>

// Package Headers
#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/ResiduePairJump.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jumping {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// forward declaration
class JumpSample;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ResiduePairJumpSetup : public BaseJumpSetup {
public:
	class JumpDef {
	public:
		JumpDef( Interval jump, Interval cut_reg ) :
			jump_( jump ), cut_reg_ ( cut_reg) { };

		Interval jump_;
		Interval cut_reg_;
	};

	std::string type_name() const {
		return "ResiduePairJumpSetup";
	}

public:

	typedef utility::vector1< JumpDef >::const_iterator const_iterator;
	typedef utility::vector1< JumpDef >::iterator iterator;


	//@brief c'stor
	ResiduePairJumpSetup() :
		total_residue_( 0 ), root_( 1 )
	{};

	//@brief c'stor
	ResiduePairJumpSetup( Size total_residue ) :
		total_residue_( total_residue )
	{};

	//@brief add a new jump to the list
	void
	add_jump( JumpDef const& jd ) {
		jumps_.push_back ( jd );
	}

	void
	add_jump( Interval const& jump, Interval const& cut_reg ) {
		add_jump( JumpDef( jump, cut_reg ) );
	}

	void
	add_jump( core::Size js, core::Size je, core::Size crs, core::Size cre ) {
		add_jump( JumpDef( Interval( js, je ), Interval( crs, cre ) ) );
	}

	void
	add_residue_pair_jump( ResiduePairJumpOP ptr ) {
		ResiduePairJumps_.push_back( ptr );
	}

	void
	set_root( core::Size root ) {
		runtime_assert( root >= 1 && root <= total_residue_ );
		root_ = root;
	};

	core::Size
	root() const {
		return root_;
	};

	core::Size
	size() const {
		return jumps_.size();
	};

	const_iterator
	begin() const {
		return jumps_.begin();
	};

	const_iterator
	end() const {
		return jumps_.end();
	};

	void
	read_file( std::string file );

	JumpSample
	create_jump_sample( ) const;

	virtual JumpSample
	clean_jumps( JumpSample const& js) const {
		std::cerr << "ERROR ERROR ERROR            unimplemented method!;";
		return js;
	}

	core::fragment::FragSetOP
	generate_jump_frags( JumpSample const&, core::kinematics::MoveMap const& mm) const;

	core::Size
	total_residue() const {
		return total_residue_;
	}
	void clear() {
		jumps_.clear();
	}

private:
	utility::vector1< JumpDef > jumps_;
	core::Size total_residue_;
	core::Size root_;
	utility::vector1 <ResiduePairJumpOP> ResiduePairJumps_;
};


} //jumping
} //protocols
#endif
