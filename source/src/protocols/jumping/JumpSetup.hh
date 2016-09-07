// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file JumpSetup
/// @brief read jump-definition file   setups fold tree an chainbreak variants
/// loop code didn't work because fold-tree to complicated ( overlapping loops )
/// @details
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jumping_JumpSetup_hh
#define INCLUDED_protocols_jumping_JumpSetup_hh


// Unit Headers
#include <protocols/jumping/JumpSetup.fwd.hh>

// Package Headers
//#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/JumpSample.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
//#include <core/scoring/constraints/ConstraintForest.fwd.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>

// Named oddly to prevent name conflicts with including .cc files (who may wish to define their own tracer).
static THREAD_LOCAL basic::Tracer hh_tr( "protocols.jumping", basic::t_info );

namespace protocols {
namespace jumping {


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief two numbers, i.e., jump start end residue, or cut-regions...
class Interval {
public:
	Interval() = default;

	Interval( core::Size start, core::Size end ) :
		start_( start ), end_( end ) {};
	core::Size start() { return start_; };
	core::Size stop() { return end_; };
	core::Size end() { return end_; };
	core::Size start_;
	core::Size end_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// forward declaration
class JumpSample;

class JumpSetup;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief virtual base class: can create a set of jumps and cuts
class BaseJumpSetup : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~BaseJumpSetup() override;
	virtual
	JumpSample
	create_jump_sample( ) const = 0;

	/// @brief returns an ordered FragSet that is compatible with the JumpSample
	/// default: generate jumps from ss-library according to JumpSample
	// if the movemap allows sampling of the down-residue but not the up-residue:
	// include a jump with torsions only for the "down" residue
	// if the movemap allows neither sampling of up or down, don't include the jump
	virtual core::fragment::FragSetOP
	generate_jump_frags( JumpSample const&, core::kinematics::MoveMap const& mm) const;

	/// @brief take from a given JumpSample only those Jumps, which could also have been created by create_jump_sample()
	virtual
	JumpSample clean_jumps( JumpSample const& ) const = 0;

	virtual std::string type_name() const = 0;

};


/// @brief output operator
std::ostream & operator <<(std::ostream & os, JumpSample const & t);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class JumpSetup : public BaseJumpSetup {
public:
	class JumpDef {
	public:
		JumpDef( Interval jump, Interval cut_reg ) :
			jump_( jump ), cut_reg_ ( cut_reg) { };

		Interval jump_;
		Interval cut_reg_;
	};

	std::string type_name() const override {
		return "JumpSetup";
	}
public:

	typedef utility::vector1< JumpDef >::const_iterator const_iterator;
	typedef utility::vector1< JumpDef >::iterator iterator;

	//@brief c'stor
	JumpSetup( Size total_residue ) :
		total_residue_( total_residue )
	{};

	//@brief add a new jump to the list
	void
	add_jump ( JumpDef const& jd ) {
		jumps_.push_back ( jd );
		jump_sample_ = JumpSample( *this );
	}

	void
	add_jump ( Interval const& jump, Interval const& cut_reg ) {
		add_jump( JumpDef( jump, cut_reg ) );
		jump_sample_ = JumpSample( *this );
	}

	void
	add_jump ( core::Size js, core::Size je, core::Size crs, core::Size cre ) {
		add_jump( JumpDef( Interval( js, je ), Interval( crs, cre ) ) );
		jump_sample_ = JumpSample( *this );
	}

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

	void clear() {
		jumps_.clear();
	}
	void
	read_file( std::string file );

	JumpSample
	create_jump_sample( ) const override
	{ return jump_sample_; }

	JumpSample
	clean_jumps( JumpSample const& js ) const override
	{
		hh_tr.Error << "ERROR: JumpSetup::clean_jumps() not implemented" << std::endl;
		return js;
	}

	void set_jump_sample( JumpSample const& jump_sample ) {
		jump_sample_ = jump_sample;
	}

	core::Size
	total_residue() const {
		return total_residue_;
	}


private:
	utility::vector1< JumpDef > jumps_;
	core::Size total_residue_;
	JumpSample jump_sample_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class JumpsFromAllPairings : public BaseJumpSetup {
public:

	JumpsFromAllPairings( const core::Size total_residue, core::scoring::dssp::PairingsList const& jumps, ObjexxFCL::FArray1D_float const& cut_probability ) : total_residue_(total_residue), jumps_(jumps), cut_prob_(cut_probability) {
	}

	JumpSample create_jump_sample() const override {
		return JumpSample(total_residue_, jumps_, cut_prob_);
	}

private:
	core::Size total_residue_;
	core::scoring::dssp::PairingsList jumps_;
	ObjexxFCL::FArray1D_float cut_prob_;
};


//class JumpsFromConstraintForest : public BaseJumpSetup {
//public:
//
// JumpsFromConstraintForest(
//  const core::Size total_residue,
//  core::scoring::constraints::ConstraintForestOP forest,
//  ObjexxFCL::FArray1D_float const& cut_probability );
//
// JumpSample create_jump_sample() const;
// ~JumpsFromConstraintForest();
//
// /// @brief Unimplemented!
// JumpSample
// clean_jumps( JumpSample const& js ) const;
//
//private:
// core::Size total_residue_;
// core::scoring::constraints::ConstraintForestOP forest_;
// ObjexxFCL::FArray1D_float cut_prob_;
//};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class JumpSelector : public BaseJumpSetup {
public:

	class JumpWeightDef {
	public:
		JumpWeightDef( Interval jump, core::Real weight ) :
			jump_( jump ), weight_( weight) { };

		Interval jump_;
		core::Real weight_;
	};

	std::string type_name() const override {
		return "JumpSelector";
	}

public:
	typedef utility::vector1< JumpWeightDef > JumpList;
	typedef JumpList::const_iterator const_iterator;
	typedef JumpList::iterator iterator;

	JumpSelector();
	JumpSelector( std::string ss );
	~JumpSelector() override;

	//@brief add a new jump to the list
	void
	add_jump ( JumpWeightDef const& jd ) {
		jumps_.push_back ( jd );
		total_weight_+=jd.weight_;
	}

	void
	add_jump ( Interval const& jump, core::Real weight ) {
		add_jump( JumpWeightDef( jump, weight ) );
	}

	void
	add_jump ( core::Size js, core::Size je, core::Real weight ) {
		add_jump( JumpWeightDef( Interval( js, je ), weight ) );
	}

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

	Interval
	select_random() const;

	JumpSample
	create_jump_sample( ) const override;

	JumpSample
	clean_jumps( JumpSample const& js ) const override {
		hh_tr.Error << "JumpSelector::clean_jumps() is NOT IMPLEMENTED." << std::endl;
		return js;
	};

	void
	set_secstruct( std::string ss ) {
		secstruct_ = ss;
	}

private:
	std::string secstruct_;
	core::pose::PoseCOP native_pose_;
	JumpList jumps_;
	core::Real total_weight_;
	core::Size min_loop_length_; // minimum length required to introduce cut
	core::Size loop_extension_; // grow loops at either side ( to consilidate things like LLL..LLL into a single loop-region
	core::Size nr_jumps_min_;
	core::Size nr_jumps_max_;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} //jumping
} //protocols
#endif
