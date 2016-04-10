// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


#ifndef INCLUDED_protocols_simple_moves_SmoothFragmentMover_HH
#define INCLUDED_protocols_simple_moves_SmoothFragmentMover_HH

// Unit Headers
#include <protocols/simple_moves/SmoothFragmentMover.fwd.hh>

// Package Headers
#include <protocols/simple_moves/FragmentMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

typedef utility::vector1< core::Real > ScoreList;

class FragmentCost : public utility::pointer::ReferenceCount {
protected:
	// Constructor protected for base class
	FragmentCost( std::string type, core::Real cutoff ) : type_( type ), cutoff_( cutoff ) {};
	virtual ~FragmentCost();
public:
	// accesor
	std::string const& type() {
		return type_;
	};

	core::Real cutoff() const {
		return cutoff_;
	};

	// the meaning of life
	virtual void score( core::fragment::Frame const&, core::pose::Pose const& pose,  ScoreList &scores ) = 0;

private:
	std::string type_;

	// SmoothFragmentMover will randomly choose from all fragments that are below cutoff_
	core::Real cutoff_;
};

class SmoothFragmentMover : virtual public ClassicFragmentMover {
public:
	typedef ClassicFragmentMover Parent;

public:
	SmoothFragmentMover(
		core::fragment::FragSetCOP fragset,
		FragmentCostOP cost );


	SmoothFragmentMover(
		core::fragment::FragSetCOP fragset,
		core::kinematics::MoveMapCOP movemap,
		FragmentCostOP cost );

	~SmoothFragmentMover();

	SmoothFragmentMoverOP shared_from_this() { return utility::pointer::dynamic_pointer_cast<SmoothFragmentMover>( Mover::shared_from_this() ); }

	virtual std::string get_name() const;

	//  void apply( core::pose::Pose & );
protected:
	SmoothFragmentMover(
		core::fragment::FragSetCOP fragset,
		core::kinematics::MoveMapCOP movemap,
		FragmentCostOP cost,
		std::string move_type );

	// frame and fragment of choice, returns false if no good fragment is found
	virtual
	bool
	choose_fragment(
		core::fragment::FrameList const&,
		core::pose::Pose const&,
		Size &frame_num,
		Size &frag_num
	) const;

	virtual
	bool
	use_ss_length_screen() const;

private:
	FragmentCostOP cost_;

	// choose randomly fragments that are below cutoff_
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// core::Real cutoff_;

};

} //simple_moves
} //protocols

#endif
