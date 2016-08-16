// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/Loop.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_Loop_HH
#define INCLUDED_core_scoring_loop_graph_Loop_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/loop_graph/Loop.fwd.hh>
#include <core/scoring/types.hh>
#include <iostream>

namespace core {
namespace scoring {
namespace loop_graph {

///  Directed edge between one position in the pose to another position.
///  Records takeoff_pos & landing_pos (in full_model numbering).
///  Also records 'domain', i.e. if there are multiple poses in a collection,
///   which pose.
class Loop: public utility::pointer::ReferenceCount {

public:

	Loop();
	Loop( Size const takeoff_pos,
		Size const landing_pos,
		Size const takeoff_domain,
		Size const landing_domain );

	//constructor
	// Undefined, commenting out to fix PyRosetta build  Loop();

	//destructor
	~Loop();

	/// @brief copy constructor
	Loop( Loop const & src );

	Loop &
	operator=( Loop const & src );

public:

	void set_takeoff_pos( core::Size const & setting ){ takeoff_pos_ = setting; }
	core::Size takeoff_pos() const{ return takeoff_pos_; }

	void set_landing_pos( core::Size const & setting ){ landing_pos_ = setting; }
	core::Size landing_pos() const{ return landing_pos_; }

	void set_takeoff_domain( core::Size const & setting ){ takeoff_domain_ = setting; }
	core::Size takeoff_domain() const{ return takeoff_domain_; }

	void set_landing_domain( core::Size const & setting ){ landing_domain_ = setting; }
	core::Size landing_domain() const{ return landing_domain_; }

	/// @brief a and b are the same atom
	friend
	inline
	bool
	operator ==(
		Loop const & a,
		Loop const & b ){
		return ( ( a.takeoff_pos_ == b.takeoff_pos_ ) &&
			( a.landing_pos_ == b.landing_pos_ ) &&
			( a.takeoff_domain_ == b.takeoff_domain_ ) &&
			( a.landing_domain_ == b.landing_domain_ ) );
	}


	/// @brief Test IO operator for debug and Python bindings
	friend
	std::ostream & operator << ( std::ostream & os, Loop const & loop){
		os << "LOOP ";
		os << "TAKEOFF:  ";
		os << loop.takeoff_pos_ << " from domain " << loop.takeoff_domain_ << "; ";
		os << "LANDING:  ";
		os << loop.landing_pos_ << " from domain " << loop.landing_domain_;
		return os;
	}

private:

	core::Size takeoff_pos_;
	core::Size landing_pos_;
	core::Size takeoff_domain_;
	core::Size landing_domain_;

};


} //loop_graph
} //scoring
} //core

#endif
