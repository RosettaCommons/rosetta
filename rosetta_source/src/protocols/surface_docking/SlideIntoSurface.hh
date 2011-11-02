// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//     under license.
// (c) The Rosetta software is developed by the contributing members of the
//     Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
//     this can be
// (c) addressed to University of Washington UW TechTransfer,
//     email: license@u.washington.edu.

/// @file SurfaceOrientMover.hh
/// @brief <add a description of the class>
/// @author Robin A Thottungal (raugust1@jhu.edu)

#ifndef INCLUDED_protocols_surface_docking_SlideIntoSurface_hh
#define INCLUDED_protocols_surface_docking_SlideIntoSurface_hh

// Unit Headers

// Package headers
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility Headers

// C++ Headers
#include <string>
#include <map>
#include <list>

#include <utility/vector1.hh>

namespace protocols {
namespace surface_docking {

class SlideIntoSurface : public moves::Mover {

public:

	SlideIntoSurface();

	SlideIntoSurface(core::Size const rb_jump);

	//destructor
	~SlideIntoSurface();

	/// virtual functions that get overloaded or
	//                           called from the inheriting classes
	void apply( core::pose::Pose & );

	virtual std::string get_name() const;

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	// which jump to use for docking
	core::Size rb_jump_;
	};

/// @brief Slides docking partners together by monitoring fa_rep.
/// @details
///		If partners are already touching, no change is made.
///		Separation will be 1A or less after calling this function.
class FaSlideIntoSurface : public moves::Mover
{
public:
	FaSlideIntoSurface();
	FaSlideIntoSurface( core::Size const rb_jump);

	//destructor
	~FaSlideIntoSurface();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	core::Size rb_jump_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real tolerance_; ///< how accurate do you want to be?

};


} // surfaceDockingProtocols
} // protocols

#endif
