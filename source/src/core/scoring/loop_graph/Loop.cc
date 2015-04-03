// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/loop_graph/Loop.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/loop_graph/Loop.hh>
namespace core {
namespace scoring {
namespace loop_graph {

	//Constructor
	Loop::Loop( Size const takeoff_pos,
							Size const landing_pos,
							Size const takeoff_domain,
							Size const landing_domain ):
		takeoff_pos_( takeoff_pos ),
		landing_pos_( landing_pos ),
		takeoff_domain_( takeoff_domain ),
		landing_domain_( landing_domain )
	{}

	//Constructor
	Loop::Loop():
		takeoff_pos_( 0 ),
		landing_pos_( 0 ),
		takeoff_domain_( 0 ),
		landing_domain_( 0 )
	{}


	//Destructor
	Loop::~Loop()
	{}

	/// @brief copy constructor
	Loop::Loop( Loop const & src ):
		ReferenceCount( src )
	{
		*this = src;
	}

	Loop &
	Loop::operator=( Loop const & src ){

		takeoff_pos_ = src.takeoff_pos_;
		landing_pos_ = src.landing_pos_;
		takeoff_domain_ = src.takeoff_domain_;
		landing_domain_ = src.landing_domain_;

		return *this;
	}

} //loop_graph
} //scoring
} //core
