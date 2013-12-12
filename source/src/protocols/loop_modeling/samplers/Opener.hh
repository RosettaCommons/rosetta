// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_loop_modeling_samplers_Opener_HH
#define INCLUDED_protocols_loop_modeling_samplers_Opener_HH

// Unit headers
#include <protocols/loop_modeling/samplers/Opener.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocols headers
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <boost/utility.hpp>

namespace protocols {
namespace loop_modeling {
namespace samplers {

using protocols::moves::Mover;
using protocols::moves::MoverOP;
using protocols::loops::Loop;
using boost::noncopyable;

class Opener : public Mover, protected noncopyable {

public:
	MoverOP setup(Loop const & loop);

protected:
	Loop loop_;

};

}
}
}

#endif

