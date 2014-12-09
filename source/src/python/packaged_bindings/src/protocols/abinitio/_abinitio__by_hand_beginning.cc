// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include "boost/python.hpp"


#include <protocols/simple_moves/FragmentMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

namespace bp = boost::python;


protocols::simple_moves::ClassicFragmentMoverOP create_ClassicFragmentMover(core::fragment::FragSet const & fragset,
																		core::kinematics::MoveMap const & movemap)
{
	return new protocols::simple_moves::ClassicFragmentMover(fragset.clone(), new core::kinematics::MoveMap(movemap));
}

protocols::simple_moves::ClassicFragmentMoverOP create_ClassicFragmentMover_CL(core::fragment::ConstantLengthFragSet const & fragset,
																		core::kinematics::MoveMap const & movemap)
{
	return new protocols::simple_moves::ClassicFragmentMover(fragset.clone(), new core::kinematics::MoveMap(movemap));
}


void __abinitio_by_hand_beginning__()
{
    bp::def("create_ClassicFragmentMover", create_ClassicFragmentMover);
    bp::def("create_ClassicFragmentMover", create_ClassicFragmentMover_CL);
}
