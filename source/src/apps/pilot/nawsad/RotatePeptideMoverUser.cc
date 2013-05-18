// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   RotatePeptideMoverUser.cpp
/// @brief  This is a pilot app for testing RotatePeptideMover.cc
/// @author Nawsad Alam

#include <devel/rotate_peptide/RotatePeptideMover.hh>
#include <core/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>


int main(int argc, char *argv[])
{
	using namespace protocols;
	using namespace protocols::jd2;

	// initialize core
	core::init(argc, argv);
	
	/////////////////////

	/////////////////////
	core::Size pep_chain = 2;
//	core::pose::Pose & my_pose
	// distribute the mover
	moves::MoverOP my_mover = new devel::rotate_peptide::RotatePeptideMover(pep_chain);
//	my_mover.apply(& pose)
	JobDistributor::get_instance()->go(my_mover);
}
