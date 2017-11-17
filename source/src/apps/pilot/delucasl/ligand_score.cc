// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/delucasl/ligand_score.cc
/// @author Sam DeLuca

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/qsar/qsarMover.hh>
#include <devel/init.hh>

int main(int argc, char* argv[])
{
    try {
	devel::init(argc,argv);


	protocols::qsar::qsarMoverOP mover = new protocols::qsar::qsarMover(40,0.25);

	mover->add_grid("atr");
	mover->add_grid("rep");
	mover->add_grid("hba");
	mover->add_grid("hbd");
	mover->add_grid("polariz");

	mover->set_chain("X");
	protocols::jd2::JobDistributor::get_instance()->go(mover);

	mover->write_all_grids("test_");

	//protocols::ligand_docking::qsar::qsarMoverOP
    } catch (utility::excn::Exception const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
    return 0;
}
