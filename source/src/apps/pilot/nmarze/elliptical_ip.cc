// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/nmarze/elliptical_ip.cc
/// @brief pilot app for elliptical initial placement
/// @author Nick Marze (nickmarze@gmail.com)


#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/exit.hh>

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <protocols/docking/EllipsoidalRandomizationMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>

static basic::Tracer TR("apps.pilot.nmarze.elliptical_ip");


/////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {

		protocols::jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		protocols::docking::DockingInitialPerturbationOP elliptical_ip = new protocols::docking::DockingInitialPerturbation( 1, true );

//		protocols::docking::EllipsoidalRandomizationMoverOP elliptical_ip = new protocols::docking::EllipsoidalRandomizationMover( 1, false );
		protocols::jd2::JobDistributor::get_instance()->go( elliptical_ip );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
				return -1;
    }
    return 0;
}
