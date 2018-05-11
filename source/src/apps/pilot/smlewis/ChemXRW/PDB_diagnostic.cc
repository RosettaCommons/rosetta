// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/smlewis/ChemXRW/PDB_diagnostic.cc
/// @brief This app is meant to report if Rosetta can successfully read a PDB, and if not, help with diagnostics on why it failed. Spiritual child of old AnchorFinder.


#include <devel/scientific_tests/PDBDiagnosticMover.hh>

#include <basic/options/option.hh>

#include <devel/init.hh>

#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

static basic::Tracer TR( "PDB_diagnostic" );

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		//force these options: don't leak memory, don't pack on read-in, and never ever design
		basic::options::option[ basic::options::OptionKeys::jd2::delete_old_poses ].value(true);
		basic::options::option[ basic::options::OptionKeys::packing::pack_missing_sidechains ].value(false);
		basic::options::option[ basic::options::OptionKeys::packing::repack_only ].value(true);

		//force no-output, otherwise we are wasting time
		bool const score_only(basic::options::option[ basic::options::OptionKeys::out::file::score_only ].user());
		if ( !score_only ) { //if not, do nothing
			//do nothing
			TR << "PDB_diagnostic requires -score_only (or you'll write god knows how many PDBs out!)" << std::endl;
		} else { //do work

			protocols::moves::MoverOP test_mover(new devel::scientific_tests::PDBDiagnosticMover);

			protocols::jd2::JobDistributor::get_instance()->go(test_mover);
		} //if score_only
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}


	TR << "*****************COMPLETED******************" << std::endl;
	return 0;
}
