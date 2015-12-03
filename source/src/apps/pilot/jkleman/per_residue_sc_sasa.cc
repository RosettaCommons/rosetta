// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    per_residue_sc_sasa.cc
/// @brief   Computes the SASA of all the sidechains in the protein
///    also computes the relative one
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/docking/membrane/MPFindInterfaceMover.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

// utility
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa/util.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.jkleman.per_residue_sc_sasa" );

//////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace utility;
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring::sasa;
		using namespace protocols::jd2;
		using namespace protocols::scoring;

		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		// cry if PDB not given
		if ( ! option[OptionKeys::in::file::s].user() ) {
			throw new utility::excn::EXCN_Msg_Exception("Please provide PDB file!");
		}

		// read in pose
		Pose pose;
		core::import_pose::pose_from_pdb( pose, option[OptionKeys::in::file::s].value_string() );

		// get residue SASA
		utility::vector1< Real > res_sasa = per_res_sc_sasa( pose );
		utility::vector1< Real > rel_res_sasa = rel_per_res_sc_sasa( pose );

		// print SASA
		for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
			TR << "PDB resn: " << pose.pdb_info()->pose2pdb( i ) << " " << pose.residue( i ).name3() << " sasa: " << res_sasa[ i ] << " rel sasa: " << rel_res_sasa[ i ] << std::endl;
		}
	}

catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
