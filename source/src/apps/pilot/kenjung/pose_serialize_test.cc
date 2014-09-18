// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ken Jung
/// @brief serialize test for pose
/// this really doesn't test a lot of things, datacache, energies, residuetypes, symmetry, etc.

#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <utility/excn/Exceptions.hh>

#include <iostream>
#include <sstream>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer trmain( "test" );

OPT_1GRP_KEY(File, m, file)

int
main( int argc, char * argv [] )
{
    try {
		// initialize core
		NEW_OPT(m::file, "file", "");

		devel::init(argc, argv);

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pose;
		using namespace core::scoring;

		core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		PoseSP tmppose(new Pose());
		core::import_pose::pose_from_pdb( *tmppose, *residue_set, option[ m::file]().name() );
		tmppose->dump_pdb(option[ m::file]().name()+".dump");

		ScoreFunctionOP scorefxn = get_score_function();
		PoseSP out;
		(*scorefxn)(*tmppose);
		std::stringstream s;
		boost::archive::binary_oarchive oa(s);
		boost::archive::binary_iarchive ia(s);
		trmain << "pose size=" << tmppose->total_residue() << std::endl;
		trmain << "pdb res num = " << tmppose->pdb_info()->number(3) << std::endl;
		trmain << "before in" << std::endl;
		oa << tmppose;
		trmain << "after in, size is " << s.str().size() <<  std::endl;
		ia >> out;
		trmain << "after out" << std::endl;
		trmain << "pose size=" << out->total_residue() << std::endl;
		trmain << "pdb res num = " << out->pdb_info()->number(3) << std::endl;
		(*scorefxn)(*out);
		trmain << "dumping" << std::endl;
		out->dump_pdb(option[ m::file]().name()+".redump");

    } catch ( utility::excn::EXCN_Base const & e ) {
			std::cerr << "caught exception " << e.msg() << std::endl;
			return -1;
    }
    return 0;
}
