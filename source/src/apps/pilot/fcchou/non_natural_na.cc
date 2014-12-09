// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/chemical/rna/util.hh>

#include <utility/io/ozstream.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <utility/vector1.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;

using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////
void test() {
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;

	ResidueTypeSetCAP rsd_set(
		core::chemical::ChemicalManager::get_instance()
		->residue_type_set( RNA ) );

	pose::Pose pose;
	pose::make_pose_from_sequence( pose, "g", *rsd_set );
	ResidueOP new_res1(
			core::conformation::ResidueFactory::create_residue(
				rsd_set->name_map( "INO" ) ) );
	pose.append_residue_by_jump( *new_res1, pose.total_residue() );

	ResidueOP new_res2(
			core::conformation::ResidueFactory::create_residue(
				rsd_set->name_map( "IGU" ) ) );
	pose.append_residue_by_jump( *new_res2, pose.total_residue() );

	ResidueOP new_res3(
			core::conformation::ResidueFactory::create_residue(
				rsd_set->name_map( "ICY" ) ) );
	pose.append_residue_by_jump( *new_res3, pose.total_residue() );

	ResidueOP new_res4(
			core::conformation::ResidueFactory::create_residue(
				rsd_set->name_map( "1AP" ) ) );
	pose.append_residue_by_jump( *new_res4, pose.total_residue() );

	pose.dump_pdb( "all.pdb" );
}
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	test();

	protocols::viewer::clear_conformation_viewers();

	exit( 0 );
}
///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
  try {

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
  return 0;
}

