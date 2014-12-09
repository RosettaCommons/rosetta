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
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergy.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>

#include <core/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/util/basic.hh>
#include <core/io/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, blah )


/////////////////////////////////////////////////////////////////////////////////
void
benzene_pair_score_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::scoring;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	// Read in pose with two benzenes. "Z" = ligand. Note need flag:
	//         -extra_res_fa BZN.params -s two_benzenes.pdb
	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	io::pdb::pose_from_pdb( pose, *rsd_set, infile );

	//This is only for graphics...
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//////////////////////////////////////////////////
	// Set up fold tree -- "chain break" between two ligands, right?
	kinematics::FoldTree f( pose.total_residue() );
	Size start( 1 ), end( 2 );
	f.new_jump( start, end, start );
	// Uh, why not?
	f.set_jump_atoms( 1," C1 "," C1 ");
	pose.fold_tree( f );

	//////////////////////////////////////////////////////////////////
	// displace in x, y, z by 2.0 A... just checking coordinate system
	//This jump should code for no translation or rotation -- two_benzenes.pdb
	// has the two benzenes perfectly superimposed.
	kinematics::Jump jump( pose.jump( 1 ) );

	jump.set_translation( Vector( 2.0, 0.0, 0.0 ) );
	pose.set_jump( 1, jump );
	pose.dump_pdb( "shift_x.pdb" );

	jump.set_translation( Vector( 0.0, 2.0, 0.0 ) );
	pose.set_jump( 1, jump );
	pose.dump_pdb( "shift_y.pdb" );

	jump.set_translation( Vector( 0.0, 0.0, 2.0 ) );
	pose.set_jump( 1, jump );
	pose.dump_pdb( "shift_z.pdb" );


	//////////////////////////////////////////////////////////////////
	// OK, how about a score function?
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	scorefxn->set_weight( fa_elec, 1.0 );

	jump.set_translation( Vector( 0.0, 3.0, 3.0 ) );
	pose.set_jump( 1, jump );
	pose.dump_pdb( "test.pdb" );
	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	//////////////////////////////////////////////////////////////////
	// Now iterate over a bunch of translations and rotations.
	using namespace core::io::silent;
	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	Size n( 0 );
	for (int i = -240; i <= 240; i++ ) {
		for (int j = 1; j <= 320; j++ ) {
			Real const x = 0.000001; //arcane thing, rosetta can't handle (0.0,0.0,0.0)
			Real const y = i / 40.0;
			Real const z = j / 40.0;

			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( 1, jump );

			(*scorefxn)( pose );

			n++;
			std::string out_tag( "S_"+lead_zero_string_of(n,4) );
			BinarySilentStruct s( pose,  out_tag);
			s.add_energy( "y", y );
			s.add_energy( "z", z );
			silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
		}
	}
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	benzene_pair_score_test();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	//Uh, options?
	NEW_OPT( blah, "blah", false );


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
