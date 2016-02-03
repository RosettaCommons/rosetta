// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file swa_monte_carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/util.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/NamedAtomID.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/farna/util.hh>
#include <protocols/viewer/viewers.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/stream_util.hh>
#include <numeric/xyz.functions.hh>

//////////////////////////////////////////////////////////
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.rhiju.block_stack_test" );

using namespace core;
using utility::tools::make_vector1;

// Turn this into a unit test soon.

///////////////////////////////////////////////////////////////
Vector
get_repl_xyz( pose::Pose const & pose,
							Size const & i, utility::vector1< std::string > const & atom_names )
{
	using namespace kinematics;
	using namespace id;

	Distance const STACK_DIST( 3.4 );
	utility::vector1< NamedAtomID > named_atom_ids;
	Size N( atom_names.size() );
	for ( Size n = 1; n <= N; n++ ) named_atom_ids.push_back( NamedAtomID( atom_names[n], i ) );

	Vector centroid;
	for ( Size n = 1; n <= N; n++ ) centroid += pose.xyz( named_atom_ids[ n ] );
	centroid /= Real( N );

	Stub stub( centroid, pose.xyz(named_atom_ids[N-2]), pose.xyz(named_atom_ids[N-1]), pose.xyz(named_atom_ids[N]) );
	for ( Size n = 1; n <= N; n++ ) {
		Vector v_local = stub.global2local( pose.xyz( named_atom_ids[ n ] ) );
		runtime_assert( std::abs( v_local.z() ) < 1.0e-2 );
	}

	// if ( i < pose.total_residue() ) {
	// 	for ( Size n = 1; n <= pose.residue( i+1 ).natoms(); n++ ) {
	// 		TR << pose.residue( i+1 ).atom_name( n ) << " --> " << stub.global2local( pose.xyz( AtomID( n, i+1 ) ) ).z() << std::endl;
	// 	}
	// }

	return ( stub.local2global( Vector( 0.0, 0.0, STACK_DIST ) ) );

}


///////////////////////////////////////////////////////////////
void
block_stack_test()
{
	using namespace core;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::import_pose;
	using namespace core::id;
  using namespace core::scoring;
  using namespace core::chemical;
  using namespace core::kinematics;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	PoseOP pose_op =	pose_from_file( "ucga_ucga_helix.pdb" , core::import_pose::PDB_file);
	Pose & pose = *pose_op;

	EnergyBaseStackList energy_stack_list = get_scored_base_stack_list( pose );
	for ( EnergyBaseStackList::const_iterator it = energy_stack_list.begin(); it != energy_stack_list.end(); ++it ) {
		BaseStack const & stack = it->second;
		TR << stack << std::endl;
	}

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_res_level_energy.wts" );

	TR << "--- START ---" << std::endl;
	(*scorefxn)( pose );
	scorefxn->show( pose );

	TR << "--- ADD to U5 ---" << std::endl;
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, 5 );
	(*scorefxn)( pose );
	scorefxn->show( pose );

	TR << "--- ADD to C6 ---" << std::endl;
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, 6 );
	(*scorefxn)( pose );
	scorefxn->show( pose );

	TR << "--- ADD to G7 ---" << std::endl;
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, 7 );
	(*scorefxn)( pose );
	scorefxn->show( pose );

	TR << "--- ADD to U8 ---" << std::endl;
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, 8 );
	scorefxn->show( pose );
	(*scorefxn)( pose );
	pose.dump_pdb( "test1.pdb" );

	// u
	pose.set_xyz( NamedAtomID( "RPB1", 5),
								get_repl_xyz( pose, 5, make_vector1( " N1 ", " C2 ", " N3 ",  " C4 ", " C5 ", " C6 " ) ) );
	// c
	pose.set_xyz( NamedAtomID( "RPB1", 6),
								get_repl_xyz( pose, 6, make_vector1( " N1 ", " C2 ", " N3 ",  " C4 ", " C5 ", " C6 " ) ) );
	// g
	pose.set_xyz( NamedAtomID( "RPB1", 7),
								get_repl_xyz( pose, 7, make_vector1( " N9 ", " C4 ", " C5 ",  " N7 ", " C8 " ) ) );
	pose.set_xyz( NamedAtomID( "RPB2", 7),
								get_repl_xyz( pose, 7, make_vector1( " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "  ) ) );
	// a
	pose.set_xyz( NamedAtomID( "RPB1", 8),
								get_repl_xyz( pose, 8, make_vector1( " N9 ", " C4 ", " C5 ",  " N7 ", " C8 " ) ) );
	pose.set_xyz( NamedAtomID( "RPB2", 8),
								get_repl_xyz( pose, 8, make_vector1( " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 "  ) ) );

	protocols::farna::print_internal_coords( pose );
	pose.dump_pdb( "test2.pdb" );


	add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, 5 );
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, 6 );
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, 7 );
	add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, 8 );
	pose.dump_pdb( "test3.pdb" );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	block_stack_test();
	protocols::viewer::clear_conformation_viewers();
  exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		core::init::init(argc, argv);
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


