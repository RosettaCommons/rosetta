// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file topo_driver.cc
/// @brief  Sandbox executable for the Environment and associated classes.
/// @author Justin Porter

#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/abinitio/abscript/RigidChunkCM.hh>

#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/pack/task/residue_selector/ResidueIndexSelector.hh>

#include <devel/init.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <numeric/xyzVector.hh>

#include <core/id/types.hh>

#include <boost/bind/bind.hpp>

#define TS_ASSERT( x ) assert( x )
#define TS_ASSERT_THROWS_NOTHING( x ) x
#define TS_ASSERT_EQUALS( x , y ) TS_ASSERT( x == y )
#define TS_ASSERT_DIFFERS( x , y ) TS_ASSERT( x != y )
#define TS_ASSERT_DELTA( x , y, d ) TS_ASSERT( std::abs( x - y ) < d )
#define TS_ASSERT_THROWS( x , y ) try{ x; } catch( y ){}
#define TS_ASSERT_LESS_THAN( x, y ) TS_ASSERT( x < y );
#define TS_TRACE( x ) std::cout << x << std::endl;
#define ANGLE_DELTA( x , y , d ) TS_ASSERT_DELTA( std::cos( x ) , std::cos( y ), d );

std::string const TEST_CHUNK_PDB = "/Users/jrporter/Rosetta/main/source/test/protocols/abinitio/abscript/ubq_frag.pdb";
core::Real const TOLERANCE = 1e-6;

void test_chunk( core::pose::Pose const pose );

core::pose::Pose make_seq( std::string const& str );

core::pose::Pose make_seq( std::string const& str ) {
  core::pose::Pose pose;

  core::pose::make_pose_from_sequence(pose, str, "fa_standard");

  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    pose.set_phi( i, -65 );
    pose.set_psi( i, -41 );
    pose.set_omega( i, 180 );
  }

  return pose;
}

void test_chunk( core::pose::Pose const ) {

		TS_TRACE( "Starting test_chunk" );

		using namespace protocols::environment;
		using namespace core::pack::task::residue_selector;
		using namespace protocols::abinitio::abscript;

		basic::datacache::DataMap datamap;

    core::pose::Pose pose = make_seq( "MQIFV" );

		std::string const SELECTOR_NAME = "sele";
		core::Size const BEGIN_CHUNK = 1;
		core::Size const END_CHUNK = 5;

		std::ostringstream os;
		os << BEGIN_CHUNK << "-" << END_CHUNK;
		datamap.add( "ResidueSelector", SELECTOR_NAME, ResidueSelectorOP( new ResidueIndexSelector( os.str() ) ) );

		std::stringstream ss;
		ss << "<RigidChunkCM name=chunk template=\"" << TEST_CHUNK_PDB << "\" selector=\"" << SELECTOR_NAME << "\" />";
		utility::tag::TagPtr tag( new utility::tag::Tag );
		tag->read( ss );

		RigidChunkCMOP rigid_chunk( new RigidChunkCM() );
		TS_ASSERT_THROWS_NOTHING( rigid_chunk->parse_my_tag( tag , datamap, protocols::filters::Filters_map(), protocols::moves::Movers_map(), core::pose::Pose() ) );
		TS_ASSERT_THROWS_NOTHING( tag->die_for_unaccessed_options() );

		EnvironmentOP env( new Environment( "env" ) );
		env->register_mover( rigid_chunk );

		core::pose::Pose ppose;
		TS_ASSERT_THROWS_NOTHING( ppose = env->start( pose ) );

		for( core::Size i = 1; i <= ppose.total_residue(); ++i ){
      if( i >= BEGIN_CHUNK && i <= END_CHUNK ){
        core::Size const seqpos = i - BEGIN_CHUNK + 1;
        if( seqpos != 1 )
          ANGLE_DELTA( ppose.phi( i ), rigid_chunk->templ().phi( seqpos ), TOLERANCE );
        if( seqpos != rigid_chunk->templ().total_residue() &&
           seqpos != ppose.total_residue() ){
          ANGLE_DELTA( ppose.psi( i ), rigid_chunk->templ().psi( seqpos ), TOLERANCE );
          ANGLE_DELTA( ppose.omega( i ) , rigid_chunk->templ().omega( seqpos ), TOLERANCE );
        }
      } else if( ( i == BEGIN_CHUNK - 1 ) ||
                ( i == END_CHUNK + 1 ) ) {
        // inconsistent contributions.
      } else {
        if( i != 1 ) ANGLE_DELTA( ppose.phi( i ), pose.phi( i ), TOLERANCE );
        if( i != pose.total_residue() ) {
          ANGLE_DELTA( ppose.psi( i ), pose.psi( i ), TOLERANCE );
          ANGLE_DELTA( ppose.omega( i ), pose.omega( i ), TOLERANCE );
        }
      }
    }
  
		core::pose::Pose end_pose;
		TS_ASSERT_THROWS_NOTHING( end_pose = env->end( ppose ) );

}

int main( int argc, char** argv ){
  devel::init( argc, argv );

  core::pose::Pose pose;
  core::pose::make_pose_from_sequence(pose, "FRMQIFVYFRIENDS", "fa_standard");

  for( core::Size i = 1; i <= pose.total_residue(); ++i ){
    pose.set_phi( i, -65 );
    pose.set_psi( i, -41 );
    pose.set_omega( i, 180 );
  }

//  pose.append_pose_by_jump( *pose.clone(), 5 );
//  core::kinematics::Jump j = pose.jump( 1 );
//  j.gaussian_move(1, 50, 0);
//  pose.set_jump(1, j);

  try {
    test_chunk( pose );
  } catch ( utility::excn::EXCN_Msg_Exception excn ){
    std::cout << excn << std::endl;
    std::exit( 1 );
  }

  std::cout << "Success!" << std::endl;

}
