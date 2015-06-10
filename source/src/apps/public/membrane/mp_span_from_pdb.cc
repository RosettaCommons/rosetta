// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		apps/pilot/membrane/mp_span_from_pdb.cc
///
/// @brief		Write span from pdb
///
/// @author		Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @note 		Last Updated: 5/18/15

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/geometry/util.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::membrane::geometry;

static thread_local basic::Tracer TR( "apps.public.membrane.mp_span_from_pdb" );

////////////////////////////////////////////////////////////////////////////////
// vector show function
template< typename T_ >
void show( utility::vector1< T_ > vector){
	for ( Size i = 1; i <= vector.size(); ++i ){
		TR << utility::to_string(vector[i]) << ", ";
	}
	TR << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

Pose read_pose() {
	
	// cry if PDB not given
	if ( ! option[OptionKeys::in::file::s].user() ){
		throw new utility::excn::EXCN_Msg_Exception("Please provide PDB file!");
	}
	
	// read in pose
	Pose pose;
	core::import_pose::pose_from_pdb( pose, option[OptionKeys::in::file::s].value_string() );
	TR.Debug << "got pose of length " << pose.total_residue() << std::endl;
	
	return pose;
	
}// read pose

////////////////////////////////////////////////////////////////////////////////

Real read_thickness() {

	// set or read in thickness
	Real thickness;
	if ( option[OptionKeys::mp::thickness].user() ){
		thickness = option[OptionKeys::mp::thickness]();
		TR << "Taking user-defined thickness: " << thickness << std::endl;
	}
	else{
		thickness = 15;
		TR << "Taking default thickness: " << thickness << std::endl;
	}
	TR.Debug << "got thickness: " << thickness << std::endl;

	return thickness;
	
}// get thickness

////////////////////////////////////////////////////////////////////////////////

void spanfile_for_each_chain( Pose & pose, Real thickness, std::string spanfile) {
	
	// split pose into chains
	utility::vector1< PoseOP > split_poses = pose.split_by_chain();
	
	// loop over chains
	for ( Size i = 1; i <= split_poses.size(); ++i ){
		
		// get pose info
		std::pair< utility::vector1< Real >, utility::vector1< Size > > split_pose_info( get_chain_and_z( *split_poses[i] ));
		utility::vector1< Real > split_z_coord( split_pose_info.first );
		utility::vector1< Size > split_chain_info( split_pose_info.second );
		utility::vector1< char > split_secstruct( get_secstruct( *split_poses[i] ) );

		// create SpanningTopology from poses
		SpanningTopologyOP topo_pose( new SpanningTopology( split_z_coord, split_chain_info, split_secstruct, thickness ) );
		
		// get filename for spanfile for each chain
		char chain( split_poses[i]->pdb_info()->chain(i) );
		std::string split_spanfile( spanfile );
		utility::trim( split_spanfile, ".span" );
		
		// output filename depends on number of chains
		if ( pose.chain( pose.total_residue() ) > 1 ){
			split_spanfile = split_spanfile + chain + ".span";
		}
		else{
			split_spanfile = split_spanfile + ".span";
		}
		
		// print SpanningTopology for poses
		topo_pose->write_spanfile( split_spanfile );
	}

}// spanfile for each chain

////////////////////////////////////////////////////////////////////////////////

void spanfile_from_pdb(){

	TR << "spanfile_from_pdb" << std::endl;

	// read input
	Pose pose = read_pose();
	Real thickness = read_thickness();
	
	// get pose info
	std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( get_chain_and_z( pose ));
	utility::vector1< Real > z_coord = pose_info.first;
	utility::vector1< Size > chain_info = pose_info.second;
	utility::vector1< char > secstruct = get_secstruct( pose );

	// for debugging
//	for ( Size i = 1; i <= secstruct.size(); ++i ) {
//		TR << "i: " << i << ", z: " << z_coord[i] << ", chain: " << chain_info[i] << ", SSE: " << secstruct[i] << std::endl;
//	}

	// create SpanningTopology from pose
	SpanningTopologyOP topo_whole_pose( new SpanningTopology( z_coord, chain_info, secstruct, thickness ) );
	
	// get filename for spanfile containing whole topology info
	std::string pdbfile( option[OptionKeys::in::file::s].value_string() );
	utility::trim( pdbfile );
	std::string spanfile( pdbfile );
	utility::trim( spanfile, ".pdb");
	spanfile = spanfile + ".span";

	// if more than one chain, write SpanningTopology for whole pose
	if ( pose.chain( pose.total_residue() ) > 1 ){
		topo_whole_pose->write_spanfile( spanfile );
	}

	// write spanfile for each individual chain
	spanfile_for_each_chain( pose, thickness, spanfile );
	
}// spanfile_from_pdb

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		// initialize option system, RNG, and all factory-registrators
		devel::init(argc, argv);

		// call my function
		spanfile_from_pdb();
		
	}
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
