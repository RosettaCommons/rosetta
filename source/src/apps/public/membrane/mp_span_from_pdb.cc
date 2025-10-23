// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/mp_span_from_pdb.cc
///
/// @brief  Write span from pdb
///
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @note   Last Updated: 5/18/15

// App headers
#include <devel/init.hh>

// Project Headers
#include <core/conformation/membrane/Span.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Package Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// utility headers

// C++ Headers
#include <string>

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::membrane;

static basic::Tracer TR( "apps.public.membrane.mp_span_from_pdb" );

////////////////////////////////////////////////////////////////////////////////
// vector show function
template< typename T_ >
void show( utility::vector1< T_ > vector){
	for ( Size i = 1; i <= vector.size(); ++i ) {
		TR << utility::to_string(vector[i]) << ", ";
	}
	TR << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

Pose read_pose() {

	// cry if PDB not given
	if ( ! option[OptionKeys::in::file::s].user() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide PDB file!");
	}

	// read in pose
	Pose pose;
	core::import_pose::pose_from_file( pose, option[OptionKeys::in::file::s].value_string() , core::import_pose::PDB_file);
	TR.Debug << "got pose of length " << pose.size() << std::endl;

	return pose;

}// read pose

////////////////////////////////////////////////////////////////////////////////

Real read_thickness() {

	// set or read in thickness
	Real thickness;
	if ( option[OptionKeys::mp::thickness].user() ) {
		thickness = option[OptionKeys::mp::thickness]();
		TR << "Taking user-defined thickness: " << thickness*2 << std::endl;
	} else {
		thickness = 15;
		TR << "Taking default thickness: " << thickness*2 << std::endl;
	}
	TR.Debug << "got thickness: " << thickness*2 << std::endl;

	return thickness;

}// get thickness

////////////////////////////////////////////////////////////////////////////////

void spanfile_for_each_chain( Pose & pose, Real thickness, std::string spanfile) {

	// split pose into chains
	utility::vector1< PoseOP > split_poses = pose.split_by_chain();

	// loop over chains
	for ( Size i = 1; i <= split_poses.size(); ++i ) {

		// get pose info
		std::pair< utility::vector1< Real >, utility::vector1< Size > > split_pose_info( get_chain_and_z( *split_poses[i] ));
		utility::vector1< Real > split_z_coord( split_pose_info.first );
		utility::vector1< Size > split_chain_info( split_pose_info.second );
		utility::vector1< char > split_secstruct( get_secstruct( *split_poses[i] ) );

		// create SpanningTopology from poses
		SpanningTopologyOP topo_pose( new SpanningTopology( split_z_coord, split_chain_info, split_secstruct, thickness ) );

		// get filename for spanfile for each chain
		std::string split_spanfile( spanfile );
		utility::trim( split_spanfile, ".span" );

		// output filename depends on number of chains
		if ( pose.chain( pose.size() ) > 1 ) {
			std::string chain( split_poses[i]->pdb_info()->chain(i) );
			split_spanfile = split_spanfile + chain + ".span";
		} else {
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
	// for ( Size i = 1; i <= secstruct.size(); ++i ) {
	//  TR << "i: " << i << ", z: " << z_coord[i] << ", chain: " << chain_info[i] << ", SSE: " << secstruct[i] << std::endl;
	// }

	// create SpanningTopology from pose
	SpanningTopologyOP topo_whole_pose( new SpanningTopology( z_coord, chain_info, secstruct, thickness ) );

	// get filename for spanfile containing whole topology info
	std::string pdbfile( option[OptionKeys::in::file::s].value_string() );
	utility::trim( pdbfile );
	std::string spanfile( pdbfile );
	utility::trim( spanfile, ".pdb");
	spanfile = spanfile + ".span";

	// if more than one chain, write SpanningTopology for whole pose
	if ( pose.chain( pose.size() ) > 1 ) {
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
catch (utility::excn::Exception const & e ) {
	e.display();
	return -1;
}
}
