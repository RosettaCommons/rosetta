// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mpdocking.cc
/// @brief   Dock two membrane proteins in the membrane
/// @details last Modified: 4/4/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/membrane/geometry/util.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>
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
using namespace core::membrane::geometry;

static basic::Tracer TR( "apps.pilot.jkleman.spanfile_from_pdb" );

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

void print_spanfile( SpanningTopologyOP topology, Size nres, std::string pdbfile, std::string spanfile ){

	TR.Debug << "printing spanfile" << std::endl;

	// print header
	utility::io::ozstream OUT;
	OUT.open( spanfile );
	OUT << "spanfile from PDB " << pdbfile << std::endl;
	OUT << topology->total_spans() << " " << nres << std::endl;
	OUT << "antiparallel" << std::endl;
	OUT << "n2c" << std::endl;
	
	// print spans
	for ( Size i = 1; i <= topology->total_spans(); ++i ){
		OUT << "\t" << topology->span(i)->start() << "\t" << topology->span(i)->end() << std::endl;
	}
	OUT.close();
	TR << "wrote " << spanfile << std::endl;
}// print spanfile

////////////////////////////////////////////////////////////////////////////////

void spanfile_from_pdb(){

	TR << "spanfile_from_pdb" << std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	// cry if PDB not given
	if ( ! option[OptionKeys::in::file::s].user() ){
		throw new utility::excn::EXCN_Msg_Exception("Please provide PDB file!");
	}
	
	// read in pose
	PoseOP pose = core::import_pose::pose_from_pdb(
							option[OptionKeys::in::file::s].value_string() );
	TR.Debug << "got pose of length " << pose->total_residue() << std::endl;
	
	// set or read in thickness
	Real thickness;
	if ( option[OptionKeys::membrane_new::thickness].user() ){
		thickness = option[OptionKeys::membrane_new::thickness]();
		TR << "Taking user-defined thickness: " << thickness << std::endl;
	}
	else{
		thickness = 12.5;
		TR << "Taking default thickness: " << thickness << std::endl;
	}
	TR.Debug << "got thickness: " << thickness << std::endl;

	// get pose info
	std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( get_chain_and_z( pose ));
	utility::vector1< Real > z_coord = pose_info.first;
	utility::vector1< Size > chain_info = pose_info.second;
	TR.Debug << "got pose info:" << std::endl;

	// create SpanningTopology from pose
	SpanningTopologyOP topo_whole_pose = new SpanningTopology( z_coord, chain_info, thickness );
	TR.Debug << "got whole pose topology: " << std::endl;
//	topo_whole_pose->show();
	
	// get filename
	std::string pdbfile( option[OptionKeys::in::file::s].value_string() );
	utility::trim( pdbfile );
	std::string spanfile( pdbfile );
	utility::trim( spanfile, ".pdb");
	spanfile = spanfile + ".span";

	// if more than one chain, print SpanningTopology for whole pose
	if ( pose->chain( pose->total_residue() ) > 1 ){
		topo_whole_pose->write_spanfile( spanfile );
	}

	// split pose into chains
	utility::vector1< PoseOP > split_poses = pose->split_by_chain();

	TR << "split pose into " << split_poses.size() << " poses." << std::endl;

	// loop over chains
	for ( Size i = 1; i <= split_poses.size(); ++i ){

		TR.Debug << "going over chains: " << i << std::endl;

		// get pose info
		std::pair< utility::vector1< Real >, utility::vector1< Size > > split_pose_info( get_chain_and_z( split_poses[i] ));
		utility::vector1< Real > split_z_coord( split_pose_info.first );
		utility::vector1< Size > split_chain_info( split_pose_info.second );

		// create SpanningTopology from poses
		SpanningTopologyOP topo_pose = new SpanningTopology( split_z_coord, split_chain_info, thickness );
		
		// get filename
		char chain( split_poses[i]->pdb_info()->chain(i) );
		std::string split_spanfile( spanfile );
		utility::trim( split_spanfile, ".span" );

		// output filename depends on number of chains
		if ( pose->chain( pose->total_residue() ) > 1 ){
			split_spanfile = split_spanfile + chain + ".span";
		}
		else{
			split_spanfile = split_spanfile + ".span";
		}
		
		// print SpanningTopology for poses
		topo_pose->write_spanfile( split_spanfile );
	}
	
}// spanfile_from_pdb

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		// initialize option system, random number generators, and all factory-registrators
		devel::init(argc, argv);
		//protocols::init(argc, argv);

		// call my function
		spanfile_from_pdb();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
