// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mp_parameters.cc
/// @brief   Checking options for MPframework
/// @details Checking options from database and user input for MPframework
///			     Last Modified: 3/24/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/MembraneUnitTestMover.hh>

// Package Headers
#include <core/types.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.pilot.jkleman.mp_parameters" );

using namespace core;

// global variables for membrane
Vector normal;
Vector center;
Real thickness;

// global variables for embedding
Vector center_start;
Real center_delta;
Real center_search_cycles;
Vector normal_start;
Real normal_angle_start;
Real normal_angle_delta;
Real normal_search_cycles;
Real chain_normal_angle_max;
Real pose_normal_angle_max;

// global parameters
bool no_interpolate_Mpair;
bool Hbond_depth_correction;

// global penalties
bool TMprojection;
bool non_helix;
bool termini;
Real wt_TMprojection;
Real wt_non_helix;
Real wt_termini;

///////////////////////////////////////////////////////////////////////////////
// read membrane file
void read_membrane (){

	// create filename
	TR << "Reading in membrane parameters from database." << std::endl;
	std::string membrane_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/membrane/membrane.mp");
	utility::io::izstream stream (membrane_file);

	// read stream
	if (! stream.good()){
		std::string membrane_file ("/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/membrane/membrane.mp");
		utility::io::izstream stream (membrane_file);	}
	if (! stream.good()){
		utility_exit_with_message("Cannot find database file " + membrane_file);
	}

	// read lines
	TR << "reading lines..." << std::endl;
	std::string line;

	while (getline(stream, line)){

		// split up line into words
		Real x, y, z;
		std::string tag;
		std::istringstream words(line);
		words >> tag >> x >> y >> z;

		// assign center
		if (tag == "center"){
				center.assign(x, y, z);
				TR << "center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;
		}
		
		// assign normal
		if (tag == "normal"){
				normal.assign(x, y, z);
				TR << "normal: " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
		}
		
		// assign thickness
		if (tag == "thickness"){
				thickness = x;
				TR << "thickness: " << thickness << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// read embeddings files
void read_embeddings (){
	
	// create filename
	TR << "Reading in embeddings parameters from database." << std::endl;
	std::string embeddings_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/membrane/embeddings.mp");
	utility::io::izstream stream (embeddings_file);
	
	// read stream
	if (! stream.good()){
		std::string embeddings_file ("/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/membrane/embeddings.mp");
		utility::io::izstream stream (embeddings_file);	}
	if (! stream.good()){
		utility_exit_with_message("Cannot find database file " + embeddings_file);
	}
	
	// read lines
	TR << "reading lines..." << std::endl;
	std::string line;
	
	while (getline(stream, line)){
		
		// split up line into words
		Real x, y, z;
		std::string tag;
		std::istringstream words(line);
		words >> tag >> x >> y >> z;
		
		// assign center_start
		if (tag == "center_start"){
			center_start.assign(x, y, z);
			TR << "center_start: " << center_start.x() << " " << center_start.y() << " " << center_start.z() << std::endl;
		}
		
		// assign center_delta
		if (tag == "center_delta"){
			center_delta = x;
			TR << "center_delta: " << center_delta << std::endl;
		}

		// assign center_search_cycles
		if (tag == "center_search_cycles"){
			center_search_cycles = x;
			TR << "center_search_cycles: " << center_search_cycles << std::endl;
		}

		// assign normal_start
		if (tag == "normal_start"){
			normal_start.assign(x, y, z);
			TR << "normal_start: " << normal_start.x() << " " << normal_start.y() << " " << normal_start.z() << std::endl;
		}

		// assign normal_angle_start
		if (tag == "normal_angle_start"){
			normal_angle_start = x;
			TR << "normal_angle_start: " << normal_angle_start << std::endl;
		}
		
		// assign normal_angle_delta
		if (tag == "normal_angle_delta"){
			normal_angle_delta = x;
			TR << "normal_angle_delta: " << normal_angle_delta << std::endl;
		}

		// assign normal_search_cycles
		if (tag == "normal_search_cycles"){
			normal_search_cycles = x;
			TR << "normal_search_cycles: " << normal_search_cycles << std::endl;
		}

		// assign chain_normal_angle_max
		if (tag == "chain_normal_angle_max"){
			chain_normal_angle_max = x;
			TR << "chain_normal_angle_max: " << chain_normal_angle_max << std::endl;
		}

		// assign pose_normal_angle_max
		if (tag == "pose_normal_angle_max"){
			pose_normal_angle_max = x;
			TR << "pose_normal_angle_max: " << pose_normal_angle_max << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// read parameter file
void read_parameters (){
	
	// create filename
	TR << "Reading in membrane parameters from database/scoring." << std::endl;
	std::string parameters_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/scoring/weights/membrane_parameters.patch");
	utility::io::izstream stream (parameters_file);
	
	// read stream
	if (! stream.good()){
		std::string parameters_file ("/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/scoring/weights/membrane_parameters.patch");
		utility::io::izstream stream (parameters_file);	}
	if (! stream.good()){
		utility_exit_with_message("Cannot find database file " + parameters_file);
	}
	
	// read lines
	TR << "reading lines..." << std::endl;
	std::string line;
	
	while (getline(stream, line)){
		
		// split up line into words
		Real x, y, z;
		std::string tag;
		std::istringstream words(line);
		words >> tag >> x >> y >> z;
		
		// assign no_interpolate_pair
		if (tag == "no_interpolate_Mpair"){
			no_interpolate_Mpair = x;
			TR << "no_interpolate_Mpair: " << no_interpolate_Mpair << std::endl;
		}

		// assign Hbond_depth_correction
		if (tag == "Hbond_depth_correction"){
			Hbond_depth_correction = x;
			TR << "Hbond_depth_correction: " << Hbond_depth_correction << std::endl;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// read penalties file
void read_penalties (){
	
	// create filename
	TR << "Reading in membrane penalties from database/scoring." << std::endl;
	std::string penalties_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/scoring/weights/membrane_penalties.patch");
	utility::io::izstream stream (penalties_file);
	
	// read stream
	if (! stream.good()){
		std::string parameters_file ("/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/scoring/weights/membrane_penalties.patch");
		utility::io::izstream stream (penalties_file);	}
	if (! stream.good()){
		utility_exit_with_message("Cannot find database file " + penalties_file);
	}
	
	// read lines
	TR << "reading lines..." << std::endl;
	std::string line;
	
	while (getline(stream, line)){
		
		// split up line into words
		Real x, y, z;
		std::string tag;
		std::istringstream words(line);
		words >> tag >> x >> y >> z;
		
		// assign TMprojection
		if (tag == "TMprojection"){
			TMprojection = x;
			wt_TMprojection = y;
			TR << "TMprojection: " << TMprojection << std::endl;
			TR << "weight for TMprojection: " << wt_TMprojection << std::endl;
		}
		
		// assign non-helix
		if (tag == "non-helix"){
			non_helix = x;
			wt_non_helix = y;
			TR << "non-helix: " << non_helix << std::endl;
			TR << "weight for non-helix: " << wt_non_helix << std::endl;
		}

		// assign termini
		if (tag == "termini"){
			termini = x;
			wt_termini = y;
			TR << "termini: " << termini << std::endl;
			TR << "weight for termini: " << wt_termini << std::endl;
		}
}
}

///////////////////////////////////////////////////////////////////////////////
// check options
void checking_options(){

	// read in membrane.mp
	read_membrane();
	
	// read in embeddings.mp
	read_embeddings();

	// read in membrane_parameters.patch
	read_parameters();
  
	// read in membrane_penalties.patch
	read_penalties();
	
	// checking options for membrane
	if (! (center.x() == 0 && center.y() == 0 && center.z() == 0)) {
		utility_exit_with_message("Center should be (0, 0, 0), but is (" + center.show());
//	utility_exit_with_message("Center should be (0, 0, 0), but is (" + string_of(center.x()) + ", " + string_of(center.y()) + ", " + string_of(center.z()) + ").");
	}
/*
	if (! (normal.x() == 0 && normal.y() == 0 && normal.z() == 1))
	utility_exit_with_message("Normal should be (0, 0, 1), but is (" + normal.x() + ", " + normal.y() + ", " + normal.z() + ").");
	}
*/
//	return NULL;
}

///////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		
		devel::init(argc, argv);
		checking_options();
//		protocols::viewer::viewer_main( checking_options() );
		
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
}

