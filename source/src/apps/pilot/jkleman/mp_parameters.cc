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

#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <numeric/xyzVector.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>
#include <string>

using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "apps.pilot.jkleman.mp_parameters" );

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
bool no_interpolate_Mpair(false);
bool Hbond_depth_correction(false);

// global penalties
bool TMprojection(false);
bool non_helix(false);
bool termini(false);
bool secstruct(false);
bool spanning(false);
Real wt_TMprojection;
Real wt_non_helix;
Real wt_termini;
Real wt_secstruct;
Real wt_spanning;

///////////////////////////////////////////////////////////////////////////////
// read membrane file
void read_membrane (){

	// create filename
	TR << "Reading in membrane parameters from database." << std::endl;
	std::string membrane_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/membrane/membrane.mp");

	if (! utility::file::file_exists(membrane_file)){
		membrane_file = "/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/membrane/membrane.mp";
    }

	// create stream
    utility::io::izstream stream (membrane_file);

    // fail if file not found
	if (stream.fail()){
		utility_exit_with_message("Cannot find database file " + membrane_file);
	}

	// read lines
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
	TR << "Reading in embedding parameters from database." << std::endl;
	std::string embedding_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/membrane/embeddings.mp");
    
	if (! utility::file::file_exists(embedding_file)){
		embedding_file = "/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/membrane/embeddings.mp";
    }
    
	// create stream
    utility::io::izstream stream (embedding_file);
    
    // fail if file not found
	if (stream.fail()){
		utility_exit_with_message("Cannot find database file " + embedding_file);
	}

	// read lines
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
	std::string parameter_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/scoring/weights/membrane_parameters.patch");
    
	if (! utility::file::file_exists(parameter_file)){
		parameter_file = "/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/scoring/weights/membrane_parameters.patch";
    }
    
	// create stream
    utility::io::izstream stream (parameter_file);
    
    // fail if file not found
	if (stream.fail()){
		utility_exit_with_message("Cannot find database file " + parameter_file);
	}
	
	// read lines
	std::string line;
	
	while (getline(stream, line)){
		
		// split up line into words
		Real x, y, z;
		std::string tag;
		std::istringstream words(line);
		words >> tag >> x >> y >> z;
		
		// assign no_interpolate_Mpair
		if (tag == "no_interpolate_Mpair"){
			x == 1 ? no_interpolate_Mpair = true : no_interpolate_Mpair = false;
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
	std::string penalty_file ("/Users/julialeman/Documents/julia/git_final/Rosetta/main/database/scoring/weights/membrane_penalties.patch");
    
	if (! utility::file::file_exists(penalty_file)){
		penalty_file = "/Users/julialeman/Postdoc/Rosetta/Rosetta/main/database/scoring/weights/membrane_penalties.patch";
    }
    
	// create stream
    utility::io::izstream stream (penalty_file);
    
    // fail if file not found
	if (stream.fail()){
		utility_exit_with_message("Cannot find database file " + penalty_file);
	}
	
	// read lines
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

		// assign secondary structure
		if (tag == "secstruct"){
			secstruct = x;
			wt_secstruct = y;
			TR << "secstruct: " << secstruct << std::endl;
			TR << "weight for secstruct: " << wt_secstruct << std::endl;
		}

        // assign spanning
		if (tag == "spanning"){
			spanning = x;
			wt_spanning = y;
			TR << "spanning: " << spanning << std::endl;
			TR << "weight for spanning: " << wt_spanning << std::endl;
		}
}
}

///////////////////////////////////////////////////////////////////////////////
// read options
void read_options(){

	// read in membrane.mp
	read_membrane();
	
	// read in embeddings.mp
	read_embeddings();

	// read in membrane_parameters.patch
	read_parameters();
  
	// read in membrane_penalties.patch
	read_penalties();
}


///////////////////////////////////////////////////////////////////////////////
// checking options
void user_overwrites_options(){

	// namespace
	using namespace basic::options;
    using namespace basic::options::OptionKeys;
	
	// read user-defined options
    // read membrane options
    if ( option[ mp::thickness ].user() ){
        thickness = option[OptionKeys::mp::thickness];
    }
    
    // read embedding options
    if ( option[ mp::center_start ].user()){
        center_start.x() = option[OptionKeys::mp::center_start]()[1];
        center_start.y() = option[OptionKeys::mp::center_start]()[2];
        center_start.z() = option[OptionKeys::mp::center_start]()[3];
    }
    
    if ( option[ mp::center_delta ].user()){
        center_delta = option[ OptionKeys::mp::center_delta ];
    }

    if ( option[ mp::center_search_cycles ].user()){
        center_search_cycles = option[ OptionKeys::mp::center_search_cycles ];
    }

    if ( option[ mp::normal_start ].user()){
        normal_start.x() = option[ OptionKeys::mp::normal_start ]()[1];
        normal_start.y() = option[ OptionKeys::mp::normal_start ]()[2];
        normal_start.z() = option[ OptionKeys::mp::normal_start ]()[3];
    }

    if ( option[ mp::normal_angle_start ].user()){
        normal_angle_start = option[ OptionKeys::mp::normal_angle_start ]();
    }

    if ( option[ mp::normal_angle_delta ].user()){
        normal_angle_delta = option[ OptionKeys::mp::normal_angle_delta ]();
    }

    if ( option[ mp::normal_search_cycles ].user()){
        normal_search_cycles = option[ OptionKeys::mp::normal_search_cycles ]();
    }

    if ( option[ mp::chain_normal_angle_max ].user()){
        chain_normal_angle_max = option[ OptionKeys::mp::chain_normal_angle_max ]();
    }

    if ( option[ mp::pose_normal_angle_max ].user()){
        pose_normal_angle_max = option[ OptionKeys::mp::pose_normal_angle_max ]();
    }

    // read scoring options
    if ( option[ mp::no_interpolate_Mpair ].user()){
        no_interpolate_Mpair = option[ OptionKeys::mp::no_interpolate_Mpair ]();
    }
    
    if ( option[ mp::Hbond_depth_correction ].user()){
        Hbond_depth_correction = option[ OptionKeys::mp::Hbond_depth_correction ]();
    }
    
    // read penalty options
    if ( option[ mp::TMprojection ].user()){
        TMprojection = option[ OptionKeys::mp::TMprojection ]();
    }

    if ( option[ mp::wt_TMprojection ].user()){
        wt_TMprojection = option[ OptionKeys::mp::wt_TMprojection ]();
    }

    if ( option[ mp::non_helix ].user()){
        non_helix = option[ OptionKeys::mp::non_helix ]();
    }

    if ( option[ mp::wt_non_helix ].user()){
        wt_non_helix = option[ OptionKeys::mp::wt_non_helix ]();
    }

    if ( option[ mp::termini ].user()){
        termini = option[ OptionKeys::mp::termini ]();
    }

    if ( option[ mp::wt_termini ].user()){
        wt_termini = option[ OptionKeys::mp::wt_termini ]();
    }
    
    if ( option[ mp::secstruct ].user() ){
        secstruct = option[ OptionKeys::mp::secstruct ]();
    }

    if ( option[ mp::wt_secstruct ].user() ){
        wt_secstruct = option[ OptionKeys::mp::wt_secstruct ]();
    }

    if ( option[ mp::spanning ].user() ){
        spanning = option[ OptionKeys::mp::spanning ]();
    }

    if ( option[ mp::wt_spanning ].user() ){
        wt_spanning = option[ OptionKeys::mp::wt_spanning ]();
    }
}
///////////////////////////////////////////////////////////////////////////////
// checking options
void check_options(){

	using namespace utility;

	// checking options for membrane
	if (center.to_string() != "(0, 0, 0)") {
		utility_exit_with_message("Center should be (0, 0, 0), but is " + center.to_string() + "!");
	}

	if (normal.to_string() != "(0, 0, 1)") {
		utility_exit_with_message("Normal should be (0, 0, 1), but is " + normal.to_string() + "!");
	}

	// checking options for embeddings
	if (center_start.to_string() != "(0, 0, 0)") {
		utility_exit_with_message("Center_start should be (0, 0, 0), but is " + to_string(center_start.to_string()) + "!");
	}

	if (center_delta < 0 || center_delta > 45) {
		utility_exit_with_message("Center_delta should be between 0 and 45 degrees, but is " + to_string(center_delta) + "!");
	}

	if (center_search_cycles < 1 || center_search_cycles > 100) {
		utility_exit_with_message("Center_search_cycles should be between 1 and 100, but is " + to_string(center_search_cycles) + "!");
	}

	if (normal_start.to_string() != "(0, 0, 1)") {
		utility_exit_with_message("Normal_start should be (0, 0, 1), but is " + normal_start.to_string() + "!");
	}

	if (normal_angle_start < 0 || normal_angle_start > 90) {
		utility_exit_with_message("Normal_angle_start should be between 0 and 90 degrees, but is " + to_string(normal_angle_start) + "!");
	}

	if (normal_angle_delta < 0 || normal_angle_delta > 45) {
		utility_exit_with_message("Normal_angle_delta should be between 0 and 45 degrees, but is " + to_string(normal_angle_delta) + "!");
	}

	if (normal_search_cycles < 1 || normal_search_cycles > 100) {
		utility_exit_with_message("Normal_search_cycles should be between 1 and 100, but is " + to_string(normal_search_cycles) + "!");
	}

	if (chain_normal_angle_max < 0 || chain_normal_angle_max > 90) {
		utility_exit_with_message("Chain_normal_angle_max should be between 0 and 90 degrees, but is " + to_string(chain_normal_angle_max) + "!");
	}

	if (pose_normal_angle_max < 0 || pose_normal_angle_max > 45) {
		utility_exit_with_message("Pose_normal_angle_max should be between 0 and 45 degrees, but is " + to_string(pose_normal_angle_max) + "!");
	}
	
	// checking parameters in scoring
	if (! no_interpolate_Mpair){
		utility_exit_with_message("No_interpolate_Mpair should be 1, but is " + to_string(no_interpolate_Mpair) + "!");
	}

	if (Hbond_depth_correction == 0){
		utility_exit_with_message("Hbond_depth_correction should be 1, but is " + to_string(Hbond_depth_correction) + "!");
	}

	// checking penalties in scoring
	if (TMprojection != 1){
		utility_exit_with_message("TMprojection should be 1, but is " + to_string(TMprojection) + "!");
	}

	if (wt_TMprojection < 1 || wt_TMprojection > 1000){
		utility_exit_with_message("Wt_TMprojection should be between 1 and 1000, but is " + to_string(wt_TMprojection) + "!");
	}

	if (non_helix != 1){
		utility_exit_with_message("Non_helix should be 1, but is " + to_string(non_helix) + "!");
	}
	
	if (wt_non_helix < 1 || wt_non_helix > 1000){
		utility_exit_with_message("Wt_non_helix should be between 1 and 1000, but is " + to_string(wt_non_helix) + "!");
	}

	if (termini != 1){
		utility_exit_with_message("Termini should be 1, but is " + to_string(termini) + "!");
	}
	
	if (wt_termini < 1 || wt_termini > 1000){
		utility_exit_with_message("Wt_termini should be between 1 and 1000, but is " + to_string(wt_termini) + "!");
	}

	if (secstruct != 1){
		utility_exit_with_message("Secstruct should be 1, but is " + to_string(secstruct) + "!");
	}
	
	if (wt_secstruct < 1 || wt_secstruct > 1000){
		utility_exit_with_message("Wt_secstruct should be between 1 and 1000, but is " + to_string(wt_secstruct) + "!");
	}

	if (spanning != 1){
		utility_exit_with_message("Spanning should be 1, but is " + to_string(spanning) + "!");
	}
	
	if (wt_spanning < 1 || wt_spanning > 1000){
		utility_exit_with_message("Wt_spanning should be between 1 and 1000, but is " + to_string(wt_spanning) + "!");
	}
}

///////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		
		devel::init(argc, argv);
		read_options();
        user_overwrites_options();
		check_options();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
}

