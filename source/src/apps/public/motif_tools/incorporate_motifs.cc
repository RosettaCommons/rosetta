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

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

//#include <core/graph/Graph.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <protocols/loops/Loops.hh>

#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/SingleMotif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/IRCollection.hh>

//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
//#include <string>
//#include <algorithm>

using namespace core;
using namespace pose;
using namespace protocols::motifs;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace basic::options::OptionKeys::motifs;

static THREAD_LOCAL basic::Tracer TR( "incorporate_motifs" );

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void
identify_targeted_positions( utility::vector1< core::Size > & pos_vec, core::pose::Pose & pose, std::string & source_file )
{

	pos_vec.clear();

	// Probe for file
	std::ifstream pos_file( source_file.c_str() );

	if ( !pos_file ) {
		TR << "Requested file " << source_file << " not found." << std::endl;
		return;
	}

	core::Size pos;
	char chain;

	pos_file >> pos;
	while ( !pos_file.eof() ) {
		pos_file >> chain;
		TR << "Reading position " << pos << " chain " << chain << " from " << source_file << std::endl;
		pos_vec.push_back( pose.pdb_info()->pdb2pose( chain, pos ) );
		pos_file >> pos;
	}
	pos_file.close();

	return;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void
read_in_flexible_regions( protocols::loops::LoopsOP & flex_regions, core::pose::Pose & pose )
{

	if ( option[ motif_flexible_loop_file].active() ) {
		std::string flex_file_name( option[ motif_flexible_loop_file ] );

		// Probe for file
		std::ifstream regions_file( flex_file_name.c_str() );

		if ( !regions_file ) {
			TR << "motif_flexible_loop_file " << flex_file_name << " not found." << std::endl;
			return;
		}

		core::Size start_pos;
		core::Size end_pos;
		char flex_chain;

		regions_file >> start_pos;
		regions_file >> end_pos;
		regions_file >> flex_chain;
		TR << "Adding flexible region from  " << start_pos << " to " << end_pos << " in chain " << flex_chain << std::endl;
		TR << "In Rosetta numbering, from  " << pose.pdb_info()->pdb2pose( flex_chain, start_pos ) <<
			" to " << pose.pdb_info()->pdb2pose( flex_chain, end_pos ) << std::endl;
		flex_regions->add_loop( pose.pdb_info()->pdb2pose( flex_chain, start_pos ), pose.pdb_info()->pdb2pose( flex_chain, end_pos ) );
		regions_file.close();
	} else {
		TR << "Error:  no specified value for motif_flexible_loop_file!" << std::endl;
	}

	return;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
find_close_motifs()
{
	using namespace conformation;
	using namespace import_pose;

	//using namespace optimization;

	Pose pose;
	std::string template_structure_file( option[ in::file::s ]()[1] );
	pose_from_pdb( pose,  template_structure_file.c_str() );

	std::string weights( "talaris2013" );
	ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( weights ) );

	MotifLibrary motif_lib;

	std::string motif_filename( option[ motifs::motif_filename ] );
	motif_lib.add_from_file( motif_filename );

	// Read in the positions to target with motifs
	utility::vector1< core::Size > target_positions;
	if ( option[ build_residue_file ].active() ) {
		std::string build_file_name( option[ build_residue_file ] );
		identify_targeted_positions( target_positions, pose, build_file_name );
	}

	// Read in the positions to trim before
	utility::vector1< core::Size > trim_positions;
	if ( option[ residue_trim_file ].active() ) {
		std::string trim_file_name( option[ residue_trim_file ] );
		identify_targeted_positions( trim_positions, pose, trim_file_name );
	}

	IRCollection inv_rotamers( pose, motif_lib, target_positions );

	// Use the Loop facility to define flexible segments
	protocols::loops::LoopsOP flexible_regions( new protocols::loops::Loops );
	//protocols::loops::Loops( option[ motif_flexible_loop_file ] );
	read_in_flexible_regions( flexible_regions, pose );

	TR << "Entering motif incorporation" << std::endl;

	inv_rotamers.incorporate_motifs( pose, flexible_regions, trim_positions );

	return;

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );

		find_close_motifs();
		TR.flush();
		return 0;
	} catch ( utility::excn::EXCN_Base const & e ) {
		TR << "caught exception" << e.msg() << std::endl;
		return -1;
	}


}
