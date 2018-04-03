// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
// /// @file burial_measure_centroid.cc
// /// @brief Determines the number of residue neighbors within a sphere around each residue by
//       measuring the distance between the target residues CEN and all other residue's CEN.
// /// @author Melanie Aprahamian (aprahamian.4@osu.edu)

#include <iostream>
#include <string>
#include <protocols/jd2/util.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/chemical/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

#include <numeric/NumericTraits.hh>

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace pose;

OPT_KEY( Real, dist_midpoint )
OPT_KEY( Real, dist_steepness )

static basic::Tracer TR( "apps.pilot.maprahamian.burial_measure_centroid" );

// Main program
int main( int argc, char * argv [] ){
	try {

		NEW_OPT( dist_midpoint, "midpoint of distance falloff", 9.0);
		NEW_OPT( dist_steepness, "exponent factor for distance falloff", 0.1);

		// Initialize Rosetta
		devel::init( argc, argv );

		// Namespaces
		using namespace core;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace pose;

		// check if input PDB file is provided from command line with option -s
		if ( !option[in::file::s].user() ) {
			// exit if no PDB file is found
			utility_exit_with_message("Input PDB file not found, please use the -s flag");
		}

		core::Real distance_fn_midpoint = option[ dist_midpoint ];
		core::Real distance_fn_steepness = option[ dist_steepness ];

		// read in and import pose
		pose::Pose pose;
		std::string pdb_file = option[in::file::s]()[1];
		core::import_pose::centroid_pose_from_pdb( pose, pdb_file, core::import_pose::PDB_file);

		// iterate over all residues to determine neighbor counts
		core::Size const num_residues (pose.total_residue());

		for ( core::Size res_count_target = 1; res_count_target <= num_residues; res_count_target++ ) {
			core::Real neighbor_count (0.0);
			for ( core::Size res_count_neighbor = 1; res_count_neighbor <= num_residues; res_count_neighbor++ ) {
				if ( pose.residue(res_count_target).seqpos() != pose.residue(res_count_neighbor).seqpos() ) {
					core::Real const distance = pose.residue(res_count_neighbor).xyz("CEN").distance(pose.residue(res_count_target).xyz("CEN"));
					neighbor_count += 1.0/(1.0 + std::exp(distance_fn_steepness*(distance-distance_fn_midpoint)));
				}
			}
			TR << pose.residue(res_count_target).type().name1() << " " << pose.residue(res_count_target).seqpos() << " " << neighbor_count << std::endl;
		}
	}
catch ( utility::excn::Exception const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
} //main

