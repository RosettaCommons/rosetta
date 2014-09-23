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
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/excn/Exceptions.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
using namespace ObjexxFCL::format;

using basic::T;
using basic::Warning;
using basic::Error;

using core::Size;
using core::Real;
using utility::vector1;

int
main( int argc, char * argv [] )
{
	try {

	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	devel::init( argc, argv );

	// setup residue types
	core::chemical::ResidueTypeSetCOP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	// read in a native pose
	core::pose::Pose pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::native ]() );
	}

	// set up score function
	core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );

	// read in constraints
	ConstraintSetOP cstset;
	std::string cstfile = core::scoring::constraints::get_cst_file_option();
	cstset = ConstraintIO::read_constraints( cstfile, new ConstraintSet, pose );

	core::scoring::EnergyMap emap;
	emap.zero();

	std::string outfile_prefix = option[ in::file::native ]();
	std::string outfile = outfile_prefix + ".cst_quality";
	if ( option[ out::file::silent ].user() ) {
		outfile = option[ out::file::silent ]();
	}
	std::ofstream output( outfile.c_str() );
	if ( ! output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
	}

	static std::string atom_name( "CA" );
	int output_width = 12;
	int precision = 4;

	bool verbose( false );
	if ( option[ out::level ]() > 300 )	verbose = true;

	if ( verbose ) {
		output
			<< A( output_width, "resi" )
			<< A( output_width, "resj" )
			<< A( output_width, "dist" )
			<< A( output_width, "score" )
			<< std::endl;
	} else {
		output
			//<< A( output_width, "atomi" )
			<< A( output_width, "resi" )
			//<< A( output_width, "atomj" )
			<< A( output_width, "resj" )
			<< A( output_width, "dist" )
			<< A( output_width, "score" )
			<< A( output_width, "dist_min" )
			<< A( output_width, "score_min" )
			<< A( output_width, "error" )
			//<< A( output_width, "kl_div" )
			<< std::endl;
	}

	using core::Size;
	using core::Real;
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		for ( Size j = i + 1; j <= pose.total_residue(); ++j ) {
			Real local_score( 0.0 );
			cstset->residue_pair_energy(
				pose.residue(i), pose.residue(j), pose, *scorefxn, emap
			);
			local_score = emap[ core::scoring::atom_pair_constraint ];
			Real const dist(
				pose.residue(i).xyz(atom_name).distance( pose.residue(j).xyz(atom_name) )
			);

			// find score_min and dist_min
			Real dist_min = dist, score_min = local_score;

			if ( local_score == 0 ) continue;

			Real kl_divergence( 0.0 );
			vector1< Real > scores, distances;
			Real stepsize = 0.1;
			Real const lower_dist( 2 );
			//Real upper_dist = std::max( dist + 2 * stepsize, 16.0 );
			Real const upper_dist( 20 );
			for ( Real r = lower_dist; r <= upper_dist; r = r + stepsize ) {
				emap.zero();
				core::conformation::Residue resi( pose.residue(i) );
				core::conformation::Residue resj( pose.residue(j) );

				resi.atom("CA").xyz( core::Vector(0,0,0) );
				resj.atom("CA").xyz( core::Vector(r,0,0) );

				cstset->residue_pair_energy( resi, resj, pose, *scorefxn, emap );
				Real const temp_score( emap[ core::scoring::atom_pair_constraint ] );
				scores.   push_back( temp_score );
				distances.push_back( r          );

				Real const prob( exp( -1 * temp_score ) );
				/// might be < 0 due to rounding errors.
				kl_divergence += std::max(
					prob * std::log( prob / dgaussian( r, 18.991, 7.353, 1.0 ) ),
					0.0
				);

				if ( temp_score < score_min ) {
					dist_min  = r;
					score_min = temp_score;
				}

				if ( std::isinf( temp_score ) ) {
					std::cout << "found an infinite value at " << r << std::endl;
				}
			} // for ( Real r = 2; r <= 16; r = r + 0.1 )

			if ( score_min == 0.0 ) {
				//std::cerr << "skipping pair (" << i << "," << j << ")" << std::endl;
				continue;
			}

			if ( verbose ) {
				for ( core::Size k = 1; k <= scores.size(); ++k ) {
					Real distance   = distances[k];
					Real score      = scores[k];

					output << I( output_width, i )
						<< I( output_width, j )
						<< I( output_width, distance )
						<< F( output_width, precision, score )
						<< std::endl;
				} // for ( k )
			} // if ( verbose )

			if ( !verbose ) {
				output
					<< I( output_width, i )
					<< I( output_width, j )
					//<< "CA " << i << " CA " << j
					<< F( output_width, precision, dist )
					<< F( output_width, precision, local_score )
					<< F( output_width, precision, dist_min )
					<< F( output_width, precision, score_min )
					<< F( output_width, precision, std::abs( dist_min - dist ) )
					//<< F( output_width, precision, kl_divergence )
					<< std::endl;
			} // for ( r )

			emap.zero();
		} // for j
	} // for i

	output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
