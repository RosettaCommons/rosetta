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

#include <protocols/viewer/viewers.hh>

// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/rigid_body_moves.hh>

// AUTO-REMOVED #include <protocols/relax_protocols.hh>


#include <core/types.hh>

// AUTO-REMOVED #include <basic/prof.hh> // profiling

// AUTO-REMOVED #include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/MixtureFunc.hh>

#include <core/scoring/constraints/util.hh>
#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Auto-header: duplicate removed #include <basic/Tracer.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
// Auto-header: duplicate removed
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/constraints/ConstraintIO.hh>

#include <apps/pilot/james/james_util.hh>

using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace ObjexxFCL::format;

// usage: constraints_matrix [-in::file::fullatom] -in::file::native pdbfile -constraints::cst_file cstfile

int
main( int argc, char * argv [] )
{
	try {

	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	devel::init( argc, argv );

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
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
	std::string outfile = outfile_prefix + ".cst_matrix";
	std::ofstream output( outfile.c_str() );
	if ( ! output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
	}

	outfile = outfile_prefix + ".cst_pair_energies";
	std::ofstream res_pair_output( outfile.c_str() );

	// iterate over residue pairs in Pose, print constraint scores if ( i < j ), print CA-CA distance
	// if ( i > j ), and print zeroes along the diagonal. Normalize each of the scores to a mean of
	// zero and a standard deviation of one.
	static std::string atom_name( "CA" );
	utility::vector1< core::Real > res_energies;
	utility::vector1< core::Real > distances;
	utility::vector1< core::Real > scores;
	ObjexxFCL::FArray2D< core::Real > matrix( pose.total_residue(), pose.total_residue(), 0.0f );
	// core::Real max_dist = 10;

	int output_width = 10;
	res_pair_output << A( output_width, "resi" ) << A( output_width, "resj" )
		<< A( output_width, "dist" ) << A( output_width, "score" )
		<< A( output_width, "dist_min" ) << A( output_width, "score_min" )
		<< std::endl;

	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		core::Real residue_energy_total = 0;
		for ( core::Size j = i+1; j <= pose.total_residue(); ++j ) {
			core::Real local_score = 0;
			cstset->residue_pair_energy( pose.residue(i), pose.residue(j), pose, *scorefxn, emap );
			local_score = emap[ core::scoring::atom_pair_constraint ];
			residue_energy_total += local_score;
			core::Real dist = pose.residue(i).xyz(atom_name).distance( pose.residue(j).xyz(atom_name) );

			scores.push_back( local_score );
			distances.push_back( dist );

			// find score_min and dist_min
			core::Real dist_min = dist, score_min = local_score;
			for ( core::Real r = 2; r <= 14; r = r + 0.1 ) {
				emap.zero();
				core::conformation::Residue resi( pose.residue(i) );
				core::conformation::Residue resj( pose.residue(j) );

				resi.atom("CA").xyz( core::Vector(0,0,0) );
				resj.atom("CA").xyz( core::Vector(r,0,0) );

				cstset->residue_pair_energy( resi, resj, pose, *scorefxn, emap );
				core::Real temp_score = emap[ core::scoring::atom_pair_constraint ];
				if ( temp_score < score_min ) {
					dist_min = r;
					score_min = temp_score;
				}
			} // for ( core::Real r = 2; r <= 30; r = r + 0.05 )

			res_pair_output << I( output_width, i )
				<< I( output_width, j )
				<< F( output_width, 3, dist )
				<< F( output_width, 3, local_score )
				<< F( output_width, 3, dist_min )
				<< F( output_width, 3, score_min )
				<< std::endl;

			//matrix( i, j ) = local_score;
			//matrix( j, i ) = std::abs( dist - dist_min );
			matrix( i, j ) = dist;
			matrix( j, i ) = dist;
			//matrix( j, i ) = -1;
			emap.zero();
		} // for j

		res_energies.push_back( residue_energy_total );
	} // for i

	// print out normalized scores
	//core::Real mean_score = mean( scores );
	//core::Real sd_score   = sd  ( scores );
	//core::Real mean_dist  = mean( distances );
	//core::Real sd_dist    = sd  ( distances );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		for ( core::Size j = 1; j <= pose.total_residue(); ++j ) {
			core::Real normalized_score = matrix( i, j );
			//if ( i < j ) {
			//	normalized_score = ( normalized_score - mean_score ) / sd_score;
			//} else if ( i > j ) {
			//	normalized_score = ( normalized_score - mean_dist	) / sd_dist;
			//}

			if ( normalized_score == 0 ) normalized_score = -1;
			//if ( i < j ) normalized_score = 10;
			output << ' ' << normalized_score;
		} // for j
		output << "\n";
	} // for i
	output.close();

	// print out the per-residue energies
	outfile = outfile_prefix + ".residue_energies";
	output.open( outfile.c_str() );
	if ( ! output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
	}

	// print out residue pair energies
	for ( core::Size i = 1; i <= res_energies.size(); ++i ) {
		output << i << " " << res_energies[i] << "\n";
	}

	output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
