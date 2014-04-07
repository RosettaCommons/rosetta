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
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#include <ObjexxFCL/string.functions.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

using core::Size;
using core::Real;
using utility::vector1;

// C++ headers
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/jobdist/Jobs.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

class ConstraintStatsMover;
typedef utility::pointer::owning_ptr< ConstraintStatsMover > ConstraintStatsMoverOP;

class ConstraintStatsMover : public protocols::moves::Mover {
public:
	ConstraintStatsMover(
		core::scoring::constraints::ConstraintSetOP cst_set,
		core::pose::Pose & native_pose
	) :
		Mover("ConstraintStats"),
		cstset_(cst_set)
	{
		calc_native_stats( native_pose );
	}

	void print_stats( std::ostream & out ) {
		using std::abs;
		using core::Size;
		using core::Real;
		using utility::vector1;

		Size const int_width( 6 );
		Size const width( 12 );
		Size const precision( 3 );
		out << A( int_width, "resi" )   << A( int_width, "resj" )
			<< A( width, "cst_min" ) << A( width, "dist_min" )
			<< A( width, "nat_sc" )  << A( width, "nat_dist" )
			<< A( width, "n_sat" ) << A( width, "n_viol" )
			<< A( width, "dist_diff" )
			<< std::endl;

		for ( Size resi = 1; resi <= cst_mins_.size(); ++resi ) {
			for ( Size resj = resi + 1; resj <= cst_mins_.size(); ++resj ) {
				if ( cst_mins_[resi][resj] == 0.0 ) continue;
				out << I( int_width, resi ) << I( int_width, resj )
					<< F( width, precision, cst_mins_ [resi][resj] )
					<< F( width, precision, dist_mins_[resi][resj] )
					<< F( width, precision, native_cst_scores_[resi][resj] )
					<< F( width, precision, native_distances_ [resi][resj] )
					<< F( width, precision, decoy_satisfactions_[resi][resj] )
					<< F( width, precision, decoy_violations_[resi][resj] )
					<< F(
						width, precision,
						abs( native_distances_ [resi][resj] - dist_mins_[resi][resj] )
					)
					<< std::endl;
			}
		}
	} // print_stats

	void calc_stats(
		core::pose::Pose & pose
	) {
		vector1< Real > dummy( pose.total_residue(), 0.0 );
		//cst_scores.resize( pose.total_residue(), dummy );
		//distances .resize( pose.total_residue(), dummy );
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		core::scoring::EnergyMap emap;
		emap.zero();
		static std::string const atom_name( "CA" );
		//static core::Real const cst_threshold_( 2.5 );
		//static core::Real const cst_threshold_( 0.69 ); // ~ a probability of 0.5

		static core::Real dist_threshold( 1.0 );

		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			for ( core::Size j = i + 1; j <= pose.total_residue(); ++j ) {
				//Real const cst_threshold_( native_cst_scores_[i][j] * 2 );
				emap.zero();
				core::conformation::Residue resi( pose.residue(i) );
				core::conformation::Residue resj( pose.residue(j) );

				//core::Real local_score = 0;
				core::Real dist = pose.residue(i).xyz(atom_name).distance(
					pose.residue(j).xyz(atom_name)
				);

				if ( std::abs( dist - native_distances_[i][j] ) < dist_threshold ) {
					decoy_satisfactions_[i][j]++;
				} else {
					decoy_violations_[i][j]++;
				}

				//cstset_->residue_pair_energy(
				//	resi, resj, pose, *scorefxn, emap
				//);
				//local_score = emap[ core::scoring::atom_pair_constraint ];
				//if ( local_score == 0 ) continue;

				//distances [i][j] = dist;
				//cst_scores[i][j] = local_score;
				//if ( local_score <= cst_threshold_ ) {
				//	decoy_satisfactions_[i][j]++;
				//} else {
				//	decoy_violations_[i][j]++;
				//}
			} // j
		} // i
	} // calc_stats

	void calc_native_stats(
		core::pose::Pose & native_pose
	) {
		using core::Real;
		using utility::vector1;
		vector1< Real > dummy( native_pose.total_residue(), 0.0 );
		cst_mins_ .resize( native_pose.total_residue(), dummy );
		dist_mins_.resize( native_pose.total_residue(), dummy );
		native_cst_scores_.resize( native_pose.total_residue(), dummy );
		native_distances_ .resize( native_pose.total_residue(), dummy );
		decoy_satisfactions_.resize( native_pose.total_residue(), dummy );
		decoy_violations_.resize( native_pose.total_residue(), dummy );

		// intentionally empty score function
		core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
		core::scoring::EnergyMap emap;
		emap.zero();
		static std::string atom_name( "CA" );
		for ( core::Size i = 1; i <= native_pose.total_residue(); ++i ) {
			for ( core::Size j = i + 1; j <= native_pose.total_residue(); ++j ) {
				emap.zero();
				core::conformation::Residue resi( native_pose.residue(i) );
				core::conformation::Residue resj( native_pose.residue(j) );

				Real local_score = 0;
				Real dist = native_pose.residue(i).xyz(atom_name).distance(
					native_pose.residue(j).xyz(atom_name)
				);
				cstset_->residue_pair_energy(
					resi, resj, native_pose, *scorefxn, emap
				);
				local_score = emap[ core::scoring::atom_pair_constraint ];

				native_distances_ [i][j] = dist;
				native_cst_scores_[i][j] = local_score;
				native_distances_ [j][i] = dist;
				native_cst_scores_[j][i] = local_score;

				// skip places where we have no constraints
				if ( local_score == 0 ) continue;
				// find score_min and dist_min for native structure
				core::Real dist_min = dist, score_min = local_score;

				vector1< Real > distances, probs;
				Real stepsize = 0.01;
				Real lower_dist = 2;
				Real upper_dist = std::max( dist + 2 * stepsize, 16.0 );
				for ( Real r = lower_dist; r <= upper_dist; r = r + stepsize ) {
					emap.zero();

					resi.atom("CA").xyz( core::Vector(0,0,0) );
					resj.atom("CA").xyz( core::Vector(r,0,0) );

					cstset_->residue_pair_energy(
						resi, resj, native_pose, *scorefxn, emap
					);
					Real score = emap[ core::scoring::atom_pair_constraint ];
					if ( score < score_min ) {
						dist_min  = r;
						score_min = score;
					}
					distances.push_back( r );
				} // for ( Real r = 2; r <= 16; r = r + 0.1 )
				cst_mins_ [i][j] = score_min;
				dist_mins_[i][j] = dist_min;
				cst_mins_ [j][i] = score_min;
				dist_mins_[j][i] = dist_min;
			} // residue j
		} // residue i
	} // calc_native_stats

	// pieces of data
	// - constraint score (per decoy)
	// - distance (per decoy)
	// - min constraint score (constraints)
	// - min distance (constraints)
	// - native constraint score (native)
	// - native distance (native)
	void apply( core::pose::Pose & pose ) {
		pose.constraint_set( cstset_ );
		calc_stats( pose );
	} // apply

private:
	// native data
	utility::vector1< utility::vector1< core::Real > > native_cst_scores_;
	utility::vector1< utility::vector1< core::Real > > native_distances_;

	// constraint data
	utility::vector1< utility::vector1< core::Real > > cst_mins_;
	utility::vector1< utility::vector1< core::Real > > dist_mins_;

	// number of decoys that "satisfy" each constraint
	utility::vector1< utility::vector1< core::Real > > decoy_satisfactions_;
	utility::vector1< utility::vector1< core::Real > > decoy_violations_;

	core::scoring::constraints::ConstraintSetOP cstset_;
};

int
main( int argc, char * argv [] )
{
	try {

	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using namespace protocols::moves;
	using namespace protocols::jobdist;
	using namespace core::import_pose::pose_stream;
	using std::string;

	devel::init( argc, argv );

	// setup residue types
	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	// read in a native pose
	core::pose::Pose native_pose;
	core::import_pose::pose_from_pdb(
		native_pose, *rsd_set, option[ in::file::native ]()
	);
	MetaPoseInputStream input = streams_from_cmd_line();

	// read in constraints
	ConstraintSetOP cstset;
	std::string cstfile = core::scoring::constraints::get_cst_file_option();
	cstset = ConstraintIO::read_constraints( cstfile, new ConstraintSet, native_pose );

	string outfile( option[ out::file::silent ]() );
	std::ofstream output( outfile.c_str() );
	if ( !output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
	}

	MoverOP mover ( new ConstraintStatsMover( cstset, native_pose ) );
	while( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		mover->apply( pose );
	}

	ConstraintStatsMoverOP downcast
		= ConstraintStatsMoverOP( static_cast< ConstraintStatsMover * > ( mover() ) );
	downcast->print_stats( output );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
