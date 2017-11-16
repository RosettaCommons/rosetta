// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/make_rot_lib/MakeRotLibMover.hh
/// @brief Implementation file for MakeRotLibMover
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )
/// @author Andy Watkins ( amw579@nyu.edu )

// unit headers
#include <protocols/make_rot_lib/MakeRotLibMover.hh>
#include <protocols/make_rot_lib/MakeRotLibJob.hh>
#include <protocols/make_rot_lib/MakeRotLibOptionsData.hh>

// protocols headers
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/simple_moves/MinMover.hh>

// core headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/id/TorsionID.hh>

#include <core/kinematics/MoveMap.hh>

// basic headers
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/make_rot_lib.OptionKeys.gen.hh>
#include <basic/options/util.hh>

// c++ headers
#include <cmath>
#include <iomanip>
#include <fstream>

namespace protocols {
namespace make_rot_lib {

static basic::Tracer TR( "protocols.make_rot_lib.MakeRotLibMover" );

/// @brief Default constructor
MakeRotLibMover::MakeRotLibMover() :
	protocols::moves::Mover("MakeRotLibMover")
{
	using namespace core::scoring;

	TR << "Creating score function..." << std::endl;

	scrfxn_ = core::scoring::get_score_function();

}

void
MakeRotLibMover::apply( core::pose::Pose & pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::make_rot_lib;


	TR << "Inside MakeRotLibMover::apply..." << std::endl;

	// get out job object
	MakeRotLibJobOP job( utility::pointer::dynamic_pointer_cast< MakeRotLibJob >( jd2::JobDistributor::get_instance()->current_job() ) );
	if ( ! job ) {
		utility_exit_with_message( "ERROR: Cannot use the MakeRotLibMover without the MakeRotLibJobInputter.");
	}

	// get the shared data from the job
	protocols::make_rot_lib::MakeRotLibOptionsDataOP mrlod( job->get_options_data() );

	// get the unique data from the job
	core::Real omg( job->get_omg() );
	bool semirotameric( job->get_semirotameric() );
	utility::vector1< core::Real > bbs( job->get_bbs() );
	utility::vector1< core::Size > bb_ids( job->get_bb_ids() );
	core::Real eps( job->get_eps() );

	KbT_ = mrlod->get_temperature();

	TR << "Working on Rotamer library for AA NAME: " << mrlod->get_name() << " OMG: " << omg;
	for ( core::Size i = 1; i <= bbs.size(); ++i ) TR << " bb" << i << ": " << bbs[ i ];
	TR << " EPS: " << eps << std::endl;

	//+-----------------------------------------+
	//|         setup          |
	//+-----------------------------------------+

	TR << "Initializing centroids..." << std::endl;
	init_centroids( mrlod->get_centroid_data(), mrlod->get_n_chi(), mrlod->get_semirotameric() );

	TR << "Initializing rotamers..." << std::endl;
	init_rotamers( mrlod->get_chi_data(), mrlod->get_n_centroids(), omg, bbs, bb_ids, eps, mrlod->get_semirotameric() );
	TR << "Minimizing rotamers..." << std::endl;
	minimize_all_rotamers( pose, bbs, bb_ids, mrlod->get_polymer_type() );

	//+-----------------------------------------+
	//|     symmetry operations       |
	//+-----------------------------------------+
	if ( option[ two_fold_symmetry_135_315 ].user() ) {
		utility::vector1< core::Size > tfs_135_315_chi( option[ two_fold_symmetry_135_315 ].value() );
		for ( core::Size i(1); i <= tfs_135_315_chi.size(); ++i ) {
			TR << "Applying two fold symmetry to minimized rotamers accross the 135/315 axis for chi:" << tfs_135_315_chi[i] << std::endl;
			make_two_fold_symmetry_135_315( rotamers_, tfs_135_315_chi[i] );
		}
	} else if ( option[ two_fold_symmetry_0_180 ].user() ) {
		utility::vector1< core::Size > tfs_0_180_chi( option[ two_fold_symmetry_0_180 ].value() );
		for ( core::Size i(1); i <= tfs_0_180_chi.size(); ++i ) {
			TR << "Applying two fold symmetry to minimized rotamers accross the 0/180 axis for chi:" << tfs_0_180_chi[i] << std::endl;
			make_two_fold_symmetry_0_180( rotamers_, tfs_0_180_chi[i] );
		}
	} else if ( option[ three_fold_symmetry_90_210_330 ].user() ) {
		utility::vector1< core::Size > tfs_90_210_330_chi( option[ three_fold_symmetry_90_210_330 ].value() );
		for ( core::Size i(1); i <= tfs_90_210_330_chi.size(); ++i ) {
			TR << "Applying three fold symmetry to minimized rotamers accross the 90/210/330 axes for chi:" << tfs_90_210_330_chi[i] << std::endl;
			make_three_fold_symmetry_90_210_330( rotamers_, tfs_90_210_330_chi[i] );
		}
	}

	//+-----------------------------------------+
	//|       clustering loop       |
	//+-----------------------------------------+

	TR << "Clustering..." << std::endl;


	if ( option[ k_medoids ].user() ) {
		seed_centroids( mrlod->get_semirotameric() );
	}

	bool clusters_changed( true );
	bool centroids_changed( true );
	core::Size num_iter(0);

	do {
		//TR << "CLUSTERS CHANGE: " << clusters_changed << " CENTROID CHANGED: " << centroids_changed << " ITER: " << num_iter << std::endl;
		//print_rot_data_vec( centroids_, TR );
		//TR << std::endl;

		++num_iter;
		calc_all_dist();
		clusters_changed = calc_rotamer_clusters();
		if ( option[ k_medoids ].user() ) {
			centroids_changed = calc_medoids();
		} else {
			centroids_changed = calc_centroids();
		}
	} while ( ( clusters_changed || centroids_changed ) && num_iter <= 499 );

	//+-----------------------------------------+
	//|        finalize         |
	//+-----------------------------------------+

	TR << "Calculating final rotamers..." << std::endl;
	calc_final_rotamers();

	TR << "Calculating final rotamer probabilities..." << std::endl;
	calc_final_rotamer_probs();

	TR << "Calculating standard deviations..." << std::endl;
	calc_standard_deviations( pose, bbs, bb_ids, mrlod->get_polymer_type() );

	if ( option[ two_fold_symmetry_135_315 ].user() ) {
		utility::vector1< core::Size > tfs_135_315_chi( option[ two_fold_symmetry_135_315 ].value() );
		for ( core::Size i(1); i <= tfs_135_315_chi.size(); ++i ) {
			TR << "Applying two fold symmetry to minimized rotamers accross the 135/315 axis for chi:" << tfs_135_315_chi[i] << std::endl;
			make_two_fold_symmetry_135_315( final_rotamers_, tfs_135_315_chi[i] );
		}
	} else if ( option[ two_fold_symmetry_0_180 ].user() ) {
		utility::vector1< core::Size > tfs_0_180_chi( option[ two_fold_symmetry_0_180 ].value() );
		for ( core::Size i(1); i <= tfs_0_180_chi.size(); ++i ) {
			TR << "Applying two fold symmetry to minimized rotamers accross the 0/180 axis for chi:" << tfs_0_180_chi[i] << std::endl;
			make_two_fold_symmetry_0_180( final_rotamers_, tfs_0_180_chi[i] );
		}
	} else if ( option[ three_fold_symmetry_90_210_330 ].user() ) {
		utility::vector1< core::Size > tfs_90_210_330_chi( option[ three_fold_symmetry_90_210_330 ].value() );
		for ( core::Size i(1); i <= tfs_90_210_330_chi.size(); ++i ) {
			TR << "Applying three fold symmetry to minimized rotamers accross the 90/210/330 axes for chi:" << tfs_90_210_330_chi[i] << std::endl;
			make_three_fold_symmetry_90_210_330( final_rotamers_, tfs_90_210_330_chi[i] );
		}
	}

	//+-----------------------------------------+
	//|        logging          |
	//+-----------------------------------------+

	// setup output files
	std::stringstream base_filename;
	base_filename << mrlod->get_name() << "_" << numeric::principal_angle_degrees( omg );
	for ( core::Size out_i = 1; out_i <= bbs.size(); ++out_i ) {
		base_filename << "_" << numeric::principal_angle_degrees( bbs[ out_i ] ); //phi )
	} // << "_" << numeric::principal_angle_degrees( psi )
	base_filename << "_" << numeric::principal_angle_degrees( eps );
	std::string log_filename( base_filename.str() + ".mrllog" ), rotlib_filename( base_filename.str() + ".rotlib" );
	std::ofstream log_out( log_filename.c_str() );
	std::ofstream rotlib_out( rotlib_filename.c_str() );

	TR << "Printing out log file to " << log_filename << "..." << std::endl;
	log_out << "ROTAMERS" << std::endl;
	print_rot_data_vec( rotamers_, bb_ids, log_out );

	log_out << "CENTROIDS" << std::endl;
	print_rot_data_vec( centroids_, bb_ids, log_out );

	log_out << "FINAL ROTAMERS" << std::endl;
	print_rot_data_vec( final_rotamers_, bb_ids, log_out );

	log_out << "AVERAGE CLUSTER CENTROID DISTANCE" << std::endl;
	print_avg_cluster_centroid_dist( log_out );

	TR << "Printing out rotamer library file to " << rotlib_filename << "..." << std::endl;
	if ( semirotameric ) {
		print_dunbrack10_rotlib( omg, bbs, bb_ids, eps, mrlod->get_polymer_type(), rotlib_out );
	} else {
		print_dunbrack02_rotlib( omg, bbs, bb_ids, eps, mrlod->get_polymer_type(), rotlib_out );

		if ( bbs_appropriate_for_definitions ( bbs ) ) {
			std::stringstream def_fn;
			def_fn << mrlod->get_name() << "_definitions.rotlib";
			std::ofstream def_os( def_fn.str().c_str() );
			print_definitions( def_os );
		}
	}
	log_out.close();
	rotlib_out.close();

}

/// @brief returns if you should print definitions from this file
bool
MakeRotLibMover::bbs_appropriate_for_definitions( utility::vector1< core::Real > bbs )
{
	// use neutral phi and psi for alpha
	// use -120 60 -120 for beta
	if ( ( bbs.size() == 3 && bbs[1] == -120 && bbs[2] == 60 && bbs[3] == -120 )
			|| ( bbs.size() == 2 && bbs[1] == -60 && bbs[2] == 60 ) ) {
		return true;

	} else if ( bbs.size() != 2 && bbs.size() != 3 ) {
		for ( Size i = 1; i <= bbs.size(); ++i ) {
			if ( bbs[ i ] != -120 ) return false;
		}
		return true;
	}
	return false;
}

/// @brief Initializes centroid arrays based on data from the MRLOptionsData
void
MakeRotLibMover::init_centroids( CentroidRotNumVecVec const & centroid_data, core::Size num_chi, bool semirotameric )
{
	using namespace core;

	// resize the centroids array to take up just enough space
	RotData rd_init( num_chi, centroid_data.size(), semirotameric );
	centroids_.resize( centroid_data.size(), rd_init );

	Size num_rotameric_chi = num_chi;
	if ( semirotameric ) num_rotameric_chi--;

	// load the centroid data fromt he MRL options data in to the centroids array
	for ( Size i(1); i <= centroid_data.size(); ++i ) {
		for ( Size j(1); j <= num_rotameric_chi; ++j ) {
			centroids_[ i ].set_inp_chi( centroid_data[ i ][ j ].angle, j );
			centroids_[ i ].set_min_chi( centroid_data[ i ][ j ].angle, j );
			centroids_[ i ].set_lib_chi_val( centroid_data[ i ][ j ].rot_num, j );
		}
	}

}

/// @brief For k-medoids, seed centroids to nearest element of rotamers_
void
MakeRotLibMover::seed_centroids( bool semirotameric )
{
	using namespace core;

	Size num_rotameric_chi = rotamers_[1].get_num_chi();
	if ( semirotameric ) num_rotameric_chi--;

	// Assign each centroid to the nearest rotamer or to no rotamer at all
	// Done WITHOUT replacement
	utility::vector1< Size > used_ids;

	for ( Size i(1); i <= centroids_.size(); ++i ) {
		Real best_distance = 360;
		Size best_id = 0;
		for ( Size j(1); j <= rotamers_.size(); ++j ) {
			Real distance = calc_dist( rotamers_[j], centroids_[i] );
			if ( distance < best_distance ) {
				bool id_clean = true;
				for ( Size id_iter = 1; id_iter <= used_ids.size(); ++id_iter ) {
					if ( used_ids[ id_iter ] == j ) {
						id_clean = false;
					}
				}
				if ( id_clean ) {
					best_distance = distance;
					best_id = j;
				}
			}
		}
		if ( best_id != 0 ) {
			centroids_[ i ] = rotamers_[ best_id ];
			used_ids.push_back( best_id );
		}
	}
}

/// @brief Initializes rotamer arrays based on data from the MRLOptionsData
void
MakeRotLibMover::init_rotamers( TorsionRangeVec const & chi_ranges, core::Size num_cluster,
	core::Real omg, utility::vector1< core::Real > bbs, utility::vector1< core::Size > /*bb_ids*/, core::Real eps, bool semirotameric )
{
	using namespace core;
	// chi_ranges is always as big as nchi so if semirotameric we go one lesser
	Size chi_max = semirotameric ? chi_ranges.size() - 1 : chi_ranges.size();

	// determin the number of chi bins for each chi angle and the size of the rotamers array
	utility::vector1< Size > total_chi_num_bins;
	total_chi_num_bins.resize( chi_max, 0 );
	Size num_rotamers( 1 );


	// TODO: check that the chi range % stepsize is zero
	for ( Size i( 1 ); i <= chi_max; ++i ) {
		total_chi_num_bins[ i ] = ( ( chi_ranges[ i ].high - chi_ranges[ i ].low ) / chi_ranges[ i ].step ) + 1;
		num_rotamers *= total_chi_num_bins[ i ];
	}

	// resize the rotamers array to take up just eneough space
	// NOTABLY we have to pass chi_ranges.size here
	RotData rd_init( chi_ranges.size(), bbs.size(), num_cluster, semirotameric );
	rotamers_.resize( num_rotamers, rd_init);

	// set all of the chi angles
	// the logic here is a little complicated but was done this way to work with any number of chi angles
	for ( Size i( 1 ); i <= chi_max; ++i ) {

		// calc all chi(i) values based on high, low, and step
		utility::vector1< Real > chi_values;
		//TR << "setting chi " << i << " from " << chi_ranges[i].low << " to " << chi_ranges[i].high << " by " << chi_ranges[i].step << std::endl;
		for ( Real k( chi_ranges[ i ].low ); k <= chi_ranges[ i ].high; k += chi_ranges[ i ].step ) {
			chi_values.push_back( k );
		}

		// iterator to keep track of where in rotamers we are
		Size rot_iter( 1 );

		// calc number of rotamer combos for larger chi numbers
		Size count_up( 1 );
		if ( i >= chi_max ) {
			count_up = 1;
		} else {
			for ( Size m( i+1 ); m <= chi_max; ++m ) {
				count_up *= total_chi_num_bins[ m ];
			}
		}

		// calc number of rotamer combos for smaller chi numbers
		Size count_down(1);
		if ( i <= 1 ) {
			count_down = 1;
		} else {
			for ( Size m( 1 ); m < i; ++m ) {
				count_down *= total_chi_num_bins[ m ];
			}
		}

		// assign chi values
		for ( Size n = 1; n <= count_down; ++n ) {
			for ( Size j = 1; j <= chi_values.size(); ++j ) {
				for ( Size l = 1; l <= count_up; ++l ) {
					//TR << "About to set input chi " << i << " for rotamer " << rot_iter << " to " << chi_values[ j ] << std::endl;
					rotamers_[ rot_iter ].set_inp_chi( chi_values[ j ], i );
					++rot_iter;
				}
			}
		}

	}

	// assign backbone torsions to rotamers array
	for ( Size i = 1; i <= rotamers_.size(); ++i ) {

		rotamers_[ i ].resize_bbs( bbs.size() );
		rotamers_[ i ].set_omg( omg );

		for ( Size bb_i = 1; bb_i <= bbs.size(); ++bb_i ) {
			rotamers_[ i ].set_bb( bb_i , bbs[ bb_i ] );
		}

		rotamers_[ i ].set_eps( eps );
	}
}

void
MakeRotLibMover::minimize_rotamer( RotData & rd, core::pose::Pose & pose, utility::vector1<core::Real> bbs, utility::vector1<core::Size> bb_ids, MakeRotLibPolymerType polymer_type )
{
	using namespace core;
	using namespace pose;
	using namespace scoring;
	using namespace chemical;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys::make_rot_lib;

	Size resi = option[ use_terminal_residues ].value() ? 2 : 1;

	// We use BB ids that count the first omega as a mainchain torsion so phi is 2 and psi is 3
	// if we have terminal residues there is no extra mainchain torsion so we decrease by one
	Size offset = option[ use_terminal_residues ].value() ? 1 : 0;

	// get number of chi
	Size const nchi( rd.get_num_chi() );

	// create filenames
	std::stringstream file_name;
	file_name << "dipeptide_" << rd.get_omg() << "_";
	for ( Size bb_i(1); bb_i <= bbs.size(); ++bb_i ) file_name << rd.get_bb(bb_i) << "_";
	file_name << rd.get_eps() << "_";
	if ( rd.get_semirotameric() ) {
		for ( Size j( 1 ); j <= nchi; ++j ) file_name << rd.get_inp_chi( j ) << "_";
	} else {
		for ( Size j( 1 ); j <= nchi-1; ++j ) file_name << rd.get_inp_chi( j ) << "_";
	}
	file_name << ".pdb";

	// set phi, psi, omg, eps based on torsion ids
	//TR << "Going to set omega to " << rd.get_omg() << std::endl;
	id::TorsionID omg_id = ( option[ use_terminal_residues ].value() ?
		id::TorsionID( resi, id::BB, 1 ) :
		id::TorsionID( 1,    id::BB, 1 ) );
	pose.set_torsion( omg_id , rd.get_omg() );

	core::Size eps_num = pose.residue( resi ).mainchain_torsions().size();
	//TR << "Going to set epsilon (torsion " << eps_num << ") to " << rd.get_eps() << std::endl;
	id::TorsionID eps_id( resi, id::BB, eps_num );
	pose.set_torsion( eps_id, rd.get_eps() );

	// bb_i is a counter that iterates through the array of backbone dihedral values stored in the RotamerData object
	// We also have an array of backbone dihedral id numbers, which comes from the OptionsData--this is used
	// to interact with the pose. These two arrays share the bb_i counter.
	for ( Size bb_i( 1 ); bb_i <= bbs.size(); ++bb_i ) {
		//TR << "About to set bb id " << bb_ids[ bb_i ] << " to " << rd.get_bb( bb_i ) << std::endl;
		id::TorsionID bb_id( resi, id::BB, bb_ids[ bb_i ] - offset );
		pose.set_torsion( bb_id, rd.get_bb( bb_i ) );
	}

	// set chi angles
	if ( rd.get_semirotameric() ) {
		for ( Size j( 1 ); j <= ( nchi - 1 ); ++j ) {
			pose.set_chi( j, resi, rd.get_inp_chi( j ) );
			TR << pose.conformation().residue( resi ).chi(j) << std::endl;
		}
	} else {
		for ( Size j( 1 ); j <= nchi; ++j ) {
			pose.set_chi( j, resi, rd.get_inp_chi( j ) );
			TR << pose.conformation().residue( resi ).chi(j) << std::endl;
		}
	}
	// setup the movemap
	kinematics::MoveMapOP mvmp( new kinematics::MoveMap );
	mvmp->set_chi( resi, true );
	if ( polymer_type == PEPTIDE ) {
		mvmp->set( omg_id, true );
		mvmp->set( eps_id, true );
	} else if ( polymer_type == PEPTOID ) {
		mvmp->set( eps_id, true );
	} else {
		utility_exit_with_message("MakeRotLib only works with PEPTIDES and PEPTOIDS.");
	}

	// setup the minimize mover to use linmin
	protocols::simple_moves::MinMover mnmvr( mvmp, scrfxn_, "linmin", 0.0001, true );

	// score the pose
	Real orig_ener( ( *scrfxn_ )( pose ) );

	// perform minimization until we are stabilized ( delta 0.001 ) or we reach 100 steps
	// should probably make these magic numbers options
	Real current_ener( orig_ener ), previous_ener;
	Size iter( 0 );
	do {
		previous_ener = current_ener;
		mnmvr.apply( pose );
		current_ener = pose.energies().total_energy();
		iter++;
		if ( iter > 100 ) break;
	} while ( current_ener + 0.001 < previous_ener );

	// score the pose
	Real min_ener( pose.energies().total_energy() );

	TR << "ORIG ENER: " << orig_ener << " MIN ENER: " << min_ener << " found in " << iter << " steps" << std::endl;

	EnergyMap em( pose.energies().residue_total_energies( resi ) );
	rd.set_inter_rep( em[ fa_rep ] );
	rd.set_inter_atr( em[ fa_atr ] );
	rd.set_solvation( em[ fa_sol ] );
	rd.set_intra_rep( em[ mm_lj_intra_rep ] );
	rd.set_intra_atr( em[ mm_lj_intra_atr ] );
	rd.set_twist( em[ mm_twist ] );

	// record min chi, energy
	for ( Size j = 1; j <= nchi; ++j ) {
		rd.set_min_chi( pose.chi( j, resi ), j );
	}

	if ( rd.get_semirotameric() ) {

		// scan rd for energies
		Real min_ener = 100000;

		MakeRotLibJobOP job( utility::pointer::dynamic_pointer_cast< MakeRotLibJob >( jd2::JobDistributor::get_instance()->current_job() ) );
		debug_assert( job );
		MakeRotLibOptionsDataOP mrlod( job->get_options_data() );

		Real step = mrlod->get_chi_data()[ nchi ].step;
		Size i  = 1;
		rd.resize_semi_vectors( ( mrlod->get_chi_data()[ nchi ].high - mrlod->get_chi_data()[ nchi ].low ) / step + 1 );

		for ( Real nchi_setting = mrlod->get_chi_data()[ nchi ].low,
				high = mrlod->get_chi_data()[ nchi ].high;
				nchi_setting <= high;
				nchi_setting += step, ++i ) {

			pose.set_chi( nchi, resi, nchi_setting );
			Real energy = ( ( * scrfxn_ ) ( pose ) );

			if ( energy < min_ener ) min_ener = energy;
			rd.set_semi_energy_dist( i, energy );

		}

		utility::vector1< Real > normalized_ener( 36, 0 );
		Real tot_prob = 0;
		for ( Size i = 1; i <= 36; ++i ) {
			normalized_ener[ i ] = rd.get_semi_energy_dist( i ) - min_ener;
			rd.set_semi_prob_dist( i, exp( ( -normalized_ener[ i ] ) / KbT_ ) );
			tot_prob += rd.get_semi_prob_dist( i );
		}
		for ( Size i = 1; i <= 36; ++i ) {
			rd.set_semi_prob_dist( i, rd.get_semi_prob_dist( i ) / tot_prob );
		}
	}

	// restore minimum chi2 from before the semirotameric scan
	pose.set_chi( nchi, resi, rd.get_min_chi( nchi ) );

	rd.set_energy( min_ener );
	rd.set_min_omg( numeric::nonnegative_principal_angle_degrees( pose.torsion( omg_id ) ) );
	rd.set_min_eps( numeric::nonnegative_principal_angle_degrees( pose.torsion( eps_id ) ) );

	// DEBUG DUMP COORDS
	//TR << "DUMPING COORDS FOR: " << file_name.str() << std::endl;
	//pose.dump_pdb( file_name.str() );
}

void
MakeRotLibMover::minimize_all_rotamers( core::pose::Pose & pose, utility::vector1< core::Real> bbs, utility::vector1< core::Size > bb_ids, MakeRotLibPolymerType polymer_type )
{
	using namespace core;

	TR << "In MakeRotLibMover::minimize_rotamers..." << std::endl;

	// iterate over rotamers array
	for ( Size i( 1 ); i <= rotamers_.size(); ++i ) {
		TR << "Working on rotamer " << i << " of " << rotamers_.size() << std::endl;
		minimize_rotamer( rotamers_[ i ], pose, bbs, bb_ids, polymer_type );
	}

}

void
MakeRotLibMover::calc_all_dist()
{
	using namespace core;

	for ( Size i( 1 ); i <= rotamers_.size(); ++i ) {
		for ( Size j( 1 ); j<= centroids_.size(); ++j ) {
			rotamers_[ i ].set_cen_dist( calc_dist( rotamers_[ i ], centroids_[ j ] ), j );
		}
	}

}

/// @brief Determins closest cluster centroid for all rotamers
bool
MakeRotLibMover::calc_rotamer_clusters()
{
	using namespace core;

	bool clust_change( false );

	// amw: record each rotamer's assigned centroid, update to nearest
	// and return if the update changed assignments
	for ( Size i( 1 ); i <= rotamers_.size(); ++i ) {
		Size oldclust( rotamers_[ i ].get_cluster_num() );
		Size clust( rotamers_[ i ].get_min_cent_dist() );
		rotamers_[ i ].set_cluster_num( clust );
		if ( oldclust != clust ) clust_change = true;
	}
	return clust_change;
}

/// @brief struct that is used in calc_centroids, values init to zero
struct running_average_pair {
	core::Real avg;
	core::Size count;
	running_average_pair() : avg( 0.0 ), count( 0 ) {}
};

bool
MakeRotLibMover::calc_centroids()
{
	using namespace core;
	using namespace utility;

	Size ncluster( centroids_.size() );
	Size nchi( rotamers_[ 1 ].get_num_chi() );

	// stores the running average for each chi of each cluster
	vector1< vector1< running_average_pair > > running_averages;

	// init running_averages to be nclusters by nchi
	//TR << "DEBUG INIT ARRAYS" << std::endl;
	running_averages.resize( ncluster );
	for ( Size i( 1 ); i <= running_averages.size(); ++i ) running_averages[ i ].resize( nchi );

	// calculate average chis for each cluster
	//TR << "DEBUG CALCULATE AVERAGES " << std::endl;
	for ( Size i( 1 ); i <= rotamers_.size(); ++i ) {
		Size cluster_num( rotamers_[ i ].get_cluster_num() );
		//TR << "GIZMO: " << i << " " << cluster_num;
		for ( Size j( 1 ); j <= nchi; ++j ) {
			//TR << " " << j << " " << rotamers_[i].get_min_chi( j ) << " " << running_averages[cluster_num][j].avg << " " << running_averages[cluster_num][j].count << " ";
			calc_running_avg( rotamers_[ i ].get_min_chi( j ), running_averages[ cluster_num ][ j ].avg, running_averages[ cluster_num ][ j ].count );
		}
		//TR << std::endl;
	}

	// assign averages to centroids
	//TR << "DEBUG ASSIGN VALUES TO CENTROIDS" << std::endl;
	// amw: examine the chis of each centroid and return if the average
	// changed for any of them
	bool centroid_change( false );
	for ( Size i( 1 ); i <= ncluster; ++i ) {
		for ( Size j( 1 ); j <= nchi; ++j ) {
			Real old_avg( centroids_[ i ].get_min_chi( j ) );
			Real new_avg( running_averages[ i ][ j ].avg );
			if ( old_avg != new_avg ) centroid_change = true;
			centroids_[ i ].set_min_chi( new_avg, j );
		}
	}

	return centroid_change;
}

bool
MakeRotLibMover::calc_medoids()
{
	using namespace core;
	using namespace utility;

	bool centroid_change( false );
	Size ncluster( centroids_.size() );

	Real score_current_config = 0;
	for ( Size k = 1; k <= rotamers_.size(); ++k ) {
		score_current_config += rotamers_[ k ].get_min_cent_dist();
	}

	for ( Size i( 1 ); i <= ncluster; ++i ) { // for each medoid

		for ( Size j( 1 ); j <= rotamers_.size(); ++j ) { // for each rotamer
			bool is_a_centroid = false;
			for ( Size k = 1; k <= ncluster; ++k ) {
				if ( rotamers_[ j ] == centroids_[ k ] ) {
					is_a_centroid = true;
				}
			}
			if ( is_a_centroid ) continue;

			RotDataVec centroids_new( centroids_ );
			RotDataVec rotamers_new( rotamers_ );
			centroids_new[ i ] = rotamers_new[ j ];

			// score new configuration with rotamers_ and centroids_new
			Real score_new_configuration = 0;
			for ( Size k = 1; k <= rotamers_new.size(); ++k ) {
				Real min_dist = 100000;
				for ( Size l = 1; l <= ncluster; ++l ) {
					Real dist = calc_dist( rotamers_new[ k ], centroids_new[ l ] );
					if ( dist < min_dist ) {
						min_dist = dist;
						rotamers_new[ k ].set_cen_dist( dist, l );
						rotamers_new[ k ].set_cluster_num( l );
					}
				}
				score_new_configuration += min_dist;
			}

			if ( score_new_configuration < score_current_config ) {
				score_current_config = score_new_configuration;
				rotamers_ = rotamers_new;
				centroids_ = centroids_new;
				centroid_change = true;

			}
		}
	}
	return centroid_change;
}

/// @brief calculates a runnig average. Tries to use the closest anglular value ( ie. the average of 350 degrees and 10 degrees should be 0 and not 180 )( ie. the average of 170 degrees and -170 degrees should be 180 and not 0 )
void
MakeRotLibMover::calc_running_avg( core::Real angle_new, core::Real & angle_old, core::Size & count )
{
	using namespace core;

	Real nnpad_angle_new( numeric::nonnegative_principal_angle_degrees( angle_new ) );
	Real pad_angle_new( numeric::principal_angle_degrees( angle_new ) );
	Real nnpad_angle_diff( fabs( angle_old - nnpad_angle_new ) );
	Real pad_angle_diff( fabs( angle_old - pad_angle_new ) );

	Real angle_blah( nnpad_angle_diff <= pad_angle_diff ? nnpad_angle_new : pad_angle_new );

	//TR << " CRA " << angle_new << " " << nnpad_angle_new << " " << pad_angle_new << " " << angle_blah << " " << angle_old << " ";

	angle_old = ( ( angle_old * count ) + angle_blah ) / (count + 1 );
	//TR << angle_old;
	count++;
}

core::Real
MakeRotLibMover::calc_dist( RotData & point1, RotData & point2 )
{
	using namespace core;

	Size nchi( point1.get_num_chi() );
	Real sd( 0 );

	for ( Size i( 1 ); i <= nchi; ++i ) sd += pow( angle_diff( point1.get_min_chi( i ), point2.get_min_chi( i ) ), 2 );

	return sqrt( sd / nchi );
}

core::Real
MakeRotLibMover::angle_diff( core::Real a1, core::Real a2 )
{
	using namespace core;

	Real t1( fabs( a1 - a2 ) );
	Real t2( 360 - t1 );
	return( t1 <= t2 ? t1 : t2 );
}

void
MakeRotLibMover::calc_final_rotamers()
{
	using namespace core;

	// resize the rotamers array to take up just enough space
	Size num_clusters( centroids_.size() );
	Size num_chi( rotamers_[ 1 ].get_num_chi() );

	// get out job object
	MakeRotLibJobOP job( utility::pointer::static_pointer_cast< MakeRotLibJob >( jd2::JobDistributor::get_instance()->current_job() ) );
	debug_assert( job );

	RotData rd_init( num_chi, rotamers_[1].get_num_bbs(), num_clusters, job->get_semirotameric() );
	final_rotamers_.resize( num_clusters, rd_init );

	//sets high energies for later minimization (probably better way to do this)
	for ( Size j( 1 ); j <= num_clusters; ++j ) final_rotamers_[ j ].set_energy( 6000 );

	//finds lowest E rotamer in each cluster
	for ( Size i( 1 ); i <= rotamers_.size(); ++i ) {
		Size cluster_num( rotamers_[ i ].get_cluster_num() );
		if ( rotamers_[ i ].get_energy() <= final_rotamers_[ cluster_num ].get_energy() ) {
			final_rotamers_[ cluster_num ] = rotamers_[ i ];
		}
	}
}

void
MakeRotLibMover::calc_final_rotamer_probs()
{
	using namespace core;
	using namespace utility;

	vector1< Real > prob_temp;
	prob_temp.resize( final_rotamers_.size(), 0 );

	vector1< Real > normalized_ener;
	normalized_ener.resize( final_rotamers_.size(), 0 );

	Real total_prob( 0 );

	// get the minimum energy
	Real min_ener( final_rotamers_[ 1 ].get_energy() ); // seed
	for ( Size i( 1 ); i <= final_rotamers_.size(); ++i ) {
		if ( final_rotamers_[ i ].get_energy() < min_ener ) {
			min_ener = final_rotamers_[ i ].get_energy();
		}
	}

	// calc the normalized energies
	for ( Size i( 1 ); i <= final_rotamers_.size(); ++i ) {
		normalized_ener[ i ] = final_rotamers_[ i ].get_energy() - min_ener;
	}

	// finds probabilities from energy
	for ( Size i( 1 ); i <= final_rotamers_.size(); ++i ) {
		prob_temp[ i ] = exp( ( -normalized_ener[ i ] ) / KbT_ ) ;
		total_prob += prob_temp[ i ];
		TR << "Cluster:\t" << i << "\t Probability=" << prob_temp[ i ] << std::endl;
	}
	TR << "Total Probability= " << total_prob << std::endl;
	TR << "Kb T value used: " << KbT_ << std::endl;

	// normalizes and sets probabilities of final_rotamers
	for ( Size i( 1 ); i <= final_rotamers_.size(); ++i ) {
		final_rotamers_[ i ].set_probability( prob_temp[ i ] / total_prob );
	}

}

void
MakeRotLibMover::calc_standard_deviations( core::pose::Pose & pose, utility::vector1<core::Real> bbs, utility::vector1<core::Size> bb_ids, MakeRotLibPolymerType polymer_type )
{
	using namespace core;

	using namespace basic::options;
	using namespace basic::options::OptionKeys::make_rot_lib;

	// was 1 in the patched case
	Size resi = option[ use_terminal_residues ].value() ? 2 : 1;

	// iterate over final rotamers
	for ( Size i( 1 ); i <= final_rotamers_.size(); ++i ) {
		Size const nchi ( final_rotamers_[ i ].get_semirotameric() ?
			final_rotamers_[ i ].get_num_chi() - 1 : final_rotamers_[ i ].get_num_chi() );

		// minimize the rotamer based on the input chi since we don't keep poses
		minimize_rotamer( final_rotamers_[ i ], pose, bbs, bb_ids, polymer_type );

		//set energy barrier to stop at
		Real cut_off_ener ( final_rotamers_[ i ].get_energy() + 0.5 );

		//test to see if getting right Energy for this rotamer
		Real test_ener ( ( *scrfxn_ )( pose ) );
		Real found_ener ( final_rotamers_[ i ].get_energy() );
		if ( std::abs( test_ener - found_ener ) >= 0.000001 ) {
			TR << "--------WARNING---------"<< std::endl;
			TR << "For Final_Rot: " << i << "\tScored E: "<< test_ener << "\tbut Min Rot E: "<< found_ener <<std::endl;
		}

		//set step size of chi
		Real chi_step ( 0.1 );

		// search around minimized chi's
		for ( Size j( 1 ); j <= nchi; ++j ) {
			Size k( 0 );
			Real new_ener( 0 );
			while ( new_ener < cut_off_ener && ( k * chi_step ) <= 30 ) {
				k = k + 1;
				pose.set_chi( j, resi, numeric::nonnegative_principal_angle_degrees( final_rotamers_[ i ].get_min_chi( j ) + k * chi_step ) );
				Real plus_ener( ( *scrfxn_ )( pose ) );
				pose.set_chi( j, resi, numeric::nonnegative_principal_angle_degrees( final_rotamers_[ i ].get_min_chi( j ) - k * chi_step ) );
				Real neg_ener( ( *scrfxn_ )( pose ) );
				// set new_ener to highest energy in either direction
				new_ener = ( plus_ener >= neg_ener ? plus_ener : neg_ener );
			}
			final_rotamers_[ i ].set_std_dev( k*chi_step, j );

			//reset pose
			pose.set_chi( j, resi, final_rotamers_[ i ].get_min_chi( j ) );
		}
	}
}


/// brief Use to symmetrize if the 2 rotamers are falling near 0 (360) and 90 (270) like the symmetric chi2 of phenylalanine
void
MakeRotLibMover::make_two_fold_symmetry_135_315( RotDataVec & rdv, core::Size chi_num )
{
	using namespace core;

	for ( Size i = 1; i <= rdv.size(); ++i ) {
		Real temp_chi( numeric::nonnegative_principal_angle_degrees( rdv[ i ].get_min_chi( chi_num ) ) );
		if      ( temp_chi >= 135 && temp_chi <= 315 ) { rdv[ i ].set_min_chi( temp_chi - 180, chi_num ); }
	}
}

/// brief Use to symmetrize if the 1 rotamer is falling near 90 (270) like the symmetric chi1 of a n-aryl peptoid
void
MakeRotLibMover::make_two_fold_symmetry_0_180( RotDataVec & rdv, core::Size chi_num )
{
	using namespace core;

	for ( Size i = 1; i <= rdv.size(); ++i ) {
		Real temp_chi( numeric::nonnegative_principal_angle_degrees( rdv[ i ].get_min_chi( chi_num ) ) );
		if      ( temp_chi >= 180 && temp_chi <= 360 ) { rdv[ i ].set_min_chi( temp_chi - 180, chi_num ); }
	}
}

/// brief Use to symmetrize if the 2 rotamers are falling near 0 and 60 like with a tert-butyl group
void
MakeRotLibMover::make_three_fold_symmetry_90_210_330( RotDataVec & rdv, core::Size chi_num )
{
	using namespace core;

	for ( Size i = 1; i <= rdv.size(); ++i ) {
		Real temp_chi( numeric::nonnegative_principal_angle_degrees( rdv[ i ].get_min_chi( chi_num ) ) );
		if  ( temp_chi >  90 && temp_chi <= 210 ) { rdv[ i ].set_min_chi( temp_chi - 120, chi_num ); }
		else if ( temp_chi > 210 && temp_chi <= 330 ) { rdv[ i ].set_min_chi( temp_chi - 240, chi_num ); }
	}
}

void
MakeRotLibMover::print_rot_data_vec( RotDataVec & rdv, utility::vector1< core::Size > bb_ids, std::ostream & os )
{
	using namespace core;

	for ( Size i(1); i <= rdv.size(); ++i ) print_rot_data( rdv[ i ], bb_ids, os );
}

void
MakeRotLibMover::print_rot_data( RotData & rd, utility::vector1< core::Size > bb_ids, std::ostream & os )
{
	using namespace std;
	using namespace core;

	core::Size const n_bb = rd.get_num_bbs();
	for ( core::Size bb_i = 1; bb_i <= n_bb; ++bb_i ) {
		os << "BB" << bb_ids[ bb_i ] <<  ": " << setw(4) << setprecision(2) << fixed << rd.get_bb( bb_i ) << " ";
	}

	os << "OMG: " << setw(8) << setprecision(4) << fixed << rd.get_min_omg() << " "
		<< "ESL: " << setw(8) << setprecision(4) << fixed << rd.get_min_eps() << " "
		<< "PRB: " << setw(4) << setprecision(4) << fixed << rd.get_probability() << " "
		<< "ENR: " << setw(8) << setprecision(4) << fixed << rd.get_energy() << " "
		<< "TWS: " << setw(8) << setprecision(4) << fixed << rd.get_twist() << " "
		<< "RAR: " << setw(8) << setprecision(4) << fixed << rd.get_intra_rep() << " "
		<< "RAA: " << setw(8) << setprecision(4) << fixed << rd.get_intra_atr() << " "
		<< "NCH: " << setw(4) << fixed << rd.get_num_chi() << " "
		<< "CLN: " << setw(4) << fixed << rd.get_cluster_num() << " ";

	Size max = rd.get_num_chi();
	if ( rd.get_semirotameric() ) {
		max--;
	}

	// print inp_chi_ vector
	os << "ICHI:";
	for ( Size i = 1; i <= max; ++i ) {
		os << setw( 7 ) << setprecision( 2 ) << fixed << rd.get_inp_chi( i ) << " ";
	}

	// print min_chi_ vector
	os << "MCHI:";
	for ( Size i = 1; i <= max; ++i ) {
		os << setw( 7 ) << setprecision( 2 ) << fixed << rd.get_min_chi( i ) << " ";
	}

	// print std_dev_ vector
	os << "STDD:";
	for ( Size i = 1; i <= max; ++i ) {
		os << setw( 7 ) << setprecision( 2 ) << fixed << rd.get_std_dev( i ) << " ";
	}

	// print cen_dst_ vector
	os << "CDST:";
	for ( Size i = 1; i <= rd.get_num_clusters(); ++i ) {
		os << setw( 7 ) << setprecision( 2 ) << fixed << rd.get_cen_dist( i ) << " ";
	}


	os << std::endl;
}


core::Real
MakeRotLibMover::print_avg_cluster_centroid_dist( std::ostream & os )
{
	using namespace core;
	using namespace utility;

	Size ncluster( centroids_.size() );

	vector1< Real > total_dist;
	vector1< Size > num_rots;
	total_dist.resize( ncluster, 0 );
	num_rots.resize( ncluster, 0 );


	for ( Size i( 1 ); i <= rotamers_.size(); ++i ) {
		Size clusternum = rotamers_[ i ].get_cluster_num();
		total_dist[ clusternum ] += rotamers_[ i ].get_cen_dist( clusternum );
		num_rots[ clusternum ]++;
	}

	//print out how many rots in each cluster
	for ( Size i( 1 ); i <= ncluster; ++i ) {
		os << "Cluster: "<< i << " contains "<< num_rots[ i ] << " of " << rotamers_.size()<< "\t or %=" << static_cast< Real >( num_rots[ i ] ) / static_cast< Real >( rotamers_.size() ) << std::endl;
	}

	os << "AVG_CLUST_CENT_DST:\t" << std::endl;
	Real avgdist( 0 );
	for ( Size i( 1 ); i <= ncluster; ++i ) {
		os << total_dist[ i ] / num_rots[ i ] << "\t";
		avgdist += ( total_dist[ i ] / num_rots[ i ] );
	}
	os << std::endl;
	os << "AVG_CENT_DST:\t";
	return avgdist / total_dist.size();
}

void
MakeRotLibMover::print_dunbrack02_rotlib( core::Real omg, utility::vector1< core::Real> bbs, utility::vector1<core::Size> /*bb_ids*/, core::Real /*eps*/, MakeRotLibPolymerType polymer_type, std::ostream & os )
{
	using namespace std;
	using namespace core;
	using namespace utility;

	// sanity check
	debug_assert( final_rotamers_.size() >= 1 && centroids_.size() >= 1 && final_rotamers_.size() == centroids_.size() );

	// get nchi and num rots
	Size nchi( final_rotamers_[ 1 ].get_num_chi() );
	Size nrot( final_rotamers_.size() );

	// get aa name
	std::string aa_name3( "UNK" );

	//counts (imaginary and I believe unused)
	Size count( 9999 );

	// sort final_rotamers based on probability adjust centroids too
	for ( Size i = 1; i <= nrot; ++i ) {
		Size max( i );

		for ( Size j = i+1; j <= nrot; ++j ) {
			if ( final_rotamers_[ j ].get_probability() > final_rotamers_[ max ].get_probability() ) {
				max = j;
			}
		}

		RotData temp_rot( final_rotamers_[ i ] );
		RotData temp_cen( centroids_[ i ] );
		final_rotamers_[ i ] = final_rotamers_[ max ];
		centroids_[ i ] = centroids_[ max ];
		final_rotamers_[ max ] = temp_rot;
		centroids_[ max ] = temp_cen;
	}

	// get "chi bin" number assignments from centroids
	vector1< vector1< Size > > chi_bin_nums; //vector1< Size > chi2_bin_nums; vector1< Size > chi3_bin_nums; vector1< Size > chi4_bin_nums;
	chi_bin_nums.resize( 4 ); // keeping to 4 chis for now!
	for ( Size i = 1; i <= 4; ++i ) {
		chi_bin_nums[ i ].resize( nrot, 0 );
		if ( i <= nchi ) {
			for ( Size j = 1; j <= nrot; ++j ) {
				chi_bin_nums[ i ][ j ] = Size( centroids_[ j ].get_lib_chi_val(i) );
			}
		}
	}

	// now print out lines
	for ( Size i = 1; i <= nrot; ++i ) {

		os << setw(3) << aa_name3 << "\t";

		// if peptoid print out the omega column
		if ( polymer_type == PEPTOID ) {
			os << setw(4) << setprecision(0) << fixed << numeric::principal_angle_degrees( omg ) << "\t";
		}

		// print other BB angles in all cases
		for ( core::Size bb_i = 1; bb_i <= bbs.size(); ++bb_i ) {
			os << setw(4) << setprecision(0) << fixed << numeric::principal_angle_degrees( bbs[ bb_i ] ) << " ";
		}

		// print count found (in this case, meaningless) and the bin numbers for these rotamers and the probability
		os << setw(4) << count << "\t\t";
		for ( Size j = 1; j <= 4; ++j ) {
			os << chi_bin_nums[ j ][ i ] << " ";
		}
		os << "   " << setw(8) << setprecision(6) << fixed << final_rotamers_[ i ].get_probability() << "\t";

		// mean chi
		for ( Size j = 1; j <= 4; ++j ) {
			if ( j <= nchi ) os << setw(6) << setprecision(1) << fixed << numeric::principal_angle_degrees( final_rotamers_[ i ].get_min_chi(j) );
			else os << setw(6) << setprecision(1) << fixed << "0.0";
			os << "\t";
		}

		os << "\t";
		// stdev chi
		for ( Size j = 1; j <= 4; ++j ) {
			if ( j <= nchi ) os << setw(6) << setprecision(1) << fixed << final_rotamers_[ i ].get_std_dev(j);
			else os << setw(6) << setprecision(1) << fixed << "0.0";
			os << "\t";
		}

		os << endl;
	}
}

void
MakeRotLibMover::print_definitions( std::ostream & os )
{
	using namespace std;
	using namespace core;
	using namespace utility;

	Size nchi( final_rotamers_[ 1 ].get_num_chi() );
	Size nrot( final_rotamers_.size() );

	vector1< vector1< Size > > chi_bin_nums; //vector1< Size > chi2_bin_nums; vector1< Size > chi3_bin_nums; vector1< Size > chi4_bin_nums;
	chi_bin_nums.resize( 4 ); // keeping to 4 chis for now!
	// keep the biggest chi bin num for n-1 chi
	vector1< Size > biggest_chi_bin_num( nchi-1, 0 );
	for ( Size i = 1; i <= 4; ++i ) {
		chi_bin_nums[ i ].resize( nrot, 0 );
		if ( i <= nchi ) {
			for ( Size j = 1; j <= nrot; ++j ) {
				chi_bin_nums[ i ][ j ] = Size( centroids_[ j ].get_lib_chi_val(i) );

				if ( j < nchi && chi_bin_nums[ i ][ j ] > biggest_chi_bin_num[ j ] ) {
					biggest_chi_bin_num[ j ] = chi_bin_nums[ i ][ j ];
				}

			}
		}
	}

	// Factor for scaling probability is this product
	Real prob_scale = 1;
	for ( Size i = 1; i <= biggest_chi_bin_num.size(); ++i ) {
		prob_scale *= biggest_chi_bin_num[ i ];
	}

	// output ind file
	// probabilities only relative for final chi
	// so scale by multiplying by the number of chi bin nums for the first two
	for ( Size i = 1; i <= nrot; ++i ) {
		for ( Size j = 1; j <= nchi; ++j ) os << setw(4) << chi_bin_nums[ j ][ i ] << " ";
		Real prob = final_rotamers_[ i ].get_probability() * prob_scale;

		os << "   " << setw(8) << setprecision(6) << fixed << prob << "\t";
		os << "   " << setw(8) << setprecision(6) << fixed << ( -std::log( prob ) ) << "\t";

		Real middle = numeric::principal_angle_degrees( final_rotamers_[ i ].get_min_chi( nchi ) );
		Real left = middle - 2 * final_rotamers_[ i ].get_std_dev( nchi );
		Real right = middle + 2 * final_rotamers_[ i ].get_std_dev( nchi );

		os << setw(8) << setprecision(3) << fixed << left;
		os << setw(8) << setprecision(3) << fixed << middle;
		os << setw(8) << setprecision(3) << fixed << right << std::endl;
	}
}

void
MakeRotLibMover::print_dunbrack10_rotlib( core::Real omg, utility::vector1< core::Real > bbs, utility::vector1<core::Size> /*bb_ids*/, core::Real /*eps*/, MakeRotLibPolymerType polymer_type, std::ostream & os )
{
	// MUST BE SEMIROTAMERIC!!!

	using namespace std;
	using namespace core;
	using namespace utility;

	// sanity check
	debug_assert( final_rotamers_.size() >= 1 && centroids_.size() >= 1 && final_rotamers_.size() == centroids_.size() );

	// get nchi and num rots
	// ASSUMES ONLY ONE CONTINUOUS ANGLE
	Size n_rotameric_chi( final_rotamers_[ 1 ].get_num_chi() - 1 );
	Size nrot( final_rotamers_.size() );

	// get aa name
	std::string aa_name3( "UNK" );

	//counts (imaginary and I believe unused)
	Size count( 9999 );

	// sort final_rotamers based on probability adjust centroids too
	for ( Size i = 1; i <= nrot; ++i ) {
		Size max( i );

		for ( Size j = i+1; j <= nrot; ++j ) {
			if ( final_rotamers_[ j ].get_probability() > final_rotamers_[ max ].get_probability() ) {
				max = j;
			}
		}

		RotData temp_rot( final_rotamers_[ i ] );
		RotData temp_cen( centroids_[ i ] );
		final_rotamers_[ i ] = final_rotamers_[ max ];
		centroids_[ i ] = centroids_[ max ];
		final_rotamers_[ max ] = temp_rot;
		centroids_[ max ] = temp_cen;
	}

	// get "chi bin" number assignments from centroids
	vector1< vector1< Size > > chi_bin_nums; //vector1< Size > chi2_bin_nums; vector1< Size > chi3_bin_nums; vector1< Size > chi4_bin_nums;
	chi_bin_nums.resize( 4 ); // keeping to 4 chis for now!
	for ( Size i = 1; i <= 4; ++i ) {
		chi_bin_nums[ i ].resize( nrot, 0 );
		if ( i <= n_rotameric_chi ) {
			for ( Size j = 1; j <= nrot; ++j ) {
				chi_bin_nums[ i ][ j ] = Size( centroids_[ j ].get_lib_chi_val(i) );
			}
		}
	}

	// now print out lines
	for ( Size i = 1; i <= nrot; ++i ) {

		os << setw(3) << aa_name3 << "\t";

		// if peptoid print out the omega column
		if ( polymer_type == PEPTOID ) {
			os << setw(4) << setprecision(0) << fixed << numeric::principal_angle_degrees( omg ) << "\t";
		}

		// print other BB angles in all cases
		for ( core::Size bb_i = 1; bb_i <= bbs.size(); ++bb_i ) {
			os << setw(4) << setprecision(0) << fixed << numeric::principal_angle_degrees( bbs[ bb_i ] ) << " ";
		}

		// print count found (in this case, meaningless) and the bin numbers for these rotamers and the probability
		os << setw(4) << count << "\t    ";
		for ( Size j = 1; j <= n_rotameric_chi; ++j ) os << chi_bin_nums[ j ][ i ] << " ";

		os << "   " << setw(8) << setprecision(6) << fixed << final_rotamers_[ i ].get_probability() << "\t";

		// mean chi
		for ( Size j = 1; j <= n_rotameric_chi; ++j ) {
			if ( j <= n_rotameric_chi ) os << setw(6) << setprecision(1) << fixed << numeric::principal_angle_degrees( final_rotamers_[ i ].get_min_chi(j) );
			else os << setw(6) << setprecision(1) << fixed << "0.0";
			os << "\t";
		}

		os << "\t";
		// stdev chi
		for ( Size j = 1; j <= n_rotameric_chi; ++j ) {
			if ( j <= n_rotameric_chi ) os << setw(6) << setprecision(1) << fixed << final_rotamers_[ i ].get_std_dev(j);
			else os << setw(6) << setprecision(1) << fixed << "0.0";
			os << "\t";
		}

		// final nonrotameric chi distribution
		for ( Size j = 1; j <= 36; ++j ) {
			os << setw(8) << setprecision(6) << final_rotamers_[ i ].get_semi_prob_dist( j ) << "  ";
		}

		os << endl;
	}
}

}//make_rot_lib
}//protocols
